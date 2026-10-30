"""Week 7 fixed-study CLI, guarded postprocessing and explicit visual sealing.

The published simulator, estimators, numerical audits and renderers are composed
unchanged. A successful run ends at ``audited_pending_visual_review``. Only a
separate seal with an actual, complete visual-review record closes the archive.
Set OPENBLAS_NUM_THREADS=1 and OMP_NUM_THREADS=1 before Python; retain the
recorded environment for readback/resumption as well as numerical execution.
"""
import argparse
import csv
import gzip
from pathlib import Path

from ..experiments import audit as ordinary_audit
from ..experiments.design import integer
from ..experiments.provenance import CORE_FILES, canonical, file_digest, utc_now
from ..experiments.storage import atomic_json, run_lock
from ..experiments.workflow import load_saved_design
from . import execution, floor_audit, plots, reporting
from .path_diagnostics import _ManifestGuard
from .selection import floor_design, followup_design


WEEK_FILES = tuple(sorted(set(plots.FILES) | {
    "ou_irregular/week7/runner.py", "docs/week7_runner.md"}))
DEPENDENCIES = set(CORE_FILES) | set(WEEK_FILES)
PAIRED_SOURCE_FIELDS = {"actual_cell_id", "seed_anchor_cell_id",
                        "seed_anchor_floor_fraction", "source_coupling_policy"}
CLOSED_STATES = {"audited_pending_visual_review", "completed"}


def design_for(kind, *, validation_replications=None):
    """Select the dated fixed study, or its separate development seed namespace."""
    if kind not in {"followup", "floor"}:
        raise ValueError("kind must be followup or floor")
    design = followup_design() if kind == "followup" else floor_design()
    if validation_replications is not None:
        integer(validation_replications, "validation replications", 1)
        design = design.validation(validation_replications)
    return execution.require_design(kind, design)


def saved_design(kind, output):
    """Reject duplicate/noncanonical JSON and every unregistered saved setting."""
    value, _ = execution._json_file(Path(output) / "resolved_design.json")
    design = load_saved_design(output)
    if canonical(value) != canonical(design.to_dict()):
        raise ValueError("saved design is not canonical")
    return execution.require_design(kind, design)


def strict_failure_header(output, kind):
    """Require the exact ordered raw schema even when the failure CSV is empty.

    Set equality would conceal duplicated/reordered column names. This new
    boundary check strengthens composition without changing the old auditor.
    """
    if kind not in {"followup", "floor"}:
        raise ValueError("kind must be followup or floor")
    fields = set(floor_audit.RAW_FIELDS)
    if kind == "followup":
        fields -= PAIRED_SOURCE_FIELDS
    with gzip.open(Path(output) / "failed_fits.csv.gz", "rt", encoding="utf-8", newline="") as stream:
        if next(csv.reader(stream), None) != sorted(fields):
            raise ValueError("failure export must have the complete exact ordered raw header")


def _root(kind, output, design):
    """Pure identity verification before any stage or lifecycle write."""
    output = plots._ordinary_path(output)
    execution.require_design(kind, design)
    execution._root_layout(output, design)
    metadata, payload = execution._json_file(output / "run_metadata.json")
    expected = {
        "schema_version": 1, "execution_kind": kind,
        "execution_contract": execution.EXECUTION_CONTRACT,
        "project": "ou-week7-paired-floor" if kind == "floor" else "ou-irregular-experiments",
        "procedure": execution._procedure(kind), "source_policy": execution.source_policy(kind, design),
        "scientific_overrides": ["paired_floor_source"] if kind == "floor" else [],
        "require_frozen": design.scope == "exploratory",
    }
    if any(key not in metadata or canonical(metadata[key]) != canonical(value)
           for key, value in expected.items()):
        raise ValueError("runner root identity differs from the fixed execution contract")
    if canonical(saved_design(kind, output).to_dict()) != canonical(design.to_dict()):
        raise ValueError("runner design differs from saved evidence")
    plots._frozen_provenance(metadata["provenance"], design, required=DEPENDENCIES,
                             exact=DEPENDENCIES, frozen=expected["require_frozen"])
    if (output / "run_metadata.json").read_bytes() != payload:
        raise ValueError("root changed during read-only runner preflight")
    return output, metadata, payload


def _saved_audit(kind, output, design):
    """Verify the original full audit and every bound input without fit replay."""
    output, metadata, payload = _root(kind, output, design)
    strict_failure_header(output, kind)
    prepared = plots._load_run_inputs(kind, output, design)
    if prepared.root_bytes != payload:
        raise ValueError("root changed during saved audit verification")
    audit, _ = execution._json_file(output / "audit.json")
    return metadata, payload, audit


def audit_study(kind, output, *, replay=True):
    """Independent repeat audit; returns a report and never writes the archive.

    Floor replay=False still reconstructs *all* paired source arrays; only the
    fixed endpoint fit replay is omitted. It is not a completion certification.
    Use verify_run for readback with no source generation or numerical replay.
    """
    if type(replay) is not bool:
        raise ValueError("replay must be Boolean")
    design = saved_design(kind, output)
    output, metadata, payload = _root(kind, output, design)
    if metadata.get("status") not in {"raw_complete_pending_audit", *CLOSED_STATES}:
        raise ValueError("incomplete or review-required evidence cannot be certified by the runner")
    strict_failure_header(output, kind)
    report = (floor_audit.audit_run(output, design, extra_files=WEEK_FILES, replay=replay)
              if kind == "floor" else ordinary_audit.audit_run(output, design, replay=replay))
    if (output / "run_metadata.json").read_bytes() != payload:
        raise ValueError("root changed during repeat audit")
    return report


def verify_run(kind, output, *, require_completed=False):
    """Read-only closed-stage verification, with zero new sources or fits."""
    if type(require_completed) is not bool:
        raise ValueError("require_completed must be Boolean")
    design = saved_design(kind, output)
    metadata, payload, audit = _saved_audit(kind, output, design)
    if metadata.get("status") not in CLOSED_STATES or (require_completed and metadata["status"] != "completed"):
        raise ValueError("run has not reached the requested audited/completed state")
    figures = plots.verify_figures(kind, output, design, extra_files=WEEK_FILES)
    if (Path(output) / "run_metadata.json").read_bytes() != payload:
        raise ValueError("root changed during closed archive verification")
    return {"metadata": metadata, "audit": audit, "figures": figures}


def postprocess(kind, output):
    """Compose tables, full audit and figures; advance only after verified success.

    Reporting and rendering own their locks. A guard spans the unlocked stage
    boundaries; the runner locks separately for audit publication and its root
    transition. On a normal exception the existing honest raw-complete state is
    retained with all partial evidence. Root drift preserves an incident and
    prevents automatic mutation; no cached root is written over it.
    """
    design = saved_design(kind, output)
    output, metadata, payload = _root(kind, output, design)
    if metadata.get("status") in CLOSED_STATES:
        return verify_run(kind, output)["audit"]
    output = execution._output_path(output)
    if metadata.get("status") != "raw_complete_pending_audit":
        raise ValueError("postprocessing requires raw-complete evidence without review-required state")
    guard = _ManifestGuard(output / "run_metadata.json", payload)
    try:
        if (output / "audit.json").exists():
            # A partial rendering already has its numerical audit. It cannot
            # enter raw execution/reporting or repeat native fits to replace it.
            _, _, report = _saved_audit(kind, output, design)
        else:
            reporting.write_tables(kind, design, output, extra_files=WEEK_FILES)
            guard.check("runner:after_reporting")
            with run_lock(output):
                guard.check("runner:before_full_audit")
                report = audit_study(kind, output, replay=True)
                guard.check("runner:before_audit_publication")
                atomic_json(output / "audit.json", report, immutable=True)
                guard.check("runner:after_audit_publication")
            _saved_audit(kind, output, design)
        guard.check("runner:before_figures")
        plots.write_figures(kind, output, design, extra_files=WEEK_FILES)
        guard.check("runner:after_figures")
        with run_lock(output):
            guard.check("runner:before_audited_transition")
            current, _, report = _saved_audit(kind, output, design)
            plots.verify_figures(kind, output, design, extra_files=WEEK_FILES)
            guard.check("runner:after_final_stage_verification")
            current.update(status="audited_pending_visual_review",
                           audit_sha256=file_digest(output / "audit.json"), analysis_finished_utc=utc_now())
            guard.write(current, "runner:audited_transition")
        return report
    except BaseException:
        if not guard.drifted:
            guard.check("runner:postprocessing_failure")
        raise


def run_study(kind, output, *, validation_replications=None, workers=1,
              partition_size=500, resume=False, progress=True):
    """Run once, or resume the existing stage without replacing saved evidence."""
    design = design_for(kind, validation_replications=validation_replications)
    integer(workers, "workers", 1)
    integer(partition_size, "partition size", 1)
    if type(resume) is not bool or type(progress) is not bool:
        raise ValueError("resume and progress must be Boolean")
    output = Path(output)
    if resume and (output / "run_metadata.json").exists():
        _, metadata, _ = _root(kind, output, design)
        status = metadata.get("status")
        if status in CLOSED_STATES:
            return verify_run(kind, output)["audit"]
        if status == "raw_complete_pending_audit":
            return postprocess(kind, output)
        if status == "implementation_review_required":
            raise ValueError("saved implementation-review state requires documented disposition, not automatic retry")
    execution.run_raw(kind, design, output, extra_files=WEEK_FILES, workers=workers,
                      partition_size=partition_size, resume=resume, progress=progress)
    return postprocess(kind, output)


def _review(review, output, figures):
    pngs = {name for name in figures["figures"] if name.endswith(".png")}
    reviewed = review.get("reviewed_pngs") if type(review) is dict else None
    if (not pngs or type(review) is not dict or review.get("all_panels_legible") is not True
            or type(reviewed) is not list or any(type(name) is not str for name in reviewed)
            or len(reviewed) != len(set(reviewed)) or set(reviewed) != pngs
            or review.get("figure_manifest_sha256") != file_digest(output / "figures/figure_manifest.json")
            or type(review.get("reviewer")) is not str or not review["reviewer"].strip()):
        raise ValueError("actual visual review must identify every current PNG exactly once and confirm legibility")


def _seal_artifacts(output):
    return {str(path.relative_to(output)): file_digest(path) for path in sorted(output.rglob("*"))
            if path.is_file() and path.name != ".execution.lock" and path != output / "run_metadata.json"}


def seal_run(kind, output, review):
    """Seal a verified archive after actual review; repeated sealing is read-only."""
    checked = verify_run(kind, output)
    output = Path(output)
    _review(review, output, checked["figures"])
    if checked["metadata"]["status"] == "completed":
        saved, _ = execution._json_file(output / "visual_review.json")
        if canonical(saved) != canonical(review):
            raise ValueError("completed visual review cannot be replaced")
        return checked["metadata"]
    output = execution._output_path(output)
    guard = _ManifestGuard(output / "run_metadata.json", (canonical(checked["metadata"]) + "\n").encode())
    with run_lock(output):
        try:
            guard.check("runner:before_visual_seal")
            checked = verify_run(kind, output)
            _review(review, output, checked["figures"])
            guard.check("runner:before_review_publication")
            atomic_json(output / "visual_review.json", review, immutable=True)
            # Root-byte protection does not protect other artifacts. Recheck
            # the saved audit and renderer bindings after review publication;
            # never bless newly changed data merely by hashing its new bytes.
            # Two equal file maps enclose this verification, so the exact map
            # written into the seal is the one whose bindings just passed.
            artifacts = _seal_artifacts(output)
            checked = verify_run(kind, output)
            _review(review, output, checked["figures"])
            saved_review, _ = execution._json_file(output / "visual_review.json")
            if canonical(saved_review) != canonical(review) or _seal_artifacts(output) != artifacts:
                raise ValueError("artifact or visual-review bytes changed during completion verification")
            metadata = checked["metadata"]
            metadata.update(status="completed", completed_utc=utc_now(),
                            visual_review_sha256=file_digest(output / "visual_review.json"))
            metadata["artifact_sha256"] = artifacts
            guard.check("runner:before_complete_transition")
            if _seal_artifacts(output) != artifacts:
                raise ValueError("archive bytes changed immediately before completion publication")
            guard.write(metadata, "runner:complete_transition")
            # A subsequent outside write invalidates this exact seal. Read it
            # back before returning success; never rewrite artifacts to hide it.
            verify_run(kind, output, require_completed=True)
            return metadata
        except BaseException:
            if not guard.drifted:
                guard.check("runner:seal_failure")
            raise


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    run = commands.add_parser("run", help="execute a fixed study; completion awaits actual figure review")
    run.add_argument("--kind", choices=("followup", "floor"), required=True)
    run.add_argument("--output", type=Path, required=True)
    run.add_argument("--validation-replications", type=int)
    run.add_argument("--workers", type=int, default=1)
    run.add_argument("--partition-size", type=int, default=500)
    run.add_argument("--resume", action="store_true")
    run.add_argument("--quiet", action="store_true")
    for name, help_text in (("postprocess", "finish tables/audit/figures without raw execution"),
                            ("verify", "read-only saved audit/figure/seal verification; no source or fit replay"),
                            ("audit", "independent numerical repeat audit; writes no archive evidence"),
                            ("seal", "seal only after actual visual inspection")):
        item = commands.add_parser(name, help=help_text)
        item.add_argument("--kind", choices=("followup", "floor"), required=True)
        item.add_argument("--run-dir", type=Path, required=True)
        if name == "audit":
            item.add_argument("--no-fit-replay", action="store_true",
                              help="floor still reconstructs all sources; result cannot certify full completion")
        if name == "verify":
            item.add_argument("--require-completed", action="store_true")
        if name == "seal":
            item.add_argument("--review", type=Path, required=True)
    args = parser.parse_args(argv)
    if args.command == "run":
        result = run_study(args.kind, args.output, validation_replications=args.validation_replications,
                           workers=args.workers, partition_size=args.partition_size,
                           resume=args.resume, progress=not args.quiet)
    elif args.command == "postprocess":
        result = postprocess(args.kind, args.run_dir)
    elif args.command == "verify":
        result = verify_run(args.kind, args.run_dir, require_completed=args.require_completed)
    elif args.command == "audit":
        result = audit_study(args.kind, args.run_dir, replay=not args.no_fit_replay)
    else:
        review = execution._loads(args.review.read_bytes())
        result = seal_run(args.kind, args.run_dir, review)
    print(canonical(result), flush=True)
    return result


if __name__ == "__main__":
    main()
