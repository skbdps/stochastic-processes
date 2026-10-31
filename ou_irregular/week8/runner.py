"""Fixed E4 orchestration and explicit, evidence-bound visual completion.

The existing E4 source, inference, durable executor, independent auditor and
saved-record renderer are composed unchanged. ``run`` stops at
``audited_pending_visual_review``. ``seal`` requires a separate record of actual
inspection of every current PNG; no automated image check supplies that judgment.
Validation is the complete seven-cell R=4 design, not a configurable smaller run.
Full E4 is the complete seven-cell R=2000 design. Set OPENBLAS_NUM_THREADS=1 and
OMP_NUM_THREADS=1 before Python and retain the pinned numerical environment.
"""
import argparse
from pathlib import Path

from ..experiments.design import integer
from ..experiments.provenance import canonical, file_digest, utc_now
from ..experiments.storage import atomic_json, run_lock
from . import audit, execution, plots, reporting
from .design import design_for, parse_design


WEEK_FILES = tuple(sorted(set(plots.FILES) | {
    "ou_irregular/week8/runner.py", "docs/week8_runner.md",
}))
CLOSED_STATES = frozenset({"audited_pending_visual_review", "completed"})


def saved_design(output):
    """Parse exactly one of the fixed complete E4 plans, never a requested subset."""
    output = execution._ordinary_path(output)
    value, _ = execution._json_file(output / "resolved_design.json")
    design = parse_design(value)
    if canonical(value) != canonical(design.to_dict()):
        raise ValueError("saved E4 design is not canonical")
    return design


def _root(output, design):
    """Pure source/runtime/identity and completed-cell-accounting verification."""
    prepared = execution.readback(design, output, extra_files=WEEK_FILES)
    if canonical(saved_design(prepared.output).to_dict()) != canonical(design.to_dict()):
        raise ValueError("runner design differs from the saved fixed E4 experiment")
    return prepared.output, prepared.metadata, prepared.expected_metadata_bytes


def _saved_audit(output, design):
    """Verify the full saved audit and its exact inputs without numerical replay."""
    output, metadata, payload = _root(output, design)
    inputs = plots.load_inputs(output, design)
    if inputs.root_bytes != payload:
        raise ValueError("root changed during saved audit verification")
    report, _ = execution._json_file(output / "audit.json")
    return metadata, payload, report


def audit_study(output, *, replay=True):
    """Read-only numerical audit; replay=False still reconstructs every source.

    A no-fit-replay result is not a full-completion certificate. Use verify_run
    for a pure saved-evidence check that generates neither sources nor fits.
    """
    if type(replay) is not bool:
        raise ValueError("replay must be Boolean")
    design = saved_design(output)
    output, metadata, payload = _root(output, design)
    if metadata.get("status") not in {"raw_complete_pending_audit", *CLOSED_STATES}:
        raise ValueError("incomplete or implementation-review-required evidence cannot be certified")
    report = audit.audit_run(output, design, extra_files=WEEK_FILES, replay=replay)
    if (output / "run_metadata.json").read_bytes() != payload:
        raise ValueError("root changed during read-only numerical audit")
    return report


def verify_run(output, *, require_completed=False):
    """Pure closed-stage readback; no sampling, fitting or archive writes."""
    if type(require_completed) is not bool:
        raise ValueError("require_completed must be Boolean")
    design = saved_design(output)
    metadata, payload, report = _saved_audit(output, design)
    if metadata.get("status") not in CLOSED_STATES or (require_completed and metadata["status"] != "completed"):
        raise ValueError("run has not reached the requested audited/completed state")
    figures = plots.verify_figures(output, design, extra_files=WEEK_FILES)
    if (Path(output) / "run_metadata.json").read_bytes() != payload:
        raise ValueError("root changed during closed-archive verification")
    return {"metadata": metadata, "audit": report, "figures": figures}


def postprocess(output):
    """Tables -> full source/numerical audit -> figures -> awaiting visual review.

    Each lower-level stage owns its lock. A root-byte guard spans stage boundaries;
    the audit publication and lifecycle transition acquire their own exclusive
    lock. Existing full audits are verified, not replaced by a convenient rerun.
    Partial figures retain their original renderer state and immutable bytes.
    """
    design = saved_design(output)
    output, metadata, payload = _root(output, design)
    if metadata.get("status") in CLOSED_STATES:
        return verify_run(output)["audit"]
    output = execution._output_path(output)
    if metadata.get("status") != "raw_complete_pending_audit":
        raise ValueError("postprocessing requires complete raw evidence without implementation-review state")
    guard = execution._ManifestGuard(output / "run_metadata.json", payload)
    try:
        if (output / "audit.json").exists():
            # Rendering interruption cannot authorize new sources or another fit
            # merely because the archive has not yet reached the next state.
            _, _, report = _saved_audit(output, design)
        else:
            reporting.write_tables(design, output, extra_files=WEEK_FILES)
            guard.check("e4_runner:after_reporting")
            with run_lock(output):
                guard.check("e4_runner:before_full_audit")
                report = audit_study(output, replay=True)
                guard.check("e4_runner:before_audit_publication")
                atomic_json(output / "audit.json", report, immutable=True)
                guard.check("e4_runner:after_audit_publication")
            _saved_audit(output, design)
        guard.check("e4_runner:before_figures")
        plots.write_figures(output, design, extra_files=WEEK_FILES)
        guard.check("e4_runner:after_figures")
        with run_lock(output):
            guard.check("e4_runner:before_audited_transition")
            current, _, report = _saved_audit(output, design)
            plots.verify_figures(output, design, extra_files=WEEK_FILES)
            guard.check("e4_runner:after_stage_verification")
            current.update(status="audited_pending_visual_review",
                           audit_sha256=file_digest(output / "audit.json"),
                           analysis_finished_utc=utc_now())
            guard.write(current, "e4_runner:audited_transition")
        verify_run(output)
        return report
    except BaseException:
        if not guard.drifted:
            guard.check("e4_runner:postprocessing_failure")
        raise


def run_study(output, *, validation=False, workers=1, partition_size=250,
              resume=False, progress=True):
    """Execute or resume exactly one fixed experiment, retaining failed outcomes.

    No seed, replication-count, estimator, threshold or cell-subset override is
    exposed. Native execution is never entered when resuming a complete raw,
    audited or sealed run. Full execution requires committed dependencies through
    the executor; the protocol additionally requires prior reviewed validation.
    """
    if type(validation) is not bool or type(resume) is not bool or type(progress) is not bool:
        raise ValueError("validation, resume and progress must be Boolean")
    integer(workers, "workers", 1)
    integer(partition_size, "partition size", 1)
    design = design_for(validation=validation)
    output = execution._ordinary_path(output)
    if resume and (output / "run_metadata.json").exists():
        _, metadata, _ = _root(output, design)
        status = metadata.get("status")
        if status in CLOSED_STATES:
            return verify_run(output)["audit"]
        if status == "raw_complete_pending_audit":
            return postprocess(output)
        if status == "implementation_review_required":
            raise ValueError("saved implementation-review state requires documented investigation, not automatic retry")
    execution.run_raw(design, output, extra_files=WEEK_FILES, workers=workers,
                      partition_size=partition_size, resume=resume, progress=progress)
    return postprocess(output)


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


def seal_run(output, review):
    """Seal verified evidence after actual inspection; repeated sealing is read-only.

    The review is an attestation by the named reviewer, not an automated visual
    assessment. The guard and repeated scientific/figure verification prevent
    hashing an altered artifact into a new completion claim after review.
    """
    checked = verify_run(output)
    output = execution._ordinary_path(output)
    _review(review, output, checked["figures"])
    if checked["metadata"]["status"] == "completed":
        saved, _ = execution._json_file(output / "visual_review.json")
        if canonical(saved) != canonical(review):
            raise ValueError("completed visual review cannot be replaced")
        return checked["metadata"]
    output = execution._output_path(output)
    guard = execution._ManifestGuard(output / "run_metadata.json", (canonical(checked["metadata"]) + "\n").encode())
    with run_lock(output):
        try:
            guard.check("e4_runner:before_visual_seal")
            checked = verify_run(output)
            _review(review, output, checked["figures"])
            guard.check("e4_runner:before_review_publication")
            atomic_json(output / "visual_review.json", review, immutable=True)
            artifacts = _seal_artifacts(output)
            checked = verify_run(output)
            _review(review, output, checked["figures"])
            saved, _ = execution._json_file(output / "visual_review.json")
            if canonical(saved) != canonical(review) or _seal_artifacts(output) != artifacts:
                raise ValueError("artifact or review bytes changed during completion verification")
            metadata = checked["metadata"]
            metadata.update(status="completed", completed_utc=utc_now(),
                            visual_review_sha256=file_digest(output / "visual_review.json"),
                            artifact_sha256=artifacts)
            guard.check("e4_runner:before_complete_transition")
            if _seal_artifacts(output) != artifacts:
                raise ValueError("archive changed immediately before completion publication")
            guard.write(metadata, "e4_runner:complete_transition")
            verify_run(output, require_completed=True)
            return metadata
        except BaseException:
            if not guard.drifted:
                guard.check("e4_runner:seal_failure")
            raise


def replot(output, destination):
    """Regenerate from saved audited tables into a separate, nonoverlapping path."""
    design = saved_design(output)
    output, _, payload = _root(output, design)
    before = plots._scientific_snapshot(output)
    result = plots.write_figures(output, design, extra_files=WEEK_FILES, destination=destination)
    if (output / "run_metadata.json").read_bytes() != payload or plots._scientific_snapshot(output) != before:
        raise ValueError("source archive changed during external replot")
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    run = commands.add_parser("run", help="fixed E4 study; completion awaits actual figure review")
    run.add_argument("--output", type=Path, required=True)
    run.add_argument("--validation", action="store_true", help="fixed seven-cell R=4 separate-seed validation")
    run.add_argument("--workers", type=int, default=1)
    run.add_argument("--partition-size", type=int, default=250)
    run.add_argument("--resume", action="store_true")
    run.add_argument("--quiet", action="store_true")
    for name, help_text in (("postprocess", "finish saved tables, full audit and figures without raw execution"),
                            ("verify", "saved audit/figure/seal readback without new sources or fits"),
                            ("audit", "independent repeat numerical audit; no archive writes"),
                            ("seal", "seal only after actual visual inspection"),
                            ("replot", "render saved tables to an external nonoverlapping destination")):
        item = commands.add_parser(name, help=help_text)
        item.add_argument("--run-dir", type=Path, required=True)
        if name == "audit": item.add_argument("--no-fit-replay", action="store_true")
        if name == "verify": item.add_argument("--require-completed", action="store_true")
        if name == "seal": item.add_argument("--review", type=Path, required=True)
        if name == "replot": item.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    if args.command == "run":
        result = run_study(args.output, validation=args.validation, workers=args.workers,
                           partition_size=args.partition_size, resume=args.resume, progress=not args.quiet)
    elif args.command == "postprocess": result = postprocess(args.run_dir)
    elif args.command == "verify": result = verify_run(args.run_dir, require_completed=args.require_completed)
    elif args.command == "audit": result = audit_study(args.run_dir, replay=not args.no_fit_replay)
    elif args.command == "replot": result = replot(args.run_dir, args.output)
    else:
        review, _ = execution._json_file(args.review)
        result = seal_run(args.run_dir, review)
    print(canonical(result), flush=True)
    return result


if __name__ == "__main__":
    main()
