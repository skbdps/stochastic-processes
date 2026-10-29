"""Post-execution tables, audit, figures and an explicit visual-review seal."""
import json
from pathlib import Path

from .analysis import write_tables
from .audit import artifact_fingerprint, audit_run, semantic_metadata_digest
from .design import Cell, RunDesign
from .provenance import file_digest, utc_now, verify
from .storage import atomic_json, run_lock


def load_saved_design(run_dir):
    data = json.loads((Path(run_dir) / "resolved_design.json").read_text())
    cells = tuple(Cell(n=c["n"], q=c["theta_mean_gap"], cv=c["cv"], theta=c["true_theta"],
                       mu=c["true_mu"], sigma=c["true_sigma"], floor_fraction=c["floor_fraction"],
                       floor_cv_threshold=c["floor_cv_threshold"]) for c in data["cells"])
    design = RunDesign(data["name"], data["seed"], data["replications"], cells, data["scope"])
    if design.to_dict() != data:
        raise ValueError("saved design is not a complete canonical design")
    return design


def analyze_and_audit(run_dir, design, figure_writer):
    output = Path(run_dir)
    with run_lock(output):
        path = output / "run_metadata.json"
        metadata = json.loads(path.read_text())
        try:
            write_tables(output, design)
            report = audit_run(output, design, replay=True)
            atomic_json(output / "audit.json", report, immutable=True)
            figure_writer(output, design)
            verify(metadata["provenance"], expected_design=design)
            metadata.update(status="audited_pending_visual_review", audit_sha256=file_digest(output / "audit.json"),
                            analysis_finished_utc=utc_now())
            atomic_json(path, metadata)
            return report
        except BaseException as exc:
            metadata.update(status="implementation_review_required", postprocessing_error_type=type(exc).__name__,
                            postprocessing_error=str(exc))
            atomic_json(path, metadata)
            raise


def verify_figures(run_dir, design):
    output = Path(run_dir)
    manifest = json.loads((output / "figures/figure_manifest.json").read_text())
    if manifest["design_sha256"] != design.digest:
        raise ValueError("figure design identity mismatch")
    tables = json.loads((output / "table_manifest.json").read_text())
    for name, digest in manifest["input_table_sha256"].items():
        if tables[name]["sha256"] != digest or file_digest(output / name) != digest:
            raise ValueError("figure inputs no longer match saved tables")
    for name, digest in manifest["figures"].items():
        if Path(name).name != name or file_digest(output / "figures" / name) != digest:
            raise ValueError("figure digest or filename mismatch")
    expected = set(manifest["figures"]) | {"figure_manifest.json"}
    if {p.name for p in (output / "figures").iterdir()} != expected:
        raise ValueError("figure artifact set differs from its manifest")
    return manifest


def seal_completed_run(run_dir, design, review):
    """Seal only after an actual visual review, not merely successful plotting.

    The review is an evidence record prepared by the researcher/reviewer after
    inspecting the PNGs. It is not an external or professor permission gate.
    """
    output = Path(run_dir)
    with run_lock(output):
        path = output / "run_metadata.json"
        metadata = json.loads(path.read_text())
        if metadata["status"] not in ("audited_pending_visual_review", "completed"):
            raise ValueError("raw work and full audit must finish before sealing")
        verify(metadata["provenance"], expected_design=design)
        audit = json.loads((output / "audit.json").read_text())
        if (file_digest(output / "audit.json") != metadata["audit_sha256"]
                or not audit["full_completion_audit"] or audit["status"] != "passed"
                or audit["run_identity_sha256"] != metadata["provenance"]["identity_sha256"]):
            raise ValueError("a full unchanged replay audit is required")
        if semantic_metadata_digest(metadata) != audit["audited_run_metadata_sha256"]:
            raise ValueError("scientific run metadata changed after its audit")
        figures = verify_figures(output, design)
        pngs = {name for name in figures["figures"] if name.endswith(".png")}
        if (not pngs or review.get("all_panels_legible") is not True
                or set(review.get("reviewed_pngs", [])) != pngs
                or review.get("figure_manifest_sha256") != file_digest(output / "figures/figure_manifest.json")
                or not isinstance(review.get("reviewer"), str) or not review["reviewer"].strip()):
            raise ValueError("visual review must identify every current PNG and confirm legibility")
        # Seal the bytes that passed the independent audit without repeating
        # its numerical work. Changes since that audit cannot acquire its seal.
        if artifact_fingerprint(output) != audit["audited_artifact_sha256"]:
            raise ValueError("scientific evidence changed after its audit")
        atomic_json(output / "visual_review.json", review, immutable=True)
        metadata.update(status="completed", completed_utc=metadata.get("completed_utc", utc_now()),
                        visual_review_sha256=file_digest(output / "visual_review.json"))
        metadata["artifact_sha256"] = {
            str(p.relative_to(output)): file_digest(p) for p in sorted(output.rglob("*"))
            if p.is_file() and p != path and p.name != ".execution.lock"}
        atomic_json(path, metadata)
        return metadata
