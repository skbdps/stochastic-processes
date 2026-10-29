import json
import shutil

import pytest

from ou_irregular.experiments.audit import AUDIT_FILES
from ou_irregular.experiments.execution import run_raw
from ou_irregular.experiments.provenance import file_digest
from ou_irregular.experiments.workflow import analyze_and_audit, load_saved_design, seal_completed_run
from ou_irregular.week5.plots import render_figures
from test_execution import small_design


def test_completion_requires_full_audit_and_current_visual_review(tmp_path):
    design = small_design()
    run_raw(design, tmp_path, extra_files=AUDIT_FILES, progress=False)
    assert load_saved_design(tmp_path) == design
    with pytest.raises(ValueError, match="audit must finish"):
        seal_completed_run(tmp_path, design, {})
    report = analyze_and_audit(tmp_path, design, render_figures)
    assert report["full_completion_audit"]
    metadata = json.loads((tmp_path / "run_metadata.json").read_text())
    assert metadata["status"] == "audited_pending_visual_review"
    with pytest.raises(ValueError, match="visual review"):
        seal_completed_run(tmp_path, design, {})
    manifest = json.loads((tmp_path / "figures/figure_manifest.json").read_text())
    # Synthetic review ONLY for a validation fixture; no claimed human review.
    review = {"reviewer": "test fixture", "all_panels_legible": True,
              "reviewed_pngs": [name for name in manifest["figures"] if name.endswith(".png")],
              "figure_manifest_sha256": file_digest(tmp_path / "figures/figure_manifest.json")}
    complete = seal_completed_run(tmp_path, design, review)
    assert complete["status"] == "completed"
    assert complete["artifact_sha256"]["parameter_summary.csv"] == file_digest(tmp_path / "parameter_summary.csv")
    assert complete["artifact_sha256"]["visual_review.json"] == file_digest(tmp_path / "visual_review.json")
    before = (tmp_path / "run_metadata.json").read_bytes()
    repeated = seal_completed_run(tmp_path, design, review)
    assert repeated["completed_utc"] == complete["completed_utc"]
    assert (tmp_path / "run_metadata.json").read_bytes() == before


@pytest.fixture(scope="module")
def audited_fixture(tmp_path_factory):
    output = tmp_path_factory.mktemp("audited_seal_fixture")
    design = small_design()
    run_raw(design, output, extra_files=AUDIT_FILES, progress=False)
    analyze_and_audit(output, design, render_figures)
    manifest = json.loads((output / "figures/figure_manifest.json").read_text())
    review = {"reviewer": "synthetic test fixture", "all_panels_legible": True,
              "reviewed_pngs": [name for name in manifest["figures"] if name.endswith(".png")],
              "figure_manifest_sha256": file_digest(output / "figures/figure_manifest.json")}
    return output, design, review


@pytest.mark.parametrize("field", ["procedure", "pipeline_errors", "scientific_overrides",
                                   "completed_cells", "code_commit"])
def test_seal_binds_audited_scientific_metadata(tmp_path, audited_fixture, field):
    source, design, review = audited_fixture
    shutil.copytree(source, tmp_path, dirs_exist_ok=True)
    path = tmp_path / "run_metadata.json"
    metadata = json.loads(path.read_text())
    if field == "procedure":
        metadata[field]["level"] = .80
    elif field == "pipeline_errors":
        metadata[field] = 1
    elif field == "scientific_overrides":
        metadata[field] = ["unreviewed_source"]
    elif field == "completed_cells":
        metadata[field] = {}
    else:
        metadata["provenance"][field] = "0"*40
    path.write_text(json.dumps(metadata))
    with pytest.raises(ValueError, match="scientific run metadata changed"):
        seal_completed_run(tmp_path, design, review)
    assert not (tmp_path / "visual_review.json").exists()


@pytest.mark.parametrize("name", ["audit.json", "figures/extra.txt"])
def test_fingerprint_exclusions_apply_only_at_run_root(tmp_path, audited_fixture, name):
    source, design, review = audited_fixture
    shutil.copytree(source, tmp_path, dirs_exist_ok=True)
    extra = tmp_path / "cells" / design.cells[0].identity / name
    extra.parent.mkdir(parents=True, exist_ok=True)
    extra.write_text("unexplained artifact")
    with pytest.raises(ValueError, match="scientific evidence changed"):
        seal_completed_run(tmp_path, design, review)


def test_seal_requires_figure_design_identity(tmp_path, audited_fixture):
    source, design, review = audited_fixture
    shutil.copytree(source, tmp_path, dirs_exist_ok=True)
    path = tmp_path / "figures/figure_manifest.json"
    manifest = json.loads(path.read_text())
    manifest["design_sha256"] = "0"*64
    path.write_text(json.dumps(manifest))
    review = {**review, "figure_manifest_sha256": file_digest(path)}
    with pytest.raises(ValueError, match="figure design identity"):
        seal_completed_run(tmp_path, design, review)


def test_seal_rejects_changed_table_manifest(tmp_path, audited_fixture):
    source, design, review = audited_fixture
    shutil.copytree(source, tmp_path, dirs_exist_ok=True)
    path = tmp_path / "table_manifest.json"
    manifest = json.loads(path.read_text())
    manifest["parameter_summary.csv"]["sha256"] = "0"*64
    path.write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="figure inputs"):
        seal_completed_run(tmp_path, design, review)


def test_resume_after_post_audit_figure_interruption_preserves_audit(tmp_path):
    design = small_design()
    run_raw(design, tmp_path, extra_files=AUDIT_FILES, progress=False)

    def interrupted_writer(*args):
        raise RuntimeError("injected interruption after numerical audit")

    with pytest.raises(RuntimeError, match="injected interruption"):
        analyze_and_audit(tmp_path, design, interrupted_writer)
    audit_before = (tmp_path / "audit.json").read_bytes()
    partitions = {str(p.relative_to(tmp_path)): p.read_bytes()
                  for p in (tmp_path / "cells").rglob("*.jsonl.gz")}
    run_raw(design, tmp_path, extra_files=AUDIT_FILES, resume=True, progress=False)
    analyze_and_audit(tmp_path, design, render_figures)
    assert (tmp_path / "audit.json").read_bytes() == audit_before
    assert all((tmp_path / name).read_bytes() == data for name, data in partitions.items())
    metadata = json.loads((tmp_path / "run_metadata.json").read_text())
    assert metadata["status"] == "audited_pending_visual_review"
    assert len(metadata["execution_attempts"]) == 2
    assert metadata["postprocessing_error_type"] == "RuntimeError"
