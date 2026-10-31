"""E4 orchestration tests plus one retained fixed R=4 native CLI invocation.

Synthetic review records exercise sealing only on disposable copies. They are
explicitly not actual visual reviews and never seal the retained native archive.
The native fixture is the complete prespecified development experiment, not E4
headline data. Its original invocation and every output are retained when an
OU_E4_RUNNER_FIXTURE_ROOT is supplied by CI.
"""
from copy import deepcopy
from dataclasses import replace
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

from ou_irregular.experiments.provenance import ROOT, canonical
from ou_irregular.week8 import audit, execution, plots, reporting, runner, source
from ou_irregular.week8.design import full_design, validation_design
from .test_audit import mutate_table


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def snapshot(root):
    root = Path(root)
    if not root.exists(): return {}
    return {str(p.relative_to(root)): ("link", str(p.readlink())) if p.is_symlink()
            else ("directory",) if p.is_dir() else ("file", sha(p))
            for p in sorted(root.rglob("*"))}


def read_json(path):
    return json.loads(Path(path).read_bytes())


def write_json(path, value):
    Path(path).write_text(canonical(value) + "\n", encoding="utf-8")


def forbidden(*args, **kwargs):
    raise AssertionError("closed/derived stage must not execute new sources or fits")


@pytest.fixture(autouse=True)
def no_science(monkeypatch):
    # The retained native CLI uses a fresh subprocess before its output is used.
    # Every in-process orchestration/corruption test must remain source-free.
    monkeypatch.setattr(execution, "run_raw", forbidden)
    monkeypatch.setattr(audit, "audit_run", forbidden)
    monkeypatch.setattr(audit, "independent_source", forbidden)
    monkeypatch.setattr(audit, "reconstruct_schedule", forbidden)
    monkeypatch.setattr(audit, "source_records", forbidden)
    monkeypatch.setattr(source.E4Source, "sample", forbidden)
    monkeypatch.setattr(source.E4Source, "gaps", forbidden)
    monkeypatch.setattr(source, "simulate_ou", forbidden)


@pytest.fixture(scope="module")
def native(tmp_path_factory):
    parent = Path(os.environ["OU_E4_RUNNER_FIXTURE_ROOT"]) if "OU_E4_RUNNER_FIXTURE_ROOT" in os.environ else tmp_path_factory.mktemp("e4_runner_native")
    parent.mkdir(parents=True, exist_ok=True)
    out = parent / "native_archive"
    assert not out.exists(), "native invocation may not overwrite prior evidence"
    command = [sys.executable, "-m", "ou_irregular.week8.runner", "run", "--validation",
               "--workers", "2", "--partition-size", "2", "--output", str(out)]
    write_json(parent / "invocation.json", {"command": command, "cwd": str(ROOT),
               "purpose": "fixed R=4 complete CLI validation; not full E4 results",
               "expected_paths": 28, "expected_fit_records": 84,
               "seed": 21261011, "replications": 4})
    with (parent / "stdout.txt").open("w") as stdout, (parent / "stderr.txt").open("w") as stderr:
        completed = subprocess.run(command, cwd=ROOT, stdout=stdout, stderr=stderr, text=True, check=False)
    write_json(parent / "exit.json", {"returncode": completed.returncode})
    assert completed.returncode == 0, (parent / "stderr.txt").read_text()
    assert read_json(out / "run_metadata.json")["status"] == "audited_pending_visual_review"
    assert not (out / "visual_review.json").exists()
    return out


@pytest.fixture
def archive(native, tmp_path):
    out = tmp_path / "copy"
    shutil.copytree(native, out)
    original_state = plots.renderer_state_path(native)
    target_state = plots.renderer_state_path(out)
    target_state.parent.mkdir()
    shutil.copyfile(original_state, target_state)
    return out


def synthetic_review(out):
    figures = read_json(out / "figures/figure_manifest.json")
    return {"reviewer": "UNIT-TEST SYNTHETIC REVIEW -- NOT ACTUAL VISUAL INSPECTION",
            "notes": ["Disposable test copy only; never use to complete retained validation or E4 data."],
            "all_panels_legible": True,
            "reviewed_pngs": sorted(name for name in figures["figures"] if name.endswith(".png")),
            "figure_manifest_sha256": sha(out / "figures/figure_manifest.json")}


def raw_copy(out, *, remove_audit=False, remove_figures=False):
    metadata = read_json(out / "run_metadata.json")
    metadata["status"] = "raw_complete_pending_audit"
    for name in ("audit_sha256", "analysis_finished_utc", "completed_utc", "visual_review_sha256", "artifact_sha256"):
        metadata.pop(name, None)
    write_json(out / "run_metadata.json", metadata)
    if (out / "visual_review.json").exists(): (out / "visual_review.json").unlink()
    if remove_audit: (out / "audit.json").unlink()
    if remove_figures:
        shutil.rmtree(out / "figures")
        shutil.rmtree(plots.renderer_state_path(out).parent)
    return metadata


def test_runner_declares_complete_previously_reviewed_dependency_closure():
    assert set(runner.WEEK_FILES) == set(plots.FILES) | {
        "ou_irregular/week8/runner.py", "docs/week8_runner.md"}
    assert all((ROOT / name).is_file() for name in runner.WEEK_FILES)


def test_native_cli_retains_all_fixed_paths_fits_audit_and_figures(native):
    before = snapshot(native)
    checked = runner.verify_run(native)
    design = runner.saved_design(native)
    assert canonical(design.to_dict()) == canonical(validation_design().to_dict())
    report = checked["audit"]
    assert report["planned_paths"] == report["source_paths_verified"] == 28
    assert report["planned_fit_records"] == 84
    assert report["cells_verified"] == 7
    assert report["replayed_source_paths"] == 14 and report["replayed_fits"] == 42
    assert report["source_array_hashes_verified"] == 56 and report["schedule_trace_hashes_verified"] == 84
    assert report["full_completion_audit"] is True and report["pipeline_errors"] == 0
    assert set(checked["figures"]["figures"]) == plots.FILENAMES
    assert checked["metadata"]["status"] == "audited_pending_visual_review"
    assert not (native / "visual_review.json").exists()
    assert snapshot(native) == before


@pytest.mark.parametrize("name,value", [("validation", 1), ("validation", None), ("resume", 1),
    ("progress", "yes"), ("workers", 0), ("workers", True), ("partition_size", 0), ("partition_size", 1.5)])
def test_invalid_controls_rejected_before_output(tmp_path, name, value):
    out = tmp_path / "missing"
    with pytest.raises(ValueError): runner.run_study(out, **{name: value})
    assert not out.exists()


@pytest.mark.parametrize("option", [["--seed", "123"], ["--replications", "1"], ["--cv", "2"],
    ["--estimator", "exact"], ["--validation-replications", "1"]])
def test_cli_has_no_scientific_override(tmp_path, option):
    with pytest.raises(SystemExit) as exc:
        runner.main(["run", "--output", str(tmp_path / "missing"), *option])
    assert exc.value.code == 2 and list(tmp_path.iterdir()) == []


def test_run_dispatch_uses_fixed_validation_plan(tmp_path, monkeypatch):
    called = []
    monkeypatch.setattr(execution, "run_raw", lambda design, out, **kw: called.append((design, out, kw)))
    monkeypatch.setattr(runner, "postprocess", lambda out: {"sentinel": "postprocessed"})
    out = tmp_path / "run"
    result = runner.run_study(out, validation=True, workers=3, partition_size=2, progress=False)
    assert result == {"sentinel": "postprocessed"} and len(called) == 1
    design, actual_out, options = called[0]
    assert design.seed == 21261011 and design.replications == 4 and len(design.cells) == 7
    assert actual_out == out and options["workers"] == 3 and options["partition_size"] == 2
    assert options["extra_files"] == runner.WEEK_FILES and options["resume"] is False


def test_full_dispatch_uses_fixed_full_plan_without_executing_it(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(execution, "run_raw", lambda design, out, **kw: calls.append(design))
    monkeypatch.setattr(runner, "postprocess", lambda out: {})
    runner.run_study(tmp_path / "full")
    assert len(calls) == 1 and canonical(calls[0].to_dict()) == canonical(full_design().to_dict())
    assert not (tmp_path / "full").exists()


@pytest.mark.parametrize("change", ["seed", "replications", "cells", "scope"])
def test_saved_design_rejects_scientific_changes(archive, change):
    path = archive / "resolved_design.json"
    value = read_json(path)
    if change == "seed": value[change] += 1
    elif change == "replications": value[change] = 3
    elif change == "cells": value[change] = value[change][:-1]
    else: value[change] = "exploratory"
    write_json(path, value)
    before = snapshot(archive)
    with pytest.raises(ValueError): runner.saved_design(archive)
    assert snapshot(archive) == before


def test_saved_design_rejects_duplicate_json_keys(archive):
    path = archive / "resolved_design.json"
    text = path.read_text(); path.write_text(text.replace('{', '{"seed":21261011,', 1))
    with pytest.raises(ValueError): runner.saved_design(archive)


def test_completed_requirement_and_boolean_types(archive):
    with pytest.raises(ValueError): runner.verify_run(archive, require_completed=True)
    with pytest.raises(ValueError): runner.verify_run(archive, require_completed=1)
    with pytest.raises(ValueError): runner.audit_study(archive, replay=1)


def test_pending_visual_resume_and_postprocess_are_byte_preserving(archive, monkeypatch):
    before = snapshot(archive)
    monkeypatch.setattr(reporting, "write_tables", forbidden)
    monkeypatch.setattr(plots, "write_figures", forbidden)
    expected = read_json(archive / "audit.json")
    assert runner.run_study(archive, validation=True, resume=True) == expected
    assert runner.postprocess(archive) == expected
    assert snapshot(archive) == before


def test_validation_archive_cannot_resume_as_full(archive):
    before = snapshot(archive)
    with pytest.raises(ValueError): runner.run_study(archive, resume=True)
    assert snapshot(archive) == before


def test_existing_output_requires_explicit_resume(archive, monkeypatch):
    # Inspect the executor's pure preflight rather than permitting native work.
    monkeypatch.setattr(execution, "run_raw", lambda design, out, **kw: execution.preflight(
        design, out, extra_files=kw["extra_files"], resume=kw["resume"]))
    before = snapshot(archive)
    with pytest.raises(FileExistsError): runner.run_study(archive, validation=True)
    assert snapshot(archive) == before


@pytest.mark.parametrize("status", ["incomplete", "running_incomplete", "implementation_review_required"])
def test_unfinished_or_review_required_state_cannot_postprocess(archive, status):
    metadata = read_json(archive / "run_metadata.json"); metadata["status"] = status
    write_json(archive / "run_metadata.json", metadata)
    before = snapshot(archive)
    with pytest.raises(ValueError): runner.postprocess(archive)
    assert snapshot(archive) == before


def test_review_required_resume_never_calls_raw(archive):
    metadata = read_json(archive / "run_metadata.json"); metadata["status"] = "implementation_review_required"
    write_json(archive / "run_metadata.json", metadata)
    before = snapshot(archive)
    with pytest.raises(ValueError, match="review"): runner.run_study(archive, validation=True, resume=True)
    assert snapshot(archive) == before


def test_raw_complete_resume_skips_raw_executor_and_reuses_saved_audit(archive, monkeypatch):
    raw_copy(archive)
    monkeypatch.setattr(reporting, "write_tables", forbidden)
    report = runner.run_study(archive, validation=True, resume=True)
    assert report == read_json(archive / "audit.json")
    assert read_json(archive / "run_metadata.json")["status"] == "audited_pending_visual_review"
    assert not (archive / "visual_review.json").exists()


def test_partial_figures_resume_without_new_science_or_tables(archive, monkeypatch):
    raw_copy(archive)
    original = {name: (archive / "figures" / name).read_bytes() for name in plots.FILENAMES}
    (archive / "figures/figure_manifest.json").unlink()
    missing = sorted(plots.FILENAMES)[0]; (archive / "figures" / missing).unlink()
    monkeypatch.setattr(reporting, "write_tables", forbidden)
    runner.postprocess(archive)
    assert all((archive / "figures" / name).read_bytes() == data for name, data in original.items())
    assert runner.verify_run(archive)["metadata"]["status"] == "audited_pending_visual_review"


def test_native_audit_failure_is_retained_without_audit_or_false_closed_state(archive, monkeypatch):
    raw_copy(archive, remove_audit=True, remove_figures=True)
    def fail(*a, **kw): raise ValueError("injected independent audit failure")
    monkeypatch.setattr(runner, "audit_study", fail)
    before = snapshot(archive)
    with pytest.raises(ValueError, match="injected"): runner.postprocess(archive)
    assert snapshot(archive) == before
    assert not (archive / "audit.json").exists() and not (archive / "figures").exists()
    assert read_json(archive / "run_metadata.json")["status"] == "raw_complete_pending_audit"


def test_nonfull_audit_cannot_advance_to_figures_or_completion(archive, monkeypatch):
    report = read_json(archive / "audit.json"); report["full_completion_audit"] = False
    raw_copy(archive, remove_audit=True, remove_figures=True)
    monkeypatch.setattr(runner, "audit_study", lambda *a, **kw: report)
    with pytest.raises(ValueError): runner.postprocess(archive)
    assert read_json(archive / "audit.json")["full_completion_audit"] is False
    assert read_json(archive / "run_metadata.json")["status"] == "raw_complete_pending_audit"
    assert not (archive / "figures").exists()


def test_root_drift_after_audit_write_is_not_overwritten(archive, monkeypatch):
    report = read_json(archive / "audit.json")
    raw_copy(archive, remove_audit=True, remove_figures=True)
    monkeypatch.setattr(runner, "audit_study", lambda *a, **kw: report)
    original = runner.atomic_json
    def inject(path, value, **kw):
        original(path, value, **kw)
        if Path(path).name == "audit.json":
            root = read_json(archive / "run_metadata.json")
            root["status"] = "incomplete"; root["external_writer_marker"] = "fault injection"
            write_json(archive / "run_metadata.json", root)
    monkeypatch.setattr(runner, "atomic_json", inject)
    with pytest.raises(execution.ManifestDriftError): runner.postprocess(archive)
    assert read_json(archive / "run_metadata.json")["external_writer_marker"] == "fault injection"
    assert read_json(archive / "run_metadata.json")["status"] == "incomplete"
    assert (archive / "audit.json").exists() and not (archive / "figures").exists()


@pytest.mark.parametrize("defect", ["missing_png", "duplicate_png", "unknown_png", "wrong_manifest", "false_legibility", "numeric_legibility", "blank_reviewer", "nonlist", "nondict"])
def test_incomplete_or_unbound_visual_review_cannot_write_or_complete(archive, defect):
    review = synthetic_review(archive)
    if defect == "missing_png": review["reviewed_pngs"].pop()
    elif defect == "duplicate_png": review["reviewed_pngs"].append(review["reviewed_pngs"][0])
    elif defect == "unknown_png": review["reviewed_pngs"][0] = "unknown.png"
    elif defect == "wrong_manifest": review["figure_manifest_sha256"] = "0" * 64
    elif defect == "false_legibility": review["all_panels_legible"] = False
    elif defect == "numeric_legibility": review["all_panels_legible"] = 1
    elif defect == "blank_reviewer": review["reviewer"] = " "
    elif defect == "nonlist": review["reviewed_pngs"] = tuple(review["reviewed_pngs"])
    else: review = None
    before = snapshot(archive)
    with pytest.raises(ValueError): runner.seal_run(archive, review)
    assert snapshot(archive) == before


def test_synthetic_seal_on_disposable_copy_and_all_closed_operations_are_pure(archive, monkeypatch):
    review = synthetic_review(archive)
    completed = runner.seal_run(archive, review)
    assert completed["status"] == "completed"
    assert completed["artifact_sha256"] == runner._seal_artifacts(archive)
    assert completed["visual_review_sha256"] == sha(archive / "visual_review.json")
    before = snapshot(archive)
    monkeypatch.setattr(reporting, "write_tables", forbidden)
    monkeypatch.setattr(plots, "write_figures", forbidden)
    assert runner.seal_run(archive, review) == completed
    assert runner.verify_run(archive, require_completed=True)["metadata"] == completed
    runner.run_study(archive, validation=True, resume=True)
    runner.postprocess(archive)
    assert snapshot(archive) == before
    altered = deepcopy(review); altered["reviewer"] += " different"
    with pytest.raises(ValueError): runner.seal_run(archive, altered)
    assert snapshot(archive) == before


@pytest.mark.parametrize("kind", ["table", "figure", "audit", "review"])
def test_artifact_mutation_after_review_publication_cannot_be_blessed(archive, monkeypatch, kind):
    review = synthetic_review(archive); original = runner.atomic_json
    def inject(path, value, **kw):
        original(path, value, **kw)
        if Path(path).name != "visual_review.json": return
        if kind == "table":
            mutate_table(archive, "parameter_summary.csv", lambda cols, rows: rows[0].update(bias="123"))
        elif kind == "figure":
            target = archive / "figures" / sorted(plots.FILENAMES)[0]
            target.write_bytes(target.read_bytes() + b"changed")
        elif kind == "audit":
            target = archive / "audit.json"; report = read_json(target)
            report["full_completion_audit"] = False; write_json(target, report)
        else:
            target = archive / "visual_review.json"; changed = dict(value)
            changed["reviewer"] = "different synthetic reviewer"; write_json(target, changed)
    monkeypatch.setattr(runner, "atomic_json", inject)
    with pytest.raises(ValueError): runner.seal_run(archive, review)
    assert read_json(archive / "run_metadata.json")["status"] == "audited_pending_visual_review"
    assert "artifact_sha256" not in read_json(archive / "run_metadata.json")


def test_artifact_drift_immediately_before_complete_transition_is_rejected(archive, monkeypatch):
    original = runner._seal_artifacts; calls = []
    def inject(out):
        calls.append(1)
        if len(calls) == 3:
            path = out / "source_summary.csv"; path.write_bytes(path.read_bytes() + b"\n")
        return original(out)
    monkeypatch.setattr(runner, "_seal_artifacts", inject)
    with pytest.raises(ValueError): runner.seal_run(archive, synthetic_review(archive))
    assert read_json(archive / "run_metadata.json")["status"] == "audited_pending_visual_review"
    assert len(calls) >= 3


def test_changed_completed_evidence_is_rejected_without_repair(archive):
    runner.seal_run(archive, synthetic_review(archive))
    path = archive / "source_summary.csv"; path.write_bytes(path.read_bytes() + b"\n")
    before = snapshot(archive)
    with pytest.raises(ValueError): runner.verify_run(archive, require_completed=True)
    with pytest.raises(ValueError): runner.run_study(archive, validation=True, resume=True)
    assert snapshot(archive) == before


def test_external_replot_of_sealed_copy_preserves_archive_and_uses_only_saved_evidence(archive, tmp_path):
    runner.seal_run(archive, synthetic_review(archive))
    before = snapshot(archive); destination = tmp_path / "external"
    result = runner.replot(archive, destination)
    assert result["new_sources_or_fits"] is False and set(result["figures"]) == plots.FILENAMES
    assert snapshot(archive) == before
    assert runner.verify_run(archive, require_completed=True)["metadata"]["status"] == "completed"


@pytest.mark.parametrize("relative", ["same", "inside", "ancestor"])
def test_replot_refuses_overlapping_destination(archive, relative):
    target = archive if relative == "same" else archive / "nested" if relative == "inside" else archive.parent
    before = snapshot(archive)
    with pytest.raises(ValueError): runner.replot(archive, target)
    assert snapshot(archive) == before


def test_numerical_audit_passes_explicit_replay_option_and_remains_read_only(archive, monkeypatch):
    seen = []
    def fake(out, design, **kw):
        seen.append((out, design, kw)); return {"sentinel": True}
    monkeypatch.setattr(audit, "audit_run", fake)
    before = snapshot(archive)
    assert runner.audit_study(archive, replay=False) == {"sentinel": True}
    assert seen[0][2] == {"extra_files": runner.WEEK_FILES, "replay": False}
    assert snapshot(archive) == before


def test_cli_routes_fixed_options_and_prints_structured_result(tmp_path, monkeypatch, capsys):
    seen = []
    monkeypatch.setattr(runner, "run_study", lambda out, **kw: seen.append((out, kw)) or {"ok": True})
    result = runner.main(["run", "--validation", "--workers", "2", "--partition-size", "2",
                          "--resume", "--quiet", "--output", str(tmp_path / "run")])
    assert result == {"ok": True} and json.loads(capsys.readouterr().out) == result
    assert seen == [(tmp_path / "run", {"validation": True, "workers": 2, "partition_size": 2, "resume": True, "progress": False})]
