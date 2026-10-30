"""Independent raw-orchestration faults using literal records and real stores.

The historical per-cell engine and CellStore remain in use. Only the point at
which it would simulate and fit a source is replaced by deterministic literal
rows. These tests generate no OU observations or scientific random streams.
"""
from copy import deepcopy
from concurrent.futures import Future
from dataclasses import replace
import gzip
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from ou_irregular.experiments import execution as historical
from ou_irregular.experiments.design import ESTIMATORS
from ou_irregular.experiments.provenance import canonical, digest_bytes
from ou_irregular.experiments.storage import CellStore, read_partition, run_lock
from ou_irregular.week7 import execution
from ou_irregular.week7.floor_source import COUPLING_POLICY, PairedFloorSource
from ou_irregular.week7.selection import floor_design, followup_design


def write_json(path, value):
    Path(path).write_text(canonical(value) + "\n")


def state(path):
    return json.loads((Path(path) / "run_metadata.json").read_text())


def snapshot(path):
    return {str(p.relative_to(path)): p.read_bytes() for p in Path(path).rglob("*") if p.is_file()}


def resign_identity(provenance):
    identity = {key: provenance[key] for key in ("design", "design_sha256", "preregistration_commit",
                                                "source_files", "environment")}
    provenance["identity_sha256"] = digest_bytes(canonical(identity).encode())


def literal_rows(cell, rep, source, services):
    """Identity-correct synthetic outcomes, with no source.sample invocation."""
    anchored = isinstance(source, PairedFloorSource)
    anchor = replace(cell, floor_fraction=1e-6) if anchored else cell
    words = tuple(int(anchor.identity[i:i + 8], 16) for i in range(0, 64, 8))
    address = [4, *words, rep]
    physical = anchor.identity if anchored and cell.cv == 0 else cell.identity
    common = {**cell.to_dict(), "rep": rep, "source_seed": source.seed,
              "source_spawn_key": canonical(address), "gap_spawn_key": canonical([*address, 0]),
              "path_spawn_key": canonical([*address, 1]), "source_valid": True,
              "source_error": "", "pipeline_error": "",
              "times_sha256": hashlib.sha256(f"literal:times:{physical}:{rep}".encode()).hexdigest(),
              "observations_sha256": hashlib.sha256(f"literal:values:{physical}:{rep}".encode()).hexdigest()}
    if anchored:
        common.update(actual_cell_id=cell.identity, seed_anchor_cell_id=anchor.identity,
                      seed_anchor_floor_fraction=1e-6, source_coupling_policy=COUPLING_POLICY)
    rows = []
    for method in ESTIMATORS:
        row = {**common, "estimator": method, "point_valid": True, "point_status": "valid_interior",
               "optimizer_success": True, "interval_valid": True, "interval_status": "valid", "nll": 1.,
               "point_diagnostics_json": '{"fixture":"literal"}', "interval_diagnostics_json": '{"fixture":"literal"}'}
        for parameter in ("theta", "mu", "sigma"):
            truth = row[f"true_{parameter}"]
            row.update({parameter: truth + .01 * rep, f"{parameter}_se": .1,
                        f"{parameter}_lower": truth - .25, f"{parameter}_upper": truth + .25})
        rows.append(row)
    return rows


@pytest.fixture
def harness(tmp_path, monkeypatch):
    calls = []

    def source_records(cell, rep, source, services):
        calls.append((cell.identity, rep, source.seed, type(source).__name__))
        return literal_rows(cell, rep, source, services)

    monkeypatch.setattr(historical, "source_records", source_records)

    def run(kind="followup", design=None, output=None, **kwargs):
        design = design or (followup_design() if kind == "followup" else floor_design()).validation(3)
        output = Path(output) if output is not None else tmp_path / "run"
        return execution.run_raw(kind, design, output, progress=False, **kwargs)

    return SimpleNamespace(root=tmp_path, output=tmp_path / "run", calls=calls,
                           run=run, source_records=source_records, monkeypatch=monkeypatch)


@pytest.mark.parametrize("kind", ["followup", "floor"])
def test_fixed_full_and_exact_validation_designs_only(kind):
    base = followup_design() if kind == "followup" else floor_design()
    execution.require_design(kind, base)
    execution.require_design(kind, base.validation(2))
    for changed in (replace(base, seed=base.seed + 1), replace(base, name=base.name + "_unregistered"),
                    replace(base, replications=base.replications - 1),
                    replace(base, cells=base.cells[:-1]), replace(base, cells=base.cells[::-1]),
                    replace(base, scope="headline"), replace(base.validation(2), seed=base.seed)):
        with pytest.raises(ValueError):
            execution.require_design(kind, changed)


@pytest.mark.parametrize("kind", ["followup", "floor"])
def test_full_run_freeze_failure_precedes_output_creation_and_source_calls(harness, monkeypatch, kind):
    seen = []

    def blocked_capture(design, extra_files=(), **kwargs):
        seen.append(kwargs.get("require_frozen"))
        raise ValueError("synthetic uncommitted dependency")

    monkeypatch.setattr(execution, "capture", blocked_capture)
    design = followup_design() if kind == "followup" else floor_design()
    with pytest.raises(ValueError, match="uncommitted"):
        harness.run(kind, design=design)
    assert seen == [True]
    assert not harness.output.exists()
    assert not harness.calls


@pytest.mark.parametrize("kind", ["followup", "floor"])
def test_raw_complete_resume_and_pure_preflight_are_read_only(harness, kind):
    design = (followup_design() if kind == "followup" else floor_design()).validation(2)
    result = harness.run(kind, design=design, partition_size=1)
    assert result["status"] == "raw_complete_pending_audit"
    assert len(harness.calls) == len(design.cells) * 2
    before = snapshot(harness.output)
    prior_calls = list(harness.calls)
    checked = execution.preflight(kind, design, harness.output, resume=True)
    assert checked.cached_completion is True
    resumed = harness.run(kind, design=design, resume=True, workers=2, partition_size=2)
    assert resumed == result
    assert snapshot(harness.output) == before
    assert harness.calls == prior_calls


@pytest.mark.parametrize("terminal", ["completed", "audited_pending_visual_review"])
def test_closed_analysis_states_reject_raw_resume_without_mutation(harness, terminal):
    harness.run(partition_size=2)
    metadata = state(harness.output)
    metadata["status"] = terminal
    write_json(harness.output / "run_metadata.json", metadata)
    before = snapshot(harness.output)
    calls = list(harness.calls)
    with pytest.raises(ValueError):
        harness.run(resume=True)
    assert snapshot(harness.output) == before and harness.calls == calls


def test_postprocessing_review_with_complete_clean_raw_data_is_cached_and_preserves_error(harness):
    harness.run(partition_size=2)
    metadata = state(harness.output)
    metadata.update(status="implementation_review_required", postprocessing_error_type="SyntheticError",
                    postprocessing_error="retain the original analysis failure")
    write_json(harness.output / "run_metadata.json", metadata)
    before = snapshot(harness.output)
    calls = list(harness.calls)
    result = harness.run(resume=True)
    assert result == metadata
    assert snapshot(harness.output) == before and harness.calls == calls


def test_cached_postprocessing_resume_preserves_partial_png_pdf_without_completion_claim(harness):
    harness.run(partition_size=2)
    metadata = state(harness.output)
    metadata.update(status="implementation_review_required", postprocessing_error_type="SyntheticPlotFailure",
                    postprocessing_error="interrupted before figure manifest creation")
    write_json(harness.output / "run_metadata.json", metadata)
    figures = harness.output / "figures"
    figures.mkdir()
    (figures / "partial_panel.png").write_bytes(b"literal partial PNG fixture; no visual-completion claim")
    (figures / "partial_panel.pdf").write_bytes(b"literal partial PDF fixture; no visual-completion claim")
    before = snapshot(harness.output)
    calls = list(harness.calls)
    checked = execution.preflight("followup", followup_design().validation(3), harness.output, resume=True)
    assert checked.cached_completion is True
    resumed = harness.run(resume=True)
    assert resumed == metadata and resumed["status"] == "implementation_review_required"
    assert not (figures / "figure_manifest.json").exists()
    assert snapshot(harness.output) == before and harness.calls == calls


def test_partial_figure_layout_still_rejects_unknown_extension_and_symbolic_link(harness):
    harness.run(partition_size=2)
    figures = harness.output / "figures"
    figures.mkdir()
    unexpected = figures / "unidentified.txt"
    unexpected.write_text("preserve unidentified partial evidence\n")
    before = snapshot(harness.output)
    calls = list(harness.calls)
    with pytest.raises(ValueError, match="(?i)(figure|artifact)"):
        harness.run(resume=True)
    assert snapshot(harness.output) == before and harness.calls == calls
    unexpected.unlink()
    target = harness.root / "literal_link_target.png"
    target.write_bytes(b"synthetic target outside the figure directory")
    (figures / "linked.png").symlink_to(target)
    before = snapshot(harness.output)
    with pytest.raises(ValueError, match="(?i)(figure|artifact)"):
        harness.run(resume=True)
    assert snapshot(harness.output) == before and harness.calls == calls


def interrupt_after_one_cell(harness, monkeypatch):
    first = followup_design().cells[0].identity

    def interrupt(cell, rep, source, services):
        if cell.identity != first:
            raise KeyboardInterrupt("synthetic interruption after one cell")
        return harness.source_records(cell, rep, source, services)

    monkeypatch.setattr(historical, "source_records", interrupt)
    with pytest.raises(KeyboardInterrupt):
        harness.run(partition_size=1)
    monkeypatch.setattr(historical, "source_records", harness.source_records)
    assert len(harness.calls) == 3


def test_incomplete_postprocessing_review_cannot_claim_cached_completion(harness, monkeypatch):
    interrupt_after_one_cell(harness, monkeypatch)
    metadata = state(harness.output)
    metadata["status"] = "implementation_review_required"
    write_json(harness.output / "run_metadata.json", metadata)
    before = snapshot(harness.output)
    with pytest.raises(ValueError):
        harness.run(resume=True)
    assert snapshot(harness.output) == before and len(harness.calls) == 3


def test_interrupted_run_resumes_only_remaining_cell_sources(harness, monkeypatch):
    interrupt_after_one_cell(harness, monkeypatch)
    assert state(harness.output)["status"] == "incomplete"
    saved = {p: data for p, data in snapshot(harness.output).items() if p.endswith(".jsonl.gz")}
    result = harness.run(resume=True, partition_size=2)
    assert result["status"] == "raw_complete_pending_audit"
    assert len(harness.calls) == 6 and len(set(harness.calls)) == 6
    assert all(snapshot(harness.output)[p] == data for p, data in saved.items())


def test_complete_cell_ahead_of_root_is_adopted_without_regeneration(harness):
    harness.run(partition_size=2)
    metadata = state(harness.output)
    metadata.update(status="incomplete", completed_cells={})
    write_json(harness.output / "run_metadata.json", metadata)
    before_parts = {p: v for p, v in snapshot(harness.output).items() if p.endswith(".jsonl.gz")}
    prior_calls = list(harness.calls)
    result = harness.run(resume=True, partition_size=1)
    assert result["status"] == "raw_complete_pending_audit"
    assert len(result["completed_cells"]) == 2
    assert harness.calls == prior_calls
    assert all(snapshot(harness.output)[p] == v for p, v in before_parts.items())


@pytest.mark.parametrize("change", ["design_bytes", "identity_digest", "dependency_omission", "environment",
                                   "procedure", "root_cell_count", "unknown_root_file", "unknown_cell_file"])
def test_resume_corruption_is_rejected_before_any_evidence_mutation(harness, change):
    harness.run(partition_size=2)
    metadata_path = harness.output / "run_metadata.json"
    metadata = state(harness.output)
    if change == "design_bytes":
        path = harness.output / "resolved_design.json"
        path.write_bytes(path.read_bytes() + b" ")
    elif change == "identity_digest":
        metadata["provenance"]["identity_sha256"] = "0" * 64
    elif change == "dependency_omission":
        metadata["provenance"]["source_files"].pop(next(iter(metadata["provenance"]["source_files"])))
        resign_identity(metadata["provenance"])
    elif change == "environment":
        metadata["provenance"]["environment"]["python"] = "0.0.0"
        resign_identity(metadata["provenance"])
    elif change == "procedure":
        metadata["procedure"]["max_condition_number"] = 1e12
    elif change == "root_cell_count":
        first = next(iter(metadata["completed_cells"]))
        metadata["completed_cells"][first]["paths"] = 999
    elif change == "unknown_root_file":
        (harness.output / "unknown.txt").write_text("unidentified evidence\n")
    else:
        first = next((harness.output / "cells").iterdir())
        (first / "unknown.txt").write_text("unidentified cell evidence\n")
    if change in {"identity_digest", "dependency_omission", "environment", "procedure", "root_cell_count"}:
        write_json(metadata_path, metadata)
    before = snapshot(harness.output)
    calls = list(harness.calls)
    with pytest.raises((ValueError, FileNotFoundError)):
        harness.run(resume=True)
    assert snapshot(harness.output) == before and harness.calls == calls


def rewrite_partition(path, mutator):
    rows = read_partition(path)
    mutator(rows)
    payload = gzip.compress(("\n".join(canonical(row) for row in rows) + "\n").encode(), compresslevel=6, mtime=0)
    Path(path).write_bytes(payload)
    return digest_bytes(payload)


def first_partition(output):
    cell_path = output / "cells" / followup_design().cells[0].identity
    manifest = json.loads((cell_path / "cell_manifest.json").read_text())
    part = manifest["partitions"][0]
    return cell_path, manifest, part, cell_path / part["filename"]


@pytest.mark.parametrize("corruption", ["nonempty", "missing", "nonstring"])
def test_pipeline_error_checks_scan_committed_parts_even_when_root_lags(harness, corruption):
    harness.run(partition_size=2)
    cell_path, manifest, part, path = first_partition(harness.output)

    def corrupt(rows):
        if corruption == "nonempty":
            rows[0].update(pipeline_error="RuntimeError: retained failure", interval_valid=False,
                           interval_status="pipeline_exception")
        elif corruption == "missing":
            rows[0].pop("pipeline_error")
        else:
            rows[0]["pipeline_error"] = None

    part["sha256"] = rewrite_partition(path, corrupt)
    write_json(cell_path / "cell_manifest.json", manifest)
    metadata = state(harness.output)
    metadata.update(status="incomplete", completed_cells={})
    write_json(harness.output / "run_metadata.json", metadata)
    before = snapshot(harness.output)
    calls = list(harness.calls)
    with pytest.raises(ValueError, match="(?i)(pipeline|error)"):
        harness.run(resume=True)
    assert snapshot(harness.output) == before and harness.calls == calls


def make_orphan(harness):
    design = followup_design().validation(5)
    harness.run(design=design, partition_size=3)
    cell_path = harness.output / "cells" / design.cells[0].identity
    manifest = json.loads((cell_path / "cell_manifest.json").read_text())
    orphan = manifest["partitions"].pop()
    assert (orphan["start"], orphan["stop"]) == (3, 5)
    manifest.update(completed_replications=3, status="running_incomplete")
    write_json(cell_path / "cell_manifest.json", manifest)
    metadata = state(harness.output)
    metadata["completed_cells"].pop(design.cells[0].identity)
    metadata["status"] = "incomplete"
    write_json(harness.output / "run_metadata.json", metadata)
    return design, cell_path / orphan["filename"]


def test_orphan_replay_preserves_original_range_when_partition_size_changes(harness):
    design, orphan = make_orphan(harness)
    before = orphan.read_bytes()
    previous = len(harness.calls)
    result = harness.run(design=design, resume=True, partition_size=1)
    assert result["status"] == "raw_complete_pending_audit"
    assert orphan.read_bytes() == before
    assert [(cell, rep) for cell, rep, *_ in harness.calls[previous:]] == [
        (design.cells[0].identity, 3), (design.cells[0].identity, 4)]
    manifest = json.loads((orphan.parent / "cell_manifest.json").read_text())
    assert manifest["recovered_orphans"][0]["method"] == "exact_byte_equal_deterministic_replay"
    assert (manifest["partitions"][-1]["start"], manifest["partitions"][-1]["stop"]) == (3, 5)


def test_byte_different_orphan_is_never_overwritten(harness):
    design, orphan = make_orphan(harness)
    with gzip.open(orphan, "rb") as handle:
        raw = handle.read()
    orphan.write_bytes(gzip.compress(raw, compresslevel=6, mtime=1))
    before = orphan.read_bytes()
    with pytest.raises(FileExistsError, match="immutable"):
        harness.run(design=design, resume=True, partition_size=1)
    assert orphan.read_bytes() == before
    assert state(harness.output)["status"] == "incomplete"


def test_pipeline_error_in_uncommitted_orphan_blocks_replay(harness):
    design, orphan = make_orphan(harness)

    def corrupt(rows):
        rows[0].update(pipeline_error="ArithmeticError: retained orphan failure", interval_valid=False,
                       interval_status="pipeline_exception")

    rewrite_partition(orphan, corrupt)
    before = snapshot(harness.output)
    calls = list(harness.calls)
    with pytest.raises(ValueError, match="(?i)(pipeline|error)"):
        harness.run(design=design, resume=True)
    assert snapshot(harness.output) == before and harness.calls == calls


def test_statistical_failure_candidates_are_retained_and_not_treated_as_pipeline_errors(harness):
    harness.run(partition_size=2)
    cell_path, manifest, part, path = first_partition(harness.output)

    def statistical_failure(rows):
        rows[0].update(point_valid=False, point_status="boundary_upper_theta", optimizer_success=False,
                       theta=1e100, interval_valid=False, interval_status="invalid_point_fit")

    part["sha256"] = rewrite_partition(path, statistical_failure)
    write_json(cell_path / "cell_manifest.json", manifest)
    before = snapshot(harness.output)
    calls = list(harness.calls)
    result = harness.run(resume=True)
    assert result["status"] == "raw_complete_pending_audit"
    assert snapshot(harness.output) == before and harness.calls == calls
    assert read_partition(path)[0]["theta"] == 1e100


def test_root_lock_excludes_concurrent_orchestrator(harness, monkeypatch):
    interrupt_after_one_cell(harness, monkeypatch)
    before = snapshot(harness.output)
    with run_lock(harness.output):
        with pytest.raises(RuntimeError, match="another process"):
            harness.run(resume=True)
    assert snapshot(harness.output) == before and len(harness.calls) == 3


def test_new_output_cannot_modify_a_preserved_archive(harness):
    archive = harness.root / "saved_archive"
    archive.mkdir()
    (archive / "run_metadata.json").write_text('{"status":"completed"}\n')
    before = snapshot(archive)
    with pytest.raises(ValueError, match="(?i)(nested|archive|saved)"):
        harness.run(output=archive / "new_run")
    assert snapshot(archive) == before and not harness.calls


def test_observed_sync_containment_rejects_mutable_output_before_creation(harness, monkeypatch):
    script = harness.root / "synthetic_sync_script.sh"
    script.write_text("synthetic presence marker; never executed\n")
    watched = harness.root / "watched_workspace"
    monkeypatch.setattr(execution, "SYNC_SCRIPT", script)
    monkeypatch.setattr(execution, "SYNC_ROOT", watched)
    with pytest.raises(ValueError, match="synchronized"):
        harness.run(output=watched / "new_run")
    assert not watched.exists() and not harness.calls


def test_output_symlink_is_rejected_before_creation(harness):
    target = harness.root / "target"
    target.mkdir()
    link = harness.root / "linked"
    link.symlink_to(target, target_is_directory=True)
    with pytest.raises(ValueError, match="symbolic"):
        harness.run(output=link / "new_run")
    assert not list(target.iterdir()) and not harness.calls


def test_fresh_pipeline_failure_is_retained_counts_all_planned_rows_and_blocks_retry(harness, monkeypatch):
    first = followup_design().cells[0].identity

    def failure(cell, rep, source, services):
        rows = harness.source_records(cell, rep, source, services)
        if cell.identity == first and rep == 0:
            rows[0].update(pipeline_error="RuntimeError: literal inference failure", interval_valid=False,
                           interval_status="pipeline_exception")
        return rows

    monkeypatch.setattr(historical, "source_records", failure)
    result = harness.run(partition_size=1)
    assert result["status"] == "implementation_review_required"
    assert result["pipeline_errors"] == 1
    assert len(harness.calls) == 6
    assert sum(value["rows"] for value in result["completed_cells"].values()) == 18
    before = snapshot(harness.output)
    with pytest.raises(ValueError, match="(?i)(pipeline|error)"):
        harness.run(resume=True)
    assert snapshot(harness.output) == before and len(harness.calls) == 6


def test_spawn_worker_exception_preserves_completed_cell_and_allows_only_unfinished_work(harness, monkeypatch):
    submitted = []
    contexts = []

    class SyntheticPool:
        def __init__(self, *, max_workers, mp_context):
            assert max_workers == 2
            contexts.append(mp_context.get_start_method())

        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

        def submit(self, function, *args):
            future = Future()
            submitted.append(future)
            if len(submitted) == 1:
                future.set_result(function(*args))
            else:
                future.set_exception(RuntimeError("synthetic worker unavailable before execution"))
            return future

    monkeypatch.setattr(execution, "ProcessPoolExecutor", SyntheticPool)
    monkeypatch.setattr(execution, "as_completed", lambda futures: list(futures))
    with pytest.raises(RuntimeError, match="synthetic worker"):
        harness.run(workers=2, partition_size=1)
    metadata = state(harness.output)
    assert contexts == ["spawn"]
    assert metadata["status"] == "incomplete"
    assert metadata["termination_type"] == "RuntimeError"
    assert "synthetic worker" in metadata["termination_reason"]
    assert len(metadata["completed_cells"]) == 1 and len(harness.calls) == 3
    result = harness.run(resume=True, workers=1, partition_size=2)
    assert result["status"] == "raw_complete_pending_audit"
    assert len(harness.calls) == 6 and len(set(harness.calls)) == 6


def test_root_drift_is_preserved_and_blocks_automatic_resume(harness, monkeypatch):
    injected = b'{"status":"unattributed_external_state","retain":"observed bytes"}\n'

    def drift(cell, rep, source, services):
        rows = harness.source_records(cell, rep, source, services)
        (harness.output / "run_metadata.json").write_bytes(injected)
        return rows

    monkeypatch.setattr(historical, "source_records", drift)
    with pytest.raises(RuntimeError, match="(?i)(manifest|drift)"):
        harness.run(partition_size=1)
    assert (harness.output / "run_metadata.json").read_bytes() == injected
    incidents = [p for p in harness.root.iterdir() if p.is_dir() and "manifest-incident-" in p.name]
    assert incidents
    assert any((p / "observed_run_metadata.json").read_bytes() == injected for p in incidents)
    expected = (incidents[0] / "expected_run_metadata.json").read_bytes()
    (harness.output / "run_metadata.json").write_bytes(expected)
    monkeypatch.setattr(historical, "source_records", harness.source_records)
    before = snapshot(harness.root)
    calls = list(harness.calls)
    with pytest.raises(ValueError, match="incident"):
        harness.run(resume=True)
    assert snapshot(harness.root) == before and harness.calls == calls
