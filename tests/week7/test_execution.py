"""Guarded orchestration contracts with literal mock records, no simulation."""
from dataclasses import replace
import json

import pytest

from ou_irregular.experiments import execution as engine
from ou_irregular.experiments.audit import semantic_metadata_digest
from ou_irregular.experiments.design import ESTIMATORS
from ou_irregular.experiments.provenance import canonical
from ou_irregular.week7 import execution
from ou_irregular.week7.selection import floor_design, followup_design


def literal_records(cell, rep, source, services):
    return [{**cell.to_dict(), "rep": rep, "estimator": method,
             "point_valid": False, "point_status": "boundary", "optimizer_success": True,
             "interval_valid": False, "interval_status": "invalid_point_fit", "nll": 0.,
             "theta": 1., "mu": 0., "sigma": .5, "pipeline_error": ""}
            for method in ESTIMATORS]


def snapshot(path):
    return {str(p.relative_to(path)): p.read_bytes() for p in path.rglob("*") if p.is_file()}


@pytest.mark.parametrize("kind,factory", [("followup", followup_design), ("floor", floor_design)])
def test_exact_design_authority_and_planned_policy_counts(kind, factory):
    full = factory()
    assert execution.require_design(kind, full) is full
    validation = full.validation(2)
    assert execution.require_design(kind, validation) is validation
    policy = execution.source_policy(kind, validation)
    assert policy["planned_source_realizations"] == len(full.cells) * 2
    assert policy["planned_fit_records"] == len(full.cells) * 6
    assert policy["planned_independent_coupling_groups"] == (8 if kind == "floor" else 4)
    assert set(policy["actual_to_seed_anchor"]) == {c.identity for c in full.cells}
    for changed in (replace(validation, cells=validation.cells[::-1]),
                    replace(validation, cells=validation.cells[:-1]),
                    replace(validation, seed=full.seed), replace(validation, name="unregistered_name"),
                    replace(full, replications=3), replace(full, scope="headline")):
        with pytest.raises(ValueError, match="design"):
            execution.require_design(kind, changed)


def test_new_preflight_is_pure_and_frozen_full_failure_creates_nothing(tmp_path, monkeypatch):
    output = tmp_path / "nested" / "new"
    design = followup_design().validation(2)
    before = snapshot(tmp_path)
    result = execution.preflight("followup", design, output)
    assert result.metadata is None and not result.cached_completion
    assert result.require_frozen is False and result.cells == {}
    assert snapshot(tmp_path) == before
    original = execution.capture
    def reject_full(*args, **kwargs):
        assert kwargs["require_frozen"] is True
        raise ValueError("synthetic unfrozen dependency refusal")
    monkeypatch.setattr(execution, "capture", reject_full)
    with pytest.raises(ValueError, match="unfrozen"):
        execution.run_raw("followup", followup_design(), output, progress=False)
    assert not output.exists() and not output.parent.exists()
    monkeypatch.setattr(execution, "capture", original)


@pytest.mark.parametrize("kind,factory", [("followup", followup_design), ("floor", floor_design)])
def test_literal_raw_completion_has_bound_policy_and_cached_noop(tmp_path, monkeypatch, kind, factory):
    monkeypatch.setattr(engine, "source_records", literal_records)
    design = factory().validation(2)
    output = tmp_path / kind
    result = execution.run_raw(kind, design, output, partition_size=1, progress=False)
    assert result["status"] == "raw_complete_pending_audit"
    assert result["pipeline_errors"] == 0
    assert len(result["completed_cells"]) == len(design.cells)
    assert result["scientific_overrides"] == (["paired_floor_source"] if kind == "floor" else [])
    if kind == "followup":
        assert result["project"] == "ou-irregular-experiments"
        assert result["procedure"] == engine.PROCEDURE
    else:
        assert result["project"] == "ou-week7-paired-floor"
        assert {k: v for k, v in result["procedure"].items() if k != "source_address"} == {
            k: v for k, v in engine.PROCEDURE.items() if k != "source_address"}
    baseline = semantic_metadata_digest(result)
    for field in ("source_policy", "execution_kind", "execution_contract", "require_frozen"):
        assert semantic_metadata_digest({**result, field: None}) != baseline
    before = snapshot(output)
    monkeypatch.setattr(execution, "_execute_cell", lambda *args: pytest.fail("cached completion called engine"))
    cached = execution.run_raw(kind, design, output, resume=True, partition_size=2, progress=False)
    assert cached == result and snapshot(output) == before


def test_known_postprocessing_manifest_and_failure_fields_survive_cached_return(tmp_path, monkeypatch):
    monkeypatch.setattr(engine, "source_records", literal_records)
    design = floor_design().validation(1)
    result = execution.run_raw("floor", design, tmp_path, progress=False)
    result.update(status="implementation_review_required", postprocessing_error_type="SyntheticError",
                  postprocessing_error="saved reporting failure")
    (tmp_path / "run_metadata.json").write_text(canonical(result) + "\n")
    (tmp_path / "epsilon_table_manifest.json").write_text('{"future_reporting_contract":true}\n')
    before = snapshot(tmp_path)
    assert execution.run_raw("floor", design, tmp_path, resume=True, progress=False) == result
    assert snapshot(tmp_path) == before


@pytest.mark.parametrize("status", ["audited_pending_visual_review", "completed", "unknown"])
def test_closed_or_unknown_lifecycle_never_reenters_raw_execution(tmp_path, monkeypatch, status):
    monkeypatch.setattr(engine, "source_records", literal_records)
    design = followup_design().validation(1)
    metadata = execution.run_raw("followup", design, tmp_path, progress=False)
    metadata["status"] = status
    (tmp_path / "run_metadata.json").write_text(canonical(metadata) + "\n")
    before = snapshot(tmp_path)
    with pytest.raises(ValueError):
        execution.run_raw("followup", design, tmp_path, resume=True, progress=False)
    assert snapshot(tmp_path) == before


def test_pipeline_error_saved_in_root_lagging_cell_blocks_before_source(tmp_path, monkeypatch):
    def failing(cell, rep, source, services):
        rows = literal_records(cell, rep, source, services)
        rows[0]["pipeline_error"] = "RuntimeError: literal preserved failure"
        return rows
    monkeypatch.setattr(engine, "source_records", failing)
    design = followup_design().validation(1)
    metadata = execution.run_raw("followup", design, tmp_path, progress=False)
    assert metadata["status"] == "implementation_review_required"
    metadata.update(status="incomplete", completed_cells={}, pipeline_errors=0)
    (tmp_path / "run_metadata.json").write_text(canonical(metadata) + "\n")
    before = snapshot(tmp_path)
    monkeypatch.setattr(engine, "source_records", lambda *args: pytest.fail("saved failure was retried"))
    with pytest.raises(execution.SavedPipelineError):
        execution.run_raw("followup", design, tmp_path, resume=True, progress=False)
    assert snapshot(tmp_path) == before


def test_sync_root_and_existing_archive_ancestor_are_blocked_before_creation(tmp_path, monkeypatch):
    marker = tmp_path / "sync_marker"
    marker.write_text("literal marker")
    watched = tmp_path / "watched"
    monkeypatch.setattr(execution, "SYNC_SCRIPT", marker)
    monkeypatch.setattr(execution, "SYNC_ROOT", watched)
    design = followup_design().validation(1)
    with pytest.raises(ValueError, match="synchronized"):
        execution.preflight("followup", design, watched / "new")
    assert not watched.exists()
    archived = tmp_path / "archive"
    archived.mkdir()
    (archived / "run_metadata.json").write_text("{}\n")
    with pytest.raises(ValueError, match="saved experiment"):
        execution.preflight("followup", design, archived / "nested")
    assert not (archived / "nested").exists()
