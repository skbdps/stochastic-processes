"""Reporting contracts on literal raw records; no observations are generated."""
from copy import deepcopy
import hashlib
import json

import pandas as pd
import pytest

from ou_irregular.experiments import execution as old_engine
from ou_irregular.experiments.design import ESTIMATORS
from ou_irregular.experiments.provenance import canonical, file_digest
from ou_irregular.week7 import execution, reporting
from ou_irregular.week7.floor_comparisons import COLUMNS
from ou_irregular.week7.floor_source import COUPLING_POLICY, PairedFloorSource, anchor_cell
from ou_irregular.week7.selection import floor_design, followup_design


def literal_records(cell, rep, source, services):
    anchor = anchor_cell(cell) if isinstance(source, PairedFloorSource) else cell
    words = [int(anchor.identity[i:i + 8], 16) for i in range(0, 64, 8)]
    key = [4, *words, rep]
    sample_label = anchor.identity + str(rep) + (str(cell.floor_fraction) if cell.cv else "")
    rows = []
    for method in ESTIMATORS:
        row = {**cell.to_dict(), "rep": rep, "estimator": method,
               "source_seed": source.seed, "source_spawn_key": json.dumps(key),
               "gap_spawn_key": json.dumps([*key, 0]), "path_spawn_key": json.dumps([*key, 1]),
               "times_sha256": hashlib.sha256(("time" + sample_label).encode()).hexdigest(),
               "observations_sha256": hashlib.sha256(("state" + sample_label).encode()).hexdigest(),
               "source_valid": True, "source_error": "", "pipeline_error": "",
               "point_valid": rep == 0, "point_status": "valid_interior" if rep == 0 else "boundary",
               "optimizer_success": True, "interval_valid": rep == 0,
               "interval_status": "valid" if rep == 0 else "invalid_point_fit", "nll": 0.,
               "realized_mean_gap": cell.q / cell.theta, "realized_cv": cell.cv,
               "realized_span": (cell.n - 1) * cell.q / cell.theta, "min_gap": .1,
               "clip_fraction": 0., "hessian_condition_number": 10.,
               "hessian_relative_difference": .001, "hessian_whitened_difference": .001,
               "standardized_score": .001}
        if isinstance(source, PairedFloorSource):
            row.update(actual_cell_id=cell.identity, seed_anchor_cell_id=anchor.identity,
                       seed_anchor_floor_fraction=1e-6, source_coupling_policy=COUPLING_POLICY)
        for parameter in ("theta", "mu", "sigma"):
            truth = cell.to_dict()["true_" + parameter]
            row.update({parameter: truth + .1 if rep == 0 else 999., parameter + "_se": .1,
                        parameter + "_lower": truth - .2, parameter + "_upper": truth + .2})
        rows.append(row)
    return rows


@pytest.fixture
def raw_run(tmp_path, monkeypatch):
    monkeypatch.setattr(old_engine, "source_records", literal_records)
    def build(kind="floor"):
        design = (floor_design() if kind == "floor" else followup_design()).validation(2)
        output = tmp_path / kind
        execution.run_raw(kind, design, output, extra_files=reporting.FILES, progress=False)
        return design, output
    return build


def snapshot(output):
    return {str(p.relative_to(output)): p.read_bytes() for p in output.rglob("*") if p.is_file()}


@pytest.mark.parametrize("kind", ["followup", "floor"])
def test_reports_exact_tables_and_preserves_root_raw_failures(raw_run, kind):
    design, output = raw_run(kind)
    before = snapshot(output)
    result = reporting.write_tables(kind, design, output)
    assert set(result) == {"standard", "epsilon"}
    assert set(result["standard"]) == set(reporting.STANDARD_FILES)
    assert result["standard"]["parameter_summary.csv"]["rows"] == 9 * len(design.cells)
    assert result["standard"]["failed_fits.csv.gz"]["rows"] == 3 * len(design.cells)
    for name, raw in before.items():
        assert (output / name).read_bytes() == raw
    if kind == "floor":
        epsilon = pd.read_csv(output / reporting.EPSILON_FILE)
        assert len(epsilon) == 72 and tuple(epsilon.columns) == COLUMNS
        assert result["epsilon"]["raw_inputs"]["source_policy"]["planned_independent_coupling_groups"] == 8
        assert result["epsilon"]["table"]["sha256"] == file_digest(output / reporting.EPSILON_FILE)
    else:
        assert result["epsilon"] is None and not (output / reporting.EPSILON_FILE).exists()
    complete = snapshot(output)
    assert reporting.write_tables(kind, design, output) == result
    assert snapshot(output) == complete


def test_preflight_is_read_only_and_input_binding_excludes_lifecycle(raw_run):
    design, output = raw_run()
    before = snapshot(output)
    prepared = reporting.preflight_reporting("floor", design, output)
    assert snapshot(output) == before
    altered = deepcopy(prepared.raw.metadata)
    altered.update(status="completed", completed_utc="different lifecycle time", audit_sha256="a" * 64)
    assert reporting.raw_input_binding(output, design, altered) == prepared.binding


@pytest.mark.parametrize("artifact", ["parameter_summary.csv", reporting.EPSILON_FILE])
def test_declared_missing_artifact_is_not_recreated(raw_run, artifact):
    design, output = raw_run()
    reporting.write_tables("floor", design, output)
    (output / artifact).unlink()
    before = snapshot(output)
    with pytest.raises((ValueError, FileNotFoundError)):
        reporting.write_tables("floor", design, output)
    assert snapshot(output) == before and not (output / artifact).exists()


@pytest.mark.parametrize("kind", ["followup", "floor"])
def test_closed_state_cannot_regenerate_tables(raw_run, kind):
    design, output = raw_run(kind)
    metadata = json.loads((output / "run_metadata.json").read_text())
    metadata["status"] = "completed"
    (output / "run_metadata.json").write_text(canonical(metadata) + "\n")
    before = snapshot(output)
    with pytest.raises(ValueError):
        reporting.write_tables(kind, design, output)
    assert snapshot(output) == before


def test_followup_rejects_epsilon_artifact_without_writing_standard_tables(raw_run):
    design, output = raw_run("followup")
    (output / reporting.EPSILON_FILE).write_text("unexpected paired-floor evidence\n")
    before = snapshot(output)
    with pytest.raises(ValueError, match="follow-up"):
        reporting.write_tables("followup", design, output)
    assert snapshot(output) == before


def test_epsilon_manifest_binds_all_raw_cells_and_exact_partition_maps(raw_run):
    design, output = raw_run()
    result = reporting.write_tables("floor", design, output)["epsilon"]
    binding = result["raw_inputs"]
    assert set(binding["cells"]) == {cell.identity for cell in design.cells}
    assert result["raw_input_sha256"] == hashlib.sha256(canonical(binding).encode()).hexdigest()
    for cid, entry in binding["cells"].items():
        assert entry["manifest_sha256"] == file_digest(output / "cells" / cid / "cell_manifest.json")
        assert sum(part["rows"] for part in entry["partitions"]) == 3 * design.replications
        for part in entry["partitions"]:
            assert part["sha256"] == file_digest(output / "cells" / cid / part["filename"])
