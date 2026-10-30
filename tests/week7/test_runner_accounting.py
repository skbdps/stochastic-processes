"""Independent accounting checks at the Week 7 runner's composition boundary.

The source rows, audit result, and rendering result are explicitly LITERAL TEST
FIXTURES. No observations, random draws, model fits, derivative calculations or
real numerical audit are performed here. Durable storage, standard/epsilon table
calculation, runner preflight, saved-audit readback and lifecycle guards are real.
These tests establish accounting/composition behavior, not scientific validity.
"""
from copy import deepcopy
import hashlib
import json
import math
from pathlib import Path
import shutil

import pandas as pd
import pytest

from ou_irregular.experiments import audit as ordinary_audit
from ou_irregular.experiments import execution as historical_execution
from ou_irregular.experiments.provenance import canonical
from ou_irregular.inference.adapters import Week3Estimator
from ou_irregular.week4.simulation import OuSampleSource
from ou_irregular.week7 import execution, floor_audit, plots, reporting, runner
from ou_irregular.week7.floor_source import PairedFloorSource

from .test_runner_reference import forbidden, literal_records, snapshot


def _write(path, value):
    Path(path).write_text(canonical(value) + "\n")


def _literal_audit_report(output, design, kind, rows):
    """Return an explicitly fabricated audit-stage result for known R=2 rows.

    The counts below are calculated from the literal six-record pattern, not
    read from the producer's summaries: five valid points, three valid
    intervals and three failures per cell. Structural hashes bind those exact
    fixtures; they do not turn them into independently audited OU observations.
    """
    metadata = json.loads((output / "run_metadata.json").read_text())
    provenance = metadata["provenance"]
    cells = len(design.cells)
    assert design.replications == 2
    endpoints = []
    for cell in design.cells:
        for rep in (0, 1):
            row = next(row for row in rows if row["cell_id"] == cell.identity
                       and row["rep"] == rep and row["estimator"] == "exact")
            endpoints.append({"cell_id": cell.identity, "replication": rep,
                              "times_sha256": row["times_sha256"],
                              "observations_sha256": row["observations_sha256"],
                              "estimators_checked": ["exact", "pfml", "euler"]})
    report = {
        "schema_version": 1, "status": "passed", "full_completion_audit": True,
        "replay_performed": True, "design_sha256": design.digest,
        "run_identity_sha256": provenance["identity_sha256"],
        "code_commit": provenance["code_commit"],
        "source_dependency_files_verified": len(provenance["source_files"]),
        "audited_artifact_sha256": ordinary_audit.artifact_fingerprint(output),
        "audited_run_metadata_sha256": ordinary_audit.semantic_metadata_digest(metadata),
        "cells_verified": cells, "partitions_verified": 2 * cells,
        "planned_paths": 2 * cells, "planned_fit_records": 6 * cells,
        "point_valid_records": 5 * cells, "interval_valid_records": 3 * cells,
        "retained_failure_records": 3 * cells, "pipeline_errors": 0,
        "replayed_source_paths": 2 * cells, "replayed_fits": 6 * cells,
        "replay_records": endpoints,
        "replay_selection": ("replications 0 and R-1 in every physical floor variant, all estimators"
                             if kind == "floor" else "replications 0 and R-1 in every cell, all estimators"),
        "tables_verified": {"parameter_summary.csv": 9 * cells,
                            "regular_comparisons.csv": 9 * cells,
                            "paired_comparisons.csv": 9 * cells,
                            "spacing_summary.csv": cells,
                            "outcome_counts.csv": 6 * cells,
                            "numerical_diagnostics.csv": 3 * cells},
        "tolerances": {
            "arithmetic_rtol": 2e-10, "arithmetic_atol": 2e-12,
            "replay_rtol": 2e-10, "replay_atol": 2e-12,
            "covariance_rtol": 2e-9,
            "covariance_atol": "2e-12 * max(1, expected matrix Frobenius norm)",
        },
        "scope": {"fixture": "literal runner accounting test; no numerical verification performed"},
    }
    if kind == "floor":
        binding = reporting.raw_input_binding(output, design, metadata)
        report.update(
            audit_contract="week7_independent_paired_floor_audit_v1",
            scientific_overrides=["paired_floor_source"], source_policy=metadata["source_policy"],
            raw_input_sha256=hashlib.sha256(canonical(binding).encode()).hexdigest(),
            planned_source_realizations=24, planned_independent_coupling_groups=8,
            source_replay_performed=True, source_paths_verified=24, source_array_hashes_verified=48,
        )
        report["tables_verified"]["epsilon_comparisons.csv"] = 72
        report["tolerances"]["source_arrays_and_diagnostics"] = "exact types, values and SHA-256 array bytes"
    else:
        report["tolerances"].update(
            span_roundoff="8*gamma_(K+8)*max(span,K*nominal_gap)+8 ulps; gamma_m=m*eps/(1-m*eps)",
            minimum_gap_roundoff="4 ulps at the realized-span scale",
        )
    return report


@pytest.fixture(scope="module", params=("followup", "floor"))
def literal_accounting_archive(request, tmp_path_factory):
    kind = request.param
    design = runner.design_for(kind, validation_replications=2)
    rows = literal_records(design, kind)
    by_source = {}
    for row in rows:
        by_source.setdefault((row["cell_id"], row["rep"]), []).append(row)
    output = tmp_path_factory.mktemp("runner_literal_accounting_" + kind) / "records"
    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(OuSampleSource, "sample", forbidden)
        patch.setattr(PairedFloorSource, "sample", forbidden)
        patch.setattr(Week3Estimator, "fit", forbidden)
        patch.setattr(historical_execution, "source_records",
                      lambda cell, rep, source, services: deepcopy(by_source[cell.identity, rep]))
        execution.run_raw(kind, design, output, extra_files=runner.WEEK_FILES,
                          workers=1, partition_size=1, progress=False)
        patch.setattr(historical_execution, "source_records", forbidden)
        patch.setattr(execution, "run_raw", forbidden)
        reporting.write_tables(kind, design, output, extra_files=runner.WEEK_FILES)
    return {"kind": kind, "design": design, "rows": rows, "output": output}


def _case_copy(case, destination):
    output = destination / "literal_study"
    shutil.copytree(case["output"], output)
    return output


def test_runner_composes_full_audit_without_changing_literal_denominators(
        literal_accounting_archive, tmp_path, monkeypatch):
    case = literal_accounting_archive
    kind, design = case["kind"], case["design"]
    output = _case_copy(case, tmp_path)
    report = _literal_audit_report(output, design, kind, case["rows"])
    calls = []

    def literal_audit(directory, supplied, **kwargs):
        assert Path(directory) == output and supplied == design
        assert kwargs["replay"] is True
        if kind == "floor":
            assert set(kwargs["extra_files"]) == set(runner.WEEK_FILES)
        else:
            assert set(kwargs) == {"replay"}
        calls.append("full_audit_stage_double")
        return deepcopy(report)

    monkeypatch.setattr(OuSampleSource, "sample", forbidden)
    monkeypatch.setattr(PairedFloorSource, "sample", forbidden)
    monkeypatch.setattr(Week3Estimator, "fit", forbidden)
    monkeypatch.setattr(execution, "run_raw", forbidden)
    monkeypatch.setattr(ordinary_audit, "audit_run", literal_audit if kind == "followup" else forbidden)
    monkeypatch.setattr(floor_audit, "audit_run", literal_audit if kind == "floor" else forbidden)
    monkeypatch.setattr(plots, "write_figures", lambda *args, **kwargs: calls.append("literal_render_stage_double"))
    monkeypatch.setattr(plots, "verify_figures", lambda *args, **kwargs: {"fixture": "no real image produced"})
    raw_before = snapshot(output / "cells")
    result = runner.postprocess(kind, output)
    assert result == report
    assert calls == ["full_audit_stage_double", "literal_render_stage_double"]
    assert snapshot(output / "cells") == raw_before
    metadata = json.loads((output / "run_metadata.json").read_text())
    assert metadata["status"] == "audited_pending_visual_review"
    assert not (output / "visual_review.json").exists()
    assert "completed_utc" not in metadata
    assert len(metadata["completed_cells"]) == (12 if kind == "floor" else 2)
    assert result["planned_fit_records"] == (72 if kind == "floor" else 12)
    assert result["point_valid_records"] == (60 if kind == "floor" else 10)
    assert result["interval_valid_records"] == (36 if kind == "floor" else 6)
    assert result["retained_failure_records"] == (36 if kind == "floor" else 6)
    # Every fixed endpoint already appears once in the fake full-audit report.
    # A cached audit must pass the REAL saved-evidence accounting readback.
    _, _, saved = runner._saved_audit(kind, output, design)
    assert saved == report


def test_runner_readback_preserves_hand_calculated_point_and_coverage_subsets(
        literal_accounting_archive, tmp_path):
    case = literal_accounting_archive
    kind, design = case["kind"], case["design"]
    output = _case_copy(case, tmp_path)
    report = _literal_audit_report(output, design, kind, case["rows"])
    _write(output / "audit.json", report)
    before = snapshot(output)
    runner._saved_audit(kind, output, design)
    assert snapshot(output) == before

    summary = pd.read_csv(output / "parameter_summary.csv", float_precision="round_trip")
    assert summary.n_planned.eq(2).all()
    assert summary.n_interval_valid.eq(1).all()
    assert summary.coverage_conditional.eq(1).all()
    assert summary.coverage_operational.eq(.5).all()
    # All rows have finite endpoints containing truth, including unavailable
    # intervals. Their exclusion from the operational numerator is essential.
    for method in ("exact", "pfml"):
        group = summary[summary.estimator == method]
        assert group.n_point_valid.eq(2).all()
        assert group.bias.to_numpy() == pytest.approx([.15] * len(group))
        assert group.rmse.to_numpy() == pytest.approx([math.sqrt(.025)] * len(group))
        assert group.mcse_bias.to_numpy() == pytest.approx([.05] * len(group))
    euler = summary[summary.estimator == "euler"]
    assert euler.n_point_valid.eq(1).all()
    assert euler.bias.to_numpy() == pytest.approx([.2] * len(euler))
    assert euler.mcse_bias.isna().all()
    mu = summary[summary.parameter == "mu"]
    assert mu.relative_bias.isna().all() and mu.relative_rmse.isna().all()

    paired = pd.read_csv(output / "paired_comparisons.csv", float_precision="round_trip")
    against_euler = paired[(paired.estimator_a == "exact") & (paired.estimator_b == "euler")]
    assert against_euler.n_common_point_valid.eq(1).all()
    assert against_euler.bias_a_common.to_numpy() == pytest.approx([.2] * len(against_euler))
    assert against_euler.bias_b_common.to_numpy() == pytest.approx([.2] * len(against_euler))
    assert against_euler.bias_difference.eq(0).all()
    # Covering indicators are (1,0) versus (0,1); differences (1,-1)
    # have sample SD sqrt(2), giving paired MCSE exactly 1 over two trials.
    assert against_euler.operational_coverage_difference.eq(0).all()
    assert against_euler.operational_coverage_difference_mcse.eq(1).all()
    assert (against_euler.operational_coverage_difference_mc_lower < -1).all()
    assert (against_euler.operational_coverage_difference_mc_upper > 1).all()
    if kind == "floor":
        epsilon = pd.read_csv(output / "epsilon_comparisons.csv", float_precision="round_trip")
        assert len(epsilon) == 72
        assert report["planned_source_realizations"] == 24
        assert report["planned_independent_coupling_groups"] == 8


@pytest.mark.parametrize("field,value", (
    ("planned_fit_records", 1),
    ("point_valid_records", 0),
    ("interval_valid_records", 0),
    ("retained_failure_records", 0),
    ("partitions_verified", 0),
    ("full_completion_audit", False),
    ("replay_performed", False),
    ("replayed_fits", 0),
))
def test_runner_rejects_inconsistent_saved_audit_counts_before_any_new_stage(
        literal_accounting_archive, tmp_path, monkeypatch, field, value):
    case = literal_accounting_archive
    output = _case_copy(case, tmp_path)
    report = _literal_audit_report(output, case["design"], case["kind"], case["rows"])
    report[field] = value
    _write(output / "audit.json", report)
    before = snapshot(output)
    for module, name in ((execution, "run_raw"), (reporting, "write_tables"),
                         (ordinary_audit, "audit_run"), (floor_audit, "audit_run"),
                         (plots, "write_figures")):
        monkeypatch.setattr(module, name, forbidden)
    with pytest.raises(ValueError, match="full saved audit|table/count"):
        runner.postprocess(case["kind"], output)
    assert snapshot(output) == before


@pytest.mark.parametrize("field,value", (
    ("planned_independent_coupling_groups", 24),
    ("source_paths_verified", 8),
    ("source_array_hashes_verified", 24),
))
def test_floor_physical_paths_are_not_miscounted_as_independent_groups(
        literal_accounting_archive, tmp_path, monkeypatch, field, value):
    case = literal_accounting_archive
    if case["kind"] != "floor":
        pytest.skip("floor-only coupling accounting; the follow-up has no paired-floor fields")
    output = _case_copy(case, tmp_path)
    report = _literal_audit_report(output, case["design"], "floor", case["rows"])
    report[field] = value
    _write(output / "audit.json", report)
    before = snapshot(output)
    monkeypatch.setattr(plots, "write_figures", forbidden)
    with pytest.raises(ValueError, match="paired-floor source override and full source scope"):
        runner.postprocess("floor", output)
    assert snapshot(output) == before
