"""Independent numerical and completion-guard checks for the supplement."""
import json
import math

import numpy as np
import pandas as pd
import pytest

from ou_irregular.experiments.provenance import file_digest
from ou_irregular.week10.cli import seal, verify_seal
from ou_irregular.week10.design import validation_design
from ou_irregular.week10.reporting import contrasts, mechanism_summary, source_summary, summaries


def synthetic_records():
    design = validation_design()
    records = []
    for cell in design.cells:
        for rep in range(4):
            error = float(rep + 1) if cell.mechanism == "mmpp" else float(1 + rep % 2)
            point = not (cell.mechanism == "mmpp" and rep == 3)
            interval = point and not (cell.mechanism == "mmpp" and rep == 2)
            for estimator in ("exact", "pfml", "euler"):
                row = {**cell.to_dict(), "rep": rep, "estimator": estimator,
                       "actual_n": 1000 + rep, "actual_transitions": 999 + rep,
                       "source_valid": True, "point_valid": point, "interval_valid": interval,
                       "optimizer_success": point, "point_status": "valid_interior" if point else "invalid",
                       "interval_status": "valid" if interval else "unavailable",
                       "nll": 1., "pipeline_error": "", "realized_span": cell.horizon,
                       "realized_mean_gap": cell.horizon / (999 + rep), "terminal_gap": .2,
                       "clip_count": 0, "min_gap": .01,
                       "times_sha256": f"time-{cell.identity}-{rep}",
                       "observations_sha256": f"observations-{cell.identity}-{rep}"}
                for parameter in ("theta", "mu", "sigma"):
                    truth = row["true_" + parameter]
                    estimate = truth + error
                    row.update({parameter: estimate, parameter + "_se": 1. if interval else np.nan,
                                parameter + "_lower": estimate - 1.5 if interval else np.nan,
                                parameter + "_upper": estimate + 1.5 if interval else np.nan})
                records.append(row)
    return design, records


def test_raw_summaries_preserve_failure_denominators_and_nominal_count_label():
    design, records = synthetic_records()
    table = summaries(records, design)
    assert len(table) == 7 * 3 * 3
    row = table[(table.mechanism == "mmpp") & (table.ratio == 5)
                & (table.estimator == "pfml") & (table.parameter == "mu")].iloc[0]
    assert row.n_planned == 4
    assert row.n_point_valid == 3
    assert row.n_interval_valid == 2
    assert row.n_covered == 1
    assert row.coverage_operational == .25
    assert row.coverage_conditional == .5
    assert row.bias == 2.
    assert row.rmse == pytest.approx(math.sqrt(14 / 3))
    assert row.mcse_bias == pytest.approx(math.sqrt(1 / 3))
    assert row.mcse_rmse == pytest.approx(np.std([1, 4, 9], ddof=1) / (2 * math.sqrt(3) * math.sqrt(14 / 3)))
    assert row.nominal_n == 1000
    assert "n" not in table
    assert row.horizon == 499.5


def test_independent_contrast_mc_uncertainty_matches_hand_calculation():
    design, records = synthetic_records()
    table = contrasts(summaries(records, design))
    assert len(table) == 3 * 3 * 3
    row = table[(table.ratio == 5) & (table.estimator == "pfml") & (table.parameter == "mu")].iloc[0]
    assert row.mmpp_n_planned == row.iid_n_planned == 4
    assert row.mmpp_n_point_valid == 3
    assert row.iid_n_point_valid == 4
    expected_bias_se = math.sqrt(np.var([1., 2., 3.], ddof=1) / 3 + np.var([1., 2., 1., 2.], ddof=1) / 4)
    expected_coverage_se = math.sqrt(.25 * .75 / 4 + .5 * .5 / 4)
    assert row.bias_difference == .5
    assert row.bias_difference_mcse == pytest.approx(expected_bias_se)
    assert row.coverage_difference == -.25
    assert row.coverage_difference_mcse == pytest.approx(expected_coverage_se)
    mse_m, mse_i = 14 / 3, 2.5
    logse = .5 * math.sqrt(np.var([1., 4., 9.], ddof=1) / (3 * mse_m**2)
                          + np.var([1., 4., 1., 4.], ddof=1) / (4 * mse_i**2))
    ratio = math.sqrt(mse_m / mse_i)
    assert row.rmse_ratio == pytest.approx(ratio)
    assert row.rmse_ratio_log_mcse == pytest.approx(logse)
    assert row.rmse_ratio_mc_lower == pytest.approx(ratio * math.exp(-1.959963984540054 * logse))
    assert row.rmse_ratio_mc_upper == pytest.approx(ratio * math.exp(1.959963984540054 * logse))


def test_missing_or_duplicate_reference_cannot_supply_a_contrast():
    design, records = synthetic_records()
    table = summaries(records, design)
    missing = table[~((table.mechanism == "iid") & (table.ratio == 5)
                      & (table.estimator == "pfml") & (table.parameter == "mu"))]
    with pytest.raises(ValueError):
        contrasts(missing)
    duplicate = table[(table.mechanism == "iid") & (table.ratio == 5)
                      & (table.estimator == "pfml") & (table.parameter == "mu")]
    with pytest.raises(ValueError):
        contrasts(pd.concat([table, duplicate]))


def test_missing_planned_fit_is_rejected_before_summary():
    design, records = synthetic_records()
    with pytest.raises(ValueError):
        summaries(records[:-1], design)


def test_source_count_is_not_tripled_by_three_estimators():
    design, records = synthetic_records()
    table = source_summary(records, design)
    assert len(table) == 7
    assert table.n_paths.eq(4).all()
    assert table.actual_n_available.eq(4).all()
    assert table.actual_n_mean.eq(1001.5).all()
    assert table.realized_span_min.eq(499.5).all()
    assert table.realized_span_max.eq(499.5).all()


def test_mechanism_prediction_calibration_retains_correlated_difference_mcse():
    design, records = synthetic_records()
    summary = summaries(records, design)
    cell = design.cells[0]
    paths = pd.DataFrame({"cell_id": [cell.identity] * 4, "mechanism": [cell.mechanism] * 4,
                          "ratio": [cell.ratio] * 4, "fixed_rho_conditional_variance": [1., 2., 3., 4.],
                          "fixed_rho_mu": [1., -2., 3., -4.], "population_working_variance": [1., 1., 1., 1.],
                          "normal_working_coverage_approximation": [.8, .9, .7, .6],
                          "point_valid": [True, True, True, False], "pfml_mu": [1.1, -1.8, 2.7, np.nan]})
    row = mechanism_summary(paths, summary).iloc[0]
    differences = np.array([0., 2., 6., 12.])
    assert row.predicted_fixed_rho_sd == pytest.approx(math.sqrt(2.5))
    assert row.empirical_fixed_rho_rmse == pytest.approx(math.sqrt(7.5))
    assert row.fixed_rho_mse_minus_predicted_variance == pytest.approx(5.)
    assert row.fixed_rho_mse_minus_predicted_variance_mcse == pytest.approx(differences.std(ddof=1) / 2)
    assert row.fitted_minus_fixed_mu_rms == pytest.approx(math.sqrt((.01 + .04 + .09) / 3))


def skeletal_package(tmp_path, *, with_figure=False):
    root = tmp_path / "package"
    (root / "evidence").mkdir(parents=True)
    (root / "figures").mkdir()
    files = {}
    if with_figure:
        for name in ("diagnostic.png", "diagnostic.pdf"):
            (root / "figures" / name).write_bytes(b"synthetic byte fixture " + name.encode())
            files[name] = file_digest(root / "figures" / name)
    (root / "figures/figure_manifest.json").write_text(json.dumps({"files": files}))
    (root / "evidence/visual_review.json").write_text(json.dumps({
        "status": "passed", "inspected_files": files,
        "figure_manifest_sha256": file_digest(root / "figures/figure_manifest.json")}))
    (root / "evidence/final_validation.json").write_text(json.dumps({"status": "passed"}))
    return root


@pytest.mark.parametrize("with_figure", [False, True])
def test_seal_cannot_complete_without_full_science_and_bound_validation(tmp_path, with_figure):
    root = skeletal_package(tmp_path, with_figure=with_figure)
    with pytest.raises((ValueError, FileNotFoundError, KeyError)):
        seal(root)
    assert not (root / "completion_seal.json").exists()


def test_changed_reviewed_figure_cannot_be_sealed(tmp_path):
    root = skeletal_package(tmp_path, with_figure=True)
    (root / "figures/diagnostic.png").write_bytes(b"changed after review")
    with pytest.raises((ValueError, FileNotFoundError, KeyError)):
        seal(root)
    assert not (root / "completion_seal.json").exists()


@pytest.mark.parametrize("mutation", ["changed", "unbound"])
def test_seal_rejects_stale_or_incomplete_final_validation_binding(tmp_path, mutation):
    root = skeletal_package(tmp_path, with_figure=True)
    for directory in ("horizon_full", "span", "analysis", "interpretation"):
        (root / directory).mkdir()
    target = root / "analysis/parameter_summary.csv"
    target.write_text("original evidence")
    binding = {str(p.relative_to(root)): file_digest(p)
               for directory in ("horizon_full", "span", "analysis", "figures", "interpretation")
               for p in (root / directory).rglob("*") if p.is_file()}
    (root / "evidence/final_validation.json").write_text(json.dumps({
        "status": "passed", "validated_files": binding, "tests": {"passed": 1, "failures": 0, "errors": 0}}))
    if mutation == "changed":
        target.write_text("changed after validation")
    else:
        (root / "analysis/added_after_validation.csv").write_text("unvalidated evidence")
    with pytest.raises(ValueError, match="validation"):
        seal(root)
    assert not (root / "completion_seal.json").exists()


def test_verify_seal_rejects_changed_or_extra_evidence(tmp_path):
    root = tmp_path / "sealed"
    root.mkdir()
    artifact = root / "evidence.txt"
    artifact.write_text("original")
    frozen = {"status": "completed", "files": {"evidence.txt": file_digest(artifact)}}
    (root / "completion_seal.json").write_text(json.dumps(frozen))
    assert verify_seal(root)["status"] == "passed"
    artifact.write_text("changed")
    with pytest.raises(ValueError):
        verify_seal(root)
    artifact.write_text("original")
    (root / "unrecorded.txt").write_text("extra")
    with pytest.raises(ValueError):
        verify_seal(root)


def test_verify_seal_does_not_ignore_nested_seal_named_evidence(tmp_path):
    root = tmp_path / "sealed"
    root.mkdir()
    artifact = root / "evidence.txt"
    artifact.write_text("original")
    frozen = {"status": "completed", "files": {"evidence.txt": file_digest(artifact)}}
    (root / "completion_seal.json").write_text(json.dumps(frozen))
    (root / "nested").mkdir()
    (root / "nested/completion_seal.json").write_text("unrecorded nested evidence")
    with pytest.raises(ValueError):
        verify_seal(root)
