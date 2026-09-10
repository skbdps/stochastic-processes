"""Hand-computable checks for denominators, Monte Carlo metrics, and controls."""

import numpy as np
import pandas as pd
import pytest

from ou_irregular.metrics import add_twin_comparisons, summarize_results
from ou_irregular.plots import plot_summary


def _results(cell_id="regular", cv=0.0, estimates=(1.0, 3.0, 100.0)):
    return pd.DataFrame([
        {"cell_id": cell_id, "theta_mean_gap": 0.2, "cv": cv, "n": 100,
         "mean_gap": 0.1, "true_theta": 2.0, "true_mu": 0.0, "true_sigma": 0.5,
         "rep": rep, "estimator": "exact", "success": rep < 2,
         "theta": estimate, "mu": (-1.0, 3.0, 100.0)[rep],
         "sigma": (0.25, 0.75, 100.0)[rep]}
        for rep, estimate in enumerate(estimates)
    ])


def _parameter(summary, parameter):
    return summary.loc[summary.parameter.eq(parameter)].iloc[0]


def test_summary_keeps_failures_and_uses_hand_computable_errors():
    summary = summarize_results(_results())
    theta = _parameter(summary, "theta")
    assert theta.n_total == 3
    assert theta.n_success == 2
    assert theta.n_failed == 1
    assert theta.failure_rate == pytest.approx(1 / 3)
    assert theta.bias == 0
    assert theta.relative_bias == 0
    assert theta.rmse == 1
    assert theta.relative_rmse == 0.5
    assert theta.mcse_bias == pytest.approx(1.0)
    assert theta.metric_scope == "conditional_on_success"
    mu = _parameter(summary, "mu")
    assert mu.bias == 1
    assert mu.rmse == pytest.approx(np.sqrt(5))
    assert np.isnan(mu.relative_bias)
    assert np.isnan(mu.relative_rmse)
    assert mu.mcse_bias == pytest.approx(2.0)
    sigma = _parameter(summary, "sigma")
    assert sigma.rmse == 0.25
    assert sigma.relative_rmse == 0.5


def test_no_intervals_means_missing_coverage_not_fake_nominal_coverage():
    summary = summarize_results(_results())
    assert summary.coverage.isna().all()
    assert summary.coverage_mcse.isna().all()
    assert summary.n_ci_valid.eq(0).all()
    assert summary.ci_status.eq("not_computed_week3").all()


def test_rmse_mcse_delta_method_and_relative_scaling_are_hand_computable():
    summary = summarize_results(_results(estimates=(3.0, 5.0, 100.0)))
    theta = _parameter(summary, "theta")
    # Successful errors are 1,3; squared errors 1,9 have sample SD sqrt(32).
    # RMSE=sqrt(5), so delta-method MCSE=sqrt(32)/(2*sqrt(2)*sqrt(5)).
    assert theta.n_success == 2
    assert theta.rmse == pytest.approx(np.sqrt(5))
    assert theta.mcse_bias == pytest.approx(1.0)
    assert theta.mcse_rmse == pytest.approx(2.0 / np.sqrt(5))
    assert theta.relative_mcse_bias == pytest.approx(0.5)
    assert theta.relative_mcse_rmse == pytest.approx(1.0 / np.sqrt(5))
    mu = _parameter(summary, "mu")
    assert mu.mcse_rmse == pytest.approx(2.0 / np.sqrt(5))
    assert np.isnan(mu.relative_mcse_bias)
    assert np.isnan(mu.relative_mcse_rmse)


def test_zero_rmse_has_undefined_delta_method_mcse():
    theta = _parameter(summarize_results(_results(estimates=(2.0, 2.0, 100.0))), "theta")
    assert theta.n_success == 2
    assert theta.rmse == 0
    assert theta.mcse_bias == 0
    assert np.isnan(theta.mcse_rmse)
    assert np.isnan(theta.relative_mcse_rmse)


def test_single_success_cannot_estimate_monte_carlo_standard_errors():
    data = _results()
    data.loc[1, "success"] = False
    summary = summarize_results(data)
    assert summary.n_success.eq(1).all()
    assert summary.rmse.notna().all()
    assert summary[["mcse_bias", "mcse_rmse", "relative_mcse_bias",
                    "relative_mcse_rmse"]].isna().all().all()


def test_actual_supplied_intervals_use_only_successful_valid_bounds():
    data = _results()
    data["theta_lower"] = [1.0, 2.5, 1.0]
    data["theta_upper"] = [2.0, 4.0, 3.0]
    theta = _parameter(summarize_results(data), "theta")
    # The first interval includes the truth exactly at its upper endpoint;
    # the failed third fit contributes neither a hit nor a denominator count.
    assert theta.n_ci_valid == 2
    assert theta.coverage == 0.5
    assert theta.coverage_mcse == pytest.approx(np.sqrt(0.25 / 2))
    assert theta.ci_status == "computed"
    assert theta.coverage_scope == "conditional_on_success_and_valid_interval"


@pytest.mark.parametrize("bad_bounds", [(np.nan, 5.0), (4.0, 3.0), (-np.inf, 5.0)])
def test_invalid_intervals_do_not_count_in_coverage_denominator(bad_bounds):
    data = _results()
    data["theta_lower"] = [1.0, bad_bounds[0], 1.0]
    data["theta_upper"] = [2.0, bad_bounds[1], 3.0]
    theta = _parameter(summarize_results(data), "theta")
    assert theta.n_ci_valid == 1
    assert theta.coverage == 1


def test_half_supplied_interval_is_an_error():
    data = _results()
    data["theta_lower"] = 0.0
    with pytest.raises(ValueError, match="supply both theta_lower and theta_upper"):
        summarize_results(data)


def test_no_valid_intervals_are_explicit():
    data = _results()
    data["theta_lower"] = np.nan
    data["theta_upper"] = np.nan
    theta = _parameter(summarize_results(data), "theta")
    assert theta.ci_status == "no_valid_intervals"
    assert theta.n_ci_valid == 0
    assert np.isnan(theta.coverage)


def test_nonfinite_estimate_makes_whole_fit_unusable():
    data = _results()
    data.loc[0, "sigma"] = np.nan
    summary = summarize_results(data)
    assert summary.n_success.eq(1).all()
    assert summary.n_failed.eq(2).all()
    assert summary.mcse_bias.isna().all()
    assert _parameter(summary, "theta").bias == 1.0


def test_all_failed_returns_missing_metrics_with_full_failure_rate():
    data = _results()
    data["success"] = False
    summary = summarize_results(data)
    assert summary.n_success.eq(0).all()
    assert summary.failure_rate.eq(1).all()
    assert summary[["bias", "rmse", "mcse_bias", "coverage"]].isna().all().all()


def test_summary_rejects_double_counted_replications_and_inconsistent_metadata():
    data = _results()
    with pytest.raises(ValueError, match="duplicate"):
        summarize_results(pd.concat([data, data.iloc[[0]]]))
    data.loc[0, "mean_gap"] = 0.7
    with pytest.raises(ValueError, match="inconsistent"):
        summarize_results(data)


def test_twin_comparisons_match_nominal_design_and_preserve_row_count():
    regular = _results()
    irregular = _results("irregular", 1.0, (2.0, 4.0, 100.0))
    summary = summarize_results(pd.concat([regular, irregular], ignore_index=True))
    compared = add_twin_comparisons(summary)
    assert len(compared) == 6
    theta = compared.loc[compared.cell_id.eq("irregular") & compared.parameter.eq("theta")].iloc[0]
    assert theta.twin_cell_id == "regular"
    assert theta.twin_bias == 0
    assert theta.excess_bias == 1
    assert theta.excess_relative_bias == 0.5
    assert theta.twin_rmse == 1
    assert theta.rmse_ratio == pytest.approx(np.sqrt(2))
    assert theta.comparison_status == "matched"
    assert np.isnan(theta.coverage_gap)
    baseline = compared.loc[compared.cell_id.eq("regular")]
    assert baseline.rmse_ratio.eq(1).all()


def test_missing_twin_is_explicit_and_different_truth_does_not_match():
    regular = _results()
    irregular = _results("irregular", 1.0)
    irregular["true_mu"] = 0.1
    compared = add_twin_comparisons(summarize_results(pd.concat([regular, irregular])))
    unmatched = compared.loc[compared.cell_id.eq("irregular")]
    assert unmatched.comparison_status.eq("missing_twin").all()
    assert unmatched.twin_cell_id.isna().all()
    assert unmatched.rmse_ratio.isna().all()


def test_zero_twin_rmse_produces_missing_ratio():
    regular = _results(estimates=(2.0, 2.0, 100.0))
    irregular = _results("irregular", 1.0)
    compared = add_twin_comparisons(summarize_results(pd.concat([regular, irregular])))
    assert compared.loc[compared.parameter.eq("theta"), "rmse_ratio"].isna().all()


def test_ambiguous_twins_are_rejected_instead_of_multiplying_rows():
    twins = pd.concat([_results(), _results("another_regular")], ignore_index=True)
    with pytest.raises(ValueError, match="multiple equidistant twins"):
        add_twin_comparisons(summarize_results(twins))


def test_twin_coverage_gap_only_uses_supplied_coverage():
    regular = _results()
    irregular = _results("irregular", 1.0)
    regular["theta_lower"], regular["theta_upper"] = 0.0, 4.0
    irregular["theta_lower"], irregular["theta_upper"] = 2.5, 4.0
    compared = add_twin_comparisons(summarize_results(pd.concat([regular, irregular])))
    theta = compared.loc[compared.cell_id.eq("irregular") & compared.parameter.eq("theta")].iloc[0]
    assert theta.coverage_gap == -1.0
    assert compared.loc[compared.parameter.eq("mu"), "coverage_gap"].isna().all()


def test_plot_summary_writes_three_nonempty_pngs(tmp_path):
    summary = summarize_results(pd.concat([_results(), _results("irregular", 1.0)]))
    paths = plot_summary(summary, tmp_path)
    assert len(paths) == 3
    for path in paths:
        assert path.exists()
        assert path.stat().st_size > 10_000
        assert path.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
