"""Hand calculations for classified exclusions and shared-path comparisons."""
import numpy as np
import pandas as pd
import pytest

from ou_irregular.metrics import (
    FIT_STATUSES, LEGACY_VALIDITY_POLICY, VALIDITY_POLICY,
    paired_estimator_comparisons, point_validity, summarize_results, validity_counts,
)


def rows(estimator="exact", errors=(1.0, 3.0), cell_id="cell", cv=0.0):
    return pd.DataFrame([
        {"cell_id": cell_id, "theta_mean_gap": 0.5, "cv": cv, "n": 250,
         "mean_gap": 0.5, "true_theta": 1.0, "true_mu": 0.0, "true_sigma": 0.5,
         "rep": rep, "estimator": estimator, "success": True, "optimizer_success": True,
         "fit_status": "valid_interior", "valid_for_point_summary": True,
         "reason": "", "solver": "fixture", "nll": 0.0,
         "theta": 1.0 + error, "mu": error, "sigma": 0.5 + abs(error)}
        for rep, error in enumerate(errors)
    ])


def invalidate(data, index, status, optimizer_success=False):
    data.loc[index, ["success", "optimizer_success"]] = optimizer_success
    data.loc[index, "fit_status"] = status
    data.loc[index, "valid_for_point_summary"] = False
    data.loc[index, "reason"] = "hand-checkable excluded fixture"


def test_status_counts_reconcile_and_boundary_success_is_excluded():
    data = rows(errors=(1.0, 3.0, 100.0, 100.0, 100.0, 100.0, 100.0))
    for index, status in enumerate(FIT_STATUSES[1:], start=2):
        invalidate(data, index, status, optimizer_success=status.startswith("boundary"))
    counts = validity_counts(data)
    assert counts["validity_policy"] == VALIDITY_POLICY
    assert (counts["n_planned"], counts["n_valid"], counts["n_excluded"],
            counts["n_optimizer_success"]) == (7, 2, 5, 4)
    assert counts["n_valid_interior"] == 2
    assert all(counts[f"n_{status}"] == 1 for status in FIT_STATUSES[1:])
    assert sum(counts[f"n_{status}"] for status in FIT_STATUSES) == 7
    assert point_validity(data).tolist() == [True, True, False, False, False, False, False]
    summary = summarize_results(data)
    assert summary.n_total.eq(7).all() and summary.n_valid.eq(2).all()
    assert summary.n_success.eq(2).all() and summary.n_optimizer_success.eq(4).all()
    assert summary.n_failed.eq(5).all() and summary.failure_rate.eq(5 / 7).all()
    assert summary.metric_scope.eq("conditional_on_point_valid").all()
    theta = summary.loc[summary.parameter.eq("theta")].iloc[0]
    assert theta.bias == 2 and theta.rmse == pytest.approx(np.sqrt(5))
    assert theta.mcse_bias == 1 and theta.mcse_rmse == pytest.approx(2 / np.sqrt(5))
    assert summary.coverage.isna().all() and summary.n_ci_valid.eq(0).all()
    # Invalid candidate 101 is preserved and never changed to true theta=1.
    assert data.loc[2:, "theta"].eq(101).all()


def test_all_invalid_cell_retains_each_row_and_missing_metrics():
    data = rows()
    invalidate(data, 0, "boundary", optimizer_success=True)
    invalidate(data, 1, "simulation_failure")
    data.loc[1, ["theta", "mu", "sigma", "nll"]] = np.nan
    summary = summarize_results(data)
    assert len(summary) == 3
    assert summary.n_planned.eq(2).all() and summary.n_valid.eq(0).all()
    assert summary.n_optimizer_success.eq(1).all() and summary.failure_rate.eq(1).all()
    assert summary[["bias", "rmse", "mcse_bias", "mcse_rmse", "coverage"]].isna().all().all()


@pytest.mark.parametrize("column,value,message", [
    ("success", False, "alias"),
    ("valid_for_point_summary", False, "agree"),
    ("fit_status", "boundary", "agree"),
    ("fit_status", "unknown", "fit_status"),
    ("sigma", 0.0, "admissible"),
    ("theta", -1.0, "admissible"),
    ("mu", np.inf, "admissible"),
    ("nll", np.inf, "admissible"),
    ("optimizer_success", "true", "Boolean"),
])
def test_contradictory_classifications_are_rejected(column, value, message):
    data = rows()
    if isinstance(value, str):
        data[column] = data[column].astype(object)
    data.loc[0, column] = value
    with pytest.raises(ValueError, match=message):
        point_validity(data)


def test_partial_or_mixed_schema_is_rejected():
    data = rows().drop(columns="optimizer_success")
    with pytest.raises(ValueError, match="all classification"):
        point_validity(data)
    legacy = rows().drop(columns=["optimizer_success", "fit_status", "valid_for_point_summary"])
    with pytest.raises(ValueError, match="Boolean"):
        point_validity(pd.concat([legacy, rows()], ignore_index=True))


def test_historical_policy_is_explicit_and_never_retroactively_reclassified():
    legacy = rows().drop(columns=["optimizer_success", "fit_status", "valid_for_point_summary"])
    # Historical policy checked finiteness, not domains. Retain that exact rule;
    # classifying this candidate under v2 requires a separately versioned refit.
    legacy.loc[0, "theta"] = -1.0
    counts = validity_counts(legacy)
    assert counts["validity_policy"] == LEGACY_VALIDITY_POLICY
    assert counts["n_valid"] == counts["n_legacy_usable"] == 2
    assert counts["n_valid_interior"] == 0
    assert summarize_results(legacy).metric_scope.eq("conditional_on_success").all()


def test_paired_differences_use_common_valid_paths_and_paired_mcse():
    exact = rows("exact", errors=(1.0, 3.0, 100.0))
    euler = rows("euler", errors=(0.0, 2.0, 4.0))
    pfml = rows("pfml", errors=(2.0, 4.0, 8.0))
    invalidate(exact, 2, "boundary", optimizer_success=True)
    result = paired_estimator_comparisons(pd.concat([exact, euler, pfml]))
    assert len(result) == 9  # Three estimator pairs, three parameters each.
    pair = result.loc[(result.estimator_a == "euler") & (result.estimator_b == "exact")
                      & result.parameter.eq("theta")].iloc[0]
    assert pair.n_planned_pairs == pair.n_recorded_pairs == 3
    assert (pair.n_valid_a, pair.n_valid_b, pair.n_common_valid) == (3, 2, 2)
    assert pair.n_excluded_pairs == 1
    assert pair.bias_a_common == 1 and pair.bias_b_common == 2
    # Shared-path differences are exactly -1,-1, despite varying marginal errors.
    assert pair.bias_difference == -1 and pair.mcse_bias_difference == 0
    # Paired squared-error differences -1,-5 have mean -3 and SE 2.
    assert pair.mse_difference == -3 and pair.mcse_mse_difference == 2
    assert pair.difference_direction == "estimator_a_minus_estimator_b"


def test_comparisons_never_pair_independent_cv_cells():
    regular = pd.concat([rows("exact"), rows("pfml")])
    irregular = pd.concat([rows("exact", (1.0, 3.0), "independent", 2.0),
                           rows("pfml", (8.0, 8.0), "independent", 2.0)])
    result = paired_estimator_comparisons(pd.concat([regular, irregular]))
    assert len(result) == 6
    assert result.n_common_valid.eq(2).all()
    assert result.loc[result.cell_id.eq("cell"), "bias_difference"].eq(0).all()
    theta = result.loc[result.cell_id.eq("independent") & result.parameter.eq("theta")].iloc[0]
    assert theta.bias_difference == -6


def test_pair_denominators_expose_missing_rows_and_all_invalid_pairs():
    first, second = rows("exact"), rows("pfml").iloc[[0]].copy()
    invalidate(first, 0, "degenerate")
    invalidate(first, 1, "numerical_failure")
    result = paired_estimator_comparisons(pd.concat([first, second]))
    assert result.n_planned_pairs.eq(2).all() and result.n_recorded_pairs.eq(1).all()
    assert result.n_common_valid.eq(0).all() and result.n_excluded_pairs.eq(2).all()
    assert result[["bias_difference", "mcse_bias_difference", "mse_difference",
                   "mcse_mse_difference"]].isna().all().all()


def test_single_common_valid_path_has_difference_but_no_mcse():
    result = paired_estimator_comparisons(pd.concat([rows("exact").iloc[[0]],
                                                    rows("pfml", (2.0,)).iloc[[0]]]))
    assert result.n_common_valid.eq(1).all()
    assert result.bias_difference.eq(-1).all()
    assert result.mcse_bias_difference.isna().all()
    assert result.mcse_mse_difference.isna().all()
