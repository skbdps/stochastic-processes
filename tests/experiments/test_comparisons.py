"""Independent hand calculations for MC contrasts, including missing outcomes."""
from dataclasses import replace

import numpy as np
import pandas as pd
import pytest

from ou_irregular.experiments.design import Cell, RunDesign, ESTIMATORS
from ou_irregular.experiments.comparisons import (
    regular_reference_comparisons, paired_estimator_comparisons,
)


def records(errors_by_cv, *, valid=None, cover=None, truth_mu=0):
    nrep = len(next(iter(errors_by_cv.values()))[0])
    design = RunDesign("fixture", 123, nrep,
                       tuple(Cell(10, .5, cv, mu=truth_mu) for cv in errors_by_cv), "validation")
    rows = []
    for cell in design.cells:
        for j, name in enumerate(ESTIMATORS):
            for rep, error in enumerate(errors_by_cv[cell.cv][j]):
                pv = valid[j][rep] if valid is not None else True
                iv = pv
                hit = cover[j][rep] if cover is not None else True
                row = {**cell.to_dict(), "rep": rep, "estimator": name, "point_valid": pv,
                       "nll": 0.0,
                       "interval_valid": iv, "interval_status": "valid" if iv else "unavailable"}
                for p, truth in zip(("theta", "mu", "sigma"), (1., truth_mu, .5)):
                    row.update({p: truth + error, f"{p}_se": .1 if iv else np.nan,
                                f"{p}_lower": truth if hit else truth + .1,
                                f"{p}_upper": truth + .2 if hit else truth + .3})
                rows.append(row)
    return design, pd.DataFrame(rows)


def pair(data, design, p="theta"):
    return paired_estimator_comparisons(data, design).query(
        "estimator_a == 'exact' and estimator_b == 'pfml' and parameter == @p").iloc[0]


def test_independent_regular_reference_hand_calculation():
    design, data = records({0: [[0, 1]] * 3, 1: [[0, 2]] * 3})
    r = regular_reference_comparisons(data, design).query("cv == 1 and estimator == 'exact' and parameter == 'theta'").iloc[0]
    assert r.excess_bias == .5
    assert r.excess_bias_mcse == pytest.approx(np.sqrt(5) / 2)
    assert r.rmse_ratio == 2
    assert r.rmse_ratio_mcse == pytest.approx(np.sqrt(2))
    assert r.comparison_type == "independent_cells"
    assert r.cell_n_point_valid == r.reference_n_point_valid == 2


def test_self_reference_is_not_independent_and_mu_ratio_is_defined():
    design, data = records({0: [[0, 1]] * 3})
    result = regular_reference_comparisons(data, design)
    assert result.comparison_type.eq("self_reference").all()
    assert result.excess_bias.eq(0).all() and result.excess_bias_mcse.eq(0).all()
    assert result.rmse_ratio.eq(1).all() and result.rmse_ratio_mcse.eq(0).all()
    assert result.query("parameter == 'mu'").relative_excess_bias.isna().all()


def test_paired_covariance_hand_calculation():
    design, data = records({0: [[1, 2, 3], [0, 1, 1], [0, 1, 2]]})
    r = pair(data, design)
    assert r.bias_a_common == 2 and r.bias_b_common == pytest.approx(2 / 3)
    assert r.bias_difference == pytest.approx(4 / 3)
    assert r.bias_difference_mcse == pytest.approx(1 / 3)
    assert r.mse_difference == 4
    assert r.mse_difference_mcse == pytest.approx(np.sqrt(13 / 3))
    assert r.rmse_ratio == pytest.approx(np.sqrt(7))
    assert r.rmse_ratio_mcse == pytest.approx(np.sqrt(3) / 4)
    assert r.squared_error_covariance_common == pytest.approx(11 / 6)
    expected_variance = ((49/3)/(4*(14/3)) + (1/3)/(4*(2/3))
                         - (11/6)/(2*np.sqrt(14/3)*np.sqrt(2/3))) / 3
    assert r.rmse_difference_mcse ** 2 == pytest.approx(expected_variance)


def test_common_valid_subset_recomputes_both_metrics_and_orders_by_rep():
    design, data = records({0: [[0, 1, 3, 999], [999, 0, 2, 100], [1, 1, 1, 1]]},
                          valid=[[True, True, True, False], [False, True, True, True], [True]*4])
    r = pair(data.sample(frac=1, random_state=4), design)
    assert r.n_planned == 4 and r.n_common_point_valid == 2
    assert r.n_point_valid_a == r.n_point_valid_b == 3
    assert r.bias_a_common == 2 and r.bias_b_common == 1
    assert r.bias_difference == 1 and r.bias_difference_mcse == 0
    assert r.rmse_a_common == pytest.approx(np.sqrt(5))
    assert r.rmse_b_common == pytest.approx(np.sqrt(2))


def test_paired_coverage_uses_every_planned_path_with_unavailable_zero():
    design, data = records({0: [[1]*4]*3},
                          valid=[[True, True, False, False], [True, False, False, False], [True]*4])
    r = pair(data, design)
    assert r.n_planned == 4 and r.n_common_point_valid == 1
    assert r.operational_coverage_difference == .25
    assert r.operational_coverage_difference_mcse == .25
    assert (r.n_both_cover, r.n_a_only_cover, r.n_b_only_cover, r.n_neither_cover) == (1, 1, 0, 2)
    assert np.isnan(r.bias_difference_mcse)


def test_identical_paired_outputs_have_zero_contrast_uncertainty():
    design, data = records({0: [[1, 2, 3]]*3}, cover=[[True, False, True]]*3)
    r = pair(data, design)
    assert r.bias_difference_mcse == r.mse_difference_mcse == r.rmse_difference_mcse == 0
    assert r.rmse_ratio_mcse == 0 and r.rmse_ratio == 1
    assert r.operational_coverage_difference_mcse == 0


@pytest.mark.parametrize("pattern", ["all_invalid", "disjoint", "one", "zero_rmse", "zero_denominator"])
def test_degenerate_or_missing_point_metrics(pattern):
    valid = [[True, True]]*3
    errors = [[0, 1]]*3
    if pattern == "all_invalid": valid = [[False, False]]*3
    if pattern == "disjoint": valid = [[True, False], [False, True], [True, True]]
    if pattern == "one": valid = [[True, False]]*3
    if pattern == "zero_rmse": errors = [[0, 0], [0, 1], [0, 1]]
    if pattern == "zero_denominator": errors = [[0, 1], [0, 0], [0, 1]]
    design, data = records({0: errors}, valid=valid)
    r = pair(data, design)
    assert np.isnan(r.rmse_ratio_mcse)
    if pattern in ("all_invalid", "disjoint"):
        assert r.n_common_point_valid == 0 and np.isnan(r.bias_difference)
    if pattern == "zero_rmse": assert r.rmse_ratio == 0
    if pattern == "zero_denominator": assert np.isnan(r.rmse_ratio)
    if pattern == "all_invalid": assert r.operational_coverage_difference == 0


def test_comparisons_reject_malformed_semantics_and_missing_groups():
    design, data = records({0: [[0, 1]]*3})
    for function in (paired_estimator_comparisons, regular_reference_comparisons):
        for kind in ("missing", "bad_bounds", "interval_without_point", "nonboolean"):
            bad = data.copy()
            if kind == "missing": bad = bad[bad.estimator != "euler"]
            elif kind == "bad_bounds": bad.loc[0, "theta_lower"] = np.inf
            elif kind == "interval_without_point": bad.loc[0, "point_valid"] = False
            else:
                bad["point_valid"] = bad.point_valid.astype(object)
                bad.loc[0, "point_valid"] = "False"
            with pytest.raises(ValueError): function(bad, design)


def test_regular_match_uses_full_identity_and_rejects_missing_reference():
    design, data = records({1: [[0, 1]]*3})
    with pytest.raises(ValueError, match="regular reference"):
        regular_reference_comparisons(data, design)


def test_negative_mu_truth_uses_magnitude_normalization():
    design, data = records({0: [[1, 2, 3], [0, 1, 1], [0, 1, 2]]}, truth_mu=-2)
    r = pair(data, design, "mu")
    assert r.relative_bias_difference == pytest.approx(2/3)
