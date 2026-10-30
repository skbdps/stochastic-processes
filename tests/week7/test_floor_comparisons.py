"""Paired-floor accounting tests using wholly fabricated records, no simulation."""
from dataclasses import replace
import hashlib
import json

import numpy as np
import pandas as pd
import pytest

from ou_irregular.experiments.design import ESTIMATORS
from ou_irregular.week7.floor_comparisons import COLUMNS, epsilon_comparisons
from ou_irregular.week7.floor_source import COUPLING_POLICY, anchor_cell
from ou_irregular.week7.selection import floor_design


def fixture(cv=None, replications=3):
    base = floor_design().validation(replications=replications)
    design = replace(base, cells=tuple(c for c in base.cells if cv is None or c.cv == cv))
    rows = []
    for cell in design.cells:
        anchor = anchor_cell(cell)
        words = [int(anchor.identity[i:i + 8], 16) for i in range(0, 64, 8)]
        for rep in range(replications):
            key = [4, *words, rep]
            # These are arbitrary syntactically valid hashes.  The test is an
            # arithmetic/accounting fixture, not a claim of actual source replay.
            label = f"literal:{cell.cv}:{rep}" + (f":{cell.floor_fraction}" if cell.cv else "")
            hashes = {name: hashlib.sha256((name + label).encode()).hexdigest()
                      for name in ("times_sha256", "observations_sha256")}
            for method in ESTIMATORS:
                row = {**cell.to_dict(), "rep": rep, "estimator": method, **hashes,
                       "actual_cell_id": cell.identity, "seed_anchor_cell_id": anchor.identity,
                       "seed_anchor_floor_fraction": 1e-6, "source_coupling_policy": COUPLING_POLICY,
                       "source_seed": design.seed, "source_spawn_key": json.dumps(key),
                       "gap_spawn_key": json.dumps([*key, 0]), "path_spawn_key": json.dumps([*key, 1]),
                       "point_valid": True, "optimizer_success": True, "point_status": "valid_interior",
                       "interval_valid": True, "interval_status": "valid", "nll": 1.,
                       "source_valid": True, "source_error": "", "pipeline_error": ""}
                scale = 1 if not cell.cv else {1e-8: 2, 1e-6: 1, 1e-4: 3}[cell.floor_fraction]
                error = .1 * (rep + 1) * scale
                for parameter in ("theta", "mu", "sigma"):
                    truth = cell.to_dict()["true_" + parameter]
                    row.update({parameter: truth + error, parameter + "_se": .1,
                                parameter + "_lower": truth - .2, parameter + "_upper": truth + .2})
                rows.append(row)
    return design, pd.DataFrame(rows)


def row(result, epsilon=1e-8, parameter="theta"):
    return result[(result.estimator == "exact") & (result.floor_fraction == epsilon)
                  & (result.parameter == parameter)].iloc[0]


def test_full_output_is72_rows_and_single_cv_is_a_complete18_row_subset():
    design, data = fixture()
    result = epsilon_comparisons(data, design)
    assert len(result) == 72
    assert tuple(result.columns) == COLUMNS and len(COLUMNS) == 71
    assert set(result.floor_fraction) == {1e-8, 1e-4}
    assert result.reference_floor_fraction.eq(1e-6).all()
    assert result.difference_direction.eq("epsilon_minus_registered_floor").all()
    assert result.ratio_direction.eq("epsilon_over_registered_floor").all()
    for cv in (0., 1.5, 2., 3.):
        subset_design = replace(design, cells=tuple(c for c in design.cells if c.cv == cv))
        subset = epsilon_comparisons(data[data.cv == cv], subset_design)
        assert len(subset) == 18
        pd.testing.assert_frame_equal(subset.reset_index(drop=True), result[result.cv == cv].reset_index(drop=True))


def test_order_and_input_bytes_do_not_change_pairing_or_values():
    design, data = fixture(1.5)
    original = data.copy(deep=True)
    forward = epsilon_comparisons(data, design)
    reverse = epsilon_comparisons(data.iloc[::-1], design)
    pd.testing.assert_frame_equal(forward, reverse)
    pd.testing.assert_frame_equal(data, original)


def test_mu_relative_metrics_remain_undefined_but_absolute_rmse_ratio_is_defined():
    design, data = fixture(2.)
    result = epsilon_comparisons(data, design)
    mu = result[result.parameter == "mu"]
    assert mu[[name for name in result if name.startswith("relative_")]].isna().all().all()
    assert row(result, parameter="mu").rmse_ratio == pytest.approx(2)
    assert np.isfinite(mu.bias_a_common).all() and np.isfinite(mu.mse_a_common).all()


def test_cv0_shared_negative_control_has_unit_ratio_and_zero_paired_uncertainty():
    design, data = fixture(0.)
    result = epsilon_comparisons(data, design)
    for name in ("bias_difference", "mse_difference", "rmse_difference", "operational_coverage_difference"):
        assert result[name].eq(0.).all()
        assert result[name + "_mcse"].eq(0.).all()
    assert result.rmse_ratio.eq(1.).all()
    assert result.rmse_ratio_mcse.eq(0.).all()
    assert result.n_both_cover.eq(design.replications).all()


def test_invalid_numerical_candidates_survive_accounting_without_entering_metrics():
    design, data = fixture(1.5)
    # Different invalid subsets; one common point remains, but all three planned
    # pairs still enter operational coverage.  Retain the invalid candidates.
    a_bad = (data.floor_fraction == 1e-8) & data.rep.eq(0)
    b_bad = (data.floor_fraction == 1e-6) & data.rep.eq(2)
    bad = a_bad | b_bad
    data.loc[bad, ["point_valid", "optimizer_success", "interval_valid"]] = False
    data.loc[bad, "point_status"] = "boundary"
    data.loc[bad, "interval_status"] = "point_invalid"
    data.loc[bad, ["theta", "mu", "sigma"]] = 999.
    before = data.copy(deep=True)
    r = row(epsilon_comparisons(data, design))
    assert r.n_planned == 3 and r.n_point_valid_a == r.n_point_valid_b == 2
    assert r.n_interval_valid_a == r.n_interval_valid_b == 2
    assert r.n_common_point_valid == r.n_common_interval_valid == 1
    assert r.bias_a_common == pytest.approx(.4)
    assert r.bias_b_common == pytest.approx(.2)
    assert np.isnan(r.bias_difference_mcse) and np.isnan(r.rmse_ratio_log_mcse)
    assert (r.n_both_cover, r.n_a_only_cover, r.n_b_only_cover, r.n_neither_cover) == (1, 1, 1, 0)
    pd.testing.assert_frame_equal(data, before)


@pytest.mark.parametrize("change", ["whole_cell", "whole_method", "one_rep", "duplicate", "missing_floor_design"])
def test_complete_cartesian_keys_and_all_three_floors_are_required(change):
    design, data = fixture(1.5)
    if change == "whole_cell":
        data = data[data.floor_fraction != 1e-8]
    elif change == "whole_method":
        data = data[data.estimator != "pfml"]
    elif change == "one_rep":
        data = data.drop(data.index[0])
    elif change == "duplicate":
        data = pd.concat([data, data.iloc[[0]]], ignore_index=True)
    else:
        design = replace(design, cells=design.cells[1:])
        data = data[data.floor_fraction != 1e-8]
    with pytest.raises(ValueError):
        epsilon_comparisons(data, design)


@pytest.mark.parametrize("change", ["hash", "estimate", "invalid_candidate", "nll", "interval", "status"])
def test_cv0_inconsistent_evidence_is_rejected(change):
    design, data = fixture(0.)
    mask = data.floor_fraction.eq(1e-8)
    if change == "hash":
        data.loc[mask, "times_sha256"] = "f" * 64
    elif change == "estimate":
        data.loc[mask, "theta"] = 1.5
    elif change == "nll":
        data.loc[mask, "nll"] = 2.
    elif change == "interval":
        data.loc[mask, "theta_lower"] = 1.2
    elif change == "status":
        data.loc[mask, "interval_valid"] = False
        data.loc[mask, "interval_status"] = "hessian_failure"
    else:
        data.loc[:, ["point_valid", "optimizer_success", "interval_valid"]] = False
        data.loc[:, "point_status"] = "boundary"
        data.loc[:, "interval_status"] = "point_invalid"
        data.loc[mask, "theta"] = 999.
    with pytest.raises(ValueError, match="negative-control"):
        epsilon_comparisons(data, design)


def test_valid_overflow_blocks_table_while_invalid_candidate_is_not_silently_dropped():
    design, data = fixture(1.5)
    data.loc[0, "theta"] = 1e308
    with pytest.raises(ValueError, match="overflow"):
        epsilon_comparisons(data, design)
    data.loc[0, ["point_valid", "optimizer_success", "interval_valid"]] = False
    data.loc[0, "point_status"] = "boundary"
    data.loc[0, "interval_status"] = "point_invalid"
    result = epsilon_comparisons(data, design)
    assert row(result).n_planned == 3 and row(result).n_point_valid_a == 2
    assert data.loc[0, "theta"] == 1e308
