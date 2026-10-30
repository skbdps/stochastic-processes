"""Independent accounting and identity checks with literal, nonrandom records."""
from copy import deepcopy
from dataclasses import replace
import hashlib
import json

import numpy as np
import pandas as pd
import pytest

from ou_irregular.experiments.design import ESTIMATORS
from ou_irregular.experiments.provenance import canonical
from ou_irregular.week7 import floor_comparisons as contrasts
from ou_irregular.week7.selection import floor_design


POLICY = "week7_paired_floor_v1_shared_gamma_stationary_x0_and_standard_normal_innovations"
PARAMETERS = ("theta", "mu", "sigma")


def records(replications=3, cvs=(0., 1.5, 2., 3.), *, all_invalid=False, errors_zero=False):
    """Canonical identities and explicit spawn keys, with no RNG construction."""
    design = floor_design().validation(replications)
    design = replace(design, cells=tuple(c for c in design.cells if c.cv in cvs))
    result = []
    for cell in design.cells:
        anchor = replace(cell, floor_fraction=1e-6)
        words = tuple(int(anchor.identity[i:i + 8], 16) for i in range(0, 64, 8))
        for rep in range(replications):
            address = [4, *words, rep]
            # CV=0 is the exact shared negative control. Irregular physical
            # arrays may differ despite correctly shared innovation addresses.
            physical_key = anchor.identity if cell.cv == 0 else cell.identity
            hashes = {name: hashlib.sha256(f"synthetic:{physical_key}:{rep}:{name}".encode()).hexdigest()
                      for name in ("times_sha256", "observations_sha256")}
            for method in ESTIMATORS:
                row = {**cell.to_dict(), "rep": rep, "estimator": method,
                       "actual_cell_id": cell.identity, "seed_anchor_cell_id": anchor.identity,
                       "seed_anchor_floor_fraction": 1e-6, "source_coupling_policy": POLICY,
                       "source_seed": design.seed,
                       "source_spawn_key": json.dumps(address, separators=(",", ":")),
                       "gap_spawn_key": json.dumps([*address, 0], separators=(",", ":")),
                       "path_spawn_key": json.dumps([*address, 1], separators=(",", ":")),
                       **hashes, "source_valid": True, "source_error": "", "pipeline_error": "",
                       "point_valid": True, "point_status": "valid_interior", "optimizer_success": True,
                       "interval_valid": True, "interval_status": "valid", "nll": 1.,
                       "point_diagnostics_json": '{"fixture":"literal","candidate_retained":true}',
                       "interval_diagnostics_json": '{"fixture":"literal"}'}
                for parameter in PARAMETERS:
                    truth = row[f"true_{parameter}"]
                    row[parameter] = truth + (0. if errors_zero else .1 + .01 * rep)
                    row[f"{parameter}_se"] = .1
                    row[f"{parameter}_lower"] = truth - .25
                    row[f"{parameter}_upper"] = truth + .25
                if all_invalid:
                    invalidate(row, huge_candidate=True)
                result.append(row)
    return design, result


def invalidate(row, *, huge_candidate=False):
    row.update(point_valid=False, point_status="boundary_upper_theta", optimizer_success=False,
               interval_valid=False, interval_status="invalid_point_fit", nll=None)
    for parameter in PARAMETERS:
        row[parameter] = 1e300 if huge_candidate else row[parameter]
        row[f"{parameter}_se"] = None
        row[f"{parameter}_lower"] = None
        row[f"{parameter}_upper"] = None


def picked(result, *, epsilon=1e-8, cv=3., method="exact", parameter="theta"):
    return result[(result.floor_fraction == epsilon) & (result.cv == cv)
                  & (result.estimator == method) & (result.parameter == parameter)].iloc[0]


def test_complete_full_grid_and_single_cv_subsets_keep_exact_output_keys():
    design, rows = records()
    result = contrasts.epsilon_comparisons(rows, design)
    assert len(result) == 72
    expected = {(cell.identity, method, parameter) for cell in design.cells if cell.floor_fraction != 1e-6
                for method in ESTIMATORS for parameter in PARAMETERS}
    assert set(result[["cell_id", "estimator", "parameter"]].itertuples(index=False, name=None)) == expected
    assert not result.cell_id.eq(result.reference_cell_id).any()
    assert result.n_planned.eq(3).all()
    assert result.point_subset.eq("intersection_of_valid_point_fits").all()
    assert result.coverage_subset.eq("all_planned_paths_unavailable_counts_as_not_covering").all()
    for cv in (0., 1.5, 2., 3.):
        subset = replace(design, cells=tuple(c for c in design.cells if c.cv == cv))
        selected = [row for row in rows if row["cv"] == cv]
        small = contrasts.epsilon_comparisons(selected, subset)
        assert len(small) == 18
        pd.testing.assert_frame_equal(small.reset_index(drop=True), result[result.cv == cv].reset_index(drop=True))


@pytest.mark.parametrize("change", ["whole_cell", "whole_method", "whole_replication", "one_fit", "duplicate",
                                   "unplanned_replication", "unplanned_method", "empty"])
def test_every_planned_cartesian_key_is_required_including_absent_groups(change):
    design, rows = records(cvs=(3.,))
    if change == "whole_cell":
        rows = [row for row in rows if row["cell_id"] != design.cells[0].identity]
    elif change == "whole_method":
        rows = [row for row in rows if row["estimator"] != "euler"]
    elif change == "whole_replication":
        rows = [row for row in rows if row["rep"] != 2]
    elif change == "one_fit":
        rows.pop()
    elif change == "duplicate":
        rows.append(deepcopy(rows[0]))
    elif change == "unplanned_replication":
        rows[0]["rep"] = 3
    elif change == "unplanned_method":
        rows[0]["estimator"] = "substituted_estimator"
    else:
        rows = []
    with pytest.raises(ValueError):
        contrasts.epsilon_comparisons(rows, design)


@pytest.mark.parametrize("change", ["missing_variant", "headline_scope", "wrong_n", "wrong_q", "wrong_cv",
                                   "wrong_floor", "wrong_truth"])
def test_design_cannot_change_the_physical_protocol_or_omit_a_floor(change):
    design, rows = records(cvs=(3.,))
    if change == "missing_variant":
        design = replace(design, cells=design.cells[:-1])
    elif change == "headline_scope":
        design = replace(design, scope="headline")
    else:
        field, value = {"wrong_n": ("n", 999), "wrong_q": ("q", .1), "wrong_cv": ("cv", 1.),
                        "wrong_floor": ("floor_fraction", 1e-7), "wrong_truth": ("theta", 2.)}[change]
        design = replace(design, cells=(replace(design.cells[0], **{field: value}), *design.cells[1:]))
    with pytest.raises(ValueError):
        contrasts.epsilon_comparisons(rows, design)


@pytest.mark.parametrize("field,bad", [
    ("actual_cell_id", "a" * 64), ("seed_anchor_cell_id", "b" * 64),
    ("seed_anchor_floor_fraction", 1e-8), ("seed_anchor_floor_fraction", "1e-6"),
    ("source_coupling_policy", "replication_index_only"), ("source_seed", 20261009),
    ("source_seed", "21261009"), ("source_seed", True),
    ("rep", True), ("rep", 0.), ("rep", "0"),
    ("n", 1000.), ("cv", "3"), ("true_mu", False), ("floor_fraction", "1e-8"),
    ("times_sha256", "bad"), ("observations_sha256", "X" * 64),
    ("point_valid", "False"), ("interval_valid", 1), ("optimizer_success", "True"),
    ("point_status", None), ("interval_status", None),
    ("point_status", "boundary_upper_theta"), ("interval_status", "invalid_point_fit"),
    ("source_valid", False), ("source_valid", 1), ("source_error", "unavailable source"),
    ("pipeline_error", None), ("theta", "1.1"), ("mu", False), ("nll", "1.0"),
    ("sigma_se", "0.1"), ("theta_lower", "0.75"), ("theta", -1.), ("sigma", 0.),
    ("nll", float("nan")), ("theta_upper", .5), ("sigma_se", 0.),
])
def test_identity_status_and_numeric_types_are_checked_before_contrasts(field, bad):
    design, rows = records(cvs=(3.,))
    rows[0][field] = bad
    with pytest.raises(ValueError):
        contrasts.epsilon_comparisons(rows, design)


@pytest.mark.parametrize("field", ["point_valid", "interval_status", "optimizer_success", "nll", "theta_se",
                                  "actual_cell_id", "source_spawn_key", "gap_spawn_key", "path_spawn_key",
                                  "times_sha256", "observations_sha256", "seed_anchor_cell_id"])
def test_required_fit_and_source_fields_cannot_disappear(field):
    design, rows = records(cvs=(3.,))
    for row in rows:
        row.pop(field)
    with pytest.raises(ValueError):
        contrasts.epsilon_comparisons(rows, design)


@pytest.mark.parametrize("change", ["namespace", "anchor_word", "replication", "gap_child", "path_child",
                                   "boolean_word", "floating_words", "extra_word", "actual_not_anchor", "not_json"])
def test_exact_anchor_spawn_keys_are_required_not_equal_replication_labels(change):
    design, rows = records(cvs=(3.,))
    row = rows[0]
    field = "source_spawn_key"
    key = json.loads(row[field])
    if change == "namespace":
        key[0] = 5
    elif change == "anchor_word":
        key[1] += 1
    elif change == "replication":
        key[-1] = 1
    elif change in ("gap_child", "path_child"):
        field = "gap_spawn_key" if change == "gap_child" else "path_spawn_key"
        key = json.loads(row[field])
        key[-1] = 1 - key[-1]
    elif change == "boolean_word":
        key[-1] = False
    elif change == "floating_words":
        key = [float(v) for v in key]
    elif change == "extra_word":
        key.append(0)
    elif change == "actual_not_anchor":
        key = [4, *(int(row["cell_id"][i:i + 8], 16) for i in range(0, 64, 8)), 0]
    row[field] = "not-json" if change == "not_json" else json.dumps(key)
    with pytest.raises(ValueError):
        contrasts.epsilon_comparisons(rows, design)


def test_every_estimator_must_share_its_variant_source_hashes():
    design, rows = records(cvs=(3.,))
    rows[1]["observations_sha256"] = "f" * 64
    with pytest.raises(ValueError, match="source"):
        contrasts.epsilon_comparisons(rows, design)


def test_irregular_variants_may_have_different_or_coincident_physical_arrays():
    design, rows = records(cvs=(3.,))
    assert len({row["observations_sha256"] for row in rows if row["rep"] == 0}) == 3
    different = contrasts.epsilon_comparisons(rows, design)
    for row in rows:
        for field in ("times_sha256", "observations_sha256"):
            row[field] = hashlib.sha256(f"coincident:{row['rep']}:{field}".encode()).hexdigest()
    coincident = contrasts.epsilon_comparisons(rows, design)
    pd.testing.assert_frame_equal(different, coincident)


@pytest.mark.parametrize("field,value", [("observations_sha256", "f" * 64), ("theta", 1.2),
                                       ("point_diagnostics_json", '{"candidate":"different"}')])
def test_cv_zero_requires_full_negative_control_equality(field, value):
    design, rows = records(cvs=(0.,))
    # Mutate all methods for one physical variant/source so the within-variant
    # source-sharing check still passes; cross-floor control equality must fail.
    for row in rows:
        if row["cell_id"] == design.cells[0].identity and row["rep"] == 0:
            row[field] = value
    with pytest.raises(ValueError, match="negative.control"):
        contrasts.epsilon_comparisons(rows, design)


def test_all_invalid_candidates_survive_denominators_without_overflow_or_input_mutation():
    design, rows = records(all_invalid=True)
    before = canonical(rows)
    result = contrasts.epsilon_comparisons(rows, design)
    assert len(result) == 72
    assert result.n_planned.eq(3).all()
    for field in ("n_point_valid_a", "n_point_valid_b", "n_interval_valid_a", "n_interval_valid_b",
                  "n_common_point_valid", "n_common_interval_valid", "n_covered_a", "n_covered_b"):
        assert result[field].eq(0).all()
    for field in ("bias_a_common", "bias_b_common", "rmse_a_common", "rmse_b_common", "bias_difference",
                  "bias_difference_mcse", "rmse_ratio", "rmse_ratio_mcse"):
        assert result[field].isna().all()
    assert result.operational_coverage_difference.eq(0).all()
    assert result.operational_coverage_difference_mcse.eq(0).all()
    assert result.n_neither_cover.eq(3).all()
    assert canonical(rows) == before
    assert all(row["theta"] == 1e300 for row in rows)


@pytest.mark.parametrize("field,bad", [("point_status", None), ("point_status", ""),
                                      ("interval_status", None), ("interval_status", "")])
def test_invalid_outcomes_still_require_explicit_status_evidence(field, bad):
    design, rows = records(cvs=(3.,), all_invalid=True)
    rows[0][field] = bad
    with pytest.raises(ValueError):
        contrasts.epsilon_comparisons(rows, design)


def test_negative_natural_scale_lower_limits_are_retained_and_counted():
    design, rows = records(cvs=(3.,))
    for row in rows:
        row["theta_lower"] = -5.
        row["sigma_lower"] = -5.
    result = contrasts.epsilon_comparisons(rows, design)
    assert result.n_covered_a.eq(3).all() and result.n_covered_b.eq(3).all()
    assert all(row["theta_lower"] == row["sigma_lower"] == -5. for row in rows)


def test_valid_interval_cannot_contradict_a_retained_pipeline_error():
    design, rows = records(cvs=(3.,))
    rows[0]["pipeline_error"] = "RuntimeError: synthetic failed inference"
    with pytest.raises(ValueError, match="(?i)(pipeline|interval|error)"):
        contrasts.epsilon_comparisons(rows, design)


def test_inference_exception_retains_valid_point_but_counts_interval_as_unavailable():
    design, rows = records(cvs=(3.,))
    rows[0].update(interval_valid=False, interval_status="pipeline_exception",
                   pipeline_error="RuntimeError: synthetic failed inference")
    for parameter in PARAMETERS:
        for suffix in ("lower", "upper", "se"):
            rows[0][f"{parameter}_{suffix}"] = None
    before = canonical(rows)
    result = picked(contrasts.epsilon_comparisons(rows, design))
    assert result.n_point_valid_a == result.n_common_point_valid == 3
    assert result.n_interval_valid_a == result.n_common_interval_valid == 2
    assert result.n_interval_valid_b == 3
    assert result.n_covered_a == 2 and result.n_covered_b == 3
    assert result.operational_coverage_difference == pytest.approx(-1 / 3)
    assert canonical(rows) == before


def test_disjoint_point_subsets_keep_own_counts_and_all_planned_coverage():
    design, rows = records(replications=2, cvs=(3.,))
    for row in rows:
        if ((row["floor_fraction"] == 1e-8 and row["rep"] == 1)
                or (row["floor_fraction"] == 1e-6 and row["rep"] == 0)):
            invalidate(row)
            # Candidate endpoints look covering, but the interval is unavailable.
            for parameter in PARAMETERS:
                truth = row[f"true_{parameter}"]
                row[f"{parameter}_lower"], row[f"{parameter}_upper"] = truth - 1., truth + 1.
    result = picked(contrasts.epsilon_comparisons(rows, design))
    assert result.n_planned == 2
    assert result.n_point_valid_a == result.n_point_valid_b == 1
    assert result.n_common_point_valid == result.n_common_interval_valid == 0
    assert np.isnan(result.bias_difference) and np.isnan(result.rmse_ratio)
    assert result.n_covered_a == result.n_covered_b == 1
    assert (result.n_both_cover, result.n_a_only_cover, result.n_b_only_cover, result.n_neither_cover) == (0, 1, 1, 0)
    assert result.operational_coverage_difference == 0
    assert result.operational_coverage_difference_mcse == 1


def test_single_replication_retains_metrics_without_fabricated_mc_uncertainty():
    design, rows = records(replications=1, cvs=(3.,))
    result = contrasts.epsilon_comparisons(rows, design)
    assert result.n_common_point_valid.eq(1).all()
    assert result.bias_difference.eq(0).all()
    assert result.rmse_ratio.eq(1).all()
    for field in ("bias_difference_mcse", "mse_difference_mcse", "rmse_difference_mcse",
                  "rmse_ratio_mcse", "rmse_ratio_log_mcse", "operational_coverage_difference_mcse"):
        assert result[field].isna().all()


def test_zero_rmse_and_zero_truth_have_explicit_undefined_metrics():
    design, rows = records(errors_zero=True)
    result = contrasts.epsilon_comparisons(rows, design)
    assert result.rmse_a_common.eq(0).all() and result.rmse_b_common.eq(0).all()
    assert result.rmse_difference.eq(0).all()
    assert result.rmse_difference_mcse.isna().all()
    assert result.rmse_ratio.isna().all() and result.rmse_ratio_mcse.isna().all()
    mu = result[result.parameter == "mu"]
    for column in ("relative_bias_a_common", "relative_bias_b_common", "relative_rmse_a_common",
                   "relative_rmse_b_common", "relative_bias_difference", "relative_bias_difference_mcse",
                   "relative_bias_difference_mc_lower", "relative_bias_difference_mc_upper"):
        assert mu[column].isna().all()


def test_valid_unrepresentable_squared_errors_block_instead_of_dropping_the_record():
    design, rows = records(cvs=(3.,))
    rows[0]["theta"] = 1e300
    before = canonical(rows)
    with pytest.raises(ValueError, match="(?i)(overflow|nonfinite)"):
        contrasts.epsilon_comparisons(rows, design)
    assert canonical(rows) == before
    assert len(rows) == 27


def test_list_frame_and_reordered_inputs_produce_same_results_without_mutation():
    design, rows = records()
    before = deepcopy(rows)
    expected = contrasts.epsilon_comparisons(rows, design)
    frame = pd.DataFrame(list(reversed(rows)))
    frame_before = frame.copy(deep=True)
    actual = contrasts.epsilon_comparisons(frame, design)
    pd.testing.assert_frame_equal(actual, expected)
    pd.testing.assert_frame_equal(frame, frame_before)
    assert canonical(rows) == canonical(before)
