"""Independent E4 reporting references; no schedules, observations, or fits.

Error/eligibility vectors and path diagnostics are literal fixtures, not study
output. Real storage, provenance verification, bindings, and writers use the
execution harness's mocked Git/provenance capture. Numeric references below
are hand calculated independently of the new reporter.
"""
from copy import deepcopy
import csv
import gzip
import hashlib
import json
import math
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from ou_irregular.experiments import execution as engine
from ou_irregular.experiments.design import ESTIMATORS
from ou_irregular.experiments.summaries import ExperimentSummarizer
from ou_irregular.week8 import execution, reporting, source
from ou_irregular.week8.design import validation_design
from .test_execution import (canonical, first_partition, harness, rewrite_partition,
                             snapshot, state, write_json)

PARAMETERS = ("theta", "mu", "sigma")
Z = 1.959963984540054
# Published source schema plus independently declared inference/raw fields.
INFERENCE_FIELDS = {
    "estimator", "point_valid", "optimizer_success", "point_status", "point_reason", "nll",
    "point_diagnostics_json", "interval_valid", "interval_status", "interval_reason",
    "interval_diagnostics_json", "covariance_type", "hessian_condition_number",
    "hessian_relative_difference", "hessian_whitened_difference", "standardized_score",
    *PARAMETERS, *(f"{p}_{s}" for p in PARAMETERS for s in ("lower", "upper", "se")),
    *(f"gradient_{i}" for i in range(3)),
    *(f"cov_{a}_{b}" for i, a in enumerate(PARAMETERS) for b in PARAMETERS[i:]),
    *(f"working_{k}_{i}_{j}" for k in ("cov", "hessian") for i in range(3) for j in range(i, 3)),
}
EXPECTED_RAW_FIELDS = (set(validation_design().cells[0].to_dict()) | set(source.SOURCE_FIELDS)
                       | INFERENCE_FIELDS | {"rep", "source_valid", "source_error", "pipeline_error"})
EXPECTED_TABLE_COUNTS = {
    "parameter_summary.csv": 63, "paired_comparisons.csv": 63,
    "regular_comparisons.csv": 54, "matched_iid_comparisons.csv": 27,
    "spacing_summary.csv": 7, "source_summary.csv": 7,
    "numerical_diagnostics.csv": 21,
}
EXPECTED_TABLES = set(EXPECTED_TABLE_COUNTS) | {"outcome_counts.csv", "failed_fits.csv.gz"}


@pytest.fixture(autouse=True)
def no_scientific_draws(monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("reporting tests must not generate sources or fit models")
    for name in ("generate_gaps", "simulate_ou"):
        monkeypatch.setattr(source, name, forbidden)
    for name in ("sample", "gaps"):
        monkeypatch.setattr(source.E4Source, name, forbidden)
    monkeypatch.setattr(np.random, "default_rng", forbidden)
    monkeypatch.setattr(engine.OuSampleSource, "sample", forbidden)
    monkeypatch.setattr(engine, "registered_wald", lambda name: SimpleNamespace(
        estimator=SimpleNamespace(fit=forbidden), infer=forbidden))


def digest_label(label):
    return hashlib.sha256(label.encode()).hexdigest()


def literal_source(cell, rep, seed=21261011):
    """Scalar fixtures only: no gap or observation array is constructed."""
    mean, regular = cell.q / cell.theta, cell.mechanism == "regular"
    floor = 1e-6 * mean if cell.cv >= 1.5 else 0.
    words = [int(cell.identity[i:i + 8], 16) for i in range(0, 64, 8)]
    address = [8, *words, rep]
    # Across-path variation exposes incorrect counting of repeated method rows.
    raw_mean = mean if regular else (1., 3., 5., 7.)[rep]
    raw_second = mean * mean if regular else (2., 10., 26., 50.)[rep]
    raw_adjacent = mean * mean if regular else (1., 4., 9., 16.)[rep]
    increment = (rep + 1) * 1e-10 if floor else 0.
    clips, counts = (rep + 1 if floor else 0), (0 if regular else 999)
    # These pure analytical calls populate constant targets, never the tested
    # empirical means/MCSEs and never a stochastic source.
    raw_theory, transformed_theory = source.moments(cell), source.moments(cell, transformed=True)
    result = {
        "family": {"mmpp": "event_stationary_mmpp", "iid": "iid_full_palm_marginal",
                   "regular": "equidistant"}[cell.mechanism],
        "target_mean_gap": mean, "requested_cv": cell.cv,
        "floor_applied": bool(floor), "floor_value": floor,
        "expected_mean_gap_after_floor": transformed_theory["mean"],
        "expected_clip_fraction": source.cdf(cell, floor) if floor else 0.,
        "clip_count": clips, "clip_fraction": clips / 999,
        "raw_below_floor_count": clips, "transformed_below_floor_count": 0,
        "transformed_at_floor_count": clips,
        "realized_mean_gap": raw_mean + increment, "realized_cv": 0. if regular else .2,
        "realized_span": 999 * (raw_mean + increment), "min_gap": mean if regular else .01,
        "raw_mean_gap": raw_mean, "raw_cv": 0. if regular else .2,
        "raw_second_moment": raw_second, "raw_adjacent_product": raw_adjacent,
        "raw_min_gap": mean if regular else .01, "raw_span": 999 * raw_mean,
        "transformed_mean_gap": raw_mean + increment, "transformed_cv": 0. if regular else .2,
        "transformed_second_moment": raw_second + increment,
        "transformed_adjacent_product": raw_adjacent + increment,
        "transformed_min_gap": mean if regular else .01,
        "transformed_span": 999 * (raw_mean + increment),
        "source_initial_phase": -1 if regular else rep % 2,
        "source_switch_count": 0 if regular else rep + 1,
        "source_max_switches_per_gap": 0 if regular else 1,
        "source_phase_reset_count": 999 if cell.mechanism == "iid" else 1 if cell.mechanism == "mmpp" else 0,
        "source_calm_start_count": counts - rep if counts else 0,
        "source_burst_start_count": rep if counts else 0,
        "source_calm_end_count": counts - rep if counts else 0,
        "source_burst_end_count": rep if counts else 0,
        "expected_raw_mean_gap": raw_theory["mean"],
        "expected_raw_second_moment": raw_theory["second_moment"],
        "expected_raw_cv": raw_theory["cv"],
        "expected_raw_adjacent_product": source.adjacent_product(cell),
        "expected_transformed_mean_gap": transformed_theory["mean"],
        "expected_transformed_second_moment": transformed_theory["second_moment"],
        "expected_transformed_cv": transformed_theory["cv"],
        "expected_transformed_adjacent_product": source.adjacent_product(cell, transformed=True),
        "source_seed": seed, "source_spawn_key": canonical(address),
        "gap_spawn_key": canonical(address + [0]), "path_spawn_key": canonical(address + [1]),
    }
    for kind in ("raw", "transformed"):
        for label, fraction in (("0p1", 0.), ("0p5", .1), ("1", .2), ("2", .3), ("5", .9)):
            result[f"{kind}_cdf_t_over_m_{label}"] = fraction
    for name in ("times_sha256", "observations_sha256", "raw_gaps_sha256",
                 "transformed_gaps_sha256", "phase_trace_sha256"):
        result[name] = digest_label(f"literal:{name}:{cell.identity}:{rep}")
    assert set(result) == set(source.SOURCE_FIELDS)
    return result


def literal_records(cell, rep, seed=21261011, *, all_valid=False):
    common = {**cell.to_dict(), **literal_source(cell, rep, seed), "rep": rep,
              "source_valid": True, "source_error": "", "pipeline_error": ""}
    result = []
    for method in ESTIMATORS:
        if cell.mechanism == "mmpp" and method == "exact":
            errors, eligible, intervals, hits = (1., 2., 3., 999.), {0, 1, 2}, {0, 1}, {0, 1}
        elif cell.mechanism == "iid" and method == "exact":
            errors, eligible, intervals, hits = (2., 999., 999., 4.), {0, 3}, {0, 3}, {0}
        elif cell.mechanism == "regular" and method == "exact":
            errors, eligible, intervals, hits = (0., 2., 4., 6.), {0, 1, 2, 3}, {0, 1, 2, 3}, {0, 1, 2}
        elif cell.mechanism == "mmpp" and method == "pfml":
            errors, eligible, intervals, hits = (2., 4., 999., 6.), {0, 1, 3}, {0, 3}, {0}
        else:
            errors, eligible, intervals, hits = (.5, 1., 1.5, 2.), {0, 1, 2, 3}, {0, 1, 2, 3}, {0, 1, 2, 3}
        if all_valid:
            errors, eligible, intervals, hits = (.5, 1., 1.5, 2.), {0, 1, 2, 3}, {0, 1, 2, 3}, {0, 1, 2, 3}
        point_valid, interval_valid = rep in eligible, rep in intervals
        row = {**common, "estimator": method, "point_valid": point_valid,
               "optimizer_success": True, "point_status": "valid_interior" if point_valid else "boundary",
               "point_reason": "literal fixture", "nll": 1., "point_diagnostics_json": '{"fixture":"literal"}',
               "interval_valid": interval_valid,
               "interval_status": "valid" if interval_valid else ("invalid_point_fit" if not point_valid else "unstable_hessian"),
               "interval_reason": "literal fixture", "interval_diagnostics_json": '{"fixture":"literal"}',
               "covariance_type": "inverse_observed_hessian_model_based",
               "hessian_condition_number": 10., "hessian_relative_difference": .001,
               "hessian_whitened_difference": .001, "standardized_score": .001}
        for parameter in PARAMETERS:
            truth = row["true_" + parameter]
            row.update({parameter: truth + errors[rep], parameter + "_se": .1,
                        parameter + "_lower": truth - .25 if rep in hits else truth + .25,
                        parameter + "_upper": truth + .25 if rep in hits else truth + .5})
        for i in range(3):
            row[f"gradient_{i}"] = 0.
            for j in range(i, 3):
                row[f"cov_{PARAMETERS[i]}_{PARAMETERS[j]}"] = .01 if i == j else 0.
                row[f"working_cov_{i}_{j}"] = .01 if i == j else 0.
                row[f"working_hessian_{i}_{j}"] = 100. if i == j else 0.
        assert set(row) == EXPECTED_RAW_FIELDS
        result.append(row)
    return result


@pytest.fixture
def records():
    return pd.DataFrame([row for cell in validation_design().cells for rep in range(4)
                         for row in literal_records(cell, rep)])


def parameter_summary(records):
    design = validation_design()
    result = ExperimentSummarizer().summarize(records, design)
    cells = {cell.identity: cell.to_dict() for cell in design.cells}
    for field in design.cells[0].to_dict():
        if field != "cell_id":
            result[field] = result.cell_id.map(lambda identity: cells[identity][field])
    return result.loc[:, sorted(result.columns)]


def select(table, **criteria):
    mask = pd.Series(True, index=table.index)
    for name, value in criteria.items():
        mask &= table[name].eq(value)
    rows = table.loc[mask]
    assert len(rows) == 1, (criteria, len(rows))
    return rows.iloc[0]


def cell_for(mechanism, ratio=1):
    return next(cell for cell in validation_design().cells
                if cell.mechanism == mechanism and (mechanism == "regular" or cell.ratio == ratio))


@pytest.fixture
def raw_run(harness, monkeypatch):
    counter = 0
    def build(*, all_valid=False, mutator=None):
        nonlocal counter
        output = harness.root / f"report-{counter}"
        counter += 1
        def saved_records(cell, rep, source_object, services):
            assert tuple(services) == ESTIMATORS
            rows = literal_records(cell, rep, source_object.seed, all_valid=all_valid)
            if mutator is not None:
                mutator(rows, cell, rep)
            return rows
        monkeypatch.setattr(engine, "source_records", saved_records)
        design = validation_design()
        harness.run(design=design, output=output, extra_files=reporting.FILES, partition_size=2)
        return design, output
    return build


def test_declared_raw_failure_schema_and_planned_counts_are_complete():
    assert set(reporting.RAW_FIELDS) == EXPECTED_RAW_FIELDS
    assert tuple(reporting.FAILURE_COLUMNS) == tuple(sorted(EXPECTED_RAW_FIELDS))
    assert set(reporting.TABLE_COLUMNS) == EXPECTED_TABLES
    assert all(tuple(columns) == tuple(sorted(columns)) for columns in reporting.TABLE_COLUMNS.values())


@pytest.mark.parametrize("kind,rows,reference_mechanism", [("regular", 54, "regular"), ("matched_iid", 27, "iid")])
def test_independent_planned_grids_and_explicit_direction(records, kind, rows, reference_mechanism):
    comparisons = reporting.independent_comparisons(parameter_summary(records), validation_design(), kind)
    assert len(comparisons) == rows
    active = [c for c in validation_design().cells if c.mechanism != "regular"
              and (kind == "regular" or c.mechanism == "mmpp")]
    expected = {(c.identity, method, p) for c in active for method in ESTIMATORS for p in PARAMETERS}
    assert set(comparisons[["cell_id", "estimator", "parameter"]].itertuples(index=False, name=None)) == expected
    assert not comparisons.duplicated(["cell_id", "estimator", "parameter"]).any()
    assert set(comparisons.comparison_type) == {"independent_cells"}
    assert set(comparisons.point_subset) == {"each_cell_own_valid_point_fits"}
    assert set(comparisons.difference_direction) == {"cell_minus_reference"}
    assert set(comparisons.ratio_direction) == {"cell_over_reference"}
    assert set(comparisons.reference_mechanism) == {reference_mechanism}
    assert not any("common" in name for name in comparisons.columns)
    for row in comparisons.itertuples(index=False):
        reference = cell_for(reference_mechanism, row.ratio)
        assert row.reference_cell_id == reference.identity
        assert row.reference_ratio == reference.ratio
        assert row.cell_id != row.reference_cell_id


@pytest.mark.parametrize("parameter", PARAMETERS)
def test_independent_own_eligibility_bias_coverage_and_positive_ratio_intervals(records, parameter):
    design, summary, a = validation_design(), parameter_summary(records), cell_for("mmpp")
    matched = select(reporting.independent_comparisons(summary, design, "matched_iid"),
                     cell_id=a.identity, estimator="exact", parameter=parameter)
    regular = select(reporting.independent_comparisons(summary, design, "regular"),
                     cell_id=a.identity, estimator="exact", parameter=parameter)
    # A: eligible errors[1,2,3]; B:[2,4] at rep0,3. Cross-cell common labels
    # would leave only rep0 and produce the wrong quantities below.
    assert (matched.cell_n_planned, matched.reference_n_planned) == (4, 4)
    assert (matched.cell_n_point_valid, matched.reference_n_point_valid) == (3, 2)
    assert (matched.cell_n_interval_valid, matched.reference_n_interval_valid) == (2, 2)
    assert (matched.cell_n_covered, matched.reference_n_covered) == (2, 1)
    assert (matched.cell_bias, matched.reference_bias) == (2., 3.)
    assert matched.cell_rmse == pytest.approx(math.sqrt(14 / 3))
    assert matched.reference_rmse == pytest.approx(math.sqrt(10))
    for row, bias_se, ratio, log_se, covdiff, condiff, condse in (
        (matched, 2 / math.sqrt(3), math.sqrt(7 / 15), math.sqrt(61) / 20, .25, .5, 1 / (2 * math.sqrt(2))),
        (regular, math.sqrt(2), 1 / math.sqrt(3), math.sqrt(7 / 48), -.25, .25, math.sqrt(3) / 8),
    ):
        assert row.excess_bias == -1
        assert row.excess_bias_mcse == pytest.approx(bias_se)
        assert row.excess_bias_mc_lower == pytest.approx(-1 - Z * bias_se)
        assert row.excess_bias_mc_upper == pytest.approx(-1 + Z * bias_se)
        assert row.rmse_ratio == pytest.approx(ratio)
        assert row.rmse_ratio_log_mcse == pytest.approx(log_se)
        assert row.rmse_ratio_mcse == pytest.approx(ratio * log_se)
        assert row.rmse_ratio_mc_lower == pytest.approx(ratio * math.exp(-Z * log_se))
        assert row.rmse_ratio_mc_upper == pytest.approx(ratio * math.exp(Z * log_se))
        assert 0 < row.rmse_ratio_mc_lower < row.rmse_ratio_mc_upper
        assert row.operational_coverage_difference == covdiff
        assert row.operational_coverage_difference_mcse == pytest.approx(math.sqrt(7) / 8)
        assert row.conditional_coverage_difference == condiff
        assert row.conditional_coverage_difference_mcse == pytest.approx(condse)
        if parameter == "mu":
            assert all(pd.isna(row["relative_excess_bias" + s]) for s in ("", "_mcse", "_mc_lower", "_mc_upper"))
        else:
            scale = 1 if parameter == "theta" else 2
            assert row.relative_excess_bias == -scale
            assert row.relative_excess_bias_mcse == pytest.approx(scale * bias_se)
    assert (regular.cell_n_point_valid, regular.reference_n_point_valid) == (3, 4)
    assert (regular.cell_n_interval_valid, regular.reference_n_interval_valid) == (2, 4)


def test_cross_cell_equal_rep_labels_never_create_pairing(records):
    design, original, altered = validation_design(), parameter_summary(records), records.copy(deep=True)
    iid = altered.mechanism.eq("iid")
    # Permute every IID source triplet together. Only each cell's own sample
    # reaches this API; no cross-cell pairing is induced by any rep label.
    altered.loc[iid, "rep"] = altered.loc[iid, "rep"].map({0: 1, 1: 0, 2: 3, 3: 2})
    changed = parameter_summary(altered)
    pd.testing.assert_frame_equal(reporting.independent_comparisons(original, design, "matched_iid"),
                                  reporting.independent_comparisons(changed, design, "matched_iid"))


@pytest.mark.parametrize("side", ["numerator", "denominator"])
def test_zero_rmse_preserves_undefined_delta_uncertainty(records, side):
    # Build a coherent zero-error sample, rather than altering a few cells in
    # an otherwise inconsistent summary row.
    changed = records.copy(deep=True)
    target = cell_for("mmpp" if side == "numerator" else "iid").identity
    mask = changed.cell_id.eq(target) & changed.estimator.eq("exact")
    changed.loc[mask, "theta"] = 1.
    row = select(reporting.independent_comparisons(parameter_summary(changed), validation_design(), "matched_iid"),
                 cell_id=cell_for("mmpp").identity, estimator="exact", parameter="theta")
    if side == "numerator":
        assert row.rmse_ratio == 0
    else:
        assert pd.isna(row.rmse_ratio)
    assert all(pd.isna(row[n]) for n in (
        "rmse_ratio_mcse", "rmse_ratio_log_mcse", "rmse_ratio_mc_lower", "rmse_ratio_mc_upper"))


@pytest.mark.parametrize("defect", ["missing_row", "duplicate_key", "extra_row", "wrong_cell_metadata", "wrong_truth", "count_overflow", "missing_column", "extra_column"])
def test_independent_summary_refuses_incomplete_or_misidentified_inputs(records, defect):
    summary = parameter_summary(records)
    if defect == "missing_row":
        summary = summary.iloc[:-1]
    elif defect == "duplicate_key":
        summary.iloc[-1] = summary.iloc[0]
    elif defect == "extra_row":
        summary = pd.concat([summary, summary.iloc[[0]]], ignore_index=True)
    elif defect == "wrong_cell_metadata":
        summary.loc[0, "mechanism"] = "unregistered"
    elif defect == "wrong_truth":
        summary.loc[0, "truth"] = 123.
    elif defect == "missing_column":
        summary = summary.drop(columns="relative_rmse")
    elif defect == "extra_column":
        summary["unregistered"] = 0
    else:
        summary.loc[0, "n_point_valid"] = 5
    with pytest.raises(ValueError):
        reporting.independent_comparisons(summary, validation_design(), "matched_iid")


def test_independent_comparison_kind_is_explicit(records):
    with pytest.raises(ValueError):
        reporting.independent_comparisons(parameter_summary(records), validation_design(), "paired")


@pytest.mark.parametrize("field,parameter,value", [
    ("true_theta", "theta", True), ("true_mu", "mu", False),
    ("truth", "theta", True), ("truth", "mu", False),
])
def test_independent_summary_rejects_boolean_numeric_metadata_or_truth(records, field, parameter, value):
    summary = parameter_summary(records)
    summary[field] = summary[field].astype(object)
    index = summary.index[summary.parameter.eq(parameter)][0]
    summary.loc[index, field] = value
    assert type(summary.loc[index, field]) is bool
    with pytest.raises(ValueError):
        reporting.independent_comparisons(summary, validation_design(), "matched_iid")


@pytest.mark.parametrize("mechanism", ["mmpp", "iid", "regular"])
def test_source_no_floor_cannot_report_coherent_positive_clipping(records, mechanism):
    changed = records.copy(deep=True)
    mask = changed.cell_id.eq(cell_for(mechanism).identity) & changed.rep.eq(0)
    assert changed.loc[mask, "floor_value"].eq(0).all()
    assert not changed.loc[mask, "floor_applied"].any()
    # All triplets and arithmetic agree; only the physical impossibility of
    # clipping a strictly positive gap at a zero floor distinguishes this.
    changed.loc[mask, ["clip_count", "raw_below_floor_count", "transformed_at_floor_count"]] = 1
    changed.loc[mask, "clip_fraction"] = 1 / 999
    with pytest.raises(ValueError):
        reporting.source_summary(changed, validation_design())



@pytest.mark.parametrize("eligible_count", [0, 1])
def test_sparse_valid_subset_keeps_counts_and_withholds_unavailable_uncertainty(records, eligible_count):
    changed = records.copy(deep=True)
    identity = cell_for("mmpp").identity
    mask = changed.cell_id.eq(identity) & changed.estimator.eq("exact") & changed.rep.ge(eligible_count)
    changed.loc[mask, "point_valid"] = False
    changed.loc[mask, "point_status"] = "boundary"
    changed.loc[mask, "interval_valid"] = False
    changed.loc[mask, "interval_status"] = "invalid_point_fit"
    row = select(reporting.independent_comparisons(parameter_summary(changed), validation_design(), "matched_iid"),
                 cell_id=identity, estimator="exact", parameter="theta")
    assert (row.cell_n_planned, row.cell_n_point_valid, row.cell_n_interval_valid) == (4, eligible_count, eligible_count)
    assert row.reference_n_point_valid == 2
    assert pd.isna(row.excess_bias_mcse) and pd.isna(row.rmse_ratio_log_mcse)
    if eligible_count:
        assert row.excess_bias == -2 and row.rmse_ratio == pytest.approx(1 / math.sqrt(10))
    else:
        assert pd.isna(row.excess_bias) and pd.isna(row.rmse_ratio)
        assert pd.isna(row.conditional_coverage_difference)
        assert row.operational_coverage_difference == -.25
        assert row.operational_coverage_difference_mcse == pytest.approx(math.sqrt(3) / 8)


def test_source_summary_uses_independent_paths_and_retains_raw_transformed_statistics(records):
    result = reporting.source_summary(records, validation_design())
    assert len(result) == 7 and not result.cell_id.duplicated().any()
    a = select(result, cell_id=cell_for("mmpp", 20).identity)
    assert (a.n_planned, a.n_unique_sources, a.n_source_valid, a.n_fit_records) == (4, 4, 4, 12)
    assert a.raw_mean_gap_n_available == 4
    assert a.raw_mean_gap_mean == 4.
    assert a.raw_mean_gap_mcse == pytest.approx(math.sqrt(5 / 3))
    assert (a.raw_mean_gap_min, a.raw_mean_gap_max) == (1., 7.)
    # Second-moment values [2,10,26,50] have mean22 and MCSEsqrt(112).
    assert a.raw_second_moment_n_available == 4
    assert a.raw_second_moment_mean == 22.
    assert a.raw_second_moment_mcse == pytest.approx(4 * math.sqrt(7))
    assert a.transformed_mean_gap_mean == pytest.approx(4.00000000025, abs=1e-14, rel=0)
    assert a.transformed_second_moment_mean == pytest.approx(22.00000000025, abs=1e-13, rel=0)
    assert a.raw_mean_gap_mean != a.transformed_mean_gap_mean
    assert a.raw_second_moment_mean != a.transformed_second_moment_mean
    assert a.expected_raw_mean_gap == pytest.approx(.5, rel=0, abs=2e-15)
    assert a.clip_count_total == 10
    assert a.clip_count_n_available == 4 and a.clip_count_mean == 2.5
    assert a.clip_count_mcse == pytest.approx(math.sqrt(5 / 12))
    assert (a.initial_phase_calm_count, a.initial_phase_burst_count, a.initial_phase_regular_count) == (2, 2, 0)
    assert a.initial_phase_burst_fraction == .5
    assert a.initial_phase_burst_mcse == pytest.approx(1 / math.sqrt(12))
    regular = select(result, cell_id=cell_for("regular").identity)
    assert regular.raw_mean_gap_mean == .5 and regular.raw_mean_gap_mcse == 0
    assert regular.source_switch_count_total == regular.clip_count_total == 0
    # Sorting rows or interleaving estimators must not determine the sample unit.
    reversed_result = reporting.source_summary(records.iloc[::-1], validation_design())
    pd.testing.assert_frame_equal(result.sort_values("cell_id").reset_index(drop=True),
                                  reversed_result.sort_values("cell_id").reset_index(drop=True))


@pytest.mark.parametrize("field", [
    "phase_trace_sha256", "raw_second_moment", "transformed_adjacent_product",
    "expected_raw_cv", "path_spawn_key", "source_calm_end_count",
    "raw_cdf_t_over_m_0p1", "source_initial_phase",
])
def test_source_triplets_require_new_and_old_diagnostics_to_agree(records, field):
    changed = records.copy(deep=True)
    old = changed.loc[1, field]
    changed.loc[1, field] = (digest_label("other fixture") if field.endswith("_sha256")
                           else "[]" if field.endswith("_key") else 1 if field == "source_initial_phase"
                           else old + .01 if "cdf" in field else old + 1)
    with pytest.raises(ValueError):
        reporting.source_summary(changed, validation_design())


@pytest.mark.parametrize("defect", [
    "missing_raw_field", "unexpected_raw_field", "missing_triplet_member", "duplicate_key",
    "wrong_seed", "wrong_spawn", "nonfinite_metric", "negative_metric", "fraction_above_one",
    "noninteger_count", "wrong_phase_total", "wrong_reset", "changing_expected_target",
])
def test_source_summary_rejects_invalid_schema_keys_addresses_or_diagnostics(records, defect):
    changed = records.copy(deep=True)
    triplet = changed.cell_id.eq(cell_for("mmpp").identity) & changed.rep.eq(0)
    if defect == "missing_raw_field":
        changed = changed.drop(columns="phase_trace_sha256")
    elif defect == "unexpected_raw_field":
        changed["unregistered_source"] = 0
    elif defect == "missing_triplet_member":
        changed = changed.iloc[1:]
    elif defect == "duplicate_key":
        changed.iloc[-1] = changed.iloc[0]
    else:
        field, value = {
            "wrong_seed": ("source_seed", 20261011),
            "wrong_spawn": ("source_spawn_key", "[]"),
            "nonfinite_metric": ("raw_second_moment", np.inf),
            "negative_metric": ("raw_second_moment", -1.),
            "fraction_above_one": ("raw_cdf_t_over_m_1", 1.01),
            "noninteger_count": ("source_switch_count", 1.5),
            "wrong_phase_total": ("source_calm_start_count", 998),
            "wrong_reset": ("source_phase_reset_count", 998),
            "changing_expected_target": ("expected_raw_mean_gap", .51),
        }[defect]
        if defect == "noninteger_count":
            changed[field] = changed[field].astype(object)
        changed.loc[triplet, field] = value
    with pytest.raises(ValueError):
        reporting.source_summary(changed, validation_design())


def test_written_tables_retain_paired_common_subset_and_operational_discordance(raw_run):
    design, output = raw_run()
    before = snapshot(output)
    report = reporting.write_tables(design, output)
    assert set(report) == {"tables", "binding"}
    assert set(report["tables"]) == EXPECTED_TABLES
    for name, count in EXPECTED_TABLE_COUNTS.items():
        assert report["tables"][name]["rows"] == count
    assert report["tables"]["failed_fits.csv.gz"]["rows"] == 18
    for name, value in before.items():
        assert snapshot(output)[name] == value
    paired = pd.read_csv(output / "paired_comparisons.csv")
    for parameter in PARAMETERS:
        row = select(paired, cell_id=cell_for("mmpp").identity,
                     estimator_a="exact", estimator_b="pfml", parameter=parameter)
        assert row.direction == "a_minus_b_and_a_over_b"
        assert row.point_subset == "intersection_of_valid_point_fits"
        assert (row.n_planned, row.n_point_valid_a, row.n_point_valid_b,
                row.n_common_point_valid, row.n_common_interval_valid) == (4, 3, 3, 2, 1)
        assert (row.bias_a_common, row.bias_b_common) == (1.5, 3.)
        assert row.bias_difference == -1.5 and row.bias_difference_mcse == .5
        assert row.mse_difference == -7.5 and row.mse_difference_mcse == 4.5
        assert row.squared_error_covariance_common == 18
        assert row.rmse_ratio == .5 and row.rmse_ratio_log_mcse == row.rmse_ratio_mcse == 0
        assert row.rmse_ratio_mc_lower == row.rmse_ratio_mc_upper == .5
        assert row.rmse_difference == pytest.approx(-math.sqrt(2.5))
        assert row.rmse_difference_mcse == pytest.approx(3 / (4 * math.sqrt(2.5)))
        assert (row.n_covered_a, row.n_covered_b, row.n_both_cover,
                row.n_a_only_cover, row.n_b_only_cover, row.n_neither_cover) == (2, 1, 1, 1, 0, 2)
        assert row.operational_coverage_difference == .25
        assert row.operational_coverage_difference_mcse == .25
        if parameter == "mu":
            assert pd.isna(row.relative_bias_difference)
    summary = pd.read_csv(output / "parameter_summary.csv")
    assert summary.n_point_valid.sum() == 216  # three parameter rows per 72 valid fits
    assert summary.n_interval_valid.sum() == 198
    outcomes = pd.read_csv(output / "outcome_counts.csv")
    assert outcomes["count"].sum() == 84
    assert outcomes.groupby(["cell_id", "estimator"])["count"].sum().eq(4).all()
    with gzip.open(output / "failed_fits.csv.gz", "rt", newline="") as stream:
        reader = csv.DictReader(stream)
        assert reader.fieldnames == sorted(EXPECTED_RAW_FIELDS)
        failures = list(reader)
    assert len(failures) == 18
    assert sum(row["point_valid"] == "True" for row in failures) == 6
    published = snapshot(output)
    assert reporting.verify_tables(design, output) == report
    assert reporting.write_tables(design, output) == report
    assert snapshot(output) == published


def test_empty_failed_export_retains_the_exact_complete_sorted_header(raw_run):
    design, output = raw_run(all_valid=True)
    result = reporting.write_tables(design, output)
    assert result["tables"]["failed_fits.csv.gz"]["rows"] == 0
    payload = gzip.decompress((output / "failed_fits.csv.gz").read_bytes()).decode()
    assert payload == ",".join(sorted(EXPECTED_RAW_FIELDS)) + "\n"
    with gzip.open(output / "failed_fits.csv.gz", "rt", newline="") as stream:
        assert list(csv.reader(stream)) == [sorted(EXPECTED_RAW_FIELDS)]


def test_pure_preflight_and_binding_cover_exact_raw_evidence_but_not_lifecycle(raw_run):
    design, output = raw_run()
    before = snapshot(output)
    prepared = reporting.preflight_reporting(design, output)
    assert snapshot(output) == before
    assert prepared.manifest is None and prepared.table_binding is None
    assert set(prepared.binding) == {"run_identity", "design_sha256", "resolved_design_sha256", "source_policy", "cells"}
    assert set(prepared.binding["cells"]) == {cell.identity for cell in design.cells}
    for cid, record in prepared.binding["cells"].items():
        directory = output / "cells" / cid
        assert record["manifest_sha256"] == hashlib.sha256((directory / "cell_manifest.json").read_bytes()).hexdigest()
        assert sum(part["rows"] for part in record["partitions"]) == 12
        for part in record["partitions"]:
            assert part["sha256"] == hashlib.sha256((directory / part["filename"]).read_bytes()).hexdigest()
    changed = deepcopy(prepared.raw.metadata)
    changed.update(status="completed", completed_utc="LIFECYCLE_ONLY", audit_sha256="a" * 64)
    assert reporting.raw_input_binding(output, design, changed) == prepared.binding
    report = reporting.write_tables(design, output)
    assert report["binding"]["raw_inputs"] == prepared.binding
    assert report["binding"]["raw_input_sha256"] == hashlib.sha256(canonical(prepared.binding).encode()).hexdigest()
    assert report["binding"]["reporting_contract"] == "week8_saved_raw_tables_v1"
    assert report["binding"]["design_sha256"] == design.digest


@pytest.mark.parametrize("artifact", ["parameter_summary.csv", "source_summary.csv", "failed_fits.csv.gz", "table_binding.json"])
def test_declared_missing_artifacts_are_not_recreated(raw_run, artifact):
    design, output = raw_run()
    reporting.write_tables(design, output)
    (output / artifact).unlink()
    before = snapshot(output)
    with pytest.raises((ValueError, FileNotFoundError)):
        reporting.write_tables(design, output)
    assert snapshot(output) == before


def test_closed_state_allows_pure_verification_but_refuses_republication(raw_run):
    design, output = raw_run()
    report = reporting.write_tables(design, output)
    metadata = state(output)
    metadata.update(status="completed", completed_utc="LITERAL_COMPLETION")
    write_json(output / "run_metadata.json", metadata)
    before = snapshot(output)
    assert reporting.verify_tables(design, output) == report
    assert reporting.raw_input_binding(output, design) == report["binding"]["raw_inputs"]
    with pytest.raises(ValueError):
        reporting.write_tables(design, output)
    assert snapshot(output) == before


@pytest.mark.parametrize("mutation", ["source_triplet", "unknown_raw_field", "invalid_source_metric"])
def test_late_cell_source_failure_cannot_partially_publish_earlier_cells(raw_run, mutation):
    def mutate(rows, cell, rep):
        # Last planned cell: a streaming writer must not publish earlier cells
        # before discovering this late evidence defect.
        if cell.mechanism == "regular" and rep == 3:
            if mutation == "source_triplet":
                rows[1]["phase_trace_sha256"] = digest_label("different literal source")
            elif mutation == "unknown_raw_field":
                rows[1]["undeclared_field"] = 1
            else:
                for row in rows:
                    row["raw_second_moment"] = -1.
    design, output = raw_run(mutator=mutate)
    before = snapshot(output)
    with pytest.raises(ValueError):
        reporting.write_tables(design, output)
    assert snapshot(output) == before
    assert not any((output / name).exists() for name in EXPECTED_TABLES)


def test_saved_raw_change_with_updated_local_hashes_cannot_rebind_existing_tables(raw_run):
    design, output = raw_run()
    reporting.write_tables(design, output)
    directory, manifest, part, path = first_partition(output)
    def change(rows):
        for row in rows:
            if row["rep"] == 0:
                row["raw_adjacent_product"] += .1
    part["sha256"] = rewrite_partition(path, change)
    write_json(directory / "cell_manifest.json", manifest)
    before = snapshot(output)
    for operation in (reporting.verify_tables, reporting.preflight_reporting, reporting.write_tables):
        with pytest.raises(ValueError):
            operation(design, output)
        assert snapshot(output) == before


def refresh_table_manifest(output, filename, rows=None):
    manifest = json.loads((output / "table_manifest.json").read_bytes())
    manifest[filename]["sha256"] = hashlib.sha256((output / filename).read_bytes()).hexdigest()
    if rows is not None:
        manifest[filename]["rows"] = rows
    write_json(output / "table_manifest.json", manifest)


def read_text_csv(path):
    """Preserve untouched tokens, avoiding incidental float-parser changes."""
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        return list(reader.fieldnames), list(reader)


def write_text_csv(path, columns, rows):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


@pytest.mark.parametrize("defect", ["duplicate_key", "unknown_column", "missing_column", "wrong_cell_metadata", "missing_row", "wrong_reference"])
def test_rehashed_table_corruption_is_rejected_by_schema_and_planned_key_checks(raw_run, defect):
    design, output = raw_run()
    reporting.write_tables(design, output)
    filename = "matched_iid_comparisons.csv" if defect == "wrong_reference" else "parameter_summary.csv"
    path = output / filename
    columns, rows = read_text_csv(path)
    if defect == "duplicate_key":
        rows[-1] = dict(rows[0])
    elif defect == "unknown_column":
        columns.append("unregistered")
        for row in rows:
            row["unregistered"] = "0"
    elif defect == "missing_column":
        columns.remove("mechanism")
        for row in rows:
            del row["mechanism"]
    elif defect == "wrong_cell_metadata":
        rows[0]["ratio"] = "123"
    elif defect == "missing_row":
        rows = rows[:-1]
    else:
        rows[0]["reference_cell_id"] = cell_for("regular").identity
    write_text_csv(path, columns, rows)
    refresh_table_manifest(output, filename, rows=len(rows))
    before = snapshot(output)
    with pytest.raises(ValueError):
        reporting.verify_tables(design, output)
    assert snapshot(output) == before


def test_empty_failed_header_cannot_drop_a_source_field_even_with_rehashed_manifest(raw_run):
    design, output = raw_run(all_valid=True)
    reporting.write_tables(design, output)
    path = output / "failed_fits.csv.gz"
    columns = sorted(EXPECTED_RAW_FIELDS - {"phase_trace_sha256"})
    path.write_bytes(gzip.compress((",".join(columns) + "\n").encode(), mtime=0))
    refresh_table_manifest(output, path.name)
    before = snapshot(output)
    with pytest.raises(ValueError):
        reporting.verify_tables(design, output)
    assert snapshot(output) == before


def test_identical_partial_publication_is_adopted_without_changing_saved_bytes(raw_run):
    design, output = raw_run()
    original = reporting.write_tables(design, output)
    expected = {name: (output / name).read_bytes() for name in EXPECTED_TABLES}
    for name in ("table_manifest.json", "table_binding.json", "regular_comparisons.csv", "outcome_counts.csv"):
        (output / name).unlink()
    before = snapshot(output)
    assert reporting.write_tables(design, output) == original
    assert all((output / name).read_bytes() == payload for name, payload in expected.items())
    assert all(snapshot(output)[name] == payload for name, payload in before.items())


def test_mismatching_late_partial_table_blocks_all_other_publication(raw_run):
    design, output = raw_run()
    reporting.write_tables(design, output)
    # Keep a structurally valid partial source table. Its changed empirical
    # mean passes shape/identity checks, so regeneration must compare bytes
    # against every candidate before publishing any other missing artifact.
    for name in EXPECTED_TABLES - {"source_summary.csv"} | {"table_manifest.json", "table_binding.json"}:
        (output / name).unlink()
    path = output / "source_summary.csv"
    columns, rows = read_text_csv(path)
    rows[0]["raw_mean_gap_mean"] = format(float(rows[0]["raw_mean_gap_mean"]) + .1, ".17g")
    write_text_csv(path, columns, rows)
    before = snapshot(output)
    assert reporting.preflight_reporting(design, output).manifest is None
    assert snapshot(output) == before
    with pytest.raises(FileExistsError, match="source_summary.csv"):
        reporting.write_tables(design, output)
    assert snapshot(output) == before
    assert not (output / "table_binding.json").exists() and not (output / "table_manifest.json").exists()


@pytest.mark.parametrize("defect", ["binding_hash", "binding_identity", "manifest_unknown_file", "root_incomplete"])
def test_reporting_preflight_failure_is_read_only(raw_run, defect):
    design, output = raw_run()
    reporting.write_tables(design, output)
    if defect.startswith("binding"):
        path = output / "table_binding.json"
        value = json.loads(path.read_bytes())
        value["raw_input_sha256" if defect == "binding_hash" else "run_identity"] = "f" * 64
        write_json(path, value)
    elif defect == "manifest_unknown_file":
        path = output / "table_manifest.json"
        value = json.loads(path.read_bytes())
        value["unregistered.csv"] = {"rows": 0, "sha256": "f" * 64}
        write_json(path, value)
    else:
        value = state(output)
        value["status"] = "incomplete"
        write_json(output / "run_metadata.json", value)
    before = snapshot(output)
    with pytest.raises(ValueError):
        reporting.preflight_reporting(design, output)
    assert snapshot(output) == before
