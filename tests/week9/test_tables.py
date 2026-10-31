"""Publication arithmetic and data-selection contracts, without new experiments."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from ou_irregular.week9.tables import (
    _validate_summary, build_claims, build_tables, classify_row, select_one,
)


ROOT = Path(__file__).resolve().parents[2]


def point(**updates):
    row = {"parameter": "theta", "truth": 1., "operational_wilson_lower": .8,
           "operational_wilson_upper": .9, "relative_bias_mc_lower": -.1,
           "relative_bias_mc_upper": .1}
    return {**row, **updates}


@pytest.mark.parametrize("row,wanted", [
    (point(), ("within_criteria", "none")),
    (point(operational_wilson_lower=.7, operational_wilson_upper=.8), ("unresolved", "uncertain")),
    (point(operational_wilson_lower=.7, operational_wilson_upper=np.nextafter(.8, 0.)), ("supported_breach", "coverage")),
    (point(relative_bias_mc_lower=.1, relative_bias_mc_upper=.12), ("unresolved", "uncertain")),
    (point(relative_bias_mc_lower=np.nextafter(.1, 1.), relative_bias_mc_upper=.12), ("supported_breach", "bias")),
    (point(relative_bias_mc_lower=-.12, relative_bias_mc_upper=-.1), ("unresolved", "uncertain")),
    (point(relative_bias_mc_lower=-.12, relative_bias_mc_upper=np.nextafter(-.1, -1.)), ("supported_breach", "bias")),
    (point(operational_wilson_lower=.7, operational_wilson_upper=.79, relative_bias_mc_lower=.11,
           relative_bias_mc_upper=.12), ("supported_breach", "bias_and_coverage")),
    (point(relative_bias_mc_lower=np.nan, relative_bias_mc_upper=np.nan), ("unresolved", "uncertain")),
    (point(operational_wilson_lower=.7, operational_wilson_upper=.79, relative_bias_mc_lower=np.nan,
           relative_bias_mc_upper=np.nan), ("supported_breach", "coverage")),
    (point(parameter="mu", truth=0.), ("not_applicable_zero_truth", "not_applicable")),
])
def test_classification_uses_unrounded_bounds_and_strict_breach(row, wanted):
    assert classify_row(row) == wanted


@pytest.mark.parametrize("updates", [
    {"parameter": "unknown"}, {"truth": 0.}, {"parameter": "mu"},
    {"operational_wilson_lower": .95, "operational_wilson_upper": .9},
    {"operational_wilson_lower": -.01}, {"operational_wilson_upper": np.nan},
    {"relative_bias_mc_lower": .1, "relative_bias_mc_upper": -.1},
])
def test_invalid_criterion_inputs_fail(updates):
    with pytest.raises(ValueError):
        classify_row(point(**updates))


def test_selector_refuses_missing_ambiguous_and_unknown_keys():
    frame = pd.DataFrame({"cell_id": ["a", "a"], "estimator": ["exact", "pfml"]})
    assert select_one(frame, cell_id="a", estimator="exact").estimator == "exact"
    for keys in ({"cell_id": "a"}, {"cell_id": "missing"}, {"not_a_key": "a"}):
        with pytest.raises(ValueError):
            select_one(frame, **keys)


@pytest.fixture(scope="module")
def evidence():
    from ou_irregular.week9.inputs import load_sources
    bundle = load_sources(root=ROOT)
    before = bundle.snapshot()
    tables = build_tables(bundle)
    claims = build_claims(bundle, tables)
    assert bundle.snapshot() == before
    return bundle, tables, claims


def test_complete_published_tables_and_class_counts(evidence):
    _, tables, _ = evidence
    sizes = {name: len(table) for name, table in tables.items()}
    assert sizes == {"parameter_summary.csv": 792, "T1_e2.csv": 360, "T1_e1.csv": 42,
                     "T1_counts.csv": 18, "classification_counts.csv": 5, "unresolved.csv": 6,
                     "availability.csv": 264, "calibration.csv": 792, "e4_matched.csv": 27,
                     "e3_floor_contrasts.csv": 72, "e4_sources.csv": 7, "f1_ecdf.csv": 126000}
    counts = tables["classification_counts.csv"].set_index("study")
    assert counts.loc["e1", ["within_criteria", "supported_breach", "unresolved", "not_applicable_zero_truth"]].tolist() == [19,22,1,21]
    assert counts.loc["e2", ["within_criteria", "supported_breach", "unresolved", "not_applicable_zero_truth"]].tolist() == [160,195,5,180]
    compact = tables["T1_counts.csv"]
    assert compact.n_cells.eq(20).all()
    assert compact[["within_criteria", "supported_breach", "unresolved"]].sum(axis=1).eq(compact.n_cells).all()


def test_unresolved_keys_and_independent_overlap_are_preserved(evidence):
    _, tables, _ = evidence
    columns = ["study", "n", "theta_mean_gap", "cv", "estimator", "parameter"]
    actual = set(tables["unresolved.csv"][columns].itertuples(index=False, name=None))
    assert actual == {("e1",1000,.5,.5,"pfml","sigma"), ("e2",250,.1,0.,"euler","theta"),
                      ("e2",250,.1,.5,"euler","theta"), ("e2",250,.1,2.,"euler","sigma"),
                      ("e2",250,.1,3.,"pfml","theta"), ("e2",1000,.5,.5,"pfml","sigma")}
    summary = tables["parameter_summary.csv"]
    overlap = summary[summary.study.isin(["e1","e2"]) & summary.n.eq(1000) & summary.theta_mean_gap.eq(.5)
                      & summary.cv.eq(.5) & summary.estimator.eq("pfml") & summary.parameter.eq("sigma")]
    assert len(overlap) == 2 and overlap.cell_id.nunique() == 1
    assert overlap.study.nunique() == 2 and overlap.source_sha256.nunique() == 2
    assert set(overlap.coverage_operational) == {.7945, .7865}
    assert not summary.duplicated(["study", "cell_id", "estimator", "parameter"]).any()


def test_map_preserves_nonmonotone_observed_classes(evidence):
    _, tables, _ = evidence
    table = tables["T1_e2.csv"]
    target = table[table.n.eq(250) & table.theta_mean_gap.eq(.1)]
    labels = {"within_criteria": "W", "supported_breach": "B", "unresolved": "U"}
    pfml = target[target.estimator.eq("pfml") & target.parameter.eq("theta")].sort_values("cv")
    euler = target[target.estimator.eq("euler") & target.parameter.eq("sigma")].sort_values("cv")
    assert [labels[x] for x in pfml.practical_classification] == ["B","B","B","W","U"]
    assert [labels[x] for x in euler.practical_classification] == ["B","B","B","U","W"]
    exact = table[table.n.eq(250) & table.theta_mean_gap.eq(.1) & table.estimator.eq("exact") & table.parameter.eq("theta")]
    assert exact.breach_driver.eq("bias").all()
    assert exact.coverage_operational.min() == .936 and exact.coverage_operational.max() == .953


def test_availability_does_not_double_count_invalid_points(evidence):
    _, tables, _ = evidence
    table = tables["availability.csv"]
    e2 = table[table.study.eq("e2")]
    assert e2.n_point_invalid.sum() == 85
    assert e2.n_point_valid_without_interval.sum() == 1
    assert e2.n_interval_unavailable.sum() == 86
    followup = table[table.study.eq("e3_followup")]
    assert followup.n_point_invalid.sum() == 421
    assert followup.n_point_valid_without_interval.sum() == 23
    assert followup.n_interval_unavailable.sum() == 444
    assert (table.n_point_invalid+table.n_point_valid_without_interval).eq(table.n_interval_unavailable).all()


def test_calibration_uses_interval_valid_sd_and_rms_se(evidence):
    _, tables, _ = evidence
    table = tables["calibration.csv"]
    row = select_one(table, study="e4", estimator="pfml", parameter="mu", mechanism="mmpp", ratio=20)
    assert row.sd_over_rms_se == row.empirical_sd_interval_valid/row.rms_reported_se
    assert row.sd_over_rms_se == pytest.approx(1.14486897, rel=1e-7)
    assert table.eligibility_subset.eq("same_interval_valid_fits_for_empirical_sd_and_reported_se").all()
    assert tables["parameter_summary.csv"].query("parameter == 'mu'")[["relative_bias", "relative_rmse"]].isna().all().all()


def test_f1_contains_every_saved_eligible_estimate_and_partition_binding(evidence):
    bundle, tables, _ = evidence
    table = tables["f1_ecdf.csv"]
    assert set(table.cv) == {0., .25, .5, 1., 1.5, 2., 3.}
    assert set(table.parameter) == {"theta", "mu", "sigma"}
    assert set(table.estimator) == {"exact", "pfml", "euler"}
    assert not table.duplicated(["study", "cell_id", "estimator", "parameter", "rep"]).any()
    for _, group in table.groupby(["cell_id", "estimator", "parameter"]):
        assert len(group) == 2000
        assert group.n_planned.eq(2000).all() and group.n_point_valid.eq(2000).all()
        assert group["rank"].tolist() == list(range(1,2001))
        assert group.estimate.is_monotonic_increasing
        assert group.ecdf.iloc[0] == 1/2000 and group.ecdf.iloc[-1] == 1.
    for path, sha in table[["source_path", "source_sha256"]].drop_duplicates().itertuples(index=False, name=None):
        assert hashlib.sha256((bundle.root/path).read_bytes()).hexdigest() == sha


def test_generated_claims_bind_saved_values_and_exclude_unsupported_scope(evidence):
    bundle, tables, claims = evidence
    ledger = {claim["claim_id"]: claim for claim in claims}
    assert len(ledger) == len(claims)
    assert ledger["e4_pfml_mu_r5"]["values"]["rmse_ratio"] == pytest.approx(1.0838740655)
    assert ledger["e4_pfml_mu_r20"]["values"]["rmse_ratio"] == pytest.approx(1.1578630861)
    assert ledger["beyond_first_order_theory"]["status"] == "unsupported_claim_excluded"
    assert ledger["general_dependent_gap_covariance"]["status"] == "unsupported_claim_excluded"
    assert ledger["e2_exact_not_uniform_lowest_rmse"]["status"] == "descriptive"
    assert ledger["e3_floor_sensitivity"]["values"]["source_realizations"] == 24000
    assert ledger["e3_floor_sensitivity"]["values"]["independent_coupling_groups"] == 8000
    assert ledger["e3_selected_replay_scope"]["values"]["evaluated_sources"] == 53
    assert ledger["e3_selected_replay_scope"]["values"]["verified_fit_replays"] == 159
    for claim in claims:
        assert claim["source_bindings"] and claim["qualifications"]
        for binding in claim["source_bindings"]:
            assert hashlib.sha256((bundle.root/binding["path"]).read_bytes()).hexdigest() == binding["sha256"]
    json.dumps(claims, allow_nan=False)


def test_claim_row_bindings_resolve_to_actual_original_csv_columns(evidence):
    bundle, _, claims = evidence
    cache = {}
    for claim in claims:
        for binding in claim["source_bindings"]:
            if not binding["path"].endswith(".csv"):
                continue
            path = binding["path"]
            if path not in cache:
                cache[path] = pd.read_csv(bundle.root/path, float_precision="round_trip")
            source = cache[path]
            assert set(binding["selector_columns"]+binding["columns"]) <= set(source.columns)
            for bound in binding["rows"]:
                row = select_one(source, **{key: bound[key] for key in binding["selector_columns"]})
                for column in binding["columns"]:
                    expected = bound[column]
                    if expected is None:
                        assert pd.isna(row[column])
                    else:
                        assert row[column] == expected


@pytest.mark.parametrize("mutation", ["classification", "relative_mu", "denominator", "duplicate"])
def test_summary_readback_detects_scientific_reporting_mutations(mutation):
    original = pd.read_csv(ROOT/"artifacts/week5_e1/parameter_summary.csv", float_precision="round_trip")
    original.insert(0,"study","e1")
    _validate_summary(original)
    changed = original.copy()
    if mutation == "classification":
        changed.loc[0,"practical_classification"] = "unresolved"
    elif mutation == "relative_mu":
        changed.loc[changed.parameter.eq("mu"),"relative_bias"] = 0.
    elif mutation == "denominator":
        changed.loc[0,"coverage_operational"] = .123
    else:
        changed = pd.concat([changed, changed.iloc[[0]]], ignore_index=True)
    with pytest.raises(ValueError):
        _validate_summary(changed)
