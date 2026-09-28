"""Week 4 workflow invariants and hand-counted coverage denominators."""
from copy import deepcopy

import numpy as np
import pytest
import yaml

from ou_irregular.week4.config import ValidationConfig, load_config
from ou_irregular.week4.simulation import OuSampleSource, array_digest, bootstrap_seed
from ou_irregular.week4.summaries import CoverageSummarizer, ThresholdPolicy, wilson_interval


@pytest.fixture
def config_data():
    # A tiny local fixture describes the same schema without running a
    # certification experiment or depending on its chosen sample size.
    return {
        "schema_version": 1, "phase": "week4_validation", "name": "test",
        "seed": 42, "truth": {"theta": 1., "mu": 0., "sigma": .5},
        "gap_floor": {"fraction": 1e-6, "cv_threshold": 1.5},
        "wald": {"level": .95, "relative_step": .001,
                 "hessian_relative_tolerance": .005, "max_whitened_difference": .01,
                 "max_condition": 1e10, "max_standardized_score": .01},
        "certification": {"q": .5, "cv": 0., "n": 30, "replications": 5, "estimator": "exact"},
        "bootstrap": {"seed": 43, "draws": 5, "level": .95, "min_valid_fraction": .95,
                      "sample_rep": 0, "estimators": ["exact", "pfml", "euler"],
                      "cells": [{"q": .5, "cv": 0., "n": 30}, {"q": .5, "cv": 2., "n": 30}]},
    }


def set_nested(value, keys, replacement):
    for key in keys[:-1]:
        value = value[key]
    value[keys[-1]] = replacement


@pytest.mark.parametrize("keys,value", [
    (("schema_version",), True), (("schema_version",), 2),
    (("phase",), "headline"), (("name",), " "),
    (("seed",), True), (("seed",), -1), (("seed",), 4.),
    (("truth", "theta"), 0.), (("truth", "sigma"), -1.), (("truth", "mu"), np.nan),
    (("gap_floor", "fraction"), 0.), (("gap_floor", "fraction"), 1.),
    (("gap_floor", "cv_threshold"), 1.),
    (("wald", "level"), 1.), (("wald", "relative_step"), 0.),
    (("wald", "max_condition"), .5),
    (("certification", "estimator"), "pfml"), (("certification", "cv"), 1.),
    (("certification", "n"), 2), (("certification", "replications"), False),
    (("bootstrap", "draws"), 1), (("bootstrap", "draws"), True),
    (("bootstrap", "sample_rep"), -1), (("bootstrap", "level"), 0.),
    (("bootstrap", "min_valid_fraction"), 0.),
    (("bootstrap", "estimators"), ["exact", "pfml", "pfml"]),
    (("bootstrap", "cells"), []),
    (("bootstrap", "cells", 0, "q"), -1.),
    (("bootstrap", "cells", 0, "cv"), True),
    (("bootstrap", "cells", 0, "n"), 3.5),
])
def test_config_rejects_invalid_values(config_data, keys, value):
    set_nested(config_data, keys, value)
    with pytest.raises(ValueError):
        ValidationConfig.from_dict(config_data)


def test_config_rejects_unknown_and_missing_keys(config_data):
    altered = deepcopy(config_data)
    altered["headlines"] = True
    with pytest.raises(ValueError):
        ValidationConfig.from_dict(altered)
    altered = deepcopy(config_data)
    del altered["wald"]["max_standardized_score"]
    with pytest.raises(ValueError):
        ValidationConfig.from_dict(altered)


@pytest.mark.parametrize("convert_numeric_types", [False, True])
def test_config_rejects_duplicate_scientific_cells(config_data, convert_numeric_types):
    first = {"q": .5, "cv": 0., "n": 30}
    duplicate = {"q": .5, "cv": 0 if convert_numeric_types else 0., "n": 30}
    config_data["bootstrap"]["cells"] = [first, duplicate]
    with pytest.raises(ValueError, match="duplicate"):
        ValidationConfig.from_dict(config_data)


def test_yaml_loader_rejects_duplicate_and_nonstring_keys(tmp_path, config_data):
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(config_data) + "seed: 99\n")
    with pytest.raises(ValueError, match="unique strings"):
        load_config(path)
    path.write_text("1: bad\n")
    with pytest.raises(ValueError, match="unique strings"):
        load_config(path)


def test_config_is_detached_from_inputs_and_digest_is_order_independent(config_data, tmp_path):
    config = ValidationConfig.from_dict(config_data)
    digest = config.digest
    equivalent = dict(reversed(list(config_data.items())))
    assert ValidationConfig.from_dict(equivalent).digest == digest
    config_data["truth"]["theta"] = 100
    returned = config.to_dict()
    returned["truth"]["theta"] = 200
    assert config.to_dict()["truth"]["theta"] == 1.
    assert config.digest == digest
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(config.to_dict()))
    assert load_config(path).digest == digest


def test_scientific_cell_identity_ignores_execution_settings(config_data):
    original = ValidationConfig.from_dict(config_data)
    spec = config_data["bootstrap"]["cells"][1]
    cell = original.scientific_cell(spec)
    config_data["bootstrap"]["draws"] = 100
    config_data["bootstrap"]["estimators"].reverse()
    config_data["certification"]["replications"] = 100
    config_data["wald"]["relative_step"] = .0001
    updated = ValidationConfig.from_dict(config_data)
    assert original.digest != updated.digest
    assert cell == updated.scientific_cell(spec)
    assert len(cell["cell_id"]) == 64
    config_data["truth"]["theta"] = 2.
    changed = ValidationConfig.from_dict(config_data).scientific_cell(spec)
    assert changed["cell_id"] != cell["cell_id"]
    assert changed["mean_gap"] == .25


def test_sample_streams_reproduce_and_observations_are_shared_across_estimators(config_data):
    config = ValidationConfig.from_dict(config_data)
    cell = config.scientific_cell(config_data["bootstrap"]["cells"][1])
    source = OuSampleSource(11)
    first, first_meta = source.sample(cell, 2)
    for name in ["euler", "exact", "pfml"]:
        # Evaluation order and estimator names must not select new datasets.
        again, metadata = source.sample({**cell, "estimator": name}, 2)
        np.testing.assert_array_equal(again.times, first.times)
        np.testing.assert_array_equal(again.values, first.values)
        assert metadata == first_meta
    assert first_meta["times_sha256"] == array_digest(first.times)
    assert first_meta["observations_sha256"] == array_digest(first.values)
    assert first_meta["gap_spawn_key"] != first_meta["path_spawn_key"]
    assert first_meta["source_spawn_key"].startswith("[4,")
    next_rep, _ = source.sample(cell, 3)
    next_seed, _ = OuSampleSource(12).sample(cell, 2)
    assert not np.array_equal(first.times, next_rep.times)
    assert not np.array_equal(first.values, next_rep.values)
    assert not np.array_equal(first.values, next_seed.values)


def test_regular_sample_keeps_spacing_identical_across_reps(config_data):
    config = ValidationConfig.from_dict(config_data)
    cell = config.scientific_cell(config_data["bootstrap"]["cells"][0])
    source = OuSampleSource(11)
    first, _ = source.sample(cell, 0)
    second, _ = source.sample(cell, 1)
    np.testing.assert_allclose(np.diff(first.times), .5)
    np.testing.assert_array_equal(first.times, second.times)
    assert not np.array_equal(first.values, second.values)


@pytest.mark.parametrize("rep", [True, np.bool_(False), -1, 1.5])
def test_sample_rejects_invalid_rep_ids(config_data, rep):
    config = ValidationConfig.from_dict(config_data)
    with pytest.raises(ValueError):
        OuSampleSource(1).sample(config.scientific_cell(config_data["bootstrap"]["cells"][0]), rep)


@pytest.mark.parametrize("seed", [True, np.bool_(False), -1, 1.5])
def test_sample_source_rejects_invalid_seed(config_data, seed):
    config = ValidationConfig.from_dict(config_data)
    with pytest.raises(ValueError):
        OuSampleSource(seed).sample(config.scientific_cell(config_data["bootstrap"]["cells"][0]), 0)


def test_bootstrap_seed_is_stable_and_distinct_by_estimator_cell_and_root(config_data):
    config = ValidationConfig.from_dict(config_data)
    first, second = [config.scientific_cell(cell)["cell_id"] for cell in config_data["bootstrap"]["cells"]]
    seed_map = {name: bootstrap_seed(11, first, name) for name in ("exact", "pfml", "euler")}
    assert seed_map == {name: bootstrap_seed(11, first, name) for name in ("euler", "exact", "pfml")}
    assert len(set(seed_map.values())) == 3
    assert all(0 <= seed < 2**64 for seed in seed_map.values())
    assert bootstrap_seed(12, first, "exact") != seed_map["exact"]
    assert bootstrap_seed(11, second, "exact") != seed_map["exact"]


@pytest.mark.parametrize("seed", [True, np.bool_(True), -1, 1.5])
def test_bootstrap_seed_rejects_invalid_root(seed):
    with pytest.raises(ValueError):
        bootstrap_seed(seed, "0" * 64, "exact")


@pytest.fixture
def records():
    rows = []
    estimates = [( .9, -.1, .4), (1.1, 0., .5), (1.2, .1, .6),
                 (np.nan, np.nan, np.nan), (np.nan, np.nan, np.nan)]
    bounds = [((.8, 1.), (-.2, 0.), (.3, .45)),
              ((1.01, 1.2), (0., .1), (.45, .55))]
    for rep, values in enumerate(estimates):
        row = {"cell_id": "test", "rep": rep, "estimator": "exact", "point_valid": rep < 3,
               "interval_valid": rep < 2, "interval_status": "valid" if rep < 2 else "unavailable",
               "true_theta": 1., "true_mu": 0., "true_sigma": .5}
        for index, p in enumerate(("theta", "mu", "sigma")):
            row[p] = values[index]
            row[f"{p}_lower"], row[f"{p}_upper"] = bounds[rep][index] if rep < 2 else (np.nan, np.nan)
        rows.append(row)
    return rows


def test_conditional_and_operational_coverage_have_hand_counted_denominators(records):
    summary = CoverageSummarizer().summarize(records, planned_replications=5).set_index("parameter")
    assert (summary.n_planned == 5).all()
    assert (summary.n_point_valid == 3).all() and (summary.n_interval_valid == 2).all()
    assert (summary.point_availability == .6).all() and (summary.interval_availability == .4).all()
    theta = summary.loc["theta"]
    assert theta.n_covered == 1 and theta.coverage_conditional == .5 and theta.coverage_operational == .2
    assert theta.bias == pytest.approx(1 / 15)
    assert theta.rmse == pytest.approx(np.sqrt(.02))
    assert theta.mcse_bias == pytest.approx(np.std([-.1, .1, .2], ddof=1) / np.sqrt(3))
    assert theta.conditional_mcse == pytest.approx(np.sqrt(.5 * .5 / 2))
    assert theta.operational_mcse == pytest.approx(np.sqrt(.2 * .8 / 5))
    mu = summary.loc["mu"]
    assert mu.n_covered == 2 and mu.coverage_conditional == 1. and mu.coverage_operational == .4
    assert np.isnan(mu.relative_bias) and mu.practical_classification == "not_applicable_zero_truth"
    assert mu.negative_lower_count == 1
    assert not summary.checkpoint_in_reference_band.any()


def test_all_invalid_fits_keep_all_planned_outcomes_without_fabricated_conditional_coverage(records):
    for row in records:
        row.update(point_valid=False, interval_valid=False, interval_status="invalid_point")
        for p in ("theta", "mu", "sigma"):
            row[p] = row[f"{p}_lower"] = row[f"{p}_upper"] = np.nan
    summary = CoverageSummarizer().summarize(records, planned_replications=5)
    assert (summary.n_planned == 5).all() and (summary.n_point_valid == 0).all()
    assert (summary.n_interval_valid == 0).all() and (summary.n_covered == 0).all()
    assert summary.coverage_conditional.isna().all() and summary.conditional_wilson_upper.isna().all()
    assert (summary.coverage_operational == 0.).all()
    assert np.isfinite(summary.operational_wilson_upper).all()
    assert summary.bias.isna().all() and summary.rmse.isna().all()
    assert not summary.checkpoint_in_reference_band.any()


def test_wilson_interval_known_values_and_empty_denominator():
    assert wilson_interval(5, 10) == pytest.approx((.236593090512564, .763406909487436))
    assert wilson_interval(0, 10)[0] == pytest.approx(0.)
    assert wilson_interval(10, 10)[1] == pytest.approx(1.)
    assert all(np.isnan(value) for value in wilson_interval(0, 0))


@pytest.mark.parametrize("successes,trials,level", [
    (1, 0, .95), (-1, 2, .95), (3, 2, .95), (0, -1, .95),
    (.5, 2, .95), (1, 2.5, .95), (True, 2, .95), (0, False, .95),
    (1, 2, 0), (1, 2, 1), (1, 2, np.nan), (1, 2, True),
])
def test_wilson_rejects_impossible_counts_and_invalid_levels(successes, trials, level):
    with pytest.raises(ValueError):
        wilson_interval(successes, trials, level)


@pytest.mark.parametrize("bounds,truth,expected", [
    ((.7, .79, -.5, .5), 1., "supported_breach"),
    ((.85, .95, .11, .2), 1., "supported_breach"),
    ((.85, .95, -.2, -.11), 1., "supported_breach"),
    ((.8, .95, -.1, .1), 1., "within_criteria"),
    ((.7, .9, -.05, .05), 1., "unresolved"),
    ((.85, .95, .08, .12), 1., "unresolved"),
    ((np.nan, np.nan, np.nan, np.nan), 1., "unresolved"),
    ((.1, .2, .2, .3), 0., "not_applicable_zero_truth"),
])
def test_threshold_classification_requires_uncertainty_bounds_not_point_estimates(bounds, truth, expected):
    assert ThresholdPolicy().classify(*bounds, truth) == expected


@pytest.mark.parametrize("mutation", ["missing", "duplicate", "replacement_id", "bool_id", "float_id"])
def test_summary_rejects_missing_duplicate_or_malformed_replication_ids(records, mutation):
    if mutation == "missing":
        records.pop()
    elif mutation == "duplicate":
        records.append(deepcopy(records[0]))
    elif mutation == "replacement_id":
        records[0]["rep"] = 5
    elif mutation == "bool_id":
        records[1]["rep"] = True
    else:
        for row in records:
            row["rep"] = float(row["rep"])
    with pytest.raises(ValueError):
        CoverageSummarizer().summarize(records, planned_replications=5)


def test_summary_rejects_boolean_planned_count(records):
    with pytest.raises(ValueError):
        CoverageSummarizer().summarize(records[:1], planned_replications=True)


@pytest.mark.parametrize("mutation", ["nonbool_point", "nonbool_interval", "interval_without_point",
                                     "nonfinite_point", "reversed_bounds", "nonfinite_bounds", "changing_truth"])
def test_summary_rejects_inconsistent_validity_claims(records, mutation):
    if mutation == "nonbool_point":
        records[0]["point_valid"] = 1
    elif mutation == "nonbool_interval":
        records[0]["interval_valid"] = "True"
    elif mutation == "interval_without_point":
        records[0]["point_valid"] = False
    elif mutation == "nonfinite_point":
        records[0]["theta"] = np.inf
    elif mutation == "reversed_bounds":
        records[0]["theta_lower"] = 2.
    elif mutation == "nonfinite_bounds":
        records[0]["theta_upper"] = np.inf
    else:
        records[1]["true_theta"] = 2.
    with pytest.raises(ValueError):
        CoverageSummarizer().summarize(records, planned_replications=5)
