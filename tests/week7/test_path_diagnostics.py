"""Selected-path math/search and exact replay tests on fixed development data.

No registered E1/E2 selected source is generated.  The one randomized fixture
uses the protocol's separate seed+1,000,000 validation namespace and n=32.
"""
from copy import deepcopy
from types import SimpleNamespace

import numpy as np
import pytest

from ou_irregular.experiments.design import Cell, ESTIMATORS, RunDesign
from ou_irregular.experiments.execution import registered_wald, source_records
from ou_irregular.inference.contracts import Sample
from ou_irregular.week4.simulation import OuSampleSource
from ou_irregular.week7 import path_diagnostics as paths


@pytest.fixture(scope="module")
def development_source():
    design = RunDesign("week7_selected_unit_validation", 20261007 + 1_000_000, 1,
                       (Cell(n=32, q=.5, cv=.5),), "validation")
    cell = design.cells[0]
    source = OuSampleSource(design.seed)
    sample, _ = source.sample(cell.to_dict(), 0)
    originals = {row["estimator"]: paths._safe(row) for row in source_records(
        cell, 0, source, {method: registered_wald(method) for method in ESTIMATORS})}
    item = {"experiment": "E1", "cell_id": cell.identity, "rep": 0,
            "reasons": ["synthetic:separate_validation_namespace"]}
    return design, sample, originals, item


def test_real_composition_replays_all_development_fields_and_retains_candidates(development_source):
    design, _, originals, item = development_source
    result = paths._one_selected(item, originals, design, "synthetic_identity")
    assert result["status"] == "evaluated"
    assert result["headline_results_replaced"] is False
    paths._payload_identity(result, item, originals, "synthetic_identity")
    assert result["original_records"] == originals
    assert set(result["original_replay"]["records"]) == set(ESTIMATORS)
    paths._compare_replay(originals, result["original_replay"]["records"])
    assert set(result["diagnostic"]["fits"]) == set(ESTIMATORS)


def test_full_replay_mismatch_is_immutable_failure_before_analysis(development_source, monkeypatch):
    design, _, originals, item = development_source
    corrupted = deepcopy(originals)
    corrupted["pfml"]["point_diagnostics_json"] = '{"invented_missing_original_field":true}'
    monkeypatch.setattr(paths, "analyze_sample", lambda *args: pytest.fail("analysis before complete replay"))
    result = paths._one_selected(item, corrupted, design, "synthetic_identity")
    assert result["status"] == "original_replay_failed"
    assert result["original_replay"]["status"] == "failed"
    assert "keys changed" in result["original_replay"]["error"]
    assert set(result["original_replay"]["attempted"]["records"]) == set(ESTIMATORS)
    assert "diagnostic" not in result
    assert result["original_records"] == corrupted
    paths._payload_identity(result, item, corrupted, "synthetic_identity")


@pytest.mark.parametrize("changed", [None, True, 7, "unavailable"])
def test_discrete_replay_fields_do_not_use_float_tolerance(changed):
    expected = {"point_valid": False, "candidate": None, "rep": 0,
                "point_diagnostics_json": '{"parameter":0.25,"accepted":false}'}
    actual = deepcopy(expected)
    actual["candidate"] = changed
    if changed is None:
        paths._compare_replay(expected, actual)
    else:
        with pytest.raises(ValueError, match="discrete"):
            paths._compare_replay(expected, actual)


def test_nested_replay_float_tolerance_is_fixed_and_all_keys_are_required():
    reference = {"point_diagnostics_json": '{"nested":{"value":1.0},"status":"boundary"}'}
    within = {"point_diagnostics_json": '{"nested":{"value":1.00000000001},"status":"boundary"}'}
    paths._compare_replay(reference, within)
    outside = {"point_diagnostics_json": '{"nested":{"value":1.000001},"status":"boundary"}'}
    with pytest.raises(ValueError, match="numerical"):
        paths._compare_replay(reference, outside)
    with pytest.raises(ValueError, match="keys"):
        paths._compare_replay(reference, {**within, "extra": None})


def test_invalid_original_points_are_skipped_not_promoted():
    sample = Sample([0., 1., 0., 2., 1.], [0., 1., 2., 3., 4.])
    rows = {method: {"estimator": method, "point_valid": False, "point_status": "boundary",
                     "interval_valid": False, "interval_status": "point_invalid",
                     "theta": 1e12, "mu": 1., "sigma": 1e6}
            for method in ESTIMATORS}
    result = paths.analyze_sample(sample, rows)
    paths._validate_diagnostic(result, rows)
    assert result["exact_profile"]["status"] == "skipped_original_point_invalid"
    for fit in result["fits"].values():
        assert fit["status"] == "skipped_original_point_invalid"
        assert "working_point" not in fit


def test_ar1_proof_scope_distinguishes_regular_from_irregular_exact():
    values = [0., 1., 0., 2., 1.]
    regular = paths._ar1_diagnostic(values, np.ones(4))
    irregular = paths._ar1_diagnostic(values, np.array([.5, 1., 1.5, 2.]))
    assert regular["phi"] == pytest.approx(-4 / 11)
    assert irregular["phi"] == regular["phi"]
    assert regular["exact_boundary_proof_applicable"] is True
    assert irregular["exact_boundary_proof_applicable"] is False
    assert regular["conclusion"] == irregular["conclusion"] == "PFML_supremum_at_phi_zero_no_finite_interior_estimate"


def test_near_regular_recognition_is_not_equal_gap_proof():
    result = paths._ar1_diagnostic([0., 1., 0., 2., 1.], [1., 1., 1. + 1e-14, 1.])
    assert result["historical_near_regular_criterion"] is True
    assert result["exact_boundary_proof_applicable"] is False


@pytest.mark.parametrize("optimizer_exception", [False, True])
def test_fixed_profile_starts_are_genuine_and_plateau_brackets_are_all_polished(monkeypatch, optimizer_exception):
    # A constant finite synthetic objective deliberately makes every interior
    # grid point a tied local-minimum bracket; this tests the difficult <= rule.
    def profile(log_theta, *args):
        return {"log_theta": np.longdouble(log_theta), "theta": np.exp(np.longdouble(log_theta)),
                "mu": np.longdouble(0), "sigma2": np.longdouble(1), "nll": np.longdouble(7)}

    starts, brackets = [], []

    def nelder(objective, x0, method, bounds, options):
        starts.append(dict(x0=x0.copy(), method=method, bounds=deepcopy(bounds), options=deepcopy(options)))
        if optimizer_exception:
            raise ArithmeticError("synthetic optimizer exception")
        return SimpleNamespace(x=x0.copy(), success=False, status=2,
                               message="synthetic fixed nonconvergence", nfev=1, nit=1)

    def bounded(objective, bounds, method, options):
        brackets.append((bounds, method, deepcopy(options)))
        if optimizer_exception:
            raise ArithmeticError("synthetic optimizer exception")
        return SimpleNamespace(x=sum(bounds) / 2, success=True, message="synthetic", nfev=1, nit=1)

    monkeypatch.setattr(paths, "profile_at", profile)
    monkeypatch.setattr(paths, "full_exact_nll", lambda *args: np.longdouble("7.00000001"))
    monkeypatch.setattr(paths, "profile_limits", lambda *args: {
        "brownian": {"alpha": np.longdouble(0), "sigma2": np.longdouble(1), "nll": np.longdouble(8)},
        "iid": {"mu": np.longdouble(0), "variance": np.longdouble(1), "nll": np.longdouble(8)}})
    monkeypatch.setattr(paths, "minimize", nelder)
    monkeypatch.setattr(paths, "minimize_scalar", bounded)
    point = np.array([11.95, .2, -.4])
    result = paths._profile_search([0., 1., 0.], [1., 1.], point)
    paths._validate_profile(result, point.tolist())
    assert len(starts) == 5
    assert len(brackets) == len(result["grid_log_theta"]) - 2
    assert np.max(np.diff(result["grid_log_theta"])) <= .05
    for offset, call, record in zip(paths.SETTINGS["start_offsets"], starts, result["local_searches"]):
        expected = np.clip(point[0] + offset, result["bounds"]["lower"], result["bounds"]["upper"])
        np.testing.assert_array_equal(call["x0"], [expected])
        assert call["method"] == "Nelder-Mead"
        assert call["options"]["xatol"] == call["options"]["fatol"] == 1e-10
        assert call["options"]["maxiter"] == 2000
        assert call["options"]["maxfev"] == 4000
        assert call["options"]["adaptive"] is False
        assert abs(np.diff(call["options"]["initial_simplex"][:, 0])[0]) == pytest.approx(.1)
        assert record["success"] is False
        assert record["candidate"]["status"] == "finite"
    assert result["disposition"] == "stored_full_fit_matches_finite_search_within_diagnostic_resolution"
    assert result["global_optimum_proved"] is False


def test_both_profile_caps_retain_unresolved_tail_evidence():
    bounds = paths._profile_bounds(np.array([1e-30, 1e30], dtype=np.longdouble))
    assert bounds["lower"] == -40 and bounds["upper"] == 40
    assert bounds["lower_cap_limited"] and bounds["upper_cap_limited"]
    assert bounds["lower_tail_reached"] is False
    assert bounds["upper_tail_reached"] is False


@pytest.mark.parametrize("values,gaps", [([1, 1, 1], [1, 1]), ([0, 1, 0], [1, 0]),
                                         ([0, 1, np.nan], [1, 1]), ([0, 1], [1])])
def test_degenerate_and_invalid_profile_data_are_not_regularized(values, gaps):
    with pytest.raises(ValueError):
        paths.profile_at(0., values, gaps)


def test_descriptor_binds_literal_selection_precision_and_fixed_numerics():
    selection = {"sources": []}
    first = paths._diagnostic_design(selection, b'{"sources":[]}\n')
    second = paths._diagnostic_design(selection, b'{ "sources": [] }\n')
    assert first.digest != second.digest
    assert first.to_dict()["selection"] == second.to_dict()["selection"]
    assert first.to_dict()["settings"] == paths.SETTINGS
    assert first.to_dict()["precision"] == paths.precision_record()
