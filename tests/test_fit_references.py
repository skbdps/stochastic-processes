"""Week 3 acceptance: fitted estimates, independent NLLs, and unit changes.

Fixture definitions and tolerances are fixed before evaluating the patched
solver. No assertion rewards closeness to simulation truth. The same paths are
used for every comparison; a failed fit is never replaced by another seed.
"""
from functools import lru_cache
import hashlib
from pathlib import Path

import numpy as np
import pytest

from ou_irregular.config import load_config
from ou_irregular.estimators import ESTIMATORS, fit_estimator
from ou_irregular.ou import simulate_ou
from ou_irregular.runner import replication_streams
from ou_irregular.spacing import generate_times
from reference_estimators import (
    NLL_ATOL, NLL_RTOL, PARAMETER_ATOL, PARAMETER_RTOL,
    ReferenceInvalid, conditional_nll, reference_fit,
)

ROOT = Path(__file__).resolve().parents[1]
# The third fixture is the exact cell/replication reported in feedback W3-01.
FIXTURES = {
    "regular_q05_rep0": (0.5, 0.0, 0),
    "irregular_q01_cv1_rep7": (0.1, 1.0, 7),
    "reported_q05_cv2_rep0": (0.5, 2.0, 0),
}
STATE_SCALES = (1e-6, 1.0, 1e6)
OFFSET_IN_SCALE_UNITS = (0.0, 100.0)
TIME_SPEEDS = (1e-3, 1.0, 1e3)


@lru_cache(maxsize=None)
def fixture_path(name):
    """Reconstruct the original smoke stream, floor policy, and stationary X0."""
    # Read the preserved run snapshot, so a later experimental configuration
    # edit cannot silently change these historical acceptance fixtures.
    config = load_config(ROOT / "artifacts/week3_smoke/config.yaml")
    q, cv, rep = FIXTURES[name]
    cell = next(cell for cell in config.cells()
                if cell["theta_mean_gap"] == q and cell["cv"] == cv and cell["n"] == 250)
    gap_rng, path_rng, _ = replication_streams(config.seed, cell["cell_id"], rep)
    times = generate_times(cell["n"], cell["mean_gap"], cell["cv"], gap_rng,
                           floor_fraction=cell["floor_fraction"],
                           floor_cv_threshold=cell["floor_cv_threshold"]).times
    observations = simulate_ou(cell["true_theta"], cell["true_mu"], cell["true_sigma"], times, path_rng)
    return observations, times


@lru_cache(maxsize=None)
def baseline_reference(name, estimator):
    x, times = fixture_path(name)
    return reference_fit(estimator, x, times)


def _normal_coordinates(parameters, original_x, original_times):
    """Unit-free coordinates for the gate; all scales come from observed data."""
    theta, mu, sigma = parameters
    state_scale = np.std(original_x)
    time_scale = np.mean(np.diff(original_times))
    return np.array([theta * time_scale,
                     (mu - np.mean(original_x)) / state_scale,
                     sigma * np.sqrt(time_scale) / state_scale])


def check_transformation(name, estimator, *, state_scale=1.0, offset_units=0.0, time_speed=1.0):
    """Return numerical evidence after asserting reference and transformation gates.

    X'=a+bX and t'=t/c imply theta'=c*theta, mu'=a+b*mu,
    sigma'=b*sqrt(c)*sigma, and NLL'=NLL+K*log(b). For the
    substantial offsets a=100*b, retained variation remains many orders above
    floating-point spacing even at b=1e-6 and b=1e6.
    """
    x, times = fixture_path(name)
    offset = offset_units * state_scale
    transformed_x = offset + state_scale * x
    transformed_times = times / time_speed
    # Guard the transformation itself: compare its round-trip data error in
    # units of the original path's observed SD, independently of fitted values.
    resolution_error = np.max(np.abs((transformed_x - offset) / state_scale - x)) / np.std(x)
    assert resolution_error < 1e-10, "fixture transformation lost floating-point resolution"

    result = fit_estimator(estimator, transformed_x, transformed_times)
    assert result["optimizer_success"], result
    assert result["fit_status"] == "valid_interior", result
    assert result["valid_for_point_summary"], result
    candidate = np.array([result[key] for key in ("theta", "mu", "sigma")])
    reference = reference_fit(estimator, transformed_x, transformed_times)
    original_reference = baseline_reference(name, estimator)

    def undo(parameters):
        theta, mu, sigma = parameters
        return np.array([theta / time_speed, (mu - offset) / state_scale,
                         sigma / (state_scale * np.sqrt(time_speed))])

    normalized_candidate = _normal_coordinates(undo(candidate), x, times)
    normalized_reference = _normal_coordinates(undo(reference.parameters), x, times)
    original_normalized_reference = _normal_coordinates(original_reference.parameters, x, times)
    for expected in (normalized_reference, original_normalized_reference):
        np.testing.assert_allclose(normalized_candidate, expected,
                                   rtol=PARAMETER_RTOL, atol=PARAMETER_ATOL)
    nll_tolerance = NLL_ATOL + NLL_RTOL * abs(reference.nll)
    independent_candidate_nll = conditional_nll(estimator, candidate, transformed_x, transformed_times)
    assert abs(result["nll"] - reference.nll) <= nll_tolerance
    assert abs(independent_candidate_nll - reference.nll) <= nll_tolerance
    assert abs(result["nll"] - independent_candidate_nll) <= nll_tolerance
    expected_nll = original_reference.nll + (len(x) - 1) * np.log(state_scale)
    assert abs(result["nll"] - expected_nll) <= NLL_ATOL + NLL_RTOL * abs(expected_nll)
    difference = np.abs(normalized_candidate - normalized_reference)
    return {
        "fixture": name, "estimator": estimator, "state_scale": state_scale,
        "offset": offset, "time_speed": time_speed,
        "fit_status": result["fit_status"], "nll": result["nll"],
        "reference_nll": reference.nll, "absolute_nll_difference": abs(result["nll"] - reference.nll),
        "allowed_nll_difference": nll_tolerance,
        "max_normalized_parameter_difference": np.max(difference),
        "max_parameter_tolerance_fraction": np.max(difference / (PARAMETER_ATOL + PARAMETER_RTOL * np.abs(normalized_reference))),
        "state_roundtrip_error_in_sd": resolution_error,
        "original_times_sha256": hashlib.sha256(times.astype("<f8").tobytes()).hexdigest(),
        "original_observations_sha256": hashlib.sha256(x.astype("<f8").tobytes()).hexdigest(),
    }


@pytest.mark.parametrize("name", FIXTURES)
@pytest.mark.parametrize("estimator", ESTIMATORS)
@pytest.mark.parametrize("state_scale", STATE_SCALES)
@pytest.mark.parametrize("offset_units", OFFSET_IN_SCALE_UNITS)
def test_fitted_state_units_against_independent_references(name, estimator, state_scale, offset_units):
    check_transformation(name, estimator, state_scale=state_scale, offset_units=offset_units)


@pytest.mark.parametrize("name", FIXTURES)
@pytest.mark.parametrize("estimator", ESTIMATORS)
@pytest.mark.parametrize("time_speed", TIME_SPEEDS)
def test_fitted_time_units_against_independent_references(name, estimator, time_speed):
    check_transformation(name, estimator, time_speed=time_speed)


@pytest.mark.parametrize("estimator", ESTIMATORS)
def test_lost_resolution_is_explicit_degeneracy(estimator):
    x, times = fixture_path("reported_q05_cv2_rep0")
    # Here variation has actually vanished from the input doubles. This is a
    # diagnosed data-resolution failure, never evidence against unit invariance.
    collapsed = x + 1e20
    assert np.unique(collapsed).size == 1
    result = fit_estimator(estimator, collapsed, times)
    assert result["fit_status"] == "degenerate"
    assert not result["valid_for_point_summary"]
    with pytest.raises(ReferenceInvalid, match="degenerate"):
        reference_fit(estimator, collapsed, times)


def test_negative_ar_slope_is_ou_boundary_but_not_euler_invalidity():
    # Reported W3-02 fixture: outside the smoke grid and never substituted for a
    # smoke path. The exact/PFML regular OU interior requires 0<phi<1; Euler's
    # coefficient 1-theta*d may be negative while theta remains admissible.
    times = np.arange(250) * 5.0
    x = simulate_ou(1.0, 0.0, 0.5, times, np.random.default_rng(1))
    for estimator in ("exact", "pfml"):
        result = fit_estimator(estimator, x, times)
        assert result["fit_status"] in ("boundary", "boundary_suspected")
        assert not result["valid_for_point_summary"]
    with pytest.raises(ReferenceInvalid, match="outside"):
        reference_fit("pfml", x, times)
    euler = fit_estimator("euler", x, times)
    reference = reference_fit("euler", x, times)
    assert euler["valid_for_point_summary"]
    assert 1 - euler["theta"] * 5 < 0
    np.testing.assert_allclose([euler[key] for key in ("theta", "mu", "sigma")],
                               reference.parameters, rtol=PARAMETER_RTOL, atol=PARAMETER_ATOL)
    assert abs(euler["nll"] - reference.nll) <= NLL_ATOL + NLL_RTOL * abs(reference.nll)


def comparison_records():
    """Re-run these fixed comparisons for a CSV evidence table, without pytest."""
    records = []
    for name in FIXTURES:
        for estimator in ESTIMATORS:
            for scale in STATE_SCALES:
                for offset in OFFSET_IN_SCALE_UNITS:
                    records.append(check_transformation(name, estimator, state_scale=scale, offset_units=offset))
            for speed in TIME_SPEEDS:
                records.append(check_transformation(name, estimator, time_speed=speed))
    return records
