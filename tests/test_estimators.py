"""Independent likelihood/closed-form checks and auditable retry failure tests."""

import numpy as np
import pytest
from scipy.optimize import OptimizeResult
from scipy.stats import norm

from ou_irregular import ou
from ou_irregular.estimators import (
    ESTIMATORS, exact_neg_loglik, pfml_neg_loglik, euler_neg_loglik, fit_estimator,
)

LIKELIHOODS = (exact_neg_loglik, pfml_neg_loglik, euler_neg_loglik)


def test_exact_likelihood_matches_independent_gaussian_transition_sum():
    theta, mu, sigma = 0.8, 1.2, 0.5
    x = np.array([0.3, 0.2, 0.7, 1.0])
    times = np.array([0.0, 0.03, 0.5, 1.7])
    gaps = np.diff(times)
    mean = mu + (x[:-1] - mu) * np.exp(-theta * gaps)
    variance = sigma ** 2 / (2 * theta) * (1 - np.exp(-2 * theta * gaps))
    expected = -norm.logpdf(x[1:], loc=mean, scale=np.sqrt(variance)).sum()
    actual = exact_neg_loglik([np.log(theta), mu, np.log(sigma)], x, times)
    assert actual == pytest.approx(expected, rel=1e-13)


def test_pfml_uses_realized_mean_and_ignores_individual_gap_order():
    x = np.array([0.3, 0.2, 0.7, 1.0])
    times = np.array([0.0, 0.03, 0.5, 1.7])
    params = np.array([0.1, 0.2, -0.7])
    regular = np.arange(x.size) * np.diff(times).mean()
    assert pfml_neg_loglik(params, x, times) == pytest.approx(exact_neg_loglik(params, x, regular))
    assert not np.isclose(pfml_neg_loglik(params, x, times), exact_neg_loglik(params, x, times))


@pytest.mark.parametrize("estimator", ESTIMATORS)
def test_regular_fit_matches_independent_ols_closed_form(estimator):
    dt = 0.2
    times = np.arange(1200) * dt
    x = ou.simulate_ou(0.8, 1.2, 0.5, times, np.random.default_rng(14))
    # OLS maximizes the unrestricted Gaussian AR(1) conditional likelihood.
    # This sample has 0<b<1, so all mapped OU/Euler parameters are interior.
    matrix = np.column_stack([np.ones(x.size - 1), x[:-1]])
    intercept, slope = np.linalg.lstsq(matrix, x[1:], rcond=None)[0]
    innovation_variance = np.mean((x[1:] - matrix @ [intercept, slope]) ** 2)
    assert 0 < slope < 1
    mu = intercept / (1 - slope)
    if estimator == "euler":
        theta = (1 - slope) / dt
        sigma = np.sqrt(innovation_variance / dt)
    else:
        theta = -np.log(slope) / dt
        sigma = np.sqrt(2 * theta * innovation_variance / (1 - slope ** 2))
    fitted = fit_estimator(estimator, x, times)
    assert fitted["success"], fitted
    np.testing.assert_allclose([fitted["theta"], fitted["mu"], fitted["sigma"]],
                               [theta, mu, sigma], rtol=2e-4, atol=2e-5)


def test_exact_and_pfml_are_identical_on_regular_schedule():
    times = np.arange(200) * 0.25
    x = ou.simulate_ou(1, 0, 0.5, times, np.random.default_rng(16))
    params = [0.2, -0.1, -0.4]
    assert exact_neg_loglik(params, x, times) == pfml_neg_loglik(params, x, times)
    exact, pfml = [fit_estimator(name, x, times) for name in ("exact", "pfml")]
    np.testing.assert_allclose([exact[k] for k in ("theta", "mu", "sigma", "nll")],
                               [pfml[k] for k in ("theta", "mu", "sigma", "nll")],
                               rtol=1e-10, atol=1e-10)


def test_subpicosecond_gaps_are_not_clamped():
    times = np.array([0.0, 1e-15, 3e-15])
    x = np.array([0.0, 1e-8, 3e-8])
    theta, mu, sigma = 1.0, 0.0, 0.5
    gaps = np.diff(times)
    mean = x[:-1] * np.exp(-theta * gaps)
    variance = sigma ** 2 / (2 * theta) * -np.expm1(-2 * theta * gaps)
    expected = -norm.logpdf(x[1:], mean, np.sqrt(variance)).sum()
    params = [np.log(theta), mu, np.log(sigma)]
    assert exact_neg_loglik(params, x, times) == pytest.approx(expected, rel=1e-13)
    assert abs(exact_neg_loglik(params, x, times) - ou.ou_neg_loglik(params, x, times)) > 1


def test_exact_variance_handles_underflow_in_theta_times_gap():
    times = np.array([0.0, 1e-50, 2e-50])
    x = np.array([0.0, 1e-25, 2e-25])
    sigma = 0.7
    expected = -norm.logpdf(x[1:], x[:-1], sigma * np.sqrt(np.diff(times))).sum()
    assert exact_neg_loglik([np.log(1e-300), 0, np.log(sigma)], x, times) == pytest.approx(expected)


def test_log_variance_avoids_intermediate_sigma_squared_overflow():
    # sigma and the transition variance are representable although sigma**2
    # overflows. The exact variance is exp(400)/2 for these very long OU gaps.
    value = exact_neg_loglik([400.0, 0.0, 400.0], [0, 0, 0], [0, 1, 2])
    assert value == pytest.approx(np.log(2 * np.pi) + 400.0 - np.log(2.0))


@pytest.mark.parametrize("likelihood", LIKELIHOODS)
@pytest.mark.parametrize("params", [[1000, 0, 0], [-1000, 0, 0], [0, 0, 1000],
                                   [0, 0, -1000], [0, np.nan, 0], [np.inf, 0, 0],
                                   [0, 0, -700], [0, 0]])
def test_invalid_or_unrepresentable_parameters_return_infinity(likelihood, params):
    assert likelihood(params, [0.1, 0.2, 0.3], [0, 1, 2]) == np.inf


@pytest.mark.parametrize("x,times", [
    ([0, 1], [0, 1]), ([0, 1, 2], [0, 1]), ([[0, 1, 2]], [0, 1, 2]),
    ([0, 1, 2], [[0, 1, 2]]), ([0, np.nan, 2], [0, 1, 2]),
    ([0, 1, 2], [0, np.inf, 2]), ([0, 1, 2], [0, 0, 1]),
    ([0, 1, 2], [0, 2, 1]), ([0, 1, 2], [-1e308, 1e308, 1.5e308]),
])
def test_invalid_data_raise_before_optimization(monkeypatch, x, times):
    def unexpected(*args, **kwargs):
        pytest.fail("optimizer must not be called for invalid observations")

    monkeypatch.setattr(ou, "fit", unexpected)
    with pytest.raises(ValueError):
        fit_estimator("exact", x, times)
    for likelihood in LIKELIHOODS:
        with pytest.raises(ValueError):
            likelihood([0, 0, 0], x, times)


@pytest.mark.parametrize("name,value", [("estimator", "naive"), ("estimator", None),
                                        ("n_starts", 0), ("n_starts", True),
                                        ("retry_starts", 2.5), ("retry_starts", 0)])
def test_invalid_controls_raise(name, value):
    kwargs = dict(estimator="exact", x=[0, 1, 2], times=[0, 1, 2])
    kwargs[name] = value
    with pytest.raises(ValueError):
        fit_estimator(**kwargs)


def _result(objective, success, starts=3):
    return dict(theta=1.0, mu=0.0, sigma=0.5, nll=objective, success=success,
                nit=12, status=0 if success else 2, message="test result",
                n_starts=starts, n_finite=starts, n_successful=starts if success else 0,
                selected_start=0, best_nll=objective, objective_gap=0.0,
                objective_tolerance=ou.OBJECTIVE_TIE_ATOL + ou.OBJECTIVE_TIE_RTOL * abs(objective))


def test_initial_success_never_retries(monkeypatch):
    calls = []

    def fake_fit(*args, n_starts, **kwargs):
        calls.append(n_starts)
        return _result(-100, True, n_starts)

    monkeypatch.setattr(ou, "fit", fake_fit)
    fitted = fit_estimator("exact", [0.1, 0.4, 0.2], [0, 1, 2])
    assert calls == [3]
    assert fitted["initial_success"] and not fitted["retried"]
    assert fitted["attempted_starts"] == 3
    assert fitted["error"] is None


def test_retry_has_five_deterministic_starts_and_does_not_change_model(monkeypatch):
    starts, objectives = [], []

    def fake_minimize(nll, start, args, **kwargs):
        starts.append(start.copy())
        objectives.append(nll(np.array([0, 0, -0.7]), *args))
        success = len(starts) > 3
        return OptimizeResult(x=np.array([0.0, 0.0, -0.7]), fun=-100.0,
                              success=success, status=0 if success else 2, nit=1,
                              message="test success" if success else "test failure")

    monkeypatch.setattr(ou, "minimize", fake_minimize)
    fitted = fit_estimator("euler", [0.1, 0.4, 0.2, 0.3], [0, 0.5, 1.5, 2])
    assert fitted["success"] and fitted["retried"] and not fitted["initial_success"]
    assert fitted["attempted_starts"] == len(starts) == 8
    assert fitted["n_starts"] == 5
    np.testing.assert_array_equal(starts[:3], starts[3:6])
    assert len({tuple(start) for start in starts[3:]}) == 5
    np.testing.assert_array_equal(objectives, np.full(8, objectives[0]))


@pytest.mark.parametrize("retry_objective,expected_success", [(-999.0, False), (-1000 + 5e-8, True), (-1001, True)])
def test_retry_respects_best_objective_instead_of_accepting_worse_convergence(monkeypatch, retry_objective, expected_success):
    outcomes = iter([_result(-1000.0, False), _result(retry_objective, True, 5)])
    monkeypatch.setattr(ou, "fit", lambda *args, **kwargs: next(outcomes))
    fitted = fit_estimator("exact", [0.1, 0.4, 0.2], [0, 1, 2])
    assert fitted["success"] == expected_success
    assert fitted["nll"] == (retry_objective if expected_success else -1000)
    assert fitted["best_nll"] == min(-1000, retry_objective)
    assert fitted["attempted_starts"] == 8


def test_optimizer_exceptions_record_nan_failure_and_are_not_dropped(monkeypatch):
    calls = []

    def fail(*args, n_starts, **kwargs):
        calls.append(n_starts)
        raise RuntimeError("no finite solution")

    monkeypatch.setattr(ou, "fit", fail)
    fitted = fit_estimator("pfml", [0.1, 0.4, 0.2], [0, 1, 2])
    assert calls == [3, 5]
    assert not fitted["success"] and fitted["retried"]
    assert np.all(np.isnan([fitted["theta"], fitted["mu"], fitted["sigma"]]))
    assert "initial fit: RuntimeError" in fitted["error"]
    assert "retry fit: RuntimeError" in fitted["error"]


def test_retry_can_recover_initial_exception_and_preserves_error(monkeypatch):
    calls = []

    def fake_fit(*args, n_starts, **kwargs):
        calls.append(n_starts)
        if len(calls) == 1:
            raise RuntimeError("initial optimizer error")
        return _result(-1000, True, n_starts)

    monkeypatch.setattr(ou, "fit", fake_fit)
    fitted = fit_estimator("exact", [0.1, 0.4, 0.2], [0, 1, 2])
    assert fitted["success"] and fitted["retried"]
    assert fitted["attempted_starts"] == 8
    assert "initial optimizer error" in fitted["error"]
