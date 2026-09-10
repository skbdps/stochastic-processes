import numpy as np
import pytest
from scipy.optimize import OptimizeResult
from ou_irregular import *
import ou_irregular.ou as ou_module

THETA, MU, SIGMA = 1.0, 0.0, 0.5


def test_one_step_law_matches_transition_moments():
    rng = np.random.default_rng(1)
    x0, d, N = 1.0, 0.3, 200_000
    m, v = transition_moments(THETA, MU, SIGMA, x0, d)
    draws = np.array([simulate_ou(THETA, MU, SIGMA, [0.0, d], rng, x0=x0)[1] for _ in range(20_000)])
    assert abs(draws.mean() - m) < 4 * np.sqrt(v / draws.size)
    assert abs(draws.var() / v - 1) < 0.03


def test_stationary_variance_and_autocorr():
    rng = np.random.default_rng(2)
    t = equidistant_times(50_000, 0.1)
    x = simulate_ou(THETA, MU, SIGMA, t, rng)
    assert abs(np.var(x) / stationary_var(THETA, SIGMA) - 1) < 0.03
    rho = np.corrcoef(x[:-1], x[1:])[0, 1]
    assert abs(rho / np.exp(-THETA * 0.1) - 1) < 0.03


def test_optimizer_matches_closed_form_on_equidistant():
    rng = np.random.default_rng(3)
    t = equidistant_times(1000, 0.1)
    x = simulate_ou(THETA, MU, SIGMA, t, rng)
    f = fit(ou_neg_loglik, x, t)
    th, mu, sg = ar1_closed_form(x, 0.1)
    assert abs(f["theta"] / th - 1) < 1e-3
    assert abs(f["mu"] - mu) < 1e-3
    assert abs(f["sigma"] / sg - 1) < 1e-3


def test_naive_equals_exact_on_equidistant():
    rng = np.random.default_rng(4)
    t = equidistant_times(500, 0.2)
    x = simulate_ou(THETA, MU, SIGMA, t, rng)
    p = np.array([0.1, 0.05, -0.6])
    assert np.isclose(naive_neg_loglik(p, x, t), ou_neg_loglik(p, x, t))


def test_gamma_gaps_mean_and_cv():
    rng = np.random.default_rng(5)
    for cv in (0.5, 1.0, 2.0):
        t, _ = gamma_gap_times(200_001, 0.5, cv, rng)
        g = np.diff(t)
        assert abs(g.mean() / 0.5 - 1) < 0.02
        assert abs(g.std() / g.mean() / cv - 1) < 0.03


def test_small_gap_exact_approaches_euler():
    # D3 small-gap limit: exact transition -> Euler moments as d -> 0
    d = 1e-4
    m, v = transition_moments(THETA, MU, SIGMA, 0.7, d)
    assert abs(m - (0.7 + THETA * (MU - 0.7) * d)) < 1e-7
    assert abs(v / (SIGMA ** 2 * d) - 1) < 1e-3


def _optimizer_result(objective, success, params=(0.0, 0.0, -0.7)):
    return OptimizeResult(x=np.asarray(params), fun=objective, success=success,
                          status=0 if success else 2, nit=12,
                          message="converged" if success else "line search failed")


def _mock_fit(monkeypatch, results):
    outcomes = iter(results)
    monkeypatch.setattr(ou_module, "minimize", lambda *args, **kwargs: next(outcomes))
    x = np.array([0.1, 0.5, 0.2, 0.4, 0.3])
    return fit(ou_neg_loglik, x, np.arange(x.size), n_starts=len(results))


def test_fit_prefers_converged_result_with_numerically_tied_objective(monkeypatch):
    result = _mock_fit(monkeypatch, [
        _optimizer_result(-1000.0, False),
        _optimizer_result(-1000.0 + 5e-8, True, params=(0.2, 0.1, -0.6)),
        _optimizer_result(-999.0, True),
    ])
    assert result["success"]
    assert result["selected_start"] == 1
    assert result["theta"] == pytest.approx(np.exp(0.2))
    assert result["status"] == 0
    assert result["message"] == "converged"
    assert result["nit"] == 12
    assert result["best_nll"] == -1000.0
    assert 0 < result["objective_gap"] <= result["objective_tolerance"]
    assert result["n_starts"] == result["n_finite"] == 3
    assert result["n_successful"] == 2


def test_fit_does_not_accept_materially_worse_converged_result(monkeypatch):
    result = _mock_fit(monkeypatch, [
        _optimizer_result(-1000.0, False),
        _optimizer_result(-999.0, True),
    ])
    assert not result["success"]
    assert result["nll"] == -1000.0
    assert result["selected_start"] == 0
    assert result["objective_gap"] == 0.0
    assert result["status"] == 2
    assert result["message"] == "line search failed"
    assert result["n_successful"] == 1


def test_fit_raises_when_all_results_are_nonfinite(monkeypatch):
    with pytest.raises(RuntimeError, match="no finite valid solution from 3 starts"):
        _mock_fit(monkeypatch, [
            _optimizer_result(np.nan, True),
            _optimizer_result(np.inf, False),
            _optimizer_result(-np.inf, False),
        ])


def test_fit_discards_nonfinite_parameters_and_objectives(monkeypatch):
    result = _mock_fit(monkeypatch, [
        _optimizer_result(np.nan, True),
        _optimizer_result(-1001.0, True, params=(np.inf, 0.0, -0.7)),
        _optimizer_result(-1000.0, True),
    ])
    assert result["success"]
    assert result["selected_start"] == 2
    assert result["n_finite"] == result["n_successful"] == 1


@pytest.mark.parametrize("n_starts", [1, 3, 5, np.int64(5)])
def test_fit_honors_requested_starts_deterministically(monkeypatch, n_starts):
    starts = []

    def mock_minimize(nll, start, **kwargs):
        starts.append(start.copy())
        return _optimizer_result(-1000.0, True)

    monkeypatch.setattr(ou_module, "minimize", mock_minimize)
    x = np.array([0.1, 0.5, 0.2, 0.4, 0.3])
    times = np.arange(x.size)
    result = fit(ou_neg_loglik, x, times, n_starts=n_starts)
    assert len(starts) == result["n_starts"] == n_starts
    assert len({tuple(start) for start in starts}) == n_starts
    initial_starts = np.asarray(starts)
    starts.clear()
    fit(ou_neg_loglik, x, times, n_starts=n_starts)
    np.testing.assert_array_equal(starts, initial_starts)


@pytest.mark.parametrize("n_starts", [0, -1, 1.5, "3", True, np.bool_(False), None])
def test_fit_rejects_invalid_start_counts_before_optimizing(monkeypatch, n_starts):
    def unexpected_minimize(*args, **kwargs):
        pytest.fail("invalid n_starts should be rejected before optimization")

    monkeypatch.setattr(ou_module, "minimize", unexpected_minimize)
    with pytest.raises(ValueError, match="n_starts must be a positive integer"):
        fit(ou_neg_loglik, np.array([0.1, 0.2, 0.3]), np.arange(3), n_starts=n_starts)
