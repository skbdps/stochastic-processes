"""Independent mathematical checks of Week 10 covariance diagnostics."""
import numpy as np
import pytest
from scipy.integrate import quad_vec
from scipy.linalg import expm, toeplitz

from ou_irregular.week8.design import E4Cell
from ou_irregular.week8.source import phase_type
from ou_irregular.week10.theory import (
    conditional_linear_variance, conditional_mean_variance,
    conditional_pfml_variance, lag_correlations, laplace_kernel,
    pfml_mean_weights, population_predictions, population_rho,
)


@pytest.mark.parametrize("ratio", [1, 5, 20])
@pytest.mark.parametrize("floor_fraction", [0, 1e-6, .15])
def test_laplace_kernel_against_direct_matrix_quadrature(ratio, floor_fraction):
    cell = E4Cell(n=10, q=.5, mechanism="mmpp", ratio=ratio,
                  floor_fraction=floor_fraction or 1e-6,
                  floor_cv_threshold=100 if floor_fraction == 0 else 1e-9)
    ph = phase_type(cell)
    f = floor_fraction * .5
    for s in [0, 1, 2.3]:
        first, _ = quad_vec(lambda t: np.exp(-s*f)*expm(-ph.A*t) @ ph.Lambda, 0, f)
        tail, _ = quad_vec(lambda t: np.exp(-s*t)*expm(-ph.A*t) @ ph.Lambda, f, np.inf)
        np.testing.assert_allclose(laplace_kernel(cell, s), first+tail, rtol=2e-12, atol=2e-13)


def test_pfml_identity_against_regression():
    x = np.array([.3, -.2, .8, .4, .7, -.1, -.3, -.4, -.1])
    intercept, slope = np.linalg.lstsq(np.column_stack([np.ones(len(x)-1), x[:-1]]), x[1:], rcond=None)[0]
    assert 0 < slope < 1
    mu = intercept / (1-slope)
    assert pfml_mean_weights(len(x), slope) @ x == pytest.approx(mu, abs=2e-16)
    assert np.mean(x[:-1]) + (x[-1]-x[0])/((len(x)-1)*(1-slope)) == pytest.approx(mu, abs=2e-16)


@pytest.mark.parametrize("times", [np.arange(12)*.5, np.array([0, 1e-7, .01, .8, .81, 1, 6])])
@pytest.mark.parametrize("rho", [.1, .7, .98])
def test_conditional_variance_against_dense_ou_covariance(times, rho):
    theta, sigma = .8, 1.2
    covariance = sigma**2/(2*theta) * np.exp(-theta*np.abs(times[:, None]-times))
    weights = pfml_mean_weights(len(times), rho)
    assert weights.sum() == pytest.approx(1)
    dense = weights @ covariance @ weights
    assert conditional_pfml_variance(times, rho, theta, sigma) == pytest.approx(dense, rel=2e-13)
    assert conditional_mean_variance(times, theta, sigma) == pytest.approx(covariance.mean(), rel=2e-13)


@pytest.mark.parametrize("mechanism,ratio", [("mmpp", 1), ("mmpp", 5), ("mmpp", 20), ("iid", 20), ("regular", 0)])
@pytest.mark.parametrize("n", [2, 7, 1000])
def test_finite_count_against_toeplitz_weight_product(mechanism, ratio, n):
    cell = E4Cell(n=1000, q=.5, mechanism=mechanism, ratio=ratio)
    result = population_predictions(cell, n)
    weights = pfml_mean_weights(n, result["population_rho"])
    covariance = .125*toeplitz(lag_correlations(cell, n-1))
    assert result["fixed_rho_mean_variance"] == pytest.approx(weights @ covariance @ weights, rel=5e-14)
    assert result["ordinary_mean_variance"] == pytest.approx(covariance.mean(), rel=5e-14)


@pytest.mark.parametrize("ratio", [1, 5, 20])
def test_iid_fixed_rho_variance_is_exact_working_value(ratio):
    cell = E4Cell(n=1000, q=.5, mechanism="iid", ratio=ratio)
    result = population_predictions(cell)
    assert result["fixed_rho_mean_variance"] == pytest.approx(result["population_working_variance"], rel=1e-14)
    assert result["normal_working_coverage_approximation"] == pytest.approx(.95, abs=1e-14)


def test_poisson_limit_and_dependence_only_beyond_one_gap():
    for r in [1, 5, 20]:
        mmpp = E4Cell(n=1000, q=.5, mechanism="mmpp", ratio=r)
        iid = E4Cell(n=1000, q=.5, mechanism="iid", ratio=r)
        assert population_rho(mmpp) == population_rho(iid)
        m, i = lag_correlations(mmpp, 20), lag_correlations(iid, 20)
        if r == 1:
            np.testing.assert_allclose(m, i, atol=1e-15)
            assert population_rho(mmpp) == pytest.approx(2/3, abs=1e-15)
        else:
            assert np.all(m[2:] > i[2:])


@pytest.mark.parametrize("rho", [0, 1, -1, np.nan, np.inf, True])
def test_invalid_slope_rejected(rho):
    with pytest.raises(ValueError):
        pfml_mean_weights(10, rho)


@pytest.mark.parametrize("times,weights", [([0, 0, 1], [1, 1, 1]), ([0, 2, 1], [1, 1, 1]), ([0, np.nan], [1, 1]), ([0, 1], [1])])
def test_invalid_schedule_rejected(times, weights):
    with pytest.raises(ValueError):
        conditional_linear_variance(times, weights)
