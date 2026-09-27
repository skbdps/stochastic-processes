"""Independent model transformations required by the experiment design.

These two cell types check units and affine state scaling without relying on
optimizer convergence or fitted Monte Carlo averages.
"""
import numpy as np
import pytest

from ou_irregular.estimators import exact_neg_loglik, pfml_neg_loglik, euler_neg_loglik
from ou_irregular.ou import simulate_ou
from ou_irregular.spacing import generate_times


@pytest.mark.parametrize("cv", [0.0, 2.0])
@pytest.mark.parametrize("objective", [exact_neg_loglik, pfml_neg_loglik, euler_neg_loglik])
def test_state_and_time_units_preserve_likelihood(cv, objective):
    times = generate_times(60, 0.3, cv, np.random.default_rng(611)).times
    theta, mu, sigma = 1.2, -0.4, 0.7
    x = simulate_ou(theta, mu, sigma, times, np.random.default_rng(612))
    parameters = [np.log(theta), mu, np.log(sigma)]
    original = objective(parameters, x, times)

    # X'=a+bX: theta'=theta, mu'=a+b*mu, sigma'=b*sigma for b>0.
    # Each conditional Gaussian density gains Jacobian 1/b, hence NLL grows
    # by (n_observations-1)*log(b); the initial density is conditioned out.
    shift, scale = 1.1, 2.4
    transformed = [np.log(theta), shift + scale * mu, np.log(scale * sigma)]
    assert objective(transformed, shift + scale * x, times) == pytest.approx(
        original + (len(x) - 1) * np.log(scale), rel=1e-10, abs=1e-9)

    # Change time units by t'=t/c: theta'=c*theta, sigma'=sqrt(c)*sigma.
    # Products theta*d and stationary variance sigma²/(2theta) stay fixed.
    speed = 3.0
    time_parameters = [np.log(speed * theta), mu, np.log(np.sqrt(speed) * sigma)]
    assert objective(time_parameters, x, times / speed) == pytest.approx(original, rel=1e-10, abs=1e-9)
    replayed = simulate_ou(speed * theta, mu, np.sqrt(speed) * sigma,
                          times / speed, np.random.default_rng(612))
    np.testing.assert_allclose(replayed, x, rtol=1e-12, atol=1e-12)
