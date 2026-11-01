"""Deterministic checks of the fitted mean's pole and working-Hessian algebra."""
import numpy as np
import pytest

from ou_irregular.estimators import fit_estimator
from ou_irregular.inference.adapters import OuLikelihood, Week3Estimator
from ou_irregular.inference.contracts import Sample
from ou_irregular.inference.wald import WaldInference
from ou_irregular.week10.theory import conditional_pfml_variance, pfml_mean_weights


def _regression(values):
    x, y = values[:-1], values[1:]
    centered = x - x.mean()
    slope = centered @ (y - y.mean()) / (centered @ centered)
    intercept = y.mean() - slope*x.mean()
    residual = y - intercept - slope*x
    return intercept, slope, residual


@pytest.mark.parametrize("prefix", [
    [0., 1., 3.], [.3, -.2, .8, .4, .7, -.1, -.3],
    [2., -1., 5., 4., -3., 2., 8.],
])
def test_cyclic_shift_proves_nonzero_pole_numerator(prefix):
    z = np.asarray(prefix)
    centered = z-z.mean()
    sxx = centered @ centered
    _, cyclic_slope, _ = _regression(np.append(z, z[0]))
    cyclic_differences = np.diff(np.append(z, z[0]))
    assert 1-cyclic_slope == pytest.approx(cyclic_differences @ cyclic_differences/(2*sxx))
    assert cyclic_slope < 1
    b = centered[-1]/sxx
    a = centered[:-1] @ z[1:]/sxx
    assert b != 0
    root = (1-a)/b
    assert root != pytest.approx(z[0])
    _, unit_slope, residual = _regression(np.append(z, root))
    assert unit_slope == pytest.approx(1)
    assert residual @ residual > 0


def test_concrete_pole_has_nonzero_residual_variance():
    intercept, slope, residual = _regression(np.array([0., 1., 3., 4.2]))
    assert intercept == pytest.approx(1.4)
    assert slope == pytest.approx(1)
    np.testing.assert_allclose(residual, [-.4, .6, -.2], atol=1e-15)
    assert residual @ residual == pytest.approx(.56)


@pytest.mark.parametrize("delta", [1e-2, 1e-3, 1e-4, 1e-5, 1e-6])
def test_unchanged_point_valid_rule_accepts_approach_to_pole(delta):
    values = np.array([0., 1., 3., 4.2-delta])
    times = np.arange(4, dtype=float)
    fit = fit_estimator("pfml", values, times)
    assert fit["valid_for_point_summary"]
    assert fit["fit_status"] == "valid_interior"
    assert fit["unrestricted_phi"] == pytest.approx(1-5*delta/14, abs=4e-16)
    assert fit["mu"] == pytest.approx(.4+3.92/delta, rel=2e-8)
    fixed_mean = pfml_mean_weights(4, .7) @ values
    assert 5.9 < fixed_mean < 6.1
    fixed_variance = conditional_pfml_variance(times, .7)
    assert 0 < fixed_variance < 1


@pytest.mark.parametrize("estimator", ["exact", "pfml", "euler"])
@pytest.mark.parametrize("offset,state_scale,time_scale", [(0., 1., 1.), (7., 2.5, .3)])
def test_observed_hessian_mean_entry_matches_joint_regression_formula(
        estimator, offset, state_scale, time_scale):
    base = np.array([.3, -.2, .8, .4, .7, -.1, -.3, -.4, -.1])
    values = offset+state_scale*base
    times = time_scale*np.arange(len(values))
    intercept, slope, residual = _regression(values)
    assert 0 < slope < 1
    mu = intercept/(1-slope)
    centered = values[:-1]-values[:-1].mean()
    k = len(values)-1
    expected = (residual @ residual/k)/(1-slope)**2 * (
        1/k+(mu-values[:-1].mean())**2/(centered @ centered))
    result = WaldInference(Week3Estimator(estimator), OuLikelihood(estimator)).infer(Sample(values, times))
    assert result.valid, result.reason
    assert result.point_fit.mu == pytest.approx(mu, rel=2e-13, abs=2e-13)
    assert result.covariance[1, 1] == pytest.approx(expected, rel=2e-5)

