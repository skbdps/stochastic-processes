"""Degeneracy and public-helper contracts from the Week 3 feedback checklist.

Fit-level unit/reference gates live in test_fit_references.py and do not call
production profiles internally. These tests cover the explicit invalid outcomes
and helper repairs, including inputs outside the original nine-cell smoke grid.
"""
from decimal import Decimal, localcontext
import numpy as np
import pytest

from ou_irregular import estimators, ou


@pytest.mark.parametrize('name', estimators.ESTIMATORS)
@pytest.mark.parametrize('x', [np.ones(10), np.r_[np.zeros(9), 1.], np.arange(10.)])
def test_degenerate_paths_have_explicit_outcomes(name, x):
    fitted = estimators.fit_estimator(name, x, np.arange(len(x), dtype=float))
    assert fitted['fit_status'] == 'degenerate', fitted
    assert not fitted['optimizer_success']
    assert not fitted['valid_for_point_summary']
    assert fitted['reason']


@pytest.mark.parametrize('x,dt', [
    ([0., 1., 2., 4., 8.], 1.), ([0., 1., 0., 1., 0.], 1.),
    (np.ones(8), 1.), ([0., 0., 0., 1.], 1.), (np.arange(5.), 1.),
    ([0., np.nan, 1.], 1.), ([0., 1.], 1.), ([[0., 1., 0.]], 1.),
    ([0., 1., 0., 0.5], 0.), ([0., 1., 0., 0.5], -1.),
    ([0., 1., 0., 0.5], np.inf), ([0., 1., 0., 0.5], True),
])
def test_ar1_helper_rejects_invalid_or_noninterior_data(x, dt):
    with pytest.raises(ValueError):
        ou.ar1_closed_form(x, dt)


def test_ar1_helper_matches_normalized_production_solution_and_preserves_units():
    t = np.arange(250) * 0.2
    x = ou.simulate_ou(1., 0., 0.5, t, np.random.default_rng(71))
    reference = ou.ar1_closed_form(x, 0.2)
    fitted = estimators.fit_estimator('exact', x, t)
    np.testing.assert_allclose(reference, [fitted[k] for k in ('theta', 'mu', 'sigma')], rtol=1e-12)
    for scale in (1e-6, 1e6):
        result = ou.ar1_closed_form(scale * (100 + x), 0.2)
        np.testing.assert_allclose([result[0], result[1] / scale - 100, result[2] / scale],
                                   reference, rtol=1e-10, atol=1e-12)


@pytest.mark.parametrize('cv', [0., 1e-8, 1e-4, 0.5, 1., 2.])
def test_pfml_raw_gamma_limit_matches_high_precision_analytic_value(cv):
    with localcontext() as ctx:
        ctx.prec = 70
        c = Decimal(str(cv))
        z = Decimal('0.5') * c * c
        ratio = (1 + z).ln() / z if z else Decimal(1)
        expected = [float(ratio), float(Decimal('0.5') * ratio.sqrt())]
    np.testing.assert_allclose(ou.pfml_limit(1., 0.5, 0.5, cv), expected, rtol=3e-16, atol=1e-16)


@pytest.mark.parametrize('args', [
    (0., .5, .5, 1.), (-1., .5, .5, 1.), (1., 0., .5, 1.),
    (1., .5, -1., 1.), (1., .5, .5, -1.), (np.nan, .5, .5, 1.),
    (1., np.inf, .5, 1.), (1., .5, np.inf, 1.), (1., .5, .5, np.nan),
    (True, .5, .5, 1.), (1e308, .5, 1e308, 1.),
    (1., .5, .5, 1e308),
])
def test_pfml_limit_rejects_invalid_and_nonrepresentable_products(args):
    with pytest.raises(ValueError):
        ou.pfml_limit(*args)


def test_irregular_exact_flat_endpoint_is_suspected_not_claimed_finite_optimum():
    # IID observations alternate strongly: high theta erases the unwanted
    # positive autoregression. Unequal gaps require profile evidence instead
    # of importing the regular-grid OLS boundary proof.
    x = np.tile([-1., 1.], 40) + np.random.default_rng(4).normal(0, .05, 80)
    t = np.cumsum(np.r_[0., np.tile([.8, 1.2], 40)[:79]])
    fitted = estimators.fit_estimator('exact', x, t)
    assert fitted['fit_status'] == 'boundary_suspected', fitted
    assert not fitted['valid_for_point_summary']
    assert np.isfinite([fitted[k] for k in ('theta', 'mu', 'sigma', 'nll')]).all()
    assert fitted['profile_logtheta_lower'] == -24
    assert fitted['profile_logtheta_upper'] >= 8
    assert fitted['profile_expansions'] >= 4


def test_scalar_optimizer_exception_is_explicit_numerical_failure(monkeypatch):
    def fail(*args, **kwargs):
        raise RuntimeError('injected scipy failure')

    monkeypatch.setattr(estimators, 'minimize_scalar', fail)
    t = np.cumsum(np.r_[0, np.random.default_rng(7).uniform(.02, .5, 249)])
    x = ou.simulate_ou(1, 0, .5, t, np.random.default_rng(16))
    fitted = estimators.fit_estimator('exact', x, t)
    assert fitted['fit_status'] == 'numerical_failure'
    assert 'injected scipy failure' in fitted['reason']
    assert not fitted['valid_for_point_summary']


def test_scalar_optimizer_exception_retains_finite_profile_candidate(monkeypatch):
    def fail(*args, **kwargs):
        raise RuntimeError('candidate retention check')

    monkeypatch.setattr(estimators, 'minimize_scalar', fail)
    t = np.cumsum(np.r_[0, np.random.default_rng(7).uniform(.02, .5, 249)])
    x = ou.simulate_ou(1, 0, .5, t, np.random.default_rng(16))
    fitted = estimators.fit_estimator('exact', x, t)
    assert fitted['fit_status'] == 'numerical_failure'
    assert np.isfinite([fitted[k] for k in ('theta', 'mu', 'sigma', 'nll')]).all()
    assert not fitted['optimizer_success'] and not fitted['valid_for_point_summary']
    assert 'candidate retention check' in fitted['error']


@pytest.mark.parametrize('state_scale,offset,time_scale', [(1., 0., 1.), (1e-6, 100., 1e3), (1e6, 100., 1e-3)])
def test_profile_best_objective_and_selection_gap_use_original_density_units(state_scale, offset, time_scale):
    t = np.cumsum(np.r_[0, np.random.default_rng(7).uniform(.02, .5, 249)])
    x = ou.simulate_ou(1, 0, .5, t, np.random.default_rng(16))
    fitted = estimators.fit_estimator('exact', state_scale * (offset + x), t / time_scale)
    assert fitted['valid_for_point_summary'], fitted
    assert abs(fitted['profile_best_nll'] + fitted['profile_selection_gap'] - fitted['nll']) <= 1e-7
