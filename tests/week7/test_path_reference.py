"""Independent analytical references for Week 7 diagnostic mathematics.

These use fixed literal data, never selected headline paths. The full Hessian
reference differentiates each Gaussian transition analytically; it does not
reuse the production likelihood or a finite-difference implementation.
"""
import math

import numpy as np
import pytest

from ou_irregular.week7.path_diagnostics import (
    central_hessian,
    full_exact_nll,
    hessian_ladder,
    profile_at,
    profile_limits,
)


VALUES = np.array([0., 1., 0., 2., 1.])
GAPS = np.array([.5, 1., 1.5, 2.])


def gls_profile(logtheta, values, gaps):
    """Independent one-column GLS solve and unprofiled Gaussian evaluation."""
    theta = np.exp(logtheta)
    phi = np.exp(-theta*gaps)
    b = -np.expm1(-theta*gaps)
    a = -np.expm1(-2*theta*gaps)/(2*theta)
    # Whiten the design and response, then solve using SVD least squares.
    design = (b/np.sqrt(a))[:, None]
    response = (values[1:]-phi*values[:-1])/np.sqrt(a)
    mu = np.linalg.lstsq(design, response, rcond=None)[0][0]
    residual = response-design[:, 0]*mu
    sigma2 = np.dot(residual, residual)/len(gaps)
    variance = sigma2*a
    raw_residual = values[1:]-(phi*values[:-1]+b*mu)
    nll = .5*np.sum(np.log(2*np.pi*variance)+raw_residual**2/variance)
    return mu, math.sqrt(sigma2), nll


def analytic_full(eta, values, gaps):
    """Full Gaussian NLL, gradient and Hessian in (log theta, mu, log sigma).

    With transition mean m, log variance r, and residual e, the Hessian is
    1/2 [r''(1-e²/v)+2 m'm'/v-2e m''/v
         +2e(r'm'+m'r')/v+(e²/v)r'r'].
    All mixed theta/mu/sigma terms are retained.
    """
    theta, mu, sigma = np.exp(eta[0]), eta[1], np.exp(eta[2])
    nll, gradient, hessian = 0., np.zeros(3), np.zeros((3, 3))
    for previous, following, gap in zip(values[:-1], values[1:], gaps):
        x = theta*gap
        phi = np.exp(-x)
        b = 1-phi
        variance = sigma**2*(1-phi**2)/(2*theta)
        residual = following-(phi*previous+b*mu)
        factor = 2*x/np.expm1(2*x)
        mean_gradient = np.array([(mu-previous)*x*phi, b, 0.])
        logvar_gradient = np.array([factor-1, 0., 2.])
        mean_hessian = np.zeros((3, 3))
        mean_hessian[0, 0] = (mu-previous)*x*phi*(1-x)
        mean_hessian[0, 1] = mean_hessian[1, 0] = x*phi
        logvar_hessian = np.zeros((3, 3))
        logvar_hessian[0, 0] = factor-factor**2*np.exp(2*x)
        u = residual**2/variance
        nll += .5*(np.log(2*np.pi*variance)+u)
        gradient += .5*(logvar_gradient*(1-u)-2*residual*mean_gradient/variance)
        hessian += .5*(logvar_hessian*(1-u)
            +2*np.outer(mean_gradient, mean_gradient)/variance
            -2*residual*mean_hessian/variance
            +2*residual*(np.outer(logvar_gradient, mean_gradient)
                         +np.outer(mean_gradient, logvar_gradient))/variance
            +u*np.outer(logvar_gradient, logvar_gradient))
    return nll, gradient, hessian


@pytest.mark.parametrize('theta', [.03, .3, 1., 3., 10.])
def test_profile_matches_independent_svd_gls(theta):
    """Both nuisance optima and the K-normalized likelihood have references."""
    expected_mu, expected_sigma, expected_nll = gls_profile(
        math.log(theta), VALUES, GAPS)
    actual = profile_at(math.log(theta), VALUES, GAPS)
    assert actual['mu'] == pytest.approx(expected_mu, rel=2e-12, abs=2e-12)
    assert actual['sigma2'] == pytest.approx(expected_sigma**2, rel=2e-12)
    assert actual['nll'] == pytest.approx(expected_nll, rel=2e-12, abs=2e-12)
    # The separately evaluated full conditional objective must agree at the
    # profiled nuisance optimum, without adding an initial stationary density.
    eta = np.array([math.log(theta), actual['mu'], .5*np.log(actual['sigma2'])],
                   dtype=np.longdouble)
    assert full_exact_nll(eta, VALUES, GAPS) == pytest.approx(
        expected_nll, rel=2e-12, abs=2e-12)


def test_limits_have_exact_rational_sufficient_statistics():
    """The theta-downarrow-zero closure retains the free Brownian drift."""
    actual = profile_limits(VALUES, GAPS)
    expected_brownian = .5*(4*(math.log(2*math.pi)+1+math.log(179/120))
                              +math.log(1.5))
    expected_iid = .5*4*(math.log(2*math.pi)+1+math.log(.5))
    assert actual['brownian']['alpha'] == pytest.approx(1/5, abs=2e-15)
    assert actual['brownian']['sigma2'] == pytest.approx(179/120, abs=2e-15)
    assert actual['brownian']['nll'] == pytest.approx(expected_brownian, abs=2e-14)
    assert actual['iid']['mu'] == pytest.approx(1., abs=2e-15)
    assert actual['iid']['variance'] == pytest.approx(.5, abs=2e-15)
    assert actual['iid']['nll'] == pytest.approx(expected_iid, abs=2e-14)
    # A zero-drift-only closure has sigma2=37/24, which is materially different.
    assert abs(actual['brownian']['sigma2']-37/24) > .04


def test_profile_approaches_both_independently_derived_limits():
    small = profile_at(-25., VALUES, GAPS)
    large = profile_at(8., VALUES, GAPS)
    expected_brownian = .5*(4*(math.log(2*math.pi)+1+math.log(179/120))
                              +math.log(1.5))
    expected_iid = .5*4*(math.log(2*math.pi)+1+math.log(.5))
    assert np.exp(np.longdouble(-25))*small['mu'] == pytest.approx(1/5, abs=2e-10)
    assert small['sigma2'] == pytest.approx(179/120, abs=2e-10)
    assert small['nll'] == pytest.approx(expected_brownian, abs=2e-10)
    assert large['mu'] == pytest.approx(1., abs=2e-14)
    assert large['sigma2']/(2*np.exp(np.longdouble(8))) == pytest.approx(.5, abs=2e-14)
    assert large['nll'] == pytest.approx(expected_iid, abs=2e-14)


def test_regular_negative_ar_slope_has_iid_boundary_infimum():
    """Here the profiled regular-gap SSR is the polynomial 2+2phi+11phi²/4.

    Its unconstrained minimum is phi=-4/11. Its derivative is positive for
    phi>=0, proving the phi=0 limit for this literal regular-gap example only.
    """
    regular_gaps = np.ones(4)
    iid_nll = .5*4*(math.log(2*math.pi)+1+math.log(.5))
    for theta in [.1, .5, 1., 2., 4.]:
        phi = math.exp(-theta)
        variance = (2+2*phi+11*phi**2/4)/4
        expected = .5*4*(math.log(2*math.pi)+1+math.log(variance))
        actual = profile_at(math.log(theta), VALUES, regular_gaps)
        assert actual['nll'] == pytest.approx(expected, abs=3e-14)
        assert actual['nll'] > iid_nll


@pytest.mark.parametrize('eta', [
    [math.log(.8), .3, math.log(1.1)],
    [math.log(1.7), -.4, math.log(.7)],
])
def test_full_extended_hessian_against_analytic_gaussian_derivatives(eta):
    values = np.array([.2, .9, -.1, 1.3, .4])
    gaps = np.array([.25, .5, 1.25, .75])
    eta = np.asarray(eta)
    expected_nll, expected_g, expected_h = analytic_full(eta, values, gaps)
    objective = lambda p: full_exact_nll(p, values, gaps)
    actual_h, actual_g = central_hessian(objective, eta, np.full(3, 1e-4),
                                          extended=True)
    assert objective(eta) == pytest.approx(expected_nll, rel=1e-13, abs=1e-13)
    np.testing.assert_allclose(actual_g, expected_g, rtol=2e-7, atol=2e-7)
    np.testing.assert_allclose(actual_h, expected_h, rtol=2e-6, atol=2e-7)
    # These nonzero mixed entries make a diagonal-only substitute fail.
    assert np.all(np.abs(expected_h[np.triu_indices(3, 1)]) > .01)


def test_profile_curvature_equals_full_hessian_schur_complement():
    values = np.array([.2, .9, -.1, 1.3, .4])
    gaps = np.array([.25, .5, 1.25, .75])
    logtheta = math.log(.8)
    mu, sigma, _ = gls_profile(logtheta, values, gaps)
    _, gradient, full = analytic_full(np.array([logtheta, mu, math.log(sigma)]),
                                     values, gaps)
    np.testing.assert_allclose(gradient[1:], 0., atol=3e-14)
    schur = full[0, 0]-full[0, 1:] @ np.linalg.solve(full[1:, 1:], full[1:, 0])
    step = np.longdouble('0.0001')
    center = profile_at(logtheta, values, gaps)['nll']
    curvature = ((profile_at(np.longdouble(logtheta)+step, values, gaps)['nll']-center)
                 +(profile_at(np.longdouble(logtheta)-step, values, gaps)['nll']-center))/step**2
    assert curvature == pytest.approx(schur, rel=2e-6, abs=2e-7)
    assert abs(curvature-full[0, 0]) > .01


def test_full_nll_observed_coordinate_jacobian():
    """Time scales parameters; the likelihood data Jacobian is K log(state_scale)."""
    normalized_values = np.array([.2, .9, -.1, 1.3, .4])
    normalized_gaps = np.array([.25, .5, 1.25, .75])
    state_scale, time_scale, offset = 3.4, 2.7, -4.2
    theta_x, mu_x, sigma_x = .8, -.3, 1.1
    eta_x = [math.log(theta_x), mu_x, math.log(sigma_x)]
    eta_u = [math.log(theta_x*time_scale), (mu_x-offset)/state_scale,
             math.log(sigma_x*math.sqrt(time_scale)/state_scale)]
    physical = full_exact_nll(eta_x, offset+state_scale*normalized_values,
                             time_scale*normalized_gaps)
    normalized = full_exact_nll(eta_u, normalized_values, normalized_gaps)
    assert physical-normalized == pytest.approx(4*math.log(state_scale), abs=3e-13)


@pytest.mark.parametrize('extended', [False, True])
def test_whitened_ladder_uses_generalized_eigenvalues_and_known_quartic_hessians(extended):
    """A nondiagonal fine Hessian discriminates the whitening orientation.

    At zero the stencil is exactly Q+diag(2*c*h²). These quartics yield a
    visible nonproportional coarse/fine difference despite the small steps.
    """
    from scipy.linalg import eigvalsh
    rotation, _ = np.linalg.qr(np.array([[1., 2., -1.], [-2., 1., 1.], [1., 3., 2.]]))
    quadratic = rotation @ np.diag([.01, 1., 10.]) @ rotation.T
    quartic = np.array([1e4, 4e5, 2e6])
    def objective(point):
        p = np.asarray(point, dtype=np.longdouble if extended else float)
        return .5*p @ quadratic @ p+np.dot(quartic, p**4)
    ladder = hessian_ladder(objective, np.zeros(3), extended=extended)
    assert [entry['relative_step'] for entry in ladder['entries']] == [
        .001, .0005, .00025, .000125]
    for entry in ladder['entries']:
        step = entry['relative_step']
        expected = quadratic+np.diag(2*quartic*step**2)
        np.testing.assert_allclose(entry['hessian'], expected, rtol=2e-12, atol=2e-12)
        assert entry['positive_definite'] is True
        assert entry['standardized_score'] == pytest.approx(0., abs=1e-15)
    for index, adjacent in enumerate(ladder['adjacent']):
        coarse = np.asarray(ladder['entries'][index]['hessian'])
        fine = np.asarray(ladder['entries'][index+1]['hessian'])
        expected_white = np.max(np.abs(eigvalsh(coarse-fine, fine)))
        expected_relative = np.linalg.norm(coarse-fine)/np.linalg.norm(fine)
        assert adjacent['whitened_spectral'] == pytest.approx(expected_white, rel=2e-12)
        assert adjacent['relative_frobenius'] == pytest.approx(expected_relative, rel=2e-12)


def test_nonpositive_and_failed_steps_remain_explicitly_unavailable():
    indefinite = lambda p: -.5*p[0]**2+.5*p[1]**2+.5*p[2]**2
    ladder = hessian_ladder(indefinite, np.zeros(3), extended=True)
    assert len(ladder['entries']) == 4
    for entry in ladder['entries']:
        assert entry['positive_definite'] is False
        assert entry['status'] == 'nonpositive_hessian'
        assert entry['standardized_score'] is None
        assert entry['theta_hessian_diagonal'] == pytest.approx(-1., abs=1e-13)
        assert entry['theta_schur_complement_at_stored_point'] == pytest.approx(-1., abs=1e-13)
        np.testing.assert_allclose(entry['eigenvalues'], [-1., 1., 1.], atol=1e-13)
    for adjacent in ladder['adjacent']:
        assert adjacent['whitened_spectral'] is None
        assert adjacent['status'] == 'fine_hessian_not_positive_definite'
    def undefined(_):
        raise ValueError('literal failed objective')
    failed = hessian_ladder(undefined, np.zeros(3), extended=True)
    assert len(failed['entries']) == 4
    assert all(entry['status'] == 'derivative_evaluation_failure'
               and entry['hessian'] is None and entry['gradient'] is None
               and entry['standardized_score'] is None for entry in failed['entries'])
    assert all(adjacent['status'] == 'derivatives_unavailable'
               and adjacent['whitened_spectral'] is None for adjacent in failed['adjacent'])
    assert failed['replacement_interval_constructed'] is False


def test_stored_point_schur_is_distinguished_from_profile_curvature():
    values = np.array([.2, .9, -.1, 1.3, .4])
    gaps = np.array([.25, .5, 1.25, .75])
    eta = np.array([math.log(.8), .3, math.log(1.1)])
    _, gradient, reference = analytic_full(eta, values, gaps)
    assert np.linalg.norm(gradient[1:]) > 1
    algebraic = (reference[0, 0]
                 -reference[0, 1:] @ np.linalg.solve(reference[1:, 1:], reference[1:, 0]))
    ladder = hessian_ladder(lambda p: full_exact_nll(p, values, gaps), eta, extended=True)
    finest = ladder['entries'][-1]
    assert finest['theta_hessian_diagonal'] == pytest.approx(reference[0, 0], abs=2e-7)
    assert finest['theta_schur_complement_at_stored_point'] == pytest.approx(algebraic, abs=2e-7)
    assert finest['schur_status'] == 'evaluated_at_stored_point'
    assert finest['schur_profile_interpretation'] == 'only_at_nuisance_stationarity_on_same_representation'
    mu, sigma, _ = gls_profile(eta[0], values, gaps)
    _, _, profiled = analytic_full(np.array([eta[0], mu, math.log(sigma)]), values, gaps)
    profile_curvature = (profiled[0, 0]
                         -profiled[0, 1:] @ np.linalg.solve(profiled[1:, 1:], profiled[1:, 0]))
    assert abs(algebraic-profile_curvature) > .3
    singular = hessian_ladder(lambda p: .5*p[0]**2+.5*p[2]**2,
                             np.zeros(3), extended=True)
    assert all(entry['theta_schur_complement_at_stored_point'] is None
               and entry['schur_status'] == 'singular_nuisance_block'
               for entry in singular['entries'])


def test_extended_stencil_retains_sub_float64_objective_differences():
    if np.finfo(np.longdouble).eps >= np.finfo(float).eps:
        pytest.skip('platform longdouble provides no additional mantissa precision')
    matrix = np.array([[1.3, .2, -.1], [.2, 2.1, .3], [-.1, .3, 3.7]],
                      dtype=np.longdouble)
    def objective(p):
        p = np.asarray(p, dtype=np.longdouble)
        return np.longdouble('1000000')+np.longdouble('.5')*p @ matrix @ p
    hessian, gradient = central_hessian(objective, np.zeros(3), np.full(3, .0005),
                                          extended=True)
    assert hessian.dtype == np.dtype(np.longdouble)
    assert gradient.dtype == np.dtype(np.longdouble)
    np.testing.assert_allclose(hessian, matrix, rtol=0, atol=2e-6)
    np.testing.assert_allclose(gradient, 0., atol=0.)


def test_fixed_tail_caps_are_reported_without_claiming_target_attainment():
    from ou_irregular.week7.path_diagnostics import _profile_bounds
    ordinary = _profile_bounds(np.ones(4))
    assert ordinary['lower'] == -30.
    assert ordinary['upper'] == 12.
    assert ordinary['lower_tail_reached'] is True
    assert ordinary['upper_tail_reached'] is True
    capped = _profile_bounds(np.array([1e-20, 1e20], dtype=np.longdouble))
    assert capped['lower'] == -40.
    assert capped['upper'] == 40.
    assert capped['lower_cap_limited'] is True
    assert capped['upper_cap_limited'] is True
    assert capped['lower_tail_reached'] is False
    assert capped['upper_tail_reached'] is False


def test_analysis_preserves_registered_ladder_and_distinct_profile_representations(monkeypatch):
    """Composition test on literal data; it does not claim these are fitted optima."""
    import ou_irregular.week7.path_diagnostics as module
    from ou_irregular.inference.adapters import OuLikelihood
    from ou_irregular.inference.contracts import Sample
    from ou_irregular.inference.coordinates import ObservationCoordinates
    sample = Sample(np.array([1.2, 2.1, .9, 3.5, .4]),
                    1e6+np.array([0., .1, .8, 1.3, 2.8]))
    coordinates = ObservationCoordinates.from_sample(sample)
    eta_x = np.array([math.log(.8), .3, math.log(1.1)])
    rows = {method: {'estimator': method, 'point_valid': True,
                    'point_status': 'valid_interior', 'interval_valid': False,
                    'interval_status': 'synthetic_not_inferred',
                    'theta': .8, 'mu': .3, 'sigma': 1.1,
                    'nll': OuLikelihood(method).nll(eta_x, sample)}
            for method in ('exact', 'pfml', 'euler')}
    original_full = module.full_exact_nll
    captured_full, captured_profile = [], []
    def recording_full(eta, values, gaps):
        captured_full.append((np.array(eta), np.array(values), np.array(gaps)))
        return original_full(eta, values, gaps)
    def minimal_ladder(objective, point, extended=False):
        # The stencil itself is independently tested above. This evaluates only
        # the composed objective to isolate which observations it receives.
        return {'synthetic_composition_probe': float(objective(np.asarray(point)))}
    def recording_profile(values, gaps, eta):
        captured_profile.append((np.array(eta), np.array(values), np.array(gaps)))
        return {'status': 'synthetic_composition_probe'}
    monkeypatch.setattr(module, 'full_exact_nll', recording_full)
    monkeypatch.setattr(module, 'hessian_ladder', minimal_ladder)
    monkeypatch.setattr(module, '_profile_search', recording_profile)
    result = module.analyze_sample(sample, rows)
    expected_eta = coordinates.to_working(np.array([.8, .3, 1.1]))
    registered_values = coordinates.normalized_sample.values.astype(np.longdouble)
    registered_gaps = np.diff(coordinates.normalized_sample.times).astype(np.longdouble)
    profile_values = ((sample.values.astype(np.longdouble)-np.longdouble(coordinates.offset))
                      /np.longdouble(coordinates.state_scale))
    profile_gaps = np.diff(sample.times.astype(np.longdouble))/np.longdouble(coordinates.time_scale)
    assert np.any(profile_gaps != registered_gaps)
    np.testing.assert_array_equal(captured_full[0][0], expected_eta)
    np.testing.assert_array_equal(captured_full[0][1], registered_values)
    np.testing.assert_array_equal(captured_full[0][2], registered_gaps)
    np.testing.assert_array_equal(captured_profile[0][0], expected_eta)
    np.testing.assert_array_equal(captured_profile[0][1], profile_values)
    np.testing.assert_array_equal(captured_profile[0][2], profile_gaps)
    assert result['coordinates']['state_offset'] == coordinates.offset
    assert result['coordinates']['source'] == 'observed_data_only'
    assert result['representation_differences']['maximum_absolute_gap_difference'] > 0
    assert result['headline_results_replaced'] is False
