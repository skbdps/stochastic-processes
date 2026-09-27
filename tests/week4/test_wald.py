"""Wald covariance, analytic curvature and unit-transform integration checks."""
from dataclasses import dataclass, replace

import numpy as np
import pytest
from scipy.stats import norm

from ou_irregular.inference.adapters import OuLikelihood, Week3Estimator
from ou_irregular.inference.contracts import PointFit, Sample
from ou_irregular.inference.coordinates import ObservationCoordinates
from ou_irregular.inference.hessian import FiniteDifferenceHessian, HessianEstimate
from ou_irregular.inference.wald import WaldInference
from ou_irregular.ou import simulate_ou
from ou_irregular.spacing import generate_times


def observed_sample(cv=0.0):
    times = generate_times(500, 0.15, cv, np.random.default_rng(840)).times
    return Sample(simulate_ou(1.1, -0.4, 0.7, times, np.random.default_rng(841)), times)


@dataclass
class FixedEstimator:
    point_fit: PointFit
    name: str = "quadratic"
    calls: int = 0

    def fit(self, sample):
        self.calls += 1
        return self.point_fit


@dataclass
class QuadraticLikelihood:
    center: np.ndarray
    precision: np.ndarray
    name: str = "quadratic"

    def nll(self, eta, sample):
        return float((eta - self.center) @ self.precision @ (eta - self.center) / 2)


@dataclass
class FixedHessian:
    matrix: np.ndarray
    gradient: np.ndarray

    def evaluate(self, objective, working_parameters):
        return HessianEstimate(self.matrix, self.gradient, True, "stable", "injected analytic result", 0.0)


def quadratic_problem(precision=None):
    sample = Sample([-1, 0.2, 0.3, 1], [0, 0.2, 0.4, 0.6])
    fit = PointFit("quadratic", 0.5, 0.4, 0.7, 0, True, True, "valid_interior", "analytic solution")
    coordinates = ObservationCoordinates.from_sample(sample)
    eta = coordinates.to_working(fit.parameters)
    precision = np.diag([0.1, 3.0, 4.0]) if precision is None else precision
    estimator = FixedEstimator(fit)
    return sample, fit, coordinates, estimator, QuadraticLikelihood(eta, precision)


def test_injected_provider_delta_covariance_and_unclipped_symmetric_interval():
    sample, fit, coordinates, estimator, likelihood = quadratic_problem()
    result = WaldInference(estimator, likelihood,
                           hessian=FixedHessian(likelihood.precision, np.zeros(3))).infer(sample)
    jacobian = np.diag([fit.theta, coordinates.state_scale, fit.sigma])
    expected = jacobian @ np.linalg.inv(likelihood.precision) @ jacobian.T
    assert result.valid, result.reason
    np.testing.assert_allclose(result.working_covariance, np.linalg.inv(likelihood.precision))
    np.testing.assert_allclose(result.covariance, expected)
    np.testing.assert_allclose(result.lower, fit.parameters - norm.ppf(0.975) * np.sqrt(np.diag(expected)))
    np.testing.assert_allclose((result.lower + result.upper) / 2, fit.parameters)
    assert result.lower[0] < 0  # no artificial truncation to the OU domain
    np.testing.assert_equal(result.point_fit.parameters, fit.parameters)
    assert result.point_fit.nll == fit.nll
    assert estimator.calls == 1


def test_supplied_point_fit_does_not_refit():
    sample, fit, _, estimator, likelihood = quadratic_problem()
    result = WaldInference(estimator, likelihood).infer(sample, point_fit=fit)
    assert result.valid, result.reason
    assert estimator.calls == 0


def test_regular_ou_hessian_matches_independent_ar1_observed_information():
    sample = observed_sample()
    fit = Week3Estimator("exact").fit(sample)
    coordinates = ObservationCoordinates.from_sample(sample)
    normalized = coordinates.normalized_sample
    eta = coordinates.to_working(fit.parameters)
    theta, mu, sigma = np.exp(eta[0]), eta[1], np.exp(eta[2])
    gap = np.diff(normalized.times).mean()
    q = theta * gap
    phi = np.exp(-q)
    variance = sigma ** 2 * (-np.expm1(-2 * q)) / (2 * theta)
    predictors = normalized.values[:-1]
    count = len(predictors)
    # Independently use AR(1) coordinates (phi, intercept, log innovation SD).
    # At OLS/MLE residuals are orthogonal to 1 and X, so cross derivatives
    # with log innovation SD vanish and its observed information is 2K.
    ar_information = np.array([[np.dot(predictors, predictors) / variance, predictors.sum() / variance, 0],
                               [predictors.sum() / variance, count / variance, 0],
                               [0, 0, 2 * count]])
    jacobian = np.array([[-q * phi, 0, 0],
                         [mu * q * phi, 1 - phi, 0],
                         [q / np.expm1(2 * q) - 0.5, 0, 1]])
    analytic_information = jacobian.T @ ar_information @ jacobian
    estimate = FiniteDifferenceHessian().evaluate(lambda p: OuLikelihood("exact").nll(p, normalized), eta)
    assert estimate.valid, estimate.reason
    np.testing.assert_allclose(estimate.hessian, analytic_information, rtol=1e-5, atol=2e-4)
    result = WaldInference(Week3Estimator("exact"), OuLikelihood("exact")).infer(sample, fit)
    assert result.valid, result.reason
    np.testing.assert_allclose(result.working_covariance, np.linalg.inv(analytic_information), rtol=1e-5, atol=2e-7)


@pytest.mark.parametrize("cv", [0, 2])
@pytest.mark.parametrize("name", ["exact", "pfml", "euler"])
def test_wald_covariance_and_intervals_transform_with_state_and_time_units(cv, name):
    sample = observed_sample(cv)
    inference = WaldInference(Week3Estimator(name), OuLikelihood(name))
    original = inference.infer(sample)
    assert original.valid, (name, cv, original.status, original.reason)
    # Simultaneously change time units and affine state units. This tests all
    # delta-method factors, including the diffusion's square-root time factor.
    shift, scale, speed = 3.7, 2.4, 3.0
    transformed = inference.infer(Sample(shift + scale * sample.values, sample.times / speed))
    assert transformed.valid, transformed.reason
    factors = np.array([speed, scale, scale * np.sqrt(speed)])
    translation = np.array([0, shift, 0])
    np.testing.assert_allclose(transformed.covariance, factors[:, None] * original.covariance * factors[None, :], rtol=8e-5, atol=1e-8)
    np.testing.assert_allclose(transformed.working_covariance, original.working_covariance, rtol=8e-5, atol=1e-8)
    np.testing.assert_allclose(transformed.lower, translation + factors * original.lower, rtol=8e-5, atol=1e-7)
    np.testing.assert_allclose(transformed.upper, translation + factors * original.upper, rtol=8e-5, atol=1e-7)


@pytest.mark.parametrize("changes", [{"valid": False}, {"optimizer_success": False},
                                        {"status": "boundary"}, {"status": "boundary_suspected"},
                                        {"status": "degenerate"}, {"theta": 0}, {"sigma": -1},
                                        {"mu": np.nan}, {"nll": np.inf}, {"status": "unknown"},
                                        {"valid": "True"}, {"optimizer_success": 1}])
def test_invalid_point_candidate_is_preserved_and_intervals_missing(changes):
    sample, fit, _, estimator, likelihood = quadratic_problem()
    candidate = replace(fit, **changes)
    result = WaldInference(estimator, likelihood).infer(sample, candidate)
    assert not result.valid
    assert result.status == "invalid_point_fit"
    np.testing.assert_equal(result.point_fit.parameters, candidate.parameters)
    assert result.point_fit.nll == candidate.nll
    assert np.isnan(result.lower).all()
    assert np.isnan(result.upper).all()
    assert np.isnan(result.covariance).all()


@pytest.mark.parametrize("matrix,gradient,status", [
    (np.diag([0, 1, 1]), np.zeros(3), "nonpositive_hessian"),
    (np.diag([-1, 1, 1]), np.zeros(3), "nonpositive_hessian"),
    (np.diag([1e-12, 1, 1]), np.zeros(3), "ill_conditioned_hessian"),
    (np.eye(3), np.array([0.02, 0, 0]), "nonstationary_point_fit"),
    (np.full((3, 3), np.nan), np.zeros(3), "invalid_hessian"),
    (np.array([[1, 0.1, 0], [0, 1, 0], [0, 0, 1]]), np.zeros(3), "invalid_hessian"),
    (np.eye(2), np.zeros(3), "invalid_hessian"),
])
def test_injected_provider_cannot_bypass_wald_validity_gates(matrix, gradient, status):
    sample, fit, _, estimator, likelihood = quadratic_problem()
    result = WaldInference(estimator, likelihood, hessian=FixedHessian(matrix, gradient)).infer(sample, fit)
    assert not result.valid
    assert result.status == status


def test_constant_observed_path_cannot_define_coordinates():
    sample, fit, _, estimator, likelihood = quadratic_problem()
    result = WaldInference(estimator, likelihood).infer(Sample(np.ones(4), sample.times), fit)
    assert not result.valid
    assert result.status == "invalid_coordinates"


def test_point_fit_from_another_estimator_is_rejected():
    sample, fit, _, estimator, likelihood = quadratic_problem()
    result = WaldInference(estimator, likelihood).infer(sample, replace(fit, estimator="other"))
    assert not result.valid
    assert result.status == "estimator_mismatch"


@pytest.mark.parametrize("field,value", [("level", 0), ("level", 1), ("level", np.nan),
                                         ("level", True), ("max_condition_number", 0.5),
                                         ("max_standardized_score", 0)])
def test_wald_settings_are_validated(field, value):
    _, _, _, estimator, likelihood = quadratic_problem()
    with pytest.raises(ValueError):
        WaldInference(estimator, likelihood, **{field: value})


def test_mismatched_estimator_likelihood_names_are_rejected():
    with pytest.raises(ValueError):
        WaldInference(Week3Estimator("exact"), OuLikelihood("pfml"))


def test_wald_result_arrays_and_diagnostics_are_read_only():
    sample, _, _, estimator, likelihood = quadratic_problem()
    result = WaldInference(estimator, likelihood).infer(sample)
    for array in (result.lower, result.upper, result.covariance, result.working_covariance):
        with pytest.raises(ValueError):
            array.flat[0] = 0
    with pytest.raises(TypeError):
        result.diagnostics["confidence_level"] = 0
    with pytest.raises(TypeError):
        result.diagnostics["hessian_diagnostics"]["relative_step"] = 0


def test_result_snapshots_mutable_point_diagnostics():
    sample, fit, _, estimator, likelihood = quadratic_problem()
    fit.diagnostics["external_list"] = [1, 2]
    result = WaldInference(estimator, likelihood).infer(sample, fit)
    fit.diagnostics["external_list"].append(3)
    assert result.point_fit.diagnostics["external_list"] == (1, 2)
