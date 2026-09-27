"""Observed-Hessian Wald inference composed from small replaceable interfaces.

This module knows no Week 3 implementation. It accepts a point estimator and
conditional likelihood, evaluates curvature in observed-data normalized log
coordinates, and maps covariance to natural units using the delta method.
The reported interval is theta_hat +/- z*SE (and similarly for mu/sigma): a
natural-scale symmetric Wald interval. Negative lower limits are deliberately
retained. Truncation or log-scale intervals would define a different procedure.

For PFML/Euler this inverse Hessian is model-based, NOT a sandwich covariance;
coverage under their misspecification must be measured rather than assumed.
"""
from dataclasses import dataclass, field, replace
from typing import Mapping

import numpy as np
from scipy.stats import norm

from .contracts import ConditionalLikelihood, PointEstimator, PointFit, Sample
from .coordinates import ObservationCoordinates
from .hessian import FiniteDifferenceHessian, HessianProvider, _immutable, _readonly


@dataclass(frozen=True)
class WaldResult:
    point_fit: PointFit
    valid: bool
    status: str
    reason: str
    lower: np.ndarray
    upper: np.ndarray
    covariance: np.ndarray
    working_covariance: np.ndarray
    diagnostics: Mapping = field(default_factory=dict)

    def __post_init__(self):
        for name in ("lower", "upper", "covariance", "working_covariance"):
            object.__setattr__(self, name, _readonly(getattr(self, name)))
        object.__setattr__(self, "point_fit", replace(
            self.point_fit, diagnostics=_immutable(self.point_fit.diagnostics)))
        object.__setattr__(self, "diagnostics", _immutable(self.diagnostics))


@dataclass(frozen=True)
class WaldInference:
    estimator: PointEstimator
    likelihood: ConditionalLikelihood
    hessian: HessianProvider = field(default_factory=FiniteDifferenceHessian)
    level: float = 0.95
    max_condition_number: float = 1e10
    max_standardized_score: float = 0.01

    def __post_init__(self):
        if self.estimator.name != self.likelihood.name:
            raise ValueError("estimator and likelihood names must agree")
        if (isinstance(self.level, (bool, np.bool_))
                or not np.isfinite(self.level) or not 0 < self.level < 1):
            raise ValueError("level must lie strictly between zero and one")
        if (isinstance(self.max_condition_number, (bool, np.bool_))
                or not np.isfinite(self.max_condition_number) or self.max_condition_number < 1):
            raise ValueError("max_condition_number must be finite and at least one")
        if (isinstance(self.max_standardized_score, (bool, np.bool_))
                or not np.isfinite(self.max_standardized_score) or self.max_standardized_score <= 0):
            raise ValueError("max_standardized_score must be positive and finite")

    @staticmethod
    def _invalid(point_fit, status, reason, diagnostics):
        # Keep the original point candidate and never replace a failed CI by
        # zero-width endpoints or silently remove it from a coverage denominator.
        return WaldResult(point_fit, False, status, reason, np.full(3, np.nan),
                          np.full(3, np.nan), np.full((3, 3), np.nan),
                          np.full((3, 3), np.nan), diagnostics)

    def infer(self, sample: Sample, point_fit: PointFit = None):
        fit = self.estimator.fit(sample) if point_fit is None else point_fit
        diagnostics = {"confidence_level": self.level,
                       "interval_scale": "natural_symmetric",
                       "conditions_on_initial_value": True,
                       "covariance_type": "inverse_observed_hessian",
                       "max_condition_number": self.max_condition_number,
                       "max_standardized_score": self.max_standardized_score}
        if fit.estimator != self.estimator.name:
            return self._invalid(fit, "estimator_mismatch", "point fit belongs to a different estimator", diagnostics)
        parameters = fit.parameters
        if (not isinstance(fit.valid, (bool, np.bool_)) or not fit.valid
                or not isinstance(fit.optimizer_success, (bool, np.bool_)) or not fit.optimizer_success
                or fit.status != "valid_interior"
                or not np.isfinite(parameters).all() or not np.isfinite(fit.nll)
                or parameters[0] <= 0 or parameters[2] <= 0):
            return self._invalid(fit, "invalid_point_fit", f"point fit is not an admissible interior optimum: {fit.reason}", diagnostics)
        try:
            coordinates = ObservationCoordinates.from_sample(sample)
            working = coordinates.to_working(parameters)
        except (ValueError, FloatingPointError, OverflowError) as exc:
            return self._invalid(fit, "invalid_coordinates", str(exc), diagnostics)
        diagnostics.update(state_offset=coordinates.offset, state_scale=coordinates.state_scale,
                           time_scale=coordinates.time_scale,
                           working_parameters=tuple(float(x) for x in working))
        objective = lambda eta: self.likelihood.nll(eta, coordinates.normalized_sample)
        estimate = self.hessian.evaluate(objective, working)
        diagnostics.update(hessian_status=estimate.status,
                           hessian_relative_difference=estimate.relative_difference,
                           hessian_diagnostics=dict(estimate.diagnostics))
        if not estimate.valid:
            return self._invalid(fit, estimate.status, estimate.reason, diagnostics)
        hessian, gradient = estimate.hessian, estimate.gradient
        if (hessian.shape != (3, 3) or gradient.shape != (3,)
                or not np.isfinite(hessian).all() or not np.isfinite(gradient).all()
                or not np.allclose(hessian, hessian.T, rtol=1e-10, atol=1e-12)):
            return self._invalid(fit, "invalid_hessian", "provider returned nonfinite, asymmetric, or incorrectly shaped derivatives", diagnostics)
        try:
            eigenvalues = np.linalg.eigvalsh(hessian)
            diagnostics["hessian_eigenvalues"] = tuple(float(x) for x in eigenvalues)
            if eigenvalues[0] <= 0:
                return self._invalid(fit, "nonpositive_hessian", "observed Hessian is not positive definite", diagnostics)
            condition = float(eigenvalues[-1] / eigenvalues[0])
            diagnostics["hessian_condition_number"] = condition
            if not np.isfinite(condition) or condition > self.max_condition_number:
                return self._invalid(fit, "ill_conditioned_hessian", "observed Hessian exceeds condition-number limit", diagnostics)
            # Newton decrement is a dimensionless local stationarity measure.
            # A small raw gradient alone depends on arbitrary parameter units.
            score = float(np.sqrt(max(0.0, gradient @ np.linalg.solve(hessian, gradient))))
            diagnostics["standardized_score"] = score
            if not np.isfinite(score) or score > self.max_standardized_score:
                return self._invalid(fit, "nonstationary_point_fit", "standardized score exceeds tolerance", diagnostics)
            working_covariance = np.linalg.solve(hessian, np.eye(3))
            jacobian = coordinates.natural_jacobian(working)
            covariance = jacobian @ working_covariance @ jacobian.T
            covariance = (covariance + covariance.T) / 2
            if not np.isfinite(covariance).all() or np.any(np.diag(covariance) <= 0):
                return self._invalid(fit, "invalid_covariance", "natural covariance is not finite and positive on its diagonal", diagnostics)
            standard_errors = np.sqrt(np.diag(covariance))
            critical = float(norm.ppf((1 + self.level) / 2))
            lower, upper = parameters - critical * standard_errors, parameters + critical * standard_errors
            if not np.isfinite(lower).all() or not np.isfinite(upper).all():
                return self._invalid(fit, "invalid_interval", "Wald endpoints are not finite", diagnostics)
        except (ValueError, np.linalg.LinAlgError, FloatingPointError, OverflowError) as exc:
            return self._invalid(fit, "numerical_failure", str(exc), diagnostics)
        diagnostics.update(standard_errors=tuple(float(x) for x in standard_errors), critical_value=critical)
        return WaldResult(fit, True, "valid", "stable positive-definite observed Hessian at an interior optimum",
                          lower, upper, covariance, working_covariance, diagnostics)
