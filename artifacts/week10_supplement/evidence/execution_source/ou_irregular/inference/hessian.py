"""Independent central differences of a conditional negative log-likelihood.

No inverse-Hessian approximation is borrowed from an optimizer. Both h and
h/2 are evaluated from the stated objective, using h_j=1e-3*max(1,abs(eta_j)).
Their relative Frobenius discrepancy and covariance-whitened discrepancy
diagnose finite-difference instability; the finer Hessian/gradient are retained.
Whitening tests weak-curvature directions that a global matrix norm can hide.
Additional curvature and stationarity gates are
applied by WaldInference, allowing this provider to be replaced independently.
"""
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Callable, Mapping, Protocol

import numpy as np


def _readonly(value):
    array = np.array(value, dtype=float, copy=True)
    array.setflags(write=False)
    return array


def _immutable(value):
    """Snapshot diagnostic containers without depending on an inference method."""
    if isinstance(value, Mapping):
        return MappingProxyType({key: _immutable(item) for key, item in value.items()})
    if isinstance(value, np.ndarray):
        result = np.array(value, copy=True)
        result.setflags(write=False)
        return result
    if isinstance(value, (tuple, list)):
        return tuple(_immutable(item) for item in value)
    if isinstance(value, set):
        return frozenset(_immutable(item) for item in value)
    return value


@dataclass(frozen=True)
class HessianEstimate:
    hessian: np.ndarray
    gradient: np.ndarray
    valid: bool
    status: str
    reason: str
    relative_difference: float = float("nan")
    diagnostics: Mapping = field(default_factory=dict)

    def __post_init__(self):
        object.__setattr__(self, "hessian", _readonly(self.hessian))
        object.__setattr__(self, "gradient", _readonly(self.gradient))
        object.__setattr__(self, "diagnostics", _immutable(self.diagnostics))


class HessianProvider(Protocol):
    def evaluate(self, objective: Callable[[np.ndarray], float],
                 working_parameters: np.ndarray) -> HessianEstimate: ...


@dataclass(frozen=True)
class FiniteDifferenceHessian:
    relative_step: float = 1e-3
    max_relative_difference: float = 0.005
    max_whitened_difference: float = 0.01

    def __post_init__(self):
        for name in ("relative_step", "max_relative_difference", "max_whitened_difference"):
            value = getattr(self, name)
            if isinstance(value, (bool, np.bool_)) or not np.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be positive and finite")

    @staticmethod
    def _central(objective, point, steps, center):
        dimension = point.size
        hessian = np.empty((dimension, dimension))
        gradient = np.empty(dimension)
        for i in range(dimension):
            ei = np.zeros(dimension)
            ei[i] = steps[i]
            plus, minus = float(objective(point + ei)), float(objective(point - ei))
            gradient[i] = (plus - minus) / (2 * steps[i])
            # Centering the differences separately reduces additive-offset
            # cancellation, but cannot remove objective evaluation rounding.
            hessian[i, i] = ((plus - center) + (minus - center)) / steps[i] ** 2
            for j in range(i):
                ej = np.zeros(dimension)
                ej[j] = steps[j]
                pp = float(objective(point + ei + ej))
                pm = float(objective(point + ei - ej))
                mp = float(objective(point - ei + ej))
                mm = float(objective(point - ei - ej))
                value = ((pp - pm) - (mp - mm)) / (4 * steps[i] * steps[j])
                hessian[i, j] = hessian[j, i] = value
        return hessian, gradient

    def evaluate(self, objective, working_parameters):
        point = np.asarray(working_parameters, dtype=float)
        if point.shape != (3,) or not np.isfinite(point).all():
            return HessianEstimate(np.full((3, 3), np.nan), np.full(3, np.nan),
                                   False, "invalid_coordinates", "three finite working coordinates required")
        steps = self.relative_step * np.maximum(1.0, np.abs(point))
        diagnostics = {"relative_step": self.relative_step,
                       "max_relative_difference": self.max_relative_difference,
                       "max_whitened_difference": self.max_whitened_difference,
                       "steps": tuple(float(x) for x in steps),
                       "selected_step_fraction": 0.5}
        try:
            center = float(objective(point))
            with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
                coarse, _ = self._central(objective, point, steps, center)
                fine, gradient = self._central(objective, point, steps / 2, center)
        except (ValueError, FloatingPointError, OverflowError, ArithmeticError) as exc:
            return HessianEstimate(np.full((3, 3), np.nan), np.full(3, np.nan),
                                   False, "nonfinite_derivatives", str(exc), diagnostics=diagnostics)
        if not (np.isfinite(center) and np.isfinite(coarse).all()
                and np.isfinite(fine).all() and np.isfinite(gradient).all()):
            return HessianEstimate(fine, gradient, False, "nonfinite_derivatives",
                                   "objective or finite differences are nonfinite", diagnostics=diagnostics)
        norm = np.linalg.norm(fine, ord="fro")
        discrepancy = np.linalg.norm(fine - coarse, ord="fro")
        relative = float(discrepancy / norm) if norm else (0.0 if discrepancy == 0 else float("inf"))
        if not np.isfinite(relative) or relative > self.max_relative_difference:
            return HessianEstimate(fine, gradient, False, "unstable_hessian",
                                   "full and half Hessian steps disagree in Frobenius norm",
                                   relative, diagnostics)
        try:
            # For H_fine=L L', normalize the discrepancy on both sides by L.
            # ||L^-1(H_coarse-H_fine)L^-T||_2 bounds relative curvature error
            # in EVERY direction, rather than only the largest eigenvalues.
            cholesky = np.linalg.cholesky(fine)
            left = np.linalg.solve(cholesky, coarse - fine)
            whitened = np.linalg.solve(cholesky, left.T).T
            whitened = (whitened + whitened.T) / 2
            discrepancy_whitened = float(np.max(np.abs(np.linalg.eigvalsh(whitened))))
        except np.linalg.LinAlgError:
            return HessianEstimate(fine, gradient, False, "nonpositive_hessian",
                                   "fine-step Hessian is not positive definite", relative, diagnostics)
        diagnostics["whitened_relative_difference"] = discrepancy_whitened
        valid = bool(np.isfinite(discrepancy_whitened)
                     and discrepancy_whitened <= self.max_whitened_difference)
        return HessianEstimate(fine, gradient, valid,
                               "stable" if valid else "unstable_hessian",
                               "full and half steps agree" if valid else "full and half Hessians disagree in a curvature direction",
                               relative, diagnostics)
