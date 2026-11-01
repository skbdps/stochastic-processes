"""Finite-design OU covariance benchmarks for the Week 10 supplement.

These functions calculate diagnostics from the known simulation law. They do
not change fitted estimators, their Hessians, or their reported intervals.
The fixed-slope mean is linear in observations; the fitted PFML mean is not.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.linalg import expm
from scipy.special import ndtr
from scipy.stats import norm

from ..week8.design import E4Cell
from ..week8.source import floor_value, phase_type


def _positive(value, name):
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be finite and positive")
    value = float(value)
    if not math.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be finite and positive")
    return value


def _integer(value, name, minimum):
    if (isinstance(value, (bool, np.bool_))
            or not isinstance(value, (int, np.integer)) or value < minimum):
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return int(value)


def _rho(value):
    value = _positive(value, "rho")
    if value >= 1:
        raise ValueError("rho must be strictly between zero and one")
    return value


def laplace_kernel(cell: E4Cell, s: float) -> np.ndarray:
    """Joint ending-phase transform E[exp(-s*used_gap); next phase].

    Both stochastic mechanisms use the same underlying arrival/switch kernel.
    The IID mechanism resets its starting phase before every gap; callers must
    not take powers of this matrix for IID multi-gap correlations.
    """
    if isinstance(s, (bool, np.bool_)) or not math.isfinite(float(s)) or s < 0:
        raise ValueError("s must be finite and nonnegative")
    if cell.mechanism == "regular":
        raise ValueError("regular schedules have no phase kernel")
    ph = phase_type(cell)
    f = floor_value(cell)
    tail = np.linalg.solve(ph.A + float(s) * np.eye(2), ph.Lambda)
    if f == 0:
        return tail
    survival = expm(-ph.A * f)
    atom = np.linalg.solve(ph.A, np.eye(2) - survival) @ ph.Lambda
    return math.exp(-float(s) * f) * (atom + survival @ tail)


def population_rho(cell: E4Cell, theta=None) -> float:
    """One-gap OU correlation under the complete used marginal law."""
    theta = _positive(cell.theta if theta is None else theta, "theta")
    if cell.mechanism == "regular":
        return math.exp(-theta * cell.q / cell.theta)
    return float(phase_type(cell).alpha @ laplace_kernel(cell, theta) @ np.ones(2))


def lag_correlations(cell: E4Cell, max_lag: int, theta=None) -> np.ndarray:
    """Stationary event-index OU correlations at lags 0,...,max_lag.

    These are unconditional fixed-count correlations. They do not describe
    a schedule conditioned on its span, its count, or a stopping boundary.
    """
    max_lag = _integer(max_lag, "max_lag", 0)
    theta = _positive(cell.theta if theta is None else theta, "theta")
    if cell.mechanism in {"iid", "regular"}:
        return population_rho(cell, theta) ** np.arange(max_lag + 1)
    kernel = laplace_kernel(cell, theta)
    phase = np.array(phase_type(cell).alpha)
    result = np.empty(max_lag + 1)
    result[0] = 1.0
    for h in range(1, max_lag + 1):
        phase = phase @ kernel
        result[h] = phase.sum()
    return result


def pfml_mean_weights(n: int, rho: float) -> np.ndarray:
    """Weights of the conditional PFML mean with its slope held fixed."""
    n = _integer(n, "n", 2)
    rho = _rho(rho)
    k = n - 1
    weights = np.full(n, 1.0 / k)
    weights[0] = -rho / (k * (1.0 - rho))
    weights[-1] = 1.0 / (k * (1.0 - rho))
    return weights


def conditional_linear_variance(times, weights, theta=1.0, sigma=0.5) -> float:
    """Exact Var(sum_i w_i X(t_i) | schedule) in O(n) time and O(n) input.

    Weights must be deterministic given the exogenous schedule. In particular,
    using a slope estimated from X does not give the fitted estimator variance.
    The stationary initial state is integrated out, as in repeated-sample
    research summaries, even though each fitted likelihood conditions on X0.
    """
    times = np.asarray(times, dtype=float)
    weights = np.asarray(weights, dtype=float)
    if (times.ndim != 1 or len(times) < 2 or weights.shape != times.shape
            or not np.isfinite(times).all() or not np.isfinite(weights).all()
            or np.any(np.diff(times) <= 0)):
        raise ValueError("times must increase strictly and weights must match finite times")
    theta, sigma = _positive(theta, "theta"), _positive(sigma, "sigma")
    # The OU innovations representation avoids cancellation for signed weights.
    # If a_i=exp(-theta*gap_i), the coefficient of innovation epsilon_i is
    # w_i + a_(i+1)*w_(i+1) + ... . Accumulate these coefficients backwards.
    correlation = np.exp(-theta * np.diff(times))
    innovation_fraction = -np.expm1(-2.0 * theta * np.diff(times))
    coefficient = float(weights[-1])
    variance_fraction = 0.0
    for i in range(len(times) - 1, 0, -1):
        variance_fraction += coefficient * coefficient * innovation_fraction[i - 1]
        coefficient = float(weights[i - 1]) + correlation[i - 1] * coefficient
    variance_fraction += coefficient * coefficient
    return float(sigma * sigma / (2.0 * theta) * variance_fraction)


def conditional_pfml_variance(times, rho, theta=1.0, sigma=0.5) -> float:
    """Conditional variance of the fixed-slope PFML mean diagnostic."""
    return conditional_linear_variance(times, pfml_mean_weights(len(times), rho), theta, sigma)


def conditional_mean_variance(times, theta=1.0, sigma=0.5) -> float:
    """Conditional variance of the ordinary average of all observed states."""
    return conditional_linear_variance(times, np.full(len(times), 1.0 / len(times)), theta, sigma)


def population_predictions(cell: E4Cell, n=None) -> dict[str, float]:
    """Finite-count fixed-slope benchmarks and leading covariance coefficients.

    Exact variance entries concern linear diagnostics, not the nonlinear
    jointly fitted PFML estimator. The Gaussian coverage entry additionally
    substitutes the population working variance and is only an approximation.
    """
    n = _integer(cell.n if n is None else n, "n", 2)
    k = n - 1
    rho = population_rho(cell)
    corr = lag_correlations(cell, k)
    stationary_variance = cell.sigma**2 / (2.0 * cell.theta)
    weights = pfml_mean_weights(n, rho)
    # Finite autocorrelations of the weights, independently checked by tests.
    lag_weights = np.empty(n)
    lag_weights[0] = np.dot(weights, weights)
    lag_weights[1:k] = (k - np.arange(1, k)) / k**2
    lag_weights[k] = -rho / (k**2 * (1.0 - rho)**2)
    fixed_variance = stationary_variance * (lag_weights[0] + 2 * np.dot(lag_weights[1:], corr[1:]))
    average_variance = stationary_variance * (n + 2 * np.dot(n - np.arange(1, n), corr[1:])) / n**2
    if cell.mechanism == "mmpp":
        kernel = laplace_kernel(cell, cell.theta)
        lag_sum = float(phase_type(cell).alpha @ kernel @ np.linalg.solve(np.eye(2) - kernel, np.ones(2)))
    else:
        lag_sum = rho / (1.0 - rho)
    correlation_factor = 1.0 + 2.0 * lag_sum
    working_variance = stationary_variance * (1.0 + rho) / (k * (1.0 - rho))
    leading_variance = stationary_variance * correlation_factor / k
    coverage = float(2 * ndtr(norm.ppf(.975) * math.sqrt(working_variance / fixed_variance)) - 1)
    return {
        "n": n, "transitions": k, "population_rho": rho,
        "stationary_variance": stationary_variance,
        "fixed_rho_mean_variance": float(fixed_variance),
        "fixed_rho_mean_sd": math.sqrt(fixed_variance),
        "ordinary_mean_variance": float(average_variance),
        "ordinary_mean_sd": math.sqrt(average_variance),
        "long_run_correlation_factor": correlation_factor,
        "leading_mean_variance": leading_variance,
        "population_working_variance": working_variance,
        "population_working_se": math.sqrt(working_variance),
        "normal_working_coverage_approximation": coverage,
        "lag1_correlation": float(corr[1]),
        "lag2_correlation": float(lag_correlations(cell, 2)[2]),
    }
