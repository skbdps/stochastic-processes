"""Reproduce the OU likelihood calculations in the expanded teaching notes.

Python 3.10+; standard library only. No external packages are needed.

The observation times must be fixed or independent of the entire OU path.
These are conditional likelihoods: the initial value is supplied, not modeled
by a parameter-dependent stationary density. Functions evaluate objectives;
this module does not claim to find a global maximum-likelihood estimate.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from collections.abc import Sequence


@dataclass(frozen=True)
class ProfileFit:
    theta: float
    mu: float
    sigma2: float
    nll: float


def validate_data(
    values: Sequence[float], gaps: Sequence[float]
) -> tuple[tuple[float, ...], tuple[float, ...]]:
    """Require n positive gaps and n+1 finite, unnoised observations."""
    x, dt = tuple(map(float, values)), tuple(map(float, gaps))
    if not dt or len(x) != len(dt) + 1:
        raise ValueError("Need n >= 1 gaps and exactly n+1 observations")
    if not all(math.isfinite(v) for v in x):
        raise ValueError("Observations must be finite")
    if not all(math.isfinite(h) and h > 0 for h in dt):
        raise ValueError("Gaps must be finite and strictly positive")
    return x, dt


def transition_coefficients(theta: float, gap: float) -> tuple[float, float, float]:
    """Return phi, b=1-phi, c such that the transition variance is sigma^2*c."""
    if not math.isfinite(theta) or theta <= 0:
        raise ValueError("theta must be finite and positive")
    if not math.isfinite(gap) or gap <= 0:
        raise ValueError("gap must be finite and positive")
    z = theta * gap
    if not math.isfinite(z) or z == 0:
        raise ArithmeticError("theta * gap is outside the supported numeric range")
    phi = math.exp(-z)
    b = -math.expm1(-z)
    # expm1 avoids cancellation at short gaps. The long-gap expression
    # avoids overflow in 2*z or 2*theta when evaluating a saturated transition.
    c = 0.5 / theta if z > 350 else gap * (-math.expm1(-2 * z) / (2 * z))
    if not math.isfinite(c) or c <= 0:
        raise ArithmeticError("Nonpositive or non-finite transition coefficient")
    return phi, b, c


def conditional_nll(
    theta: float, mu: float, sigma: float,
    values: Sequence[float], gaps: Sequence[float]
) -> float:
    """Evaluate the exact conditional OU negative log-likelihood.

    The recorded predecessor, not the previous predicted mean, is used in
    each transition. All parameter-dependent variance terms are retained.
    """
    x, dt = validate_data(values, gaps)
    if not math.isfinite(mu) or not math.isfinite(sigma) or sigma <= 0:
        raise ValueError("mu must be finite; sigma must be finite and positive")
    terms = []
    for previous, observed, h in zip(x[:-1], x[1:], dt):
        phi, b, c = transition_coefficients(theta, h)
        mean = phi * previous + b * mu
        variance = (sigma * sigma) * c
        if not math.isfinite(variance) or variance <= 0:
            raise ArithmeticError("Transition variance is outside numeric range")
        standardized = (observed - mean) / math.sqrt(variance)
        term = 0.5 * (math.log(2 * math.pi) + math.log(variance)
                      + standardized * standardized)
        if not math.isfinite(term):
            raise ArithmeticError("Non-finite likelihood contribution")
        terms.append(term)
    return math.fsum(terms)


def profile_at_theta(
    theta: float, values: Sequence[float], gaps: Sequence[float]
) -> ProfileFit:
    """Analytically fit free mu and positive sigma^2 at a FIXED theta.

    Return the one-dimensional profile objective plus the associated level
    and variance. A zero residual sum is degenerate, not an interior MLE.
    """
    x, dt = validate_data(values, gaps)
    coefficients = [transition_coefficients(theta, h) for h in dt]
    b = [item[1] for item in coefficients]
    c = [item[2] for item in coefficients]
    d = [next_value - item[0] * previous
         for previous, next_value, item in zip(x[:-1], x[1:], coefficients)]
    denominator = math.fsum(bi * bi / ci for bi, ci in zip(b, c))
    numerator = math.fsum(bi * di / ci for bi, di, ci in zip(b, d, c))
    if not math.isfinite(denominator) or denominator <= 0:
        raise ArithmeticError("Degenerate profile weights")
    mu = numerator / denominator
    squared_error = math.fsum((di - bi * mu) ** 2 / ci
                             for bi, di, ci in zip(b, d, c))
    sigma2 = squared_error / len(dt)
    if not math.isfinite(mu) or not math.isfinite(sigma2) or sigma2 <= 0:
        raise ValueError("Degenerate or non-finite profiled fit")
    nll = 0.5 * (len(dt) * (math.log(2 * math.pi) + 1 + math.log(sigma2))
                 + math.fsum(math.log(ci) for ci in c))
    if not math.isfinite(nll):
        raise ArithmeticError("Non-finite profile objective")
    return ProfileFit(theta, mu, sigma2, nll)


def main() -> None:
    values = (29.0, 27.5, 26.0, 25.8)
    times = (0.0, 0.2, 1.0, 1.3)
    gaps = tuple(right - left for left, right in zip(times[:-1], times[1:]))
    print("Candidate evaluations, NOT a global MLE search")
    print("theta   fixed(mu=25,sigma=2) NLL     profiled mu    profiled sigma^2    profile NLL")
    for theta in (0.5, 1.0, 2.0):
        nll = conditional_nll(theta, 25.0, 2.0, values, gaps)
        fit = profile_at_theta(theta, values, gaps)
        print(f"{theta:4.1f}    {nll:24.8f}  {fit.mu:14.8f}  {fit.sigma2:18.8f}  {fit.nll:13.8f}")
        # The direct and profiled formulas must describe the same objective
        # when evaluated at identical parameters and data.
        direct = conditional_nll(theta, fit.mu, math.sqrt(fit.sigma2), values, gaps)
        if not math.isclose(direct, fit.nll, rel_tol=1e-10, abs_tol=1e-10):
            raise AssertionError("Profile and direct objective disagree")


if __name__ == "__main__":
    main()
