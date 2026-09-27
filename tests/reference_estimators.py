"""Independent conditional-likelihood checks used only by acceptance tests.

These routines deliberately do not import a production objective, fitter,
normalizer, or optimizer helper. Exact OU profiling uses the stationary variance
parameter and roots of its analytic envelope score; the production fitter uses
bounded minimization of another profile parameterization. A finite bracketing
search is a numerical reference on these fixtures, not a proof of a global
maximum for every irregular path.

The reference engineering gates were fixed before validating the patch:
parameters rtol=1e-4, atol=1e-6 in the original fixture's normalized coordinates;
original-scale NLL absolute difference <= 1e-7 + 1e-9*abs(reference NLL).
"""
from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq

PARAMETER_RTOL = 1e-4
PARAMETER_ATOL = 1e-6
NLL_ATOL = 1e-7
NLL_RTOL = 1e-9


class ReferenceInvalid(ValueError):
    """No admissible, nondegenerate interior reference was established."""


@dataclass(frozen=True)
class ReferenceFit:
    theta: float
    mu: float
    sigma: float
    nll: float

    @property
    def parameters(self):
        return np.array([self.theta, self.mu, self.sigma])


def _data(x, times):
    x, times = np.asarray(x, dtype=float), np.asarray(times, dtype=float)
    if (x.ndim != 1 or times.ndim != 1 or x.size != times.size or x.size < 3
            or not np.isfinite(x).all() or not np.isfinite(times).all()):
        raise ReferenceInvalid("finite matching vectors with at least three observations required")
    gaps = np.diff(times)
    if not np.isfinite(gaps).all() or np.any(gaps <= 0):
        raise ReferenceInvalid("strictly positive finite gaps required")
    # First-observation / max-deviation coordinates are independent of the
    # production mean/RMS normalization. No simulation truth enters a fit.
    location = x[0]
    scale = np.max(np.abs(x - location))
    if not np.isfinite(scale) or scale <= 0:
        raise ReferenceInvalid("degenerate state variation")
    time_scale = np.median(gaps)
    return x, gaps, (x - location) / scale, gaps / time_scale, location, scale, time_scale


def conditional_nll(estimator, parameters, x, times):
    """Original conditional NLL, written independently of production code.

    Each of the K=n-1 transition densities is counted; there is no initial
    stationary-density term. PFML uses the realized arithmetic mean gap.
    """
    theta, mu, sigma = np.asarray(parameters, dtype=float)
    x, times = np.asarray(x, dtype=float), np.asarray(times, dtype=float)
    gaps = np.diff(times)
    if estimator == "pfml":
        gaps = np.repeat(np.mean(gaps), gaps.size)
    if estimator == "euler":
        mean = x[:-1] + theta * gaps * (mu - x[:-1])
        variance = sigma ** 2 * gaps
    elif estimator in ("exact", "pfml"):
        increment = theta * gaps
        mean = x[:-1] + (-np.expm1(-increment)) * (mu - x[:-1])
        variance = sigma ** 2 * (-np.expm1(-2 * increment)) / (2 * theta)
    else:
        raise ValueError(estimator)
    return float(0.5 * np.sum(np.log(2 * np.pi * variance)
                              + (x[1:] - mean) ** 2 / variance))


def _finish(estimator, normalized, x, times, location, scale, time_scale):
    theta, mu, sigma = normalized
    parameters = np.array([theta / time_scale, location + scale * mu,
                           scale * sigma / np.sqrt(time_scale)])
    if (not np.isfinite(parameters).all() or parameters[0] <= 0
            or parameters[2] <= 0):
        raise ReferenceInvalid("nonpositive or nonfinite reference parameters")
    return ReferenceFit(*parameters, conditional_nll(estimator, parameters, x, times))


def _ols_reference(y, gaps):
    """Unweighted AR(1) OLS followed by the admissible realized-gap OU map."""
    design = np.column_stack((np.ones(y.size - 1), y[:-1]))
    coefficients, _, rank, _ = np.linalg.lstsq(design, y[1:], rcond=None)
    intercept, phi = coefficients
    if rank != 2:
        raise ReferenceInvalid("degenerate predictor variation")
    if not 0 < phi < 1:
        raise ReferenceInvalid("AR slope is outside the finite interior OU domain")
    residual = y[1:] - design @ coefficients
    variance = np.dot(residual, residual) / residual.size
    if variance <= np.finfo(float).eps ** 2:
        raise ReferenceInvalid("nonpositive or numerically degenerate residual variance")
    theta = -np.log(phi) / np.mean(gaps)
    return theta, intercept / (1 - phi), np.sqrt(2 * theta * variance / (1 - phi ** 2))


def _euler_reference(y, gaps):
    """WLS: Δy/sqrt(d)=beta0*sqrt(d)-theta*y_previous*sqrt(d)."""
    root_gap = np.sqrt(gaps)
    design = np.column_stack((root_gap, -y[:-1] * root_gap))
    response = np.diff(y) / root_gap
    coefficients, _, rank, _ = np.linalg.lstsq(design, response, rcond=None)
    intercept, theta = coefficients
    if rank != 2 or theta <= 0:
        raise ReferenceInvalid("degenerate or nonpositive Euler mean-reversion coefficient")
    residual = response - design @ coefficients
    variance = np.dot(residual, residual) / residual.size
    if variance <= np.finfo(float).eps ** 2:
        raise ReferenceInvalid("nonpositive or numerically degenerate residual variance")
    # There is intentionally no requirement that 1-theta*d be positive.
    return theta, intercept / theta, np.sqrt(variance)


def _exact_reference(y, gaps):
    """Profile μ and stationary variance and solve the analytic profile score.

    Set phi_i=exp(-theta*d_i), a_i=1-phi_i, w_i=1-phi_i²,
    z_i=y_i-phi_i*y_previous. Then μ=sum(a*z/w)/sum(a²/w),
    s²=mean((z-a*μ)²/w). Envelope differentiation holds fitted μ,s²
    fixed: sum[w'/(2w)*(1-e²/(s²w)) + e*z_centered'/(s²w)].
    The returned diffusion variance is 2*theta*s².
    """
    count = y.size - 1

    def profile(log_theta):
        theta = np.exp(log_theta)
        phi = np.exp(-theta * gaps)
        a = -np.expm1(-theta * gaps)
        w = -np.expm1(-2 * theta * gaps)
        z = y[1:] - phi * y[:-1]
        denominator = np.sum(a * a / w)
        mu = np.sum(a * z / w) / denominator
        residual = z - a * mu
        stationary_variance = np.sum(residual * residual / w) / count
        if stationary_variance <= 0 or not np.isfinite(stationary_variance):
            return np.inf, np.nan, mu, stationary_variance
        nll = 0.5 * (count * (np.log(2 * np.pi * stationary_variance) + 1)
                     + np.sum(np.log(w)))
        w_prime = 2 * gaps * phi ** 2
        score = np.sum(0.5 * w_prime / w * (1 - residual ** 2 / (stationary_variance * w))
                       + residual * gaps * phi * (y[:-1] - mu) / (stationary_variance * w))
        return nll, theta * score, mu, stationary_variance

    # This broad, pre-fixed dimensionless range brackets the deterministic
    # acceptance fixtures. Explicitly reject unresolved endpoints; do not label
    # a finite search endpoint as a proven finite MLE.
    grid = np.linspace(-20, 20, 321)
    values = [profile(point) for point in grid]
    candidates = []
    for left, right, left_value, right_value in zip(grid[:-1], grid[1:], values[:-1], values[1:]):
        if left_value[1] < 0 < right_value[1]:
            optimum = brentq(lambda value: profile(value)[1], left, right,
                             xtol=1e-13, rtol=1e-14)
            candidates.append((profile(optimum)[0], optimum))
    if not candidates:
        raise ReferenceInvalid("no interior profile minimum bracketed")
    nll, log_theta = min(candidates)
    if nll >= min(values[0][0], values[-1][0]):
        raise ReferenceInvalid("interior candidate does not improve the profile endpoints")
    _, _, mu, stationary_variance = profile(log_theta)
    theta = np.exp(log_theta)
    return theta, mu, np.sqrt(2 * theta * stationary_variance)


def reference_fit(estimator, x, times):
    """Fit using only independent reference formulas and scalar score roots."""
    original_x, _, y, gaps, location, scale, time_scale = _data(x, times)
    if estimator == "exact":
        normalized = _exact_reference(y, gaps)
    elif estimator == "pfml":
        normalized = _ols_reference(y, gaps)
    elif estimator == "euler":
        normalized = _euler_reference(y, gaps)
    else:
        raise ValueError(estimator)
    return _finish(estimator, normalized, original_x, times, location, scale, time_scale)
