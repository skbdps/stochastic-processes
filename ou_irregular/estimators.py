"""Strict Week 3 likelihoods and a shared, auditable optimizer retry policy.

All three objectives condition on X_0: the initial stationary density is not
included. Exact and PFML use the exact OU transition; PFML replaces each
observed gap by the realized arithmetic mean. Euler is a Gaussian contrast
with Euler moments fitted to EXACTLY simulated data, not an Euler simulator.
Week 2's public functions in ou.py remain unchanged for reproduction.
"""

import numpy as np

from . import ou

ESTIMATORS = ("exact", "pfml", "euler")
_LOG_2 = np.log(2.0)
_LOG_2PI = np.log(2.0 * np.pi)


def _validated_data(x, times):
    try:
        x = np.asarray(x, dtype=float)
        times = np.asarray(times, dtype=float)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("x and times must be finite one-dimensional arrays") from exc
    if x.ndim != 1 or times.ndim != 1 or x.size != times.size or x.size < 3:
        raise ValueError("x and times must be one-dimensional, matching, with at least 3 observations")
    if np.any(~np.isfinite(x)) or np.any(~np.isfinite(times)):
        raise ValueError("x and times must contain only finite values")
    with np.errstate(over="ignore", invalid="ignore"):
        gaps = np.diff(times)
    if np.any(~np.isfinite(gaps)) or np.any(gaps <= 0):
        raise ValueError("times must have finite, strictly positive gaps")
    return x, times, gaps


def _strict_nll(params, x, gaps, estimator):
    """Evaluate a Gaussian transition objective without changing any data gap.

    Log coordinates enforce theta,sigma>0. Impossible floating-point parameter
    proposals receive +inf rather than a clipped theta, variance, or gap. The
    optimizer therefore always sees the stated mathematical model.
    """
    try:
        params = np.asarray(params, dtype=float)
    except (TypeError, ValueError, OverflowError):
        return float("inf")
    if params.shape != (3,) or np.any(~np.isfinite(params)):
        return float("inf")
    log_theta, mu, log_sigma = params
    with np.errstate(over="ignore", under="ignore", invalid="ignore", divide="ignore"):
        theta, sigma = np.exp(log_theta), np.exp(log_sigma)
        if not (np.isfinite(theta) and theta > 0 and np.isfinite(sigma) and sigma > 0):
            return float("inf")
        if estimator == "euler":
            # Euler: m_i=x_{i-1}+theta*(mu-x_{i-1})*d_i, v_i=sigma^2*d_i.
            mean = x[:-1] + theta * gaps * (mu - x[:-1])
            log_variance = 2.0 * log_sigma + np.log(gaps)
        else:
            # Exact: phi_i=e^(-theta*d_i), v_i=(sigma^2/(2theta))*(1-phi_i^2).
            # expm1 retains precision when theta*d_i is close to zero.
            td = theta * gaps
            one_minus_phi = -np.expm1(-td)
            mean = np.exp(-td) * x[:-1] + one_minus_phi * mu
            u = 2.0 * td
            small = u < 1e-4
            log_variance = np.empty_like(gaps)
            # v=sigma^2*d*[(-expm1(-u))/u], with the ratio's limit 1 at u=0.
            # This branch handles theta*d underflow without losing valid gaps.
            ratio = np.ones_like(u[small])
            np.divide(-np.expm1(-u[small]), u[small], out=ratio, where=u[small] != 0)
            log_variance[small] = 2.0 * log_sigma + np.log(gaps[small]) + np.log(ratio)
            # Log variance avoids premature overflow of sigma^2 or 2*theta.
            log_variance[~small] = (2.0 * log_sigma - _LOG_2 - log_theta
                                    + np.log(-np.expm1(-u[~small])))
        variance = np.exp(log_variance)
        if (np.any(~np.isfinite(mean)) or np.any(~np.isfinite(variance))
                or np.any(variance <= 0)):
            return float("inf")
        standardized = (x[1:] - mean) / np.sqrt(variance)
        value = 0.5 * np.sum(_LOG_2PI + log_variance + standardized ** 2)
    return float(value) if np.isfinite(value) else float("inf")


def exact_neg_loglik(params, x, times):
    """Exact OU conditional NLL at the actual observed positive gaps."""
    x, _, gaps = _validated_data(x, times)
    return _strict_nll(params, x, gaps, "exact")


def pfml_neg_loglik(params, x, times):
    """Pretend-fixed-mesh NLL: replace gaps by their REALIZED arithmetic mean."""
    x, _, gaps = _validated_data(x, times)
    return _strict_nll(params, x, np.full_like(gaps, gaps.mean()), "pfml")


def euler_neg_loglik(params, x, times):
    """Euler Gaussian contrast conditional on X_0 at actual observed gaps."""
    x, _, gaps = _validated_data(x, times)
    return _strict_nll(params, x, gaps, "euler")


def _start_count(value, name):
    if (isinstance(value, (bool, np.bool_))
            or not isinstance(value, (int, np.integer)) or value < 1):
        raise ValueError(f"{name} must be a positive integer")
    return int(value)


def _failed_result(n_starts, message):
    return {
        "theta": float("nan"), "mu": float("nan"), "sigma": float("nan"),
        "nll": float("inf"), "success": False, "nit": 0, "status": -1,
        "message": message, "n_starts": n_starts, "n_finite": 0,
        "n_successful": 0, "selected_start": -1, "best_nll": float("inf"),
        "objective_gap": float("nan"), "objective_tolerance": float("nan"),
    }


def _select_across_attempts(initial, retry):
    """Preserve ou.fit's objective-aware convergence rule across retry calls."""
    finite = [result for result in (initial, retry) if np.isfinite(result["nll"])]
    if not finite:
        return dict(retry)
    best_nll = min(float(result["best_nll"]) for result in finite)
    tolerance = ou.OBJECTIVE_TIE_ATOL + ou.OBJECTIVE_TIE_RTOL * abs(best_nll)
    eligible = [result for result in finite
                if result["success"] and result["nll"] - best_nll <= tolerance]
    selected = dict(min(eligible or finite, key=lambda result: result["nll"]))
    selected.update(best_nll=best_nll, objective_gap=selected["nll"] - best_nll,
                    objective_tolerance=tolerance)
    return selected


def fit_estimator(estimator, x, times, n_starts=3, retry_starts=5):
    """Fit one model, retrying deterministically only after initial failure.

    Default policy: one 3-start ou.fit call; if it is unsuccessful (including
    an exception), one 5-start ou.fit call from the SAME deterministic schedule.
    This is 8 attempted starts, including the repeated first three. No random
    seed or model changes occur on retry. A converged retry is accepted only if
    its NLL is within ou.fit's numerical tolerance of the best finite objective
    across both calls; a materially better failed optimum remains a failure.

    Invalid data, estimator names, or start counts raise before optimization.
    Optimizer exceptions are recorded in ``error`` and cannot silently remove
    a replicate. If neither call yields a finite fit, estimates are NaN.
    ``attempted_starts`` counts all requested starts; ``n_starts``, ``n_finite``,
    ``n_successful`` and ``selected_start`` describe the selected ou.fit call.
    ``best_nll`` and its tolerance/gap cover both calls after a retry.
    """
    if not isinstance(estimator, str) or estimator not in ESTIMATORS:
        raise ValueError(f"estimator must be one of {ESTIMATORS}")
    n_starts = _start_count(n_starts, "n_starts")
    retry_starts = _start_count(retry_starts, "retry_starts")
    x, times, gaps = _validated_data(x, times)
    if estimator == "pfml":
        gaps = np.full_like(gaps, gaps.mean())

    # Validate/calculate gaps once, so optimization does not repeatedly allocate
    # differences or accidentally call Week 2's gap-clamping likelihoods.
    def objective(params, _x, _times):
        return _strict_nll(params, x, gaps, estimator)

    errors = []

    def attempt(starts, label):
        try:
            return ou.fit(objective, x, times, n_starts=starts)
        except Exception as exc:
            message = f"{label} fit: {type(exc).__name__}: {exc}"
            errors.append(message)
            return _failed_result(starts, message)

    initial = attempt(n_starts, "initial")
    initial_success = bool(initial["success"])
    retried = not initial_success
    selected = dict(initial)
    attempted_starts = n_starts
    if retried:
        retry = attempt(retry_starts, "retry")
        attempted_starts += retry_starts
        selected = _select_across_attempts(initial, retry)
    selected.update(estimator=estimator, initial_success=initial_success,
                    retried=retried, attempted_starts=attempted_starts,
                    error="; ".join(errors) if errors else None)
    return selected
