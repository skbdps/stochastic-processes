"""Strict Week 3 likelihoods and data-normalized, auditable fitting.

All three objectives condition on X_0: the initial stationary density is not
included. Exact and PFML use the exact OU transition; PFML replaces each
observed gap by the realized arithmetic mean. Euler is a Gaussian contrast
with Euler moments fitted to EXACTLY simulated data, not an Euler simulator.
Week 2's optimizer, likelihoods and simulator remain unchanged for reproduction.
"""

import numpy as np
from scipy.optimize import minimize_scalar

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


# Validity is deliberately separate from a numerical solver's termination flag.
# An AR(1) boundary can be solved successfully without identifying a finite OU
# theta. Downstream point summaries must use valid_for_point_summary.
FIT_STATUSES = ("valid_interior", "boundary", "boundary_suspected", "degenerate", "numerical_failure")
FIT_VALIDITY_VERSION = "week3-normalized-v1"
_RESIDUAL_TOL = 64 * np.finfo(float).eps ** 2


class _DegenerateData(ValueError):
    pass


def _normalization(x, gaps):
    """Use only observed data; never use simulation truth or an absolute floor.

    Solve for y=(x-a)/b and u=t/m, where a is the observed mean, b the RMS
    centered state, and m the realized mean gap. Extended intermediates prevent
    avoidable overflow while constructing the dimensionless problem. If the
    original state resolution is already poor, normalization cannot recover it.
    """
    extended = x.astype(np.longdouble)
    a = float(np.mean(extended))
    centered = extended - np.longdouble(a)
    amplitude = np.max(np.abs(centered))
    if not np.isfinite(amplitude) or amplitude == 0:
        raise _DegenerateData("constant path: diffusion variance is not positive")
    b_extended = amplitude * np.sqrt(np.mean((centered / amplitude) ** 2))
    if not 0 < b_extended <= np.finfo(float).max:
        raise _DegenerateData("state scale is not representable")
    b = float(b_extended)
    if b == 0 or not np.isfinite(b):
        raise _DegenerateData("state scale is not representable")
    # A representable offset is allowed, but a path resolved by only a few
    # floating-point steps cannot support the declared 1e-4 parameter tolerance.
    resolution = np.max(np.abs(np.spacing(np.abs(x))))
    if resolution / b > 1e-6:
        raise _DegenerateData("insufficient floating-point state resolution relative to variation")
    y = np.asarray(centered / b, dtype=float)
    m = float(np.mean(gaps.astype(np.longdouble)))
    g = gaps / m
    if not np.isfinite(g).all() or np.any(g <= 0):
        raise _DegenerateData("normalized gaps are not representable")
    predictor = y[:-1] - np.mean(y[:-1])
    if np.dot(predictor, predictor) <= _RESIDUAL_TOL:
        raise _DegenerateData("insufficient predictor variation")
    return y, g, a, b, m


def _variance(residual):
    value = float(np.mean(np.square(residual)))
    if not np.isfinite(value) or value <= _RESIDUAL_TOL:
        raise _DegenerateData("nonpositive or numerically zero residual variance")
    return value


def _candidate(theta, mu, sigma, *, solver, fit_status="valid_interior", reason="interior solution", **diagnostics):
    return dict(theta=float(theta), mu=float(mu), sigma=float(sigma),
                solver=solver, fit_status=fit_status, reason=reason,
                optimizer_success=True, nit=0, **diagnostics)


def _ar1_solution(y, gaps):
    """OLS for regular exact/PFML; project only a REPORTED boundary candidate.

    With an intercept, profiled residual SSE is a convex quadratic in phi.
    Therefore an unrestricted slope outside (0,1) proves that the constrained
    likelihood supremum is on its closure. The finite candidate below is merely
    an approach to that boundary and is never eligible for point summaries.
    """
    X, Y = y[:-1], y[1:]
    xc, yc = X - X.mean(), Y - Y.mean()
    slope = float(np.dot(xc, yc) / np.dot(xc, xc))
    _variance(yc - slope * xc)  # reject a noiseless regression, including outside the model
    status, reason = "valid_interior", "admissible AR(1) OLS slope"
    if slope <= 0:
        phi = 1e-12
        status, reason = "boundary", "AR(1) slope <= 0: theta tends to infinity (phi=0 boundary)"
    elif slope >= 1:
        phi = float(np.exp(-1e-10))
        status, reason = "boundary", "AR(1) slope >= 1: theta tends to zero (phi=1 boundary)"
    else:
        phi = slope
    intercept = float(Y.mean() - phi * X.mean())
    variance = _variance(Y - intercept - phi * X)
    theta = -np.log(phi) / float(np.mean(gaps))
    mu = intercept / (1 - phi)
    sigma = np.sqrt(2 * theta * variance / (-np.expm1(2 * np.log(phi))))
    return _candidate(theta, mu, sigma, solver="normalized_ar1_ols",
                      fit_status=status, reason=reason, unrestricted_phi=slope,
                      candidate_phi=phi)


def _euler_solution(y, gaps):
    """WLS for the Euler Gaussian contrast, with sigma²=SSE/K.

    Regress dx/sqrt(d) on sqrt(d), -x_previous*sqrt(d). The coefficients
    are theta*mu and theta. Negative 1-theta*d is permitted by this model.
    Only theta<=0 violates its parameter space; a finite positive candidate is
    retained in that case, with the intercept reprofiled at the boundary.
    """
    root_gap = np.sqrt(gaps)
    response = np.diff(y) / root_gap
    matrix = np.column_stack((root_gap, -y[:-1] * root_gap))
    coeff, _, rank, _ = np.linalg.lstsq(matrix, response, rcond=None)
    if rank < 2:
        raise _DegenerateData("rank-deficient Euler weighted regression")
    _variance(response - matrix @ coeff)
    drift, theta = map(float, coeff)
    status, reason = "valid_interior", "interior Euler weighted-regression solution"
    unrestricted_theta = theta
    if theta <= 0:
        theta = 1e-10
        drift = float(np.sum(np.diff(y) + theta * y[:-1] * gaps) / np.sum(gaps))
        status, reason = "boundary", "Euler unrestricted theta <= 0: theta=0 boundary"
    residual = response - (drift - theta * y[:-1]) * root_gap
    return _candidate(theta, drift / theta, np.sqrt(_variance(residual)),
                      solver="normalized_euler_wls", fit_status=status,
                      reason=reason, unrestricted_theta=unrestricted_theta)


def _profile_at(log_theta, y, gaps):
    """Profile drift alpha=theta*mu and diffusion variance at fixed theta.

    Parameterizing alpha avoids dividing tiny differences by theta when the
    lower search tail approaches Brownian motion with drift. All calculations
    are the original strict conditional objective at the actual normalized gaps.
    """
    theta = float(np.exp(log_theta))
    td = theta * gaps
    h = -np.expm1(-td)
    A = h / theta
    w = -np.expm1(-2 * td) / (2 * theta)
    response = y[1:] - np.exp(-td) * y[:-1]
    alpha = np.sum(A * response / w) / np.sum(A * A / w)
    residual = response - alpha * A
    variance = float(np.mean(residual * residual / w))
    if not np.isfinite(variance) or variance <= _RESIDUAL_TOL:
        return float("inf"), (theta, float(alpha / theta), float("nan"))
    nll = 0.5 * (len(gaps) * (_LOG_2PI + 1 + np.log(variance)) + np.log(w).sum())
    return float(nll), (theta, float(alpha / theta), float(np.sqrt(variance)))


def _exact_solution(y, gaps):
    """Expand and inspect the scalar profile; do not claim a global proof.

    Search bounds begin at log(theta*mean_gap) in [-8,8], then expand to
    -24 (evaluating both endpoints at each expansion) and high enough that the SHORTEST observed gap has phi<=exp(-40), capped
    at log(theta)=40. Every sampled interior local minimum is polished.
    Compare the selected candidate with both tail endpoints and analytic
    Brownian-drift/iid limits. An unresolved endpoint or a limit tied within
    objective tolerance is boundary_suspected, not a finite identified optimum.
    This finite grid cannot prove absence of a narrow missed local minimum.
    """
    target_upper = float(min(40.0, max(8.0, np.log(40.0) - np.log(np.min(gaps)))))
    lower, upper = -8.0, 8.0
    expansion_history = []
    while True:
        expansion_history.append((lower, upper, _profile_at(lower, y, gaps)[0],
                                  _profile_at(upper, y, gaps)[0]))
        if lower <= -24 and upper >= target_upper:
            break
        lower, upper = max(-24.0, lower - 4), min(target_upper, upper + 4)
    grid = np.linspace(lower, upper, int(np.ceil((upper - lower) / 0.5)) + 1)
    values = np.array([_profile_at(point, y, gaps)[0] for point in grid])
    finite = np.isfinite(values)
    if not finite.any():
        raise _DegenerateData("no positive residual variance on the exact profile")
    candidates = [(float(values[i]), float(grid[i]), False, 0) for i in np.flatnonzero(finite)]
    optimizer_failures = 0
    errors = []
    for index in range(1, len(grid) - 1):
        if values[index] <= min(values[index - 1], values[index + 1]) and np.isfinite(values[index]):
            try:
                result = minimize_scalar(lambda z: _profile_at(z, y, gaps)[0],
                                         bounds=(grid[index - 1], grid[index + 1]),
                                         method="bounded", options={"xatol": 1e-12, "maxiter": 300})
                if np.isfinite(result.fun):
                    candidates.append((float(result.fun), float(result.x), bool(result.success), int(result.nfev)))
                optimizer_failures += not bool(result.success)
            except Exception as exc:
                # Keep the best available finite candidate even when SciPy
                # fails; unsuccessful refinement cannot certify an interior fit.
                errors.append(f"{type(exc).__name__}: {exc}")
                optimizer_failures += 1
    minimum = min(candidates, key=lambda item: item[0])
    tie = 1e-12 + 1e-12 * abs(minimum[0])
    polished = [item for item in candidates if item[2] and item[0] <= minimum[0] + tie]
    value, point, converged, evaluations = min(polished, key=lambda item: item[0]) if polished else minimum
    _, parameters = _profile_at(point, y, gaps)
    dx = np.diff(y)
    drift = float(np.sum(dx) / np.sum(gaps))
    brownian_variance = _variance((dx - drift * gaps) / np.sqrt(gaps))
    lower_limit = 0.5 * (len(gaps) * (_LOG_2PI + 1 + np.log(brownian_variance)) + np.log(gaps).sum())
    upper_limit = 0.5 * len(gaps) * (_LOG_2PI + 1 + np.log(_variance(y[1:] - y[1:].mean())))
    tolerance = 1e-7 + 1e-9 * abs(value)
    resolved = (lower + 1e-5 < point < upper - 1e-5
                and min(values[0], values[-1], lower_limit, upper_limit) > value + tolerance)
    status = "valid_interior" if resolved else "boundary_suspected"
    reason = ("interior profile minimum separated from expanded endpoints and limiting objectives"
              if resolved else "expanded profile endpoint or limiting objective unresolved; no finite optimum certified")
    # A boundary scan can finish successfully without an interior polishing
    # step. An unresolved/interrupted interior refinement cannot do so.
    if not resolved and optimizer_failures == 0:
        converged = True
    if errors or not converged:
        status, reason = "numerical_failure", "selected scalar profile refinement did not converge"
        if errors:
            reason += ": " + "; ".join(errors)
        converged = False
    result = _candidate(*parameters, solver="normalized_exact_profile", fit_status=status, reason=reason,
                       profile_logtheta_lower=lower, profile_logtheta_upper=upper,
                       profile_expansions=len(expansion_history) - 1,
                       profile_search_ranges=";".join(f"{lo:g},{hi:g}" for lo, hi, _, _ in expansion_history),
                       profile_best_nll=float(minimum[0]), profile_selection_gap=float(value - minimum[0]),
                       profile_lower_nll=float(values[0]), profile_upper_nll=float(values[-1]),
                       profile_lower_limit_nll=float(lower_limit), profile_upper_limit_nll=float(upper_limit),
                       profile_grid_points=len(grid), profile_local_failures=int(optimizer_failures))
    result.update(optimizer_success=converged, nit=evaluations,
                  error="; ".join(errors) if errors else None)
    return result


def fit_estimator(estimator, x, times, n_starts=3, retry_starts=5):
    """Fit the strict conditional objective in data-derived dimensionless units.

    Week 3 validation replaces the scale-sensitive joint optimizer with analytic
    AR(1)/Euler regressions and a one-dimensional exact profile. Start-count
    controls remain validated for YAML/backward compatibility, but no fictitious
    starts or retries are reported. ``success`` aliases ``optimizer_success``;
    only ``valid_for_point_summary`` authorizes inclusion in point metrics.

    For y=(x-a)/b and u=t/m, map theta_x=theta_y/m, mu_x=a+b*mu_y,
    sigma_x=b*sigma_y/sqrt(m), NLL_x=NLL_y+K*log(b). Constant or poorly
    resolved data return a degenerate outcome; invalid array/timestamp contracts
    still raise ValueError. Solver exceptions return numerical_failure and do
    not remove a planned result. Boundary candidates are retained, never truth
    imputed. Week 2's ou.fit remains available with its original retry behavior.
    """
    if not isinstance(estimator, str) or estimator not in ESTIMATORS:
        raise ValueError(f"estimator must be one of {ESTIMATORS}")
    n_starts = _start_count(n_starts, "n_starts")
    retry_starts = _start_count(retry_starts, "retry_starts")
    x, _, gaps = _validated_data(x, times)
    base = dict(estimator=estimator, theta=float("nan"), mu=float("nan"), sigma=float("nan"),
                nll=float("inf"), success=False, optimizer_success=False, fit_status="numerical_failure",
                valid_for_point_summary=False, reason="fit not completed", solver="not_started",
                validity_version=FIT_VALIDITY_VERSION, nit=0, status=-1, message="", error=None,
                initial_success=False, retried=False, attempted_starts=0, n_starts=0,
                requested_n_starts=n_starts, requested_retry_starts=retry_starts,
                n_finite=0, n_successful=0, selected_start=-1, best_nll=float("inf"),
                objective_gap=float("nan"), objective_tolerance=float("nan"))
    try:
        y, g, a, b, m = _normalization(x, gaps)
        base.update(state_center=a, state_scale=b, time_scale=m)
        solver_gaps = np.ones_like(g) if estimator == "pfml" else g
        # Near-equidistance is a numerical assumption: arange/cumulative time
        # coordinates accumulate O(K*eps) gap subtraction error. Treat deviations
        # within 128*K*eps of mean-normalized one as regular for the OLS shortcut;
        # evaluate the returned exact NLL at the ACTUAL gaps below. The fit-level
        # reference tests also check original-gap likelihoods. High-CV data never
        # enter this shortcut; no gap is clamped or replaced in the exact objective.
        regular = np.max(np.abs(g - 1)) <= 128 * np.finfo(float).eps * len(g)
        if estimator == "pfml" or (estimator == "exact" and regular):
            result = _ar1_solution(y, solver_gaps)
        elif estimator == "euler":
            result = _euler_solution(y, solver_gaps)
        else:
            result = _exact_solution(y, solver_gaps)
        base.update(result)
        normalized_nll = _strict_nll([np.log(result["theta"]), result["mu"], np.log(result["sigma"])],
                                     y, solver_gaps, estimator)
        shift = len(g) * np.log(b)
        with np.errstate(over="ignore", under="ignore", invalid="ignore"):
            theta = result["theta"] / m
            mu = float(np.longdouble(a) + np.longdouble(b) * result["mu"])
            sigma = b * (result["sigma"] / np.sqrt(m))
        base.update(theta=float(theta), mu=mu, sigma=float(sigma), nll=float(normalized_nll + shift),
                    normalized_nll=normalized_nll, nll_state_jacobian=float(shift))
        # Profile diagnostics use the same original-scale density units as nll.
        for name in ("profile_lower_nll", "profile_upper_nll", "profile_lower_limit_nll", "profile_upper_limit_nll", "profile_best_nll"):
            if name in base:
                base[name] += shift
        if not (np.isfinite([theta, mu, sigma, base["nll"]]).all() and theta > 0 and sigma > 0):
            base.update(fit_status="numerical_failure", optimizer_success=False,
                        reason="mapped parameters or objective are not representable")
        base.update(success=bool(base["optimizer_success"]),
                    valid_for_point_summary=bool(base["optimizer_success"] and base["fit_status"] == "valid_interior"),
                    initial_success=bool(base["optimizer_success"]), n_finite=int(np.isfinite(base["nll"])),
                    n_successful=int(base["optimizer_success"]), best_nll=base["nll"], objective_gap=0.0,
                    objective_tolerance=1e-7 + 1e-9 * abs(base["nll"]))
    except _DegenerateData as exc:
        base.update(fit_status="degenerate", reason=str(exc), error=f"degenerate: {exc}")
    except Exception as exc:
        base.update(fit_status="numerical_failure", reason=f"{type(exc).__name__}: {exc}",
                    error=f"solver: {type(exc).__name__}: {exc}", optimizer_success=False,
                    success=False, valid_for_point_summary=False)
    base.update(status=0 if base["optimizer_success"] else -1, message=base["reason"])
    return base
