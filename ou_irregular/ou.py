"""
Exact OU simulation + three estimators (exact MLE, naive/PFML, Euler contrast).

Model:  dX_t = theta (mu - X_t) dt + sigma dW_t,  theta > 0.
Exact transition over a gap d (companion doc D3):
    X_{t+d} | X_t = x  ~  N( mu + (x - mu) e^{-theta d},  sigma^2/(2 theta) (1 - e^{-2 theta d}) )
Under gaps d_i this is a Gaussian AR(1) with phi_i = e^{-theta d_i} (D4); the exact
log-likelihood conditions on x_0 and sums the transition log-densities (D5).
"""
import numpy as np
from scipy.optimize import minimize

GAP_FLOOR = 1e-12  # numerical floor on gaps inside likelihoods (never on data generation)
OBJECTIVE_TIE_ATOL = 1e-8
OBJECTIVE_TIE_RTOL = 1e-10


# ----------------------------------------------------------------------------- model
def stationary_var(theta, sigma):
    return sigma ** 2 / (2.0 * theta)


def transition_moments(theta, mu, sigma, x_prev, d):
    """Exact conditional mean and variance of X_{t+d} given X_t = x_prev."""
    phi = np.exp(-theta * d)
    mean = mu + (x_prev - mu) * phi
    var = stationary_var(theta, sigma) * (-np.expm1(-2.0 * theta * d))  # = (sigma^2/2theta)(1 - e^{-2 theta d})
    return mean, var


def simulate_ou(theta, mu, sigma, times, rng, x0=None):
    """Exact simulation: recurse the transition density along the observation times.
    If x0 is None, draw X_0 from the stationary law N(mu, sigma^2 / 2 theta)."""
    times = np.asarray(times, dtype=float)
    d = np.diff(times)
    if np.any(d <= 0):
        raise ValueError("times must be strictly increasing")
    n = times.size
    x = np.empty(n)
    x[0] = rng.normal(mu, np.sqrt(stationary_var(theta, sigma))) if x0 is None else float(x0)
    phi = np.exp(-theta * d)
    sd = np.sqrt(stationary_var(theta, sigma) * (-np.expm1(-2.0 * theta * d)))
    z = rng.standard_normal(n - 1)
    for i in range(1, n):
        x[i] = mu + (x[i - 1] - mu) * phi[i - 1] + sd[i - 1] * z[i - 1]
    return x


# ------------------------------------------------------------------- spacing families
def equidistant_times(n, dt):
    """n observation times 0, dt, 2dt, ... (n-1 gaps)."""
    return np.arange(n, dtype=float) * dt


def exponential_gap_times(n, mean_gap, rng):
    """n observation times with i.i.d. exponential gaps (CV = 1); first time = 0."""
    gaps = rng.exponential(mean_gap, size=n - 1)
    return np.concatenate(([0.0], np.cumsum(gaps)))


def gamma_gap_times(n, mean_gap, cv, rng, floor_frac=None):
    """n observation times with i.i.d. Gamma gaps: shape k = CV^-2, scale s = mean_gap * CV^2
    (so E[gap] = mean_gap, CV as requested). cv = 0 returns equidistant times.
    floor_frac: optional pre-registered gap floor, gap <- max(gap, floor_frac * mean_gap).
    Returns (times, clip_fraction)."""
    if cv == 0:
        return equidistant_times(n, mean_gap), 0.0
    k, s = cv ** -2, mean_gap * cv ** 2
    gaps = rng.gamma(k, s, size=n - 1)
    clip = 0.0
    if floor_frac is not None:
        floor = floor_frac * mean_gap
        clip = float(np.mean(gaps < floor))
        gaps = np.maximum(gaps, floor)
    return np.concatenate(([0.0], np.cumsum(gaps))), clip


# ---------------------------------------------------------------------- likelihoods
def _unpack(params):
    log_theta, mu, log_sigma = params
    return np.exp(log_theta), mu, np.exp(log_sigma)


def ou_neg_loglik(params, x, times, gaps=None):
    """Exact negative log-likelihood, conditional on x_0.  params = (log theta, mu, log sigma)."""
    theta, mu, sigma = _unpack(params)
    d = np.diff(times) if gaps is None else gaps
    d = np.maximum(d, GAP_FLOOR)
    phi = np.exp(-theta * d)
    m = mu + (x[:-1] - mu) * phi
    v = (sigma ** 2 / (2.0 * theta)) * (-np.expm1(-2.0 * theta * d))
    r = x[1:] - m
    return 0.5 * np.sum(np.log(2.0 * np.pi * v) + r ** 2 / v)


def naive_neg_loglik(params, x, times):
    """Naive / PFML estimator (Aït-Sahalia & Mykland 2003): the exact likelihood with every
    gap replaced by the mean gap — timestamps enter only through their mean."""
    d = np.diff(times)
    return ou_neg_loglik(params, x, times, gaps=np.full(d.size, d.mean()))


def euler_neg_loglik(params, x, times):
    """Euler / Gaussian-contrast quasi-likelihood (Florens-Zmirou 1989):
    mean x + theta (mu - x) d, variance sigma^2 d.  (Week-3 preview.)"""
    theta, mu, sigma = _unpack(params)
    d = np.maximum(np.diff(times), GAP_FLOOR)
    m = x[:-1] + theta * (mu - x[:-1]) * d
    v = sigma ** 2 * d
    r = x[1:] - m
    return 0.5 * np.sum(np.log(2.0 * np.pi * v) + r ** 2 / v)


# -------------------------------------------------------------------------- fitting
def _moment_start(x, times):
    d = np.diff(times)
    rho = np.corrcoef(x[:-1], x[1:])[0, 1]
    rho = float(np.clip(rho, 0.02, 0.995))
    theta0 = -np.log(rho) / d.mean()
    sigma0 = np.sqrt(2.0 * theta0 * max(np.var(x), 1e-12))
    return np.log(theta0), float(np.mean(x)), np.log(sigma0)


def fit(nll, x, times, n_starts=3):
    """Fit from a positive integer number of deterministic L-BFGS-B starts.

    The first start uses moments; subsequent starts alternate lower and higher
    theta/sigma perturbations, with increasing magnitude in each pair. Prefer a
    converged result only when its objective is within
    ``1e-8 + 1e-10 * abs(best_nll)`` of the best finite objective. A materially
    better failed result remains marked unsuccessful for the caller to handle.
    Raise RuntimeError if no start returns finite, valid parameters and objective.
    Diagnostics describe the selected result; ``selected_start`` is zero-based.
    """
    if (isinstance(n_starts, (bool, np.bool_))
            or not isinstance(n_starts, (int, np.integer)) or n_starts < 1):
        raise ValueError("n_starts must be a positive integer")
    n_starts = int(n_starts)
    lt0, mu0, ls0 = _moment_start(x, times)
    starts = [(lt0, mu0, ls0)]
    for index in range(1, n_starts):
        level = (index + 1) // 2
        theta_factor, sigma_factor = (0.5, 0.7) if index % 2 else (2.0, 1.4)
        starts.append((lt0 + level * np.log(theta_factor), mu0,
                       ls0 + level * np.log(sigma_factor)))

    candidates = []
    for index, s in enumerate(starts):
        res = minimize(nll, np.asarray(s, float), args=(x, times), method="L-BFGS-B",
                       options={"ftol": 1e-11, "gtol": 1e-6, "maxiter": 500})
        with np.errstate(over="ignore", under="ignore", invalid="ignore"):
            params = _unpack(res.x)
        if (np.isfinite(res.fun) and np.all(np.isfinite(params))
                and params[0] > 0 and params[2] > 0):
            candidates.append((res, params, index))
    if not candidates:
        raise RuntimeError(f"OU fit failed: no finite valid solution from {n_starts} starts")

    minimum = min(candidates, key=lambda candidate: candidate[0].fun)
    best_nll = float(minimum[0].fun)
    tolerance = OBJECTIVE_TIE_ATOL + OBJECTIVE_TIE_RTOL * abs(best_nll)
    converged = [candidate for candidate in candidates if candidate[0].success]
    eligible = [candidate for candidate in converged
                if candidate[0].fun - best_nll <= tolerance]
    selected = min(eligible, key=lambda candidate: candidate[0].fun) if eligible else minimum
    best, (theta, mu, sigma), selected_start = selected
    return {"theta": theta, "mu": mu, "sigma": sigma, "nll": float(best.fun),
            "success": bool(best.success), "nit": int(best.nit),
            "status": int(best.status), "message": str(best.message),
            "n_starts": n_starts, "n_finite": len(candidates),
            "n_successful": len(converged), "selected_start": selected_start,
            "best_nll": best_nll, "objective_gap": float(best.fun - best_nll),
            "objective_tolerance": tolerance}


# ------------------------------------------------------------ analytic cross-checks
def ar1_closed_form(x, dt):
    """Closed-form conditional MLE on EQUIDISTANT data: the exact likelihood is a Gaussian AR(1)
    regression of x_i on x_{i-1} with intercept, so the MLE is OLS + the map
    theta = -log(phi)/dt, mu = a/(1-phi), sigma^2 = 2 theta v/(1 - phi^2), v = mean sq. residual."""
    X, Y = x[:-1], x[1:]
    b = np.cov(X, Y, ddof=0)[0, 1] / np.var(X)
    a = Y.mean() - b * X.mean()
    resid = Y - a - b * X
    v = np.mean(resid ** 2)
    theta = -np.log(b) / dt
    mu = a / (1.0 - b)
    sigma2 = 2.0 * theta * v / (1.0 - b ** 2)
    return theta, mu, np.sqrt(sigma2)


def pfml_limit(theta, sigma, mean_gap, cv=1.0):
    """Large-n limit of the naive/PFML estimator under i.i.d. Gamma gaps independent of the path.
    PFML's phi estimate converges to Corr(X_i, X_{i-1}) = E[e^{-theta gap}] (Gamma MGF at -theta);
    its stationary-variance estimate is consistent, so sigma_hat^2 / (2 theta_hat) = sigma^2/(2 theta)."""
    if cv == 0:
        return theta, sigma
    k, s = cv ** -2, mean_gap * cv ** 2
    phi_bar = (1.0 + theta * s) ** (-k)          # E[exp(-theta * Gamma(k, s))]
    theta_lim = -np.log(phi_bar) / mean_gap
    sigma_lim = sigma * np.sqrt(theta_lim / theta)
    return theta_lim, sigma_lim


def euler_limit_equidistant(theta, sigma, dt):
    """Large-n limit of the Euler contrast on equidistant data with gap dt:
    theta_hat -> (1 - e^{-theta dt})/dt,  sigma_hat^2 -> (sigma^2 / (2 theta dt)) (1 - e^{-2 theta dt})."""
    theta_lim = -np.expm1(-theta * dt) / dt
    sigma_lim = np.sqrt((sigma ** 2 / (2.0 * theta * dt)) * (-np.expm1(-2.0 * theta * dt)))
    return theta_lim, sigma_lim
