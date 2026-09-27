"""Week 3 observation schedules and diagnostics for the actual sampled gaps.

For CV=c>0, D ~ Gamma(k=c**-2, s=mean_gap*c**2), with shape/scale
parameterization. Thus E[D]=k*s=mean_gap and SD(D)/E[D]=1/sqrt(k)=c.
CV=0 is the deterministic limiting case, and CV=1 is exponential sampling.
The optional high-CV floor is censoring D below a positive threshold, not a
new Gamma distribution. It changes both the mean and CV and is never followed
by rescaling. Diagnostics make that change visible in every experiment.
"""

from dataclasses import dataclass

import numpy as np
from scipy.special import gammainc, gammaincc


@dataclass
class SpacingSample:
    """An observation schedule and JSON-compatible spacing diagnostics."""

    times: np.ndarray
    diagnostics: dict


def _positive_finite(value, name):
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be finite and positive")
    try:
        value = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be finite and positive") from exc
    if not np.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be finite and positive")
    return value


def _moments(gaps):
    # Scale before averaging/squaring so representable gaps do not overflow
    # merely because their unscaled sum or second moment would overflow.
    maximum = float(np.max(gaps))
    if maximum == 0:
        return 0.0, float("nan")
    mean = maximum * float(np.mean(gaps / maximum))
    cv = float(np.std(gaps / maximum, ddof=0) / np.mean(gaps / maximum))
    return mean, cv


def generate_times(n, mean_gap, cv, rng, floor_fraction=1e-6,
                   floor_cv_threshold=1.5):
    """Generate ``n`` times starting at zero, with ``n-1`` independent gaps.

    A floor ``a=floor_fraction*mean_gap`` is applied ONLY when
    ``cv >= floor_cv_threshold``. In that regime epsilon must be positive;
    ``None`` is accepted only in the lower-CV regime where no floor is used.
    The transformation is literally ``max(D,a)``; sample means are not fixed.

    For positive CV, the expected floored mean is
    E[max(D,a)] = a*P(k,a/s) + mean_gap*Q(k+1,a/s),
    where P and Q are the regularized lower/upper incomplete Gamma functions.
    This follows by splitting at a and using d*f_k(d)=k*s*f_{k+1}(d).
    ``expected_clip_fraction`` is P(D<a)=P(k,a/s), not a fitted quantity.

    Empirical CVs use population SD (ddof=0) of this finite schedule. The
    returned realized statistics use timestamp differences. If cumulative
    floating-point addition loses a positive gap, raise instead of secretly
    increasing a gap or changing the experimental spacing distribution.
    """
    if (isinstance(n, (bool, np.bool_))
            or not isinstance(n, (int, np.integer)) or n < 3):
        raise ValueError("n must be an integer at least 3")
    n = int(n)
    mean_gap = _positive_finite(mean_gap, "mean_gap")
    floor_cv_threshold = _positive_finite(floor_cv_threshold, "floor_cv_threshold")
    if isinstance(cv, (bool, np.bool_)):
        raise ValueError("cv must be finite and nonnegative")
    try:
        cv = float(cv)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("cv must be finite and nonnegative") from exc
    if not np.isfinite(cv) or cv < 0:
        raise ValueError("cv must be finite and nonnegative")
    if floor_fraction is not None:
        floor_fraction = _positive_finite(floor_fraction, "floor_fraction")
    floored = cv >= floor_cv_threshold
    if floored and floor_fraction is None:
        raise ValueError("a positive floor_fraction is required in the high-CV regime")

    family = "equidistant" if cv == 0 else "exponential" if cv == 1 else "gamma"
    floor_value = 0.0
    expected_mean = mean_gap
    expected_clip = 0.0
    clip_count = 0
    if cv == 0:
        # Multiplication makes the intended equidistant schedule explicit.
        # Record its exact design moments rather than round-off in diff(times).
        with np.errstate(over="ignore", invalid="ignore"):
            times = np.arange(n, dtype=float) * mean_gap
        raw_mean = realized_mean = mean_gap
        raw_cv = realized_cv = 0.0
    else:
        with np.errstate(over="ignore", under="ignore", invalid="ignore", divide="ignore"):
            shape = float(np.float64(cv) ** -2)
            scale = float(np.float64(mean_gap) * np.float64(cv) ** 2)
        if not (np.isfinite(shape) and shape > 0 and np.isfinite(scale) and scale > 0):
            raise ValueError("Gamma shape and scale must be representable and positive")
        raw = np.asarray(rng.gamma(shape, scale, size=n - 1), dtype=float)
        if raw.shape != (n - 1,) or np.any(~np.isfinite(raw)) or np.any(raw < 0):
            raise ValueError("generated Gamma gaps must be finite and nonnegative")
        raw_mean, raw_cv = _moments(raw)
        gaps = raw
        if floored:
            floor_value = floor_fraction * mean_gap
            if not np.isfinite(floor_value) or floor_value <= 0:
                raise ValueError("the high-CV gap floor must be representable and positive")
            clip_count = int(np.count_nonzero(raw < floor_value))
            gaps = np.maximum(raw, floor_value)
            with np.errstate(over="ignore", divide="ignore"):
                argument = np.float64(floor_value) / scale
            expected_clip = float(gammainc(shape, argument))
            expected_mean = float(floor_value * expected_clip
                                  + mean_gap * gammaincc(shape + 1.0, argument))
        with np.errstate(over="ignore", invalid="ignore"):
            times = np.concatenate(([0.0], np.cumsum(gaps)))

    with np.errstate(over="ignore", invalid="ignore"):
        actual_gaps = np.diff(times)
    if (np.any(~np.isfinite(times)) or np.any(~np.isfinite(actual_gaps))
            or np.any(actual_gaps <= 0)):
        raise ValueError("times must be finite and strictly increasing; a gap may be lost to floating-point accumulation")
    if cv > 0:
        realized_mean, realized_cv = _moments(actual_gaps)

    return SpacingSample(times=times, diagnostics={
        "family": family,
        "target_mean_gap": mean_gap,
        "requested_cv": cv,
        "floor_applied": bool(floored),
        "floor_value": float(floor_value),
        "clip_count": clip_count,
        "clip_fraction": clip_count / (n - 1),
        "raw_mean_gap": raw_mean,
        "realized_mean_gap": realized_mean,
        "raw_cv": raw_cv,
        "realized_cv": realized_cv,
        "realized_span": float(times[-1] - times[0]),
        "min_gap": mean_gap if cv == 0 else float(np.min(actual_gaps)),
        "expected_mean_gap_after_floor": expected_mean,
        "expected_clip_fraction": expected_clip,
    })
