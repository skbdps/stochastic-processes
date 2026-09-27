"""Conditional, model-based parametric bootstrap behind inference ports.

The bootstrap fixes the observed timestamps AND observed initial value. For
each independent random stream it simulates the exact OU transition model at
the supplied estimator's fitted parameters, then calls that same estimator.
No population truth, Week 3 runner, or particular optimizer enters this module.

For Exact MLE this is an ordinary plug-in parametric bootstrap. For PFML and
Euler it is a *diagnostic* model-based plug-in bootstrap: fitting an approximate
or misspecified objective does not make its plug-in parameters the true OU
parameters. Its percentile intervals therefore have no claimed robust coverage.
Invalid draws are counted, never replaced by additional simulations.
"""
from dataclasses import dataclass, replace
from math import ceil
from types import MappingProxyType
from typing import Mapping

import numpy as np

from .contracts import PathSimulator, PointEstimator, PointFit, Sample


def _immutable(value):
    """Snapshot diagnostics so the frozen result has no mutable containers."""
    if isinstance(value, Mapping):
        return MappingProxyType({key: _immutable(item) for key, item in value.items()})
    if isinstance(value, np.ndarray):
        result = np.array(value, copy=True)
        result.setflags(write=False)
        return result
    if isinstance(value, (list, tuple)):
        return tuple(_immutable(item) for item in value)
    if isinstance(value, set):
        return frozenset(_immutable(item) for item in value)
    return value


@dataclass(frozen=True)
class BootstrapResult:
    """Auditable bootstrap outcome, including every attempted draw.

    Parameter-vector order is ``theta, mu, sigma``. ``planned_count`` remains
    the configured B even when an invalid original point fit prevents all
    simulations; ``attempted_count`` then equals zero. The reported sample
    standard deviations describe the eligible draws, while intervals are
    withheld unless the preconfigured valid-fraction gate is satisfied.
    """

    point_fit: PointFit
    valid: bool
    status: str
    reason: str
    lower: np.ndarray
    upper: np.ndarray
    standard_errors: np.ndarray
    draw_records: tuple
    planned_count: int
    valid_count: int
    diagnostics: Mapping

    def __post_init__(self):
        for name in ("lower", "upper", "standard_errors"):
            values = np.array(getattr(self, name), dtype=float, copy=True)
            if values.shape != (3,):
                raise ValueError(f"{name} must have shape (3,)")
            values.setflags(write=False)
            object.__setattr__(self, name, values)
        object.__setattr__(self, "point_fit", replace(
            self.point_fit, diagnostics=_immutable(self.point_fit.diagnostics)))
        object.__setattr__(self, "draw_records", tuple(_immutable(row) for row in self.draw_records))
        object.__setattr__(self, "diagnostics", _immutable(self.diagnostics))

    @property
    def attempted_count(self):
        return len(self.draw_records)

    @property
    def invalid_count(self):
        """Failed/invalid attempted draws; skipped draws are not failures."""
        return self.attempted_count - self.valid_count


def _positive_integer(value, label, *, allow_zero=False):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)):
        raise ValueError(f"{label} must be an integer (not bool)")
    value = int(value)
    if value < (0 if allow_zero else 1):
        raise ValueError(f"{label} must be {'nonnegative' if allow_zero else 'positive'}")
    return value


def _probability(value, label, *, closed_upper=False):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, float, np.number)):
        raise ValueError(f"{label} must be a finite probability")
    try:
        value = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{label} must be a finite probability") from exc
    if not np.isfinite(value) or value <= 0 or (value > 1 if closed_upper else value >= 1):
        raise ValueError(f"{label} must be in (0, {'1]' if closed_upper else '1)'}")
    return value


def _eligible(fit):
    """Point validity is necessary but cannot override a boundary diagnosis."""
    return (isinstance(fit, PointFit)
            and isinstance(fit.valid, (bool, np.bool_)) and bool(fit.valid)
            and isinstance(fit.optimizer_success, (bool, np.bool_)) and bool(fit.optimizer_success)
            and fit.status == "valid_interior"
            and np.isfinite(fit.parameters).all() and np.isfinite(fit.nll)
            and fit.theta > 0 and fit.sigma > 0)


class ParametricBootstrap:
    """Simulate/refit service whose dependencies implement small protocols.

    ``draws`` is the total number of planned draws, not a target number of
    successes. At least two eligible draws and ``ceil(min_valid_fraction*B)``
    are needed for intervals. Percentiles use NumPy's explicit ``linear``
    quantile convention; standard errors are sample SDs with denominator M-1,
    where M is the eligible-draw count. Monte Carlo precision is not guaranteed
    merely by passing the valid-fraction gate.
    """

    def __init__(self, estimator: PointEstimator, simulator: PathSimulator,
                 draws=500, level=.95, min_valid_fraction=.95):
        if (not isinstance(getattr(estimator, "name", None), str)
                or not estimator.name or not callable(getattr(estimator, "fit", None))):
            raise ValueError("estimator must provide a nonempty name and callable fit")
        if not callable(getattr(simulator, "simulate", None)):
            raise ValueError("simulator must provide callable simulate")
        self.estimator = estimator
        self.simulator = simulator
        self.draws = _positive_integer(draws, "draws")
        self.level = _probability(level, "level")
        self.min_valid_fraction = _probability(min_valid_fraction, "min_valid_fraction", closed_upper=True)

    def run(self, sample: Sample, point_fit: PointFit | None = None, *, seed=0):
        """Return intervals and raw draw records without altering the sample.

        A fresh SeedSequence is constructed on every call. Its child spawn key
        is the draw index, so a fixed seed reproduces each draw and extending B
        preserves its earlier streams. An ineligible original fit is reported
        as ``skipped_invalid_point`` with zero attempted draws.
        """
        if not isinstance(sample, Sample):
            raise ValueError("sample must be a validated Sample")
        seed = _positive_integer(seed, "seed", allow_zero=True)
        if point_fit is None:
            try:
                point_fit = self.estimator.fit(sample)
            except Exception as exc:
                # A dependency failure is visible as a skipped original fit;
                # KeyboardInterrupt/SystemExit remain interruptible.
                point_fit = PointFit(self.estimator.name, np.nan, np.nan, np.nan,
                                     np.inf, False, False, "numerical_failure",
                                     f"{type(exc).__name__}: {exc}")
        if not isinstance(point_fit, PointFit):
            raise ValueError("estimator must return PointFit")
        if point_fit.estimator != self.estimator.name:
            raise ValueError("point_fit estimator must match the bootstrap estimator")

        required = max(2, ceil(self.min_valid_fraction * self.draws))
        interpretation = ("conditional_exact_ou_parametric_bootstrap"
                          if self.estimator.name == "exact"
                          else "diagnostic_model_based_plugin_bootstrap")
        diagnostics = {
            "method": "percentile", "quantile_method": "linear", "level": self.level,
            "min_valid_fraction": self.min_valid_fraction, "required_valid_count": required,
            "seed_entropy": seed, "rng": "PCG64", "condition_on_times": True,
            "condition_on_x0": True, "x0": float(sample.values[0]),
            "interpretation": interpretation,
            "coverage_guarantee": False,
            "standard_error_definition": "sample_sd_of_valid_draws_ddof_1",
        }
        missing = np.full(3, np.nan)
        if not _eligible(point_fit):
            diagnostics.update(valid_fraction=np.nan, attempted_count=0,
                               invalid_attempted_count=0)
            return BootstrapResult(point_fit, False, "skipped_invalid_point",
                                   f"Original point fit is not a valid interior fit: {point_fit.status}; {point_fit.reason}",
                                   missing, missing, missing, (), self.draws, 0, diagnostics)

        records = []
        estimates = []
        streams = np.random.SeedSequence(seed).spawn(self.draws)
        for draw, stream in enumerate(streams):
            record = {"draw": draw, "seed_entropy": seed,
                      "seed_spawn_key": tuple(stream.spawn_key),
                      "theta": np.nan, "mu": np.nan, "sigma": np.nan, "nll": np.nan,
                      "point_valid": False, "optimizer_success": False,
                      "status": "numerical_failure", "reason": "draw did not finish",
                      "dependency_exception": False, "exception_type": ""}
            try:
                generated = self.simulator.simulate(
                    point_fit.parameters.copy(), sample.times.copy(),
                    np.random.Generator(np.random.PCG64(stream)), x0=float(sample.values[0]))
                # Enforce conditioning at the port boundary: an injected
                # simulator must not silently redraw x0 or sampling times.
                if (not isinstance(generated, Sample)
                        or not np.array_equal(generated.times, sample.times)
                        or generated.values[0] != sample.values[0]):
                    raise ValueError("simulator must preserve all observed timestamps and original x0")
                fitted = self.estimator.fit(generated)
                if not isinstance(fitted, PointFit) or fitted.estimator != self.estimator.name:
                    raise ValueError("bootstrap refit must return PointFit for the same estimator")
                eligible = bool(_eligible(fitted))
                record.update(theta=float(fitted.theta), mu=float(fitted.mu), sigma=float(fitted.sigma),
                              nll=float(fitted.nll), point_valid=eligible,
                              optimizer_success=bool(fitted.optimizer_success),
                              status=fitted.status, reason=fitted.reason)
                if not eligible and fitted.status == "valid_interior":
                    record.update(status="invalid_point_contract",
                                  reason="Interior fit has invalid parameters, objective, or success/validity flags")
                if eligible:
                    estimates.append(fitted.parameters)
            except Exception as exc:
                # Preserve this planned slot, its seed, and the actual error.
                # Never continue drawing until a preferred success count.
                record.update(status="numerical_failure", reason=f"{type(exc).__name__}: {exc}",
                              dependency_exception=True, exception_type=type(exc).__name__)
            records.append(record)

        count = len(estimates)
        se = missing.copy()
        lower, upper = missing.copy(), missing.copy()
        enough = count >= required
        if count >= 2:
            estimates = np.asarray(estimates, dtype=float)
            se = np.std(estimates, axis=0, ddof=1)
            if enough:
                alpha = 1.0 - self.level
                lower, upper = np.quantile(estimates, [alpha / 2, 1 - alpha / 2],
                                           axis=0, method="linear")
        diagnostics.update(valid_fraction=count / self.draws,
                           attempted_count=len(records), invalid_attempted_count=len(records) - count)
        reason = ("Percentile intervals computed from eligible bootstrap fits"
                  if enough else f"Only {count} of {self.draws} draws were valid; require at least {required}")
        if self.estimator.name != "exact":
            reason += "; model-based plug-in diagnostic under possible misspecification, not a robust interval"
        return BootstrapResult(point_fit, enough, "valid" if enough else "insufficient_valid_draws",
                               reason, lower, upper, se, tuple(records), self.draws, count, diagnostics)
