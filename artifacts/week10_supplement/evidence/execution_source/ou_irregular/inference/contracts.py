"""Small value objects and ports separating fitting, likelihood and inference.

Core inference consumes these contracts instead of importing Week 3's runner,
configuration, records, or solver internals. Concrete adapters live separately.
The data/parameter order is always theta, mu, sigma; working coordinates are
log(theta), mu, log(sigma). Model truths are absent from all inference inputs.
"""
from dataclasses import dataclass, field
from typing import Mapping, Protocol

import numpy as np


@dataclass(frozen=True)
class Sample:
    values: np.ndarray
    times: np.ndarray

    def __post_init__(self):
        values = np.array(self.values, dtype=float, copy=True)
        times = np.array(self.times, dtype=float, copy=True)
        if (values.ndim != 1 or times.ndim != 1 or values.shape != times.shape
                or values.size < 3 or not np.isfinite(values).all()
                or not np.isfinite(times).all()):
            raise ValueError("finite matching one-dimensional arrays with at least three observations required")
        with np.errstate(over="ignore", invalid="ignore"):
            gaps = np.diff(times)
        if not np.isfinite(gaps).all() or np.any(gaps <= 0):
            raise ValueError("times must have finite strictly positive gaps")
        values.setflags(write=False)
        times.setflags(write=False)
        object.__setattr__(self, "values", values)
        object.__setattr__(self, "times", times)


@dataclass(frozen=True)
class PointFit:
    estimator: str
    theta: float
    mu: float
    sigma: float
    nll: float
    valid: bool
    optimizer_success: bool
    status: str
    reason: str
    diagnostics: Mapping = field(default_factory=dict)

    @property
    def parameters(self):
        return np.array([self.theta, self.mu, self.sigma], dtype=float)


class PointEstimator(Protocol):
    name: str

    def fit(self, sample: Sample) -> PointFit: ...


class ConditionalLikelihood(Protocol):
    name: str

    def nll(self, working_parameters: np.ndarray, sample: Sample) -> float: ...


class PathSimulator(Protocol):
    def simulate(self, parameters: np.ndarray, times: np.ndarray,
                 rng: np.random.Generator, *, x0: float) -> Sample: ...
