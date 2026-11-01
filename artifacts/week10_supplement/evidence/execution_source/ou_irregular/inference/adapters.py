"""Explicit one-way adapters to the stable Week 3 public APIs.

Only this integration boundary imports the existing fitter/objectives/simulator.
Week 3 never imports inference and does not change its validity or output schema.
New likelihood classes can implement the ports without changing the core.
"""
from dataclasses import dataclass

from ..estimators import (ESTIMATORS, fit_estimator, exact_neg_loglik,
                          pfml_neg_loglik, euler_neg_loglik)
from ..ou import simulate_ou
from .contracts import Sample, PointFit


@dataclass(frozen=True)
class Week3Estimator:
    name: str = "exact"

    def __post_init__(self):
        if self.name not in ESTIMATORS:
            raise ValueError(f"estimator must be one of {ESTIMATORS}")

    def fit(self, sample):
        result = fit_estimator(self.name, sample.values, sample.times)
        return PointFit(self.name, result["theta"], result["mu"], result["sigma"],
                        result["nll"], result["valid_for_point_summary"],
                        result["optimizer_success"], result["fit_status"],
                        result["reason"], dict(result))


@dataclass(frozen=True)
class OuLikelihood:
    name: str = "exact"

    def __post_init__(self):
        if self.name not in ESTIMATORS:
            raise ValueError(f"likelihood must be one of {ESTIMATORS}")

    def nll(self, working_parameters, sample):
        objective = {"exact": exact_neg_loglik, "pfml": pfml_neg_loglik,
                     "euler": euler_neg_loglik}[self.name]
        return objective(working_parameters, sample.values, sample.times)


class ExactOuSimulator:
    def simulate(self, parameters, times, rng, *, x0):
        theta, mu, sigma = parameters
        return Sample(simulate_ou(theta, mu, sigma, times, rng, x0=x0), times)
