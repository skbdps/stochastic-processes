"""Immutable experiment identities, independent of execution and inference.

Cell identity deliberately reproduces the preregistered Week 4 hash. Execution
order, worker count and storage partitions must never enter a scientific seed.
"""
from dataclasses import asdict, dataclass, replace
import hashlib
import math

from ..week4.config import canonical_json


ESTIMATORS = ("exact", "pfml", "euler")


def integer(value, name, minimum):
    if type(value) is not int or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")


@dataclass(frozen=True)
class Cell:
    n: int
    q: float
    cv: float
    theta: float = 1.0
    mu: float = 0.0
    sigma: float = 0.5
    floor_fraction: float = 1e-6
    floor_cv_threshold: float = 1.5

    def __post_init__(self):
        integer(self.n, "n", 3)
        for name in ("q", "cv", "theta", "mu", "sigma", "floor_fraction", "floor_cv_threshold"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
                raise ValueError(f"{name} must be finite numeric")
            if name != "mu" and (value < 0 if name == "cv" else value <= 0):
                raise ValueError(f"{name} outside its domain")
            object.__setattr__(self, name, float(value))
        if self.floor_fraction >= 1:
            raise ValueError("floor_fraction must be below one")

    def to_dict(self):
        cell = {"n": self.n, "theta_mean_gap": self.q, "mean_gap": self.q / self.theta,
                "cv": self.cv, "true_theta": self.theta, "true_mu": self.mu,
                "true_sigma": self.sigma, "floor_fraction": self.floor_fraction,
                "floor_cv_threshold": self.floor_cv_threshold}
        return {"cell_id": hashlib.sha256(canonical_json(cell).encode()).hexdigest(), **cell}

    @property
    def identity(self):
        return self.to_dict()["cell_id"]


@dataclass(frozen=True)
class RunDesign:
    name: str
    seed: int
    replications: int
    cells: tuple[Cell, ...]
    scope: str = "headline"

    def __post_init__(self):
        integer(self.seed, "seed", 0)
        integer(self.replications, "replications", 1)
        if not isinstance(self.name, str) or not self.name.strip():
            raise ValueError("name must be nonempty")
        if self.scope not in ("headline", "validation", "exploratory"):
            raise ValueError("unknown scientific scope")
        object.__setattr__(self, "cells", tuple(self.cells))
        if not self.cells or not all(isinstance(c, Cell) for c in self.cells):
            raise ValueError("cells must contain scientific Cell objects")
        if len({c.identity for c in self.cells}) != len(self.cells):
            raise ValueError("duplicate scientific cells")

    def to_dict(self):
        return {"name": self.name, "seed": self.seed, "replications": self.replications,
                "scope": self.scope, "estimators": list(ESTIMATORS),
                "cells": [c.to_dict() for c in self.cells]}

    @property
    def digest(self):
        return hashlib.sha256(canonical_json(self.to_dict()).encode()).hexdigest()

    def validation(self, replications=3):
        """Use a separate stream; development runs cannot be labelled headlines."""
        return replace(self, name=self.name + "_validation", scope="validation",
                       seed=self.seed + 1_000_000, replications=replications)
