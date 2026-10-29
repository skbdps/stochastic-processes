"""E2's exact registered 60-cell design, independent of Week 5's entry point."""
from pathlib import Path

import yaml

from ..config import UniqueKeyLoader
from ..experiments.design import Cell, RunDesign
from ..week4.config import canonical_json


REGISTERED = {"experiment": "E2", "seed": 20261002, "replications": 2000,
              "truth": {"theta": 1.0, "mu": 0.0, "sigma": 0.5},
              "q": [0.1, 0.5, 1.0, 2.0], "cv": [0.0, 0.5, 1.0, 2.0, 3.0],
              "n": [250, 1000, 4000], "floor_fraction": 1e-6,
              "floor_cv_threshold": 1.5}


def load_config(path):
    data = yaml.load(Path(path).read_text(), Loader=UniqueKeyLoader)
    if canonical_json(data) != canonical_json(REGISTERED):
        raise ValueError("E2 configuration differs from the frozen registration")
    return RunDesign("week6_e2", data["seed"], data["replications"], tuple(
        Cell(n=n, q=q, cv=cv, **data["truth"], floor_fraction=data["floor_fraction"],
             floor_cv_threshold=data["floor_cv_threshold"])
        for n in data["n"] for q in data["q"] for cv in data["cv"]))
