"""E1's frozen design; changes require an explicitly separate experiment."""
from pathlib import Path

import yaml

from ..config import UniqueKeyLoader
from ..experiments.design import Cell, RunDesign
from ..week4.config import canonical_json


REGISTERED = {"experiment": "E1", "seed": 20261001, "replications": 2000,
              "truth": {"theta": 1.0, "mu": 0.0, "sigma": 0.5},
              "q": [0.5], "cv": [0.0, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0], "n": [1000],
              "floor_fraction": 1e-6, "floor_cv_threshold": 1.5}


def load_config(path):
    data = yaml.load(Path(path).read_text(), Loader=UniqueKeyLoader)
    # Canonical JSON distinguishes booleans from numbers and rejects NaN.
    if canonical_json(data) != canonical_json(REGISTERED):
        raise ValueError("E1 configuration differs from the frozen registration")
    return RunDesign("week5_e1", data["seed"], data["replications"], tuple(
        Cell(n=n, q=q, cv=cv, **data["truth"], floor_fraction=data["floor_fraction"],
             floor_cv_threshold=data["floor_cv_threshold"])
        for n in data["n"] for q in data["q"] for cv in data["cv"]))
