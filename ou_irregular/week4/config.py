"""Strict Week 4 configuration without extending Week 3's smoke schema."""
from copy import deepcopy
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path

import numpy as np
import yaml


def canonical_json(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


class _UniqueLoader(yaml.SafeLoader):
    pass


def _mapping(loader, node, deep=False):
    result = {}
    for key_node, value_node in node.value:
        key = loader.construct_object(key_node, deep=deep)
        if not isinstance(key, str) or key in result:
            raise ValueError("configuration keys must be unique strings")
        result[key] = loader.construct_object(value_node, deep=deep)
    return result


_UniqueLoader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, _mapping)


def _keys(value, expected, label):
    if not isinstance(value, dict) or set(value) != set(expected):
        raise ValueError(f"{label} requires exactly {sorted(expected)}")


def _number(value, label, lower=None, upper=None, inclusive=False):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not np.isfinite(value):
        raise ValueError(f"{label} must be finite and numeric")
    if lower is not None and (value < lower if inclusive else value <= lower):
        raise ValueError(f"{label} below its permitted range")
    if upper is not None and (value > upper if inclusive else value >= upper):
        raise ValueError(f"{label} above its permitted range")


def _integer(value, label, minimum):
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError(f"{label} must be an integer >= {minimum}")


def _cell(value):
    _keys(value, ("q", "cv", "n"), "cell")
    _number(value["q"], "q", lower=0)
    _number(value["cv"], "cv", lower=0, inclusive=True)
    _integer(value["n"], "n", 3)


@dataclass(frozen=True)
class ValidationConfig:
    _data: dict

    @classmethod
    def from_dict(cls, source):
        data = deepcopy(source)
        _keys(data, ("schema_version", "phase", "name", "seed", "truth", "gap_floor",
                     "wald", "certification", "bootstrap"), "Week 4 config")
        if type(data["schema_version"]) is not int or data["schema_version"] != 1:
            raise ValueError("Week 4 schema_version must be 1")
        if data["phase"] != "week4_validation":
            raise ValueError("This runner is restricted to Week 4 validation, not headline experiments")
        if not isinstance(data["name"], str) or not data["name"].strip():
            raise ValueError("name must be nonempty")
        _integer(data["seed"], "seed", 0)
        _keys(data["truth"], ("theta", "mu", "sigma"), "truth")
        for name, value in data["truth"].items():
            _number(value, name, lower=None if name == "mu" else 0)
        _keys(data["gap_floor"], ("fraction", "cv_threshold"), "gap_floor")
        _number(data["gap_floor"]["fraction"], "floor fraction", lower=0, upper=1)
        if data["gap_floor"]["cv_threshold"] != 1.5:
            raise ValueError("floor CV threshold must retain the plan's 1.5")
        _keys(data["wald"], ("level", "relative_step", "hessian_relative_tolerance",
                             "max_whitened_difference", "max_condition", "max_standardized_score"), "wald")
        for name, value in data["wald"].items():
            _number(value, name, lower=0, upper=1 if name == "level" else None)
        if data["wald"]["max_condition"] < 1:
            raise ValueError("max_condition must be at least one")
        if data["wald"]["level"] != .95:
            raise ValueError("the Week 4 checkpoint evaluates 95% Wald intervals")
        cert = data["certification"]
        _keys(cert, ("q", "cv", "n", "replications", "estimator"), "certification")
        _cell({name: cert[name] for name in ("q", "cv", "n")})
        _integer(cert["replications"], "replications", 1)
        if cert["estimator"] != "exact" or cert["cv"] != 0:
            raise ValueError("certification uses exact estimation with regular spacing")
        boot = data["bootstrap"]
        _keys(boot, ("seed", "draws", "level", "min_valid_fraction", "sample_rep", "estimators", "cells"), "bootstrap")
        _integer(boot["seed"], "bootstrap seed", 0)
        _integer(boot["draws"], "bootstrap draws", 2)
        _integer(boot["sample_rep"], "sample_rep", 0)
        _number(boot["level"], "bootstrap level", lower=0, upper=1)
        if boot["level"] != .95:
            raise ValueError("the Week 4 diagnostics use 95% bootstrap intervals")
        _number(boot["min_valid_fraction"], "min_valid_fraction", lower=0, upper=1, inclusive=True)
        if boot["min_valid_fraction"] == 0:
            raise ValueError("min_valid_fraction must be positive")
        if not isinstance(boot["estimators"], list) or sorted(boot["estimators"]) != ["euler", "exact", "pfml"]:
            raise ValueError("bootstrap must list each of exact, pfml, euler once")
        if not isinstance(boot["cells"], list) or not boot["cells"]:
            raise ValueError("bootstrap cells must be nonempty")
        for cell in boot["cells"]:
            _cell(cell)
        if len({(float(cell["q"]), float(cell["cv"]), int(cell["n"]))
                for cell in boot["cells"]}) != len(boot["cells"]):
            raise ValueError("duplicate bootstrap cells")
        return cls(data)

    def to_dict(self):
        return deepcopy(self._data)

    @property
    def digest(self):
        return hashlib.sha256(canonical_json(self._data).encode()).hexdigest()

    def scientific_cell(self, specification):
        """Same scientific identity convention, independent Week 4 stream namespace."""
        truth, floor = self._data["truth"], self._data["gap_floor"]
        cell = {"n": int(specification["n"]), "theta_mean_gap": float(specification["q"]),
                "mean_gap": float(specification["q"] / truth["theta"]),
                "cv": float(specification["cv"]),
                **{f"true_{key}": float(value) for key, value in truth.items()},
                "floor_fraction": float(floor["fraction"]), "floor_cv_threshold": float(floor["cv_threshold"])}
        return {"cell_id": hashlib.sha256(canonical_json(cell).encode()).hexdigest(), **cell}


def load_config(path):
    with Path(path).open() as handle:
        return ValidationConfig.from_dict(yaml.load(handle, Loader=_UniqueLoader))
