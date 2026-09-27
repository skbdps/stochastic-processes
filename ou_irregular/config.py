"""Strict, serializable Week 3 experiment configuration.

The YAML grid uses dimensionless coarseness theta*mean_gap. Dividing by the
true theta gives the nominal mean gap in the simulator's time units. Floors
act on generated gaps without rescaling, so this is a nominal design value,
not a promise about every realized path's average or total span.
"""
from dataclasses import asdict, dataclass
import hashlib
import itertools
import json
import math
from pathlib import Path

import yaml


class UniqueKeyLoader(yaml.SafeLoader):
    """Reject duplicate settings instead of silently keeping the last YAML value."""


def _mapping(loader, node, deep=False):
    result = {}
    for key_node, value_node in node.value:
        key = loader.construct_object(key_node, deep=deep)
        if not isinstance(key, str) or key in result:
            raise ValueError(f"YAML settings need unique string keys: {key!r}")
        result[key] = loader.construct_object(value_node, deep=deep)
    return result


UniqueKeyLoader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, _mapping)


def canonical_json(value):
    """Stable encoding used for cell identity, configuration snapshots and replay."""
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _keys(value, expected, label):
    if not isinstance(value, dict) or set(value) != set(expected):
        raise ValueError(f"{label} must have exactly these keys: {sorted(expected)}")


def _number(value, label, *, minimum=None, strict=False):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"{label} must be a finite number")
    value = float(value)
    if minimum is not None and (value <= minimum if strict else value < minimum):
        raise ValueError(f"{label} is outside its allowed domain")
    return value


def _integer(value, label, minimum):
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError(f"{label} must be an integer >= {minimum}")
    return value


def _grid(values, label, integer=False, minimum=0, strict=False):
    if not isinstance(values, list) or not values:
        raise ValueError(f"{label} must be a nonempty list")
    checked = [(_integer(x, label, minimum) if integer else
                _number(x, label, minimum=minimum, strict=strict)) for x in values]
    if len(set(checked)) != len(checked):
        raise ValueError(f"{label} contains duplicate cells")
    return tuple(sorted(checked))


@dataclass(frozen=True)
class ExperimentConfig:
    schema_version: int
    name: str
    phase: str
    seed: int
    replications: int
    truth: dict
    grid: dict
    estimators: tuple
    gap_floor: dict
    optimizer: dict

    def to_dict(self):
        # Round-trip through JSON converts immutable tuples to serializable lists.
        return json.loads(canonical_json(asdict(self)))

    @classmethod
    def from_dict(cls, data):
        _keys(data, cls.__dataclass_fields__, "experiment")
        version = _integer(data["schema_version"], "schema_version", 1)
        if version != 1 or data["phase"] != "smoke":
            raise ValueError("This Week 3 harness accepts schema_version=1, phase=smoke only; headline runs require Week 4 validation")
        if not isinstance(data["name"], str) or not data["name"].strip():
            raise ValueError("name must be nonempty")
        _keys(data["truth"], ("theta", "mu", "sigma"), "truth")
        truth = {key: _number(data["truth"][key], key, minimum=None if key == "mu" else 0,
                              strict=key != "mu") for key in ("theta", "mu", "sigma")}
        _keys(data["grid"], ("theta_mean_gap", "cv", "n"), "grid")
        grid = {"theta_mean_gap": _grid(data["grid"]["theta_mean_gap"], "theta_mean_gap", strict=True),
                "cv": _grid(data["grid"]["cv"], "cv"),
                "n": _grid(data["grid"]["n"], "n", integer=True, minimum=3)}
        if 0.0 not in grid["cv"]:
            raise ValueError("Include CV=0 for each irregular cell's equidistant comparison")
        if len(grid["theta_mean_gap"]) > 3 or len(grid["n"]) > 3:
            raise ValueError("Week 3 smoke plots support at most 3 coarseness values and 3 sample sizes")
        estimators = data["estimators"]
        if not isinstance(estimators, list) or set(estimators) != {"exact", "pfml", "euler"} or len(estimators) != 3:
            raise ValueError("Week 3 requires each of exact, pfml and euler once")
        _keys(data["gap_floor"], ("fraction", "cv_threshold"), "gap_floor")
        floor = {"fraction": _number(data["gap_floor"]["fraction"], "floor fraction", minimum=0, strict=True),
                 "cv_threshold": _number(data["gap_floor"]["cv_threshold"], "floor CV threshold", minimum=0, strict=True)}
        if floor["fraction"] >= 1 or floor["cv_threshold"] != 1.5:
            raise ValueError("Use 0 < floor fraction < 1 and the plan's CV threshold of 1.5")
        _keys(data["optimizer"], ("n_starts", "retry_starts"), "optimizer")
        optimizer = {key: _integer(data["optimizer"][key], key, 1) for key in ("n_starts", "retry_starts")}
        if optimizer["retry_starts"] < optimizer["n_starts"]:
            raise ValueError("retry_starts must be at least n_starts")
        return cls(version, data["name"], "smoke", _integer(data["seed"], "seed", 0),
                   _integer(data["replications"], "replications", 1), truth, grid,
                   ("exact", "pfml", "euler"), floor, optimizer)

    def cells(self):
        """Stable cell IDs exclude grid ordering, replication count and optimizers.

        Thus adding a cell, increasing replications or changing estimator order
        does not change existing simulated paths. All data-generating parameters,
        including floor settings, are part of the SHA-256 identity.
        """
        for n, coarseness, cv in itertools.product(self.grid["n"], self.grid["theta_mean_gap"], self.grid["cv"]):
            cell = {"n": n, "theta_mean_gap": coarseness, "mean_gap": coarseness / self.truth["theta"],
                    "cv": cv, **{f"true_{key}": value for key, value in self.truth.items()},
                    "floor_fraction": self.gap_floor["fraction"],
                    "floor_cv_threshold": self.gap_floor["cv_threshold"]}
            identity = hashlib.sha256(canonical_json(cell).encode()).hexdigest()
            yield {"cell_id": identity, **cell}


def load_config(path):
    """Read safe YAML, reject typo/duplicate keys, and normalize before hashing."""
    with Path(path).open() as handle:
        return ExperimentConfig.from_dict(yaml.load(handle, Loader=UniqueKeyLoader))
