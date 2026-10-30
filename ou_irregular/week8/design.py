"""Frozen E4 scientific identities, with no source generation or inference.

Public API:
``E4Cell`` is a frozen Cell subclass; ``cv`` is derived, not a free input.
``raw_cv_squared`` evaluates the protocol's raw Palm-gap calibration.
``full_design``, ``validation_design`` and ``design_for`` select the two fixed
fitted-run plans. ``validate_design`` checks their exact root contract.
``parse_cell`` checks a complete canonical cell mapping, including its hash;
it permits structurally valid cells outside the fixed execution grid.
``parse_design`` accepts a mapping or canonical JSON text/bytes and accepts
only the complete full or validation root plan. ``source_policy`` declares
addressing/independence/floor semantics; ``reference_map`` supplies explicit
independent regular and matched-IID references.

The unchanged RunDesign still permits temporary cell/partition views needed
by CellStore. Those views never pass validate_design or become resolved E4
root plans. Source-only validation is a separate harness, not a fitted run.
"""
from collections.abc import Mapping
from dataclasses import dataclass, field
import hashlib
import json
import math

from ..experiments.design import Cell, ESTIMATORS, RunDesign
from ..week4.config import canonical_json


NAME = "week8_e4_clustered_sampling"
FULL_SEED = 20261011
FULL_REPLICATIONS = 2000
VALIDATION_SEED = 21261011
VALIDATION_REPLICATIONS = 4
SOURCE_VALIDATION_SEED = 22261011
SOURCE_VALIDATION_REPLICATIONS = 256
SOURCE_NAMESPACE = 8
SOURCE_POLICY_VERSION = "week8_exogenous_palm_mmpp_v1"
FLOOR_POLICY_VERSION = "week8_posthoc_raw_palm_cv_floor_v1"
RATIOS = (1.0, 5.0, 20.0)
MECHANISMS = ("mmpp", "iid", "regular")
SOURCE_ADDRESS = "SeedSequence(seed,spawn_key=(8,*eight_cell_id_words,rep)).spawn(2)"
FILES = (
    "ou_irregular/week8/__init__.py", "ou_irregular/week8/design.py",
    "docs/week8_gate.md", "docs/week8_exploratory_protocol.md",
    "docs/week8_mathematics.md",
)

_CALIBRATION_FIELDS = (
    "calm_sojourn_mean", "burst_sojourn_mean", "calm_switch_rate", "burst_switch_rate",
    "calm_arrival_rate", "burst_arrival_rate", "calm_time_probability", "burst_time_probability",
    "calm_event_probability", "burst_event_probability", "raw_cv_squared",
)
_CELL_FIELDS = frozenset({
    "cell_id", "n", "theta_mean_gap", "mean_gap", "cv", "true_theta", "true_mu", "true_sigma",
    "floor_fraction", "floor_cv_threshold", "mechanism", "ratio", "source_namespace",
    "source_policy_version", "floor_policy_version", *_CALIBRATION_FIELDS,
})
_ROOT_FIELDS = frozenset({"name", "seed", "replications", "scope", "estimators", "cells"})


def _ratio(value, *, minimum):
    if type(value) not in (int, float):
        raise ValueError("ratio must be a finite real number, not a Boolean or string")
    try:
        result = float(value)
    except (ValueError, OverflowError) as exc:
        raise ValueError("ratio must be finite") from exc
    if not math.isfinite(result) or result < minimum:
        raise ValueError(f"ratio must be finite and >= {minimum}")
    return result


def raw_cv_squared(r):
    """Raw event-gap CV squared for r>=1, calm/burst sojourns 10m/2m.

    The equivalent rational form below avoids squaring a potentially large
    finite ratio: (53r²+110r+125)/(3r²+210r+75). The regular r=0 sentinel
    is deliberately outside this function's domain.
    """
    inverse = 1.0 / _ratio(r, minimum=1.0)
    return ((53.0 + 110.0 * inverse + 125.0 * inverse * inverse)
            / (3.0 + 210.0 * inverse + 75.0 * inverse * inverse))


@dataclass(frozen=True, kw_only=True)
class E4Cell(Cell):
    """Flat mechanism-aware identity; n/q/truth/floor remain ordinary Cell fields.

    CV and all calibration fields are derived from mechanism, ratio and nominal
    mean gap. Regular cells have explicit zero calibration sentinels and never
    evaluate an MMPP expression at ratio zero. Constructor inputs normalize to
    the same scalar representation as Cell; parse_cell requires that already
    normalized representation and rejects contradictory derived metadata.
    """

    cv: float = field(init=False)
    mechanism: str
    ratio: float

    def __post_init__(self):
        if type(self.mechanism) is not str or self.mechanism not in MECHANISMS:
            raise ValueError("mechanism must be mmpp, iid or regular")
        ratio = _ratio(self.ratio, minimum=0.0 if self.mechanism == "regular" else 1.0)
        if self.mechanism == "regular":
            if ratio != 0.0:
                raise ValueError("regular cells require the ratio=0 sentinel")
            ratio, cv = 0.0, 0.0
        else:
            cv = math.sqrt(raw_cv_squared(ratio))
        object.__setattr__(self, "ratio", ratio)
        object.__setattr__(self, "cv", cv)
        super().__post_init__()
        mean_gap = self.q / self.theta
        if not math.isfinite(mean_gap) or mean_gap <= 0:
            raise ValueError("nominal mean gap must be finite and positive")
        try:
            calibration = self._calibration()
        except (OverflowError, ZeroDivisionError) as exc:
            raise ValueError("E4 calibration is outside finite arithmetic") from exc
        if any(not math.isfinite(v) or (v <= 0 if self.mechanism != "regular" else v != 0)
               for v in calibration.values()):
            raise ValueError("E4 calibration must have finite admissible scalars")

    def _calibration(self):
        if self.mechanism == "regular":
            return {name: 0.0 for name in _CALIBRATION_FIELDS}
        m, r = self.q / self.theta, self.ratio
        return {
            "calm_sojourn_mean": 10.0 * m, "burst_sojourn_mean": 2.0 * m,
            "calm_switch_rate": 1.0 / (10.0 * m), "burst_switch_rate": 1.0 / (2.0 * m),
            "calm_arrival_rate": (6.0 / (r + 5.0)) / m,
            "burst_arrival_rate": (6.0 * (r / (r + 5.0))) / m,
            "calm_time_probability": 5.0 / 6.0, "burst_time_probability": 1.0 / 6.0,
            "calm_event_probability": 5.0 / (r + 5.0), "burst_event_probability": r / (r + 5.0),
            "raw_cv_squared": raw_cv_squared(r),
        }

    def to_dict(self):
        common = super().to_dict()
        del common["cell_id"]
        scientific = {
            **common, "mechanism": self.mechanism, "ratio": self.ratio,
            "source_namespace": SOURCE_NAMESPACE, "source_policy_version": SOURCE_POLICY_VERSION,
            "floor_policy_version": FLOOR_POLICY_VERSION, **self._calibration(),
        }
        identity = hashlib.sha256(canonical_json(scientific).encode()).hexdigest()
        return {"cell_id": identity, **scientific}


def full_design():
    """The complete ordered exploratory seven-cell, R=2000 plan."""
    cells = tuple(E4Cell(n=1000, q=.5, mechanism=mechanism, ratio=ratio)
                  for mechanism in ("mmpp", "iid") for ratio in RATIOS)
    cells += (E4Cell(n=1000, q=.5, mechanism="regular", ratio=0.0),)
    return RunDesign(NAME, FULL_SEED, FULL_REPLICATIONS, cells, "exploratory")


def validation_design():
    """Complete seven-cell development fit plan, fixed R=4 and separate seed."""
    return full_design().validation(VALIDATION_REPLICATIONS)


def design_for(*, validation=False):
    if type(validation) is not bool:
        raise ValueError("validation must be Boolean")
    return validation_design() if validation else full_design()


def _mapping(value, fields, label):
    if (not isinstance(value, Mapping) or any(type(key) is not str for key in value)
            or set(value) != fields):
        raise ValueError(f"{label} requires its exact canonical field set")
    return dict(value)


def _same(actual, expected, label):
    """Exact scalar/container types as well as canonical values; bool is not int."""
    if type(actual) is not type(expected):
        raise ValueError(f"noncanonical type in {label}")
    if isinstance(expected, dict):
        if set(actual) != set(expected):
            raise ValueError(f"noncanonical fields in {label}")
        for key, value in expected.items():
            _same(actual[key], value, f"{label}.{key}")
    elif isinstance(expected, list):
        if len(actual) != len(expected):
            raise ValueError(f"noncanonical length in {label}")
        for index, (a, b) in enumerate(zip(actual, expected)):
            _same(a, b, f"{label}[{index}]")
    elif canonical_json(actual) != canonical_json(expected):
        raise ValueError(f"noncanonical value in {label}")


def parse_cell(mapping):
    """Validate every supplied field/ID, permitting structurally valid off-grid cells."""
    value = _mapping(mapping, _CELL_FIELDS, "E4 cell")
    cell = E4Cell(n=value["n"], q=value["theta_mean_gap"], theta=value["true_theta"],
                  mu=value["true_mu"], sigma=value["true_sigma"],
                  floor_fraction=value["floor_fraction"], floor_cv_threshold=value["floor_cv_threshold"],
                  mechanism=value["mechanism"], ratio=value["ratio"])
    _same(value, cell.to_dict(), "E4 cell")
    return cell


def validate_design(design):
    """Accept only the exact complete full/validation E4 root, including order."""
    if type(design) is not RunDesign or any(type(cell) is not E4Cell for cell in design.cells):
        raise ValueError("E4 root requires RunDesign with E4Cell objects")
    expected = validation_design() if design.scope == "validation" else full_design()
    _same(design.to_dict(), expected.to_dict(), "E4 root design")
    return design


def _unique(pairs):
    value = {}
    for key, item in pairs:
        if key in value:
            raise ValueError("duplicate JSON object key")
        value[key] = item
    return value


def _nonfinite(value):
    raise ValueError("nonfinite JSON constant: " + value)


def parse_design(payload):
    """Strict root parser, with canonical JSON optionally ending in one newline.

    Mapping input must already have exact canonical scalar/container types.
    Serialized input additionally rejects duplicate keys, nonfinite constants,
    formatting aliases and surplus whitespace. No file is read or written.
    """
    if isinstance(payload, (str, bytes)):
        try:
            text = payload.decode("utf-8") if isinstance(payload, bytes) else payload
            value = json.loads(text, object_pairs_hook=_unique, parse_constant=_nonfinite)
        except (UnicodeError, json.JSONDecodeError) as exc:
            raise ValueError("invalid E4 design JSON") from exc
        if text not in (canonical_json(value), canonical_json(value) + "\n"):
            raise ValueError("E4 design JSON must use canonical serialization")
    else:
        value = payload
    value = _mapping(value, _ROOT_FIELDS, "E4 root design")
    if type(value["name"]) is not str or type(value["scope"]) is not str:
        raise ValueError("E4 root name and scope must be canonical strings")
    if type(value["seed"]) is not int or type(value["replications"]) is not int:
        raise ValueError("E4 root seed and replications must be canonical integers")
    _same(value["estimators"], list(ESTIMATORS), "E4 estimator order")
    if type(value["cells"]) is not list:
        raise ValueError("E4 root cells must be a canonical list")
    cells = tuple(parse_cell(item) for item in value["cells"])
    design = RunDesign(value["name"], value["seed"], value["replications"], cells, value["scope"])
    _same(value, design.to_dict(), "E4 root design")
    return validate_design(design)


def reference_map(design):
    """Six independent regular references and three MMPP-to-matched-IID references.

    There is no regular self-comparison and no cross-cell paired subset.
    Values are scientific IDs, never surrogate cells made by replacing CV.
    """
    validate_design(design)
    regular = next(cell for cell in design.cells if cell.mechanism == "regular")
    iid = {cell.ratio: cell.identity for cell in design.cells if cell.mechanism == "iid"}
    return {
        "regular": {cell.identity: regular.identity for cell in design.cells if cell.mechanism != "regular"},
        "matched_iid": {cell.identity: iid[cell.ratio] for cell in design.cells if cell.mechanism == "mmpp"},
    }


def source_policy(design):
    """Versioned source contract and planned accounting, never completion claims."""
    validate_design(design)
    sources = len(design.cells) * design.replications
    return {
        "name": SOURCE_POLICY_VERSION, "source_namespace": SOURCE_NAMESPACE,
        "root_seed": design.seed, "source_address": SOURCE_ADDRESS,
        "stream_children": {"schedule_and_phase": 0, "stationary_x0_and_ou_innovations": 1},
        "actual_to_seed_anchor": {cell.identity: cell.identity for cell in design.cells},
        "planned_source_realizations": sources, "planned_independent_source_paths": sources,
        "planned_fit_records": sources * len(ESTIMATORS),
        "independent_cells_and_replications": True, "estimators_share_source_arrays": True,
        "schedule_independent_of_ou": True, "ou_initialization": "independent_stationary",
        "state_order": ["calm", "burst"], "time_origin": "event_at_zero",
        "mechanism_initialization": {
            "mmpp": "palm_initial_phase_then_retain_arrival_phase",
            "iid": "independent_palm_phase_reset_before_each_gap",
            "regular": "constant_nominal_mean_gap",
        },
        "floor_policy_version": FLOOR_POLICY_VERSION,
        "floor_trigger": "nominal_raw_event_gap_cv_ge_cell_threshold",
        "floor_transform": "max_raw_gap_with_cell_fraction_times_nominal_mean",
        "floor_after_raw_gap_generation": True, "advance_modulator_during_floor_padding": False,
        "renormalize_floored_gaps": False,
    }
