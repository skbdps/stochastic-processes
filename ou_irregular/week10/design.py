"""Fixed-horizon supplementary identities, distinct from the completed E4 study."""
from dataclasses import dataclass
import hashlib
import json
import math

from ..experiments.design import RunDesign
from ..week4.config import canonical_json
from ..week8.design import E4Cell, RATIOS

NAME = "week10_fixed_horizon_supplement"
FULL_SEED = 20261101
VALIDATION_SEED = 21261101
FULL_REPLICATIONS = 2000
VALIDATION_REPLICATIONS = 4
SOURCE_NAMESPACE = 10
POLICY = "week10_palm_fixed_horizon_endpoints_v1"
ENDPOINT_POLICY = "event_at_zero_append_exact_horizon_without_flooring_terminal_gap"


@dataclass(frozen=True, kw_only=True)
class HorizonCell(E4Cell):
    """n is the nominal E4 reference, never the realized observation count."""

    horizon: float = 499.5

    def __post_init__(self):
        super().__post_init__()
        if (isinstance(self.horizon, bool) or not isinstance(self.horizon, (int, float))
                or not math.isfinite(self.horizon) or self.horizon <= 0):
            raise ValueError("horizon must be finite and positive")
        object.__setattr__(self, "horizon", float(self.horizon))

    def to_dict(self):
        data = super().to_dict()
        del data["cell_id"]
        data.update(horizon=self.horizon, nominal_observations=self.n,
                    observation_count_policy="random_with_fixed_endpoints",
                    source_namespace=SOURCE_NAMESPACE, source_policy_version=POLICY,
                    endpoint_policy=ENDPOINT_POLICY)
        return {"cell_id": hashlib.sha256(canonical_json(data).encode()).hexdigest(), **data}

    def parent_cell(self):
        """The unchanged uncensored E4 parent-gap law, without a seed identity."""
        return E4Cell(n=self.n, q=self.q, theta=self.theta, mu=self.mu, sigma=self.sigma,
                      mechanism=self.mechanism, ratio=self.ratio,
                      floor_fraction=self.floor_fraction, floor_cv_threshold=self.floor_cv_threshold)


def full_design():
    cells = tuple(HorizonCell(n=1000, q=.5, mechanism=mechanism, ratio=ratio)
                  for mechanism in ("mmpp", "iid") for ratio in RATIOS)
    cells += (HorizonCell(n=1000, q=.5, mechanism="regular", ratio=0.),)
    return RunDesign(NAME, FULL_SEED, FULL_REPLICATIONS, cells, "exploratory")


def validation_design():
    return full_design().validation(VALIDATION_REPLICATIONS)


def validate_design(design):
    expected = validation_design() if getattr(design, "scope", None) == "validation" else full_design()
    if (type(design) is not RunDesign or any(type(c) is not HorizonCell for c in design.cells)
            or canonical_json(design.to_dict()) != canonical_json(expected.to_dict())):
        raise ValueError("only the complete fixed full or validation horizon design is executable")
    return design


def parse_cell(mapping):
    if not isinstance(mapping, dict):
        raise ValueError("horizon cell requires a canonical mapping")
    try:
        cell = HorizonCell(n=mapping["n"], q=mapping["theta_mean_gap"], theta=mapping["true_theta"],
                           mu=mapping["true_mu"], sigma=mapping["true_sigma"],
                           mechanism=mapping["mechanism"], ratio=mapping["ratio"],
                           floor_fraction=mapping["floor_fraction"],
                           floor_cv_threshold=mapping["floor_cv_threshold"], horizon=mapping["horizon"])
    except (KeyError, TypeError) as exc:
        raise ValueError("incomplete horizon cell mapping") from exc
    if canonical_json(mapping) != canonical_json(cell.to_dict()):
        raise ValueError("noncanonical or inconsistent horizon cell identity")
    return cell


def parse_design(payload):
    value = json.loads(payload) if isinstance(payload, (str, bytes)) else payload
    if not isinstance(value, dict):
        raise ValueError("horizon design requires a mapping")
    try:
        design = RunDesign(value["name"], value["seed"], value["replications"],
                           tuple(parse_cell(c) for c in value["cells"]), value["scope"])
    except (KeyError, TypeError) as exc:
        raise ValueError("incomplete horizon design") from exc
    if canonical_json(value) != canonical_json(design.to_dict()):
        raise ValueError("noncanonical horizon design")
    return validate_design(design)


def source_policy(design):
    validate_design(design)
    return {"version": POLICY, "namespace": SOURCE_NAMESPACE,
            "address": "SeedSequence(seed,spawn_key=(10,*eight_cell_id_words,rep)).spawn(2)",
            "children": {"schedule": 0, "stationary_initial_value_and_ou_innovations": 1},
            "origin": "event_stationary_palm_event_at_zero", "endpoint": ENDPOINT_POLICY,
            "parent_gap_law": "unchanged_E4_full_marginal_with_same_floor",
            "finite_window_observed_gaps_are_not_unconditionally_parent_distributed": True,
            "n_field_is_nominal_reference": True, "actual_n_is_per_source": True,
            "planned_sources": len(design.cells) * design.replications,
            "planned_fits": len(design.cells) * design.replications * 3}
