"""Focused contract checks; no schedules, OU samples or fits are generated."""
from copy import deepcopy
from dataclasses import FrozenInstanceError, replace
from fractions import Fraction
import hashlib
import json
import math

import pytest

from ou_irregular.experiments.design import Cell, RunDesign
from ou_irregular.experiments.storage import CellStore
from ou_irregular.week8 import design as e4


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def test_fixed_plan_counts_order_and_separate_validation_namespace():
    full, validation = e4.full_design(), e4.validation_design()
    assert type(full) is RunDesign
    assert (full.name, full.seed, full.replications, full.scope) == (
        "week8_e4_clustered_sampling", 20261011, 2000, "exploratory")
    assert [(c.mechanism, c.ratio) for c in full.cells] == [
        ("mmpp", 1.), ("mmpp", 5.), ("mmpp", 20.),
        ("iid", 1.), ("iid", 5.), ("iid", 20.), ("regular", 0.)]
    assert len({c.identity for c in full.cells}) == 7
    assert all((c.n, c.q, c.theta, c.mu, c.sigma) == (1000, .5, 1., 0., .5) for c in full.cells)
    assert all((c.floor_fraction, c.floor_cv_threshold) == (1e-6, 1.5) for c in full.cells)
    assert len(full.cells) * full.replications == 14000
    assert (validation.name, validation.seed, validation.replications, validation.scope) == (
        "week8_e4_clustered_sampling_validation", 21261011, 4, "validation")
    assert validation.cells == full.cells
    assert validation.digest != full.digest
    assert e4.design_for() == full and e4.design_for(validation=True) == validation
    assert e4.validate_design(full) is full
    assert e4.validate_design(validation) is validation
    with pytest.raises(FrozenInstanceError):
        full.cells[0].ratio = 5.


@pytest.mark.parametrize("r,cv2", [(1, Fraction(1)), (5, Fraction(5, 3)), (20, Fraction(941, 219))])
def test_independent_exact_calibration_values(r, cv2):
    cell = e4.E4Cell(n=1000, q=.5, mechanism="mmpp", ratio=r)
    data = cell.to_dict()
    m = Fraction(1, 2)
    calm = Fraction(6, r + 5) / m
    burst = r * calm
    assert e4.raw_cv_squared(r) == pytest.approx(float(cv2), rel=2e-15)
    assert cell.cv == pytest.approx(math.sqrt(float(cv2)), rel=2e-15)
    assert data["raw_cv_squared"] == pytest.approx(float(cv2), rel=2e-15)
    assert data["calm_sojourn_mean"] == 5.
    assert data["burst_sojourn_mean"] == 1.
    assert data["calm_switch_rate"] == .2 and data["burst_switch_rate"] == 1.
    assert data["calm_arrival_rate"] == pytest.approx(float(calm), rel=2e-15)
    assert data["burst_arrival_rate"] == pytest.approx(float(burst), rel=2e-15)
    assert data["calm_time_probability"] == float(Fraction(5, 6))
    assert data["burst_time_probability"] == float(Fraction(1, 6))
    assert data["calm_event_probability"] == float(Fraction(5, r + 5))
    assert data["burst_event_probability"] == float(Fraction(r, r + 5))
    average_rate = (data["calm_time_probability"] * data["calm_arrival_rate"]
                    + data["burst_time_probability"] * data["burst_arrival_rate"])
    assert average_rate == pytest.approx(2., rel=2e-15)
    iid = replace(cell, mechanism="iid")
    assert iid.identity != cell.identity
    assert {k: v for k, v in data.items() if k not in {"cell_id", "mechanism"}} == {
        k: v for k, v in iid.to_dict().items() if k not in {"cell_id", "mechanism"}}


def test_floor_trigger_and_regular_sentinels_without_mmpp_evaluation(monkeypatch):
    cells = e4.full_design().cells
    assert [(c.mechanism, c.ratio) for c in cells if c.cv >= c.floor_cv_threshold] == [
        ("mmpp", 20.), ("iid", 20.)]
    monkeypatch.setattr(e4, "raw_cv_squared", lambda value: pytest.fail("regular evaluated MMPP formula"))
    regular = e4.E4Cell(n=1000, q=.5, mechanism="regular", ratio=0)
    saved = regular.to_dict()
    assert regular.ratio == regular.cv == saved["raw_cv_squared"] == 0.
    assert all(saved[k] == 0. for k in saved if k.startswith(("calm_", "burst_")))
    assert e4.parse_cell(saved) == regular


@pytest.mark.parametrize("value", [True, False, None, "5", 0, -1, .5, float("nan"), float("inf"), 10**1000])
def test_invalid_raw_ratio_rejected(value):
    with pytest.raises(ValueError):
        e4.raw_cv_squared(value)


def test_finite_large_ratio_calibration_does_not_overflow():
    assert e4.raw_cv_squared(1e300) == pytest.approx(53 / 3)


@pytest.mark.parametrize("changes", [
    {"mechanism": "gamma"}, {"mechanism": True}, {"ratio": True}, {"ratio": 0},
    {"mechanism": "regular", "ratio": 1}, {"q": 0}, {"n": True},
    {"floor_fraction": 0}, {"floor_cv_threshold": False}, {"sigma": float("nan")},
])
def test_bad_cell_configuration_rejected(changes):
    arguments = {"n": 1000, "q": .5, "mechanism": "mmpp", "ratio": 5.}
    arguments.update(changes)
    with pytest.raises(ValueError):
        e4.E4Cell(**arguments)


@pytest.mark.parametrize("changes", [
    {"n": 500}, {"q": 1.}, {"theta": 2.}, {"mu": .1}, {"sigma": .7},
    {"mechanism": "iid"}, {"ratio": 20.}, {"floor_fraction": 1e-4}, {"floor_cv_threshold": 2.5},
])
def test_full_scientific_identity_changes_and_off_grid_cells_parse(changes):
    original = e4.full_design().cells[1]
    cell = replace(original, **changes)
    saved = cell.to_dict()
    assert cell.identity != original.identity
    assert e4.parse_cell(saved) == cell
    body = {k: v for k, v in saved.items() if k != "cell_id"}
    assert cell.identity == hashlib.sha256(canonical(body).encode()).hexdigest()
    assert cell.identity != Cell(cell.n, cell.q, cell.cv, cell.theta, cell.mu, cell.sigma,
                                 cell.floor_fraction, cell.floor_cv_threshold).identity
    assert all(type(v) in (str, int, float) for v in saved.values())


@pytest.mark.parametrize("field,value", [
    ("cell_id", "0" * 64), ("source_namespace", 4), ("source_namespace", 8.),
    ("source_policy_version", "registered_ou_source_v1"), ("floor_policy_version", "different"),
    ("ratio", 5), ("n", 1000.), ("true_theta", True), ("calm_sojourn_mean", 10.),
    ("burst_arrival_rate", 5.), ("calm_event_probability", .4), ("cv", 1.),
    ("raw_cv_squared", float("nan")), ("mean_gap", .6),
])
def test_cell_parser_rejects_noncanonical_or_contradictory_metadata(field, value):
    saved = e4.full_design().cells[1].to_dict()
    saved[field] = value
    with pytest.raises(ValueError):
        e4.parse_cell(saved)


def test_parser_does_not_ignore_surplus_missing_or_rehashed_derived_fields():
    saved = e4.full_design().cells[1].to_dict()
    for modified in ({**saved, "undeclared": 1}, {k: v for k, v in saved.items() if k != "calm_switch_rate"}):
        with pytest.raises(ValueError):
            e4.parse_cell(modified)
    saved["calm_switch_rate"] = .3
    saved["cell_id"] = hashlib.sha256(canonical({k: v for k, v in saved.items() if k != "cell_id"}).encode()).hexdigest()
    with pytest.raises(ValueError):
        e4.parse_cell(saved)


@pytest.mark.parametrize("validation", [False, True])
def test_root_parser_canonical_roundtrip(validation):
    design = e4.design_for(validation=validation)
    payload = design.to_dict()
    before = deepcopy(payload)
    for value in (payload, canonical(payload), (canonical(payload) + "\n").encode()):
        parsed = e4.parse_design(value)
        assert parsed == design and parsed.digest == design.digest
        assert all(type(cell) is e4.E4Cell for cell in parsed.cells)
    assert payload == before


@pytest.mark.parametrize("field,value", [
    ("name", "other"), ("seed", 20261012), ("seed", True), ("replications", 4),
    ("scope", "headline"), ("estimators", ["exact", "euler", "pfml"]),
    ("estimators", ("exact", "pfml", "euler")), ("extra", "unplanned"),
    ("scope", []), ("name", None), ("seed", 20261011.), ("cells", ()),
])
def test_root_parser_rejects_changed_or_noncanonical_plan(field, value):
    saved = e4.full_design().to_dict()
    saved[field] = value
    with pytest.raises(ValueError):
        e4.parse_design(saved)


def test_root_plan_rejects_wrong_order_subsets_duplicate_cells_and_source_only_plan():
    full = e4.full_design()
    wrong = [replace(full, cells=full.cells[::-1]), replace(full, cells=full.cells[:1]),
             full.validation(3), replace(full.validation(256), seed=22261011, cells=full.cells[:6])]
    for design in wrong:
        with pytest.raises(ValueError):
            e4.validate_design(design)
        with pytest.raises(ValueError):
            e4.parse_design(design.to_dict())
    saved = full.to_dict()
    saved["cells"][1] = saved["cells"][0]
    with pytest.raises(ValueError):
        e4.parse_design(saved)
    ordinary = RunDesign(full.name, full.seed, full.replications,
                         tuple(Cell(c.n, c.q, c.cv) for c in full.cells[:3]), full.scope)
    with pytest.raises(ValueError):
        e4.validate_design(ordinary)


def test_serialized_parser_rejects_duplicate_keys_nonfinite_and_format_aliases():
    text = canonical(e4.full_design().to_dict())
    bad = [text + "\n\n", " " + text, json.dumps(e4.full_design().to_dict(), indent=2),
           text[:-1] + ',"seed":20261011}', text.replace('"ratio":1.0', '"ratio":NaN', 1),
           text.replace('"ratio":1.0', '"ratio":Infinity', 1), b"\xff", "[]"]
    for value in bad:
        with pytest.raises(ValueError):
            e4.parse_design(value)


@pytest.mark.parametrize("value", [0, 1, None, "yes"])
def test_validation_selector_does_not_coerce_truthy_values(value):
    with pytest.raises(ValueError):
        e4.design_for(validation=value)


@pytest.mark.parametrize("validation", [False, True])
def test_source_policy_and_explicit_independent_reference_keys(validation):
    design = e4.design_for(validation=validation)
    policy, references = e4.source_policy(design), e4.reference_map(design)
    count = 28 if validation else 14000
    assert policy["name"] == "week8_exogenous_palm_mmpp_v1"
    assert policy["source_namespace"] == 8 and policy["root_seed"] == design.seed
    assert policy["planned_source_realizations"] == policy["planned_independent_source_paths"] == count
    assert policy["planned_fit_records"] == 3 * count
    assert policy["actual_to_seed_anchor"] == {c.identity: c.identity for c in design.cells}
    assert policy["independent_cells_and_replications"] and policy["estimators_share_source_arrays"]
    assert not policy["advance_modulator_during_floor_padding"] and not policy["renormalize_floored_gaps"]
    assert references["regular"] == {c.identity: design.cells[-1].identity for c in design.cells[:-1]}
    assert references["matched_iid"] == {design.cells[i].identity: design.cells[i + 3].identity for i in range(3)}
    assert len(references["regular"]) == 6 and len(references["matched_iid"]) == 3
    policy["stream_children"]["schedule_and_phase"] = 99
    references["regular"].clear()
    assert e4.source_policy(design)["stream_children"]["schedule_and_phase"] == 0
    assert len(e4.reference_map(design)["regular"]) == 6


def test_cellstore_partition_views_preserve_new_metadata_without_becoming_root_plans(tmp_path):
    design = e4.validation_design()
    cell = design.cells[2]
    partition_view = replace(design, cells=(cell,), replications=2)
    assert partition_view.cells == (cell,)
    with pytest.raises(ValueError):
        e4.validate_design(partition_view)

    # Literal ineligible candidates exercise only storage/metadata; these are
    # not generated sources or fitted scientific results.
    def rows(start, stop):
        return [{**cell.to_dict(), "rep": rep, "estimator": method,
                 "point_valid": False, "optimizer_success": True, "point_status": "boundary",
                 "interval_valid": False, "interval_status": "invalid_point_fit",
                 "theta": 1., "mu": 0., "sigma": .5}
                for rep in range(start, stop) for method in ("exact", "pfml", "euler")]

    store = CellStore(tmp_path, cell, design, "literal-design-contract")
    store.append(rows(0, 2), 2)
    resumed = CellStore(tmp_path, cell, design, "literal-design-contract", resume=True)
    resumed.append(rows(2, 4), 4)
    assert resumed.verify() and resumed.state["status"] == "complete"
    records = list(resumed.records())
    assert len(records) == 12
    assert all(all(row[k] == v for k, v in cell.to_dict().items()) for row in records)
    assert e4.validate_design(design) is design
