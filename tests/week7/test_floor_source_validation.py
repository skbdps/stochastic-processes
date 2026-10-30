"""Independent identity/domain checks for the exploratory floor-source adapter.

The source-construction spy uses literal observations. The sole unmocked source
check uses the separate development seed, never the full exploratory stream.
Direct Gamma/OU numerical references belong to test_floor_source_reference.py.
"""
from copy import deepcopy
from dataclasses import FrozenInstanceError, replace
from types import MappingProxyType, SimpleNamespace

import numpy as np
import pytest

from ou_irregular.experiments.design import Cell, ESTIMATORS
from ou_irregular.experiments.execution import failed_result, source_records
from ou_irregular.inference.contracts import PointFit, Sample
from ou_irregular.week4.simulation import array_digest
from ou_irregular.week7 import floor_source as floor


POLICY = "week7_paired_floor_v1_shared_gamma_stationary_x0_and_standard_normal_innovations"
VALIDATION_SEED = 21261009
CELLS = tuple(Cell(n=1000, q=.5, cv=cv, floor_fraction=epsilon)
              for cv in (0, 1.5, 2, 3) for epsilon in (1e-8, 1e-6, 1e-4))


@pytest.mark.parametrize("cell", CELLS)
def test_all_twelve_canonical_actual_cells_are_accepted_without_mutation(cell):
    supplied = cell.to_dict()
    before = deepcopy(supplied)
    actual = floor.validate_cell(MappingProxyType(supplied))
    anchor = floor.anchor_cell(supplied)
    assert actual == cell
    assert anchor == Cell(n=1000, q=.5, cv=cell.cv, floor_fraction=1e-6)
    assert floor.anchor_cell(cell) == anchor
    assert supplied == before


@pytest.mark.parametrize("field,bad", [
    ("n", True), ("n", 1000.), ("n", "1000"), ("n", np.int64(1000)),
    ("theta_mean_gap", True), ("theta_mean_gap", ".5"), ("mean_gap", .50000000001),
    ("mean_gap", True), ("mean_gap", np.nan), ("mean_gap", float("inf")),
    ("cv", 3), ("cv", True), ("cv", "3"), ("cv", float("nan")),
    ("true_theta", 1), ("true_theta", False), ("true_theta", 0.),
    ("true_mu", 0), ("true_mu", False), ("true_mu", -0.0), ("true_mu", float("inf")),
    ("true_sigma", 0.), ("true_sigma", "0.5"),
    ("floor_fraction", 0.), ("floor_fraction", True), ("floor_fraction", "1e-8"),
    ("floor_cv_threshold", 1.4), ("floor_cv_threshold", True),
    ("cell_id", "a" * 64), ("cell_id", "ABCDEF" * 10 + "ABCD"), ("cell_id", None),
])
def test_noncanonical_or_inconsistent_actual_mapping_is_rejected(field, bad):
    value = CELLS[-3].to_dict()
    value[field] = bad
    with pytest.raises(ValueError):
        floor.validate_cell(value)
    with pytest.raises(ValueError):
        floor.anchor_cell(value)


@pytest.mark.parametrize("mode", ["missing", "extra", "list", "cell_object"])
def test_actual_mapping_requires_exact_declared_fields(mode):
    value = CELLS[0].to_dict()
    if mode == "missing":
        value.pop("mean_gap")
    elif mode == "extra":
        value["seed_anchor_cell_id"] = value["cell_id"]
    elif mode == "list":
        value = list(value.items())
    else:
        value = CELLS[0]
    with pytest.raises(ValueError):
        floor.validate_cell(value)


@pytest.mark.parametrize("cell", [
    Cell(n=999, q=.5, cv=3, floor_fraction=1e-8),
    Cell(n=1000, q=.1, cv=3, floor_fraction=1e-8),
    Cell(n=1000, q=.5, cv=1, floor_fraction=1e-8),
    Cell(n=1000, q=.5, cv=3, floor_fraction=1e-7),
    Cell(n=1000, q=.5, cv=3, floor_fraction=1e-8, theta=2.),
    Cell(n=1000, q=.5, cv=3, floor_fraction=1e-8, mu=.1),
    Cell(n=1000, q=.5, cv=3, floor_fraction=1e-8, sigma=1.),
    Cell(n=1000, q=.5, cv=3, floor_fraction=1e-8, floor_cv_threshold=2.),
])
def test_other_canonical_scientific_cells_do_not_silently_enter_this_protocol(cell):
    with pytest.raises(ValueError):
        floor.validate_cell(cell.to_dict())
    with pytest.raises(ValueError):
        floor.anchor_cell(cell)


@pytest.fixture
def delegate_spy(monkeypatch):
    calls = []
    sample = Sample(values=[.2, .1, -.3, .4], times=[0., .2, .4, .6])
    diagnostics = {"source_seed": VALIDATION_SEED, "source_spawn_key": "[4,1,0]",
                   "gap_spawn_key": "[4,1,0,0]", "path_spawn_key": "[4,1,0,1]",
                   "times_sha256": array_digest(sample.times),
                   "observations_sha256": array_digest(sample.values),
                   "floor_value": 5e-9, "clip_count": 1, "sentinel": "unchanged_delegate_fields"}

    class SpySource:
        def __init__(self, seed):
            self.seed = seed

        def sample(self, mapping, rep):
            calls.append((mapping, deepcopy(mapping), rep, self.seed))
            return sample, diagnostics

    monkeypatch.setattr(floor, "OuSampleSource", SpySource)
    return SimpleNamespace(calls=calls, sample=sample, diagnostics=diagnostics)


@pytest.mark.parametrize("cell", CELLS)
def test_only_the_scientific_address_is_overridden(cell, delegate_spy):
    actual = cell.to_dict()
    before = deepcopy(actual)
    delegate_before = deepcopy(delegate_spy.diagnostics)
    sample, diagnostics = floor.PairedFloorSource(VALIDATION_SEED).sample(actual, 1)
    assert sample is delegate_spy.sample
    assert len(delegate_spy.calls) == 1
    mapping, snapshot, rep, seed = delegate_spy.calls[0]
    expected = dict(before, cell_id=Cell(n=1000, q=.5, cv=cell.cv, floor_fraction=1e-6).identity)
    assert mapping is not actual
    assert mapping == snapshot == expected
    assert mapping["floor_fraction"] == cell.floor_fraction
    assert rep == 1 and seed == VALIDATION_SEED
    assert actual == before
    assert delegate_spy.diagnostics == delegate_before
    assert diagnostics is not delegate_spy.diagnostics
    assert {key: diagnostics[key] for key in delegate_before} == delegate_before
    assert "cell_id" not in diagnostics
    assert diagnostics["actual_cell_id"] == cell.identity
    assert diagnostics["seed_anchor_cell_id"] == expected["cell_id"]
    assert diagnostics["seed_anchor_floor_fraction"] == 1e-6
    assert diagnostics["source_coupling_policy"] == POLICY


@pytest.mark.parametrize("bad_rep", [True, -1, .5, "0", None])
def test_invalid_replication_is_rejected_by_existing_source_before_random_addressing(bad_rep, monkeypatch):
    def forbidden_random_address(*args, **kwargs):
        raise AssertionError("invalid replication reached random addressing")
    monkeypatch.setattr(np.random, "SeedSequence", forbidden_random_address)
    with pytest.raises(ValueError):
        floor.PairedFloorSource(VALIDATION_SEED).sample(CELLS[0].to_dict(), bad_rep)


@pytest.mark.parametrize("bad_seed", [True, -1, .5, "0", None])
def test_invalid_root_seed_is_rejected(bad_seed):
    with pytest.raises(ValueError):
        floor.PairedFloorSource(bad_seed)


def test_invalid_cell_is_rejected_before_delegation(delegate_spy):
    cell = CELLS[0].to_dict()
    cell["cell_id"] = "f" * 64
    with pytest.raises(ValueError):
        floor.PairedFloorSource(VALIDATION_SEED).sample(cell, 0)
    assert delegate_spy.calls == []


def test_adapter_does_not_retry_or_replace_delegate_failure(monkeypatch):
    calls = []

    class Source:
        def __init__(self, seed):
            self.seed = seed

        def sample(self, cell, rep):
            calls.append((deepcopy(cell), rep, self.seed))
            raise ArithmeticError("fixed synthetic source failure")

    monkeypatch.setattr(floor, "OuSampleSource", Source)
    with pytest.raises(ArithmeticError, match="fixed synthetic source failure"):
        floor.PairedFloorSource(VALIDATION_SEED).sample(CELLS[-1].to_dict(), 0)
    assert len(calls) == 1
    assert calls[0][0]["floor_fraction"] == 1e-4


@pytest.mark.parametrize("key", ["cell_id", "actual_cell_id", "seed_anchor_cell_id",
                               "seed_anchor_floor_fraction", "source_coupling_policy",
                               "floor_fraction", "n", "true_theta", "rep"])
def test_upstream_diagnostics_cannot_overwrite_actual_physical_or_address_identity(key, delegate_spy):
    delegate_spy.diagnostics[key] = "synthetic_colliding_value"
    with pytest.raises(ValueError, match="collid"):
        floor.PairedFloorSource(VALIDATION_SEED).sample(CELLS[-3].to_dict(), 0)


def synthetic_services():
    class Service:
        def __init__(self, name):
            self.name = name
            self.estimator = SimpleNamespace(fit=lambda sample: PointFit(
                name, np.nan, np.nan, np.nan, np.nan, False, False,
                "synthetic_invalid_point", "identity integration fixture"))

        def infer(self, sample, point_fit):
            return failed_result(self.name, "synthetic_no_interval", "identity integration fixture", point_fit)

    return {name: Service(name) for name in ESTIMATORS}


def test_composed_source_records_preserve_actual_cell_identity_and_all_three_rows(delegate_spy):
    actual = CELLS[-3]
    rows = source_records(actual, 0, floor.PairedFloorSource(VALIDATION_SEED), synthetic_services())
    assert len(rows) == 3
    assert {row["estimator"] for row in rows} == set(ESTIMATORS)
    for row in rows:
        assert row["cell_id"] == row["actual_cell_id"] == actual.identity
        assert row["seed_anchor_cell_id"] != actual.identity
        assert row["seed_anchor_cell_id"] == replace(actual, floor_fraction=1e-6).identity
        assert row["floor_fraction"] == 1e-8
        assert row["source_coupling_policy"] == POLICY
        assert row["source_valid"] is True
        assert row["pipeline_error"] == ""
        assert row["times_sha256"] == array_digest(delegate_spy.sample.times)
        assert row["observations_sha256"] == array_digest(delegate_spy.sample.values)
    assert len(delegate_spy.calls) == 1


def test_source_instances_are_frozen_and_sampling_does_not_change_global_rng():
    source = floor.PairedFloorSource(VALIDATION_SEED)
    before = np.random.get_state()
    with pytest.raises((FrozenInstanceError, AttributeError)):
        source.seed = VALIDATION_SEED + 1
    sample, diagnostics = source.sample(CELLS[-3].to_dict(), 0)
    after = np.random.get_state()
    assert before[0] == after[0]
    np.testing.assert_array_equal(before[1], after[1])
    assert before[2:] == after[2:]
    assert sample.values.shape == sample.times.shape == (1000,)
    assert not sample.values.flags.writeable and not sample.times.flags.writeable
    assert diagnostics["source_seed"] == VALIDATION_SEED
