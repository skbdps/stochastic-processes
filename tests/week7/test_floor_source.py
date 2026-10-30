"""Paired-source integration and identity tests, with validation seed only."""
from copy import deepcopy
from dataclasses import FrozenInstanceError
import json
from types import MappingProxyType

import pytest

from ou_irregular.experiments.design import ESTIMATORS
from ou_irregular.experiments.execution import registered_wald, source_records
from ou_irregular.experiments.provenance import canonical
from ou_irregular.week4.simulation import OuSampleSource
from ou_irregular.week7.floor_source import (
    ANCHOR_FLOOR_FRACTION, COUPLING_POLICY, PairedFloorSource, anchor_cell, validate_cell)
from ou_irregular.week7.selection import floor_design


VALIDATION_SEED = 20261009 + 1_000_000


@pytest.mark.parametrize("cell", floor_design().cells, ids=lambda cell: f"CV={cell.cv}-epsilon={cell.floor_fraction}")
def test_frozen_actual_grid_and_anchor_identity(cell):
    actual = validate_cell(cell.to_dict())
    assert actual == cell
    assert actual.to_dict() == cell.to_dict()
    anchor = anchor_cell(cell)
    assert anchor == anchor_cell(MappingProxyType(cell.to_dict()))
    assert anchor.floor_fraction == ANCHOR_FLOOR_FRACTION
    actual_fields = {k: v for k, v in cell.to_dict().items() if k not in {"cell_id", "floor_fraction"}}
    anchor_fields = {k: v for k, v in anchor.to_dict().items() if k not in {"cell_id", "floor_fraction"}}
    assert actual_fields == anchor_fields


def test_delegate_receives_only_address_override_without_mutating_input(monkeypatch):
    cell = next(c for c in floor_design().cells if c.cv == 2 and c.floor_fraction == 1e-4)
    supplied = cell.to_dict()
    before = deepcopy(supplied)
    calls = []
    sentinel_sample = object()
    returned_diagnostics = {"clip_count": 7, "floor_value": cell.floor_fraction * cell.q / cell.theta}

    def sample(source, addressed, rep):
        calls.append((source.seed, deepcopy(addressed), rep))
        return sentinel_sample, returned_diagnostics

    monkeypatch.setattr(OuSampleSource, "sample", sample)
    source = PairedFloorSource(VALIDATION_SEED)
    result, diagnostics = source.sample(MappingProxyType(supplied), 3)
    assert result is sentinel_sample
    assert supplied == before
    assert len(calls) == 1
    seed, addressed, rep = calls[0]
    assert seed == VALIDATION_SEED and rep == 3
    assert {key for key in before if before[key] != addressed[key]} == {"cell_id"}
    assert addressed["cell_id"] == anchor_cell(cell).identity
    assert addressed["floor_fraction"] == 1e-4
    assert diagnostics["actual_cell_id"] == cell.identity
    assert diagnostics["seed_anchor_cell_id"] == anchor_cell(cell).identity
    assert diagnostics["seed_anchor_floor_fraction"] == 1e-6
    assert diagnostics["source_coupling_policy"] == COUPLING_POLICY
    assert "cell_id" not in diagnostics
    assert returned_diagnostics == {"clip_count": 7, "floor_value": cell.floor_fraction * cell.q / cell.theta}
    with pytest.raises(FrozenInstanceError):
        source.seed = 0


@pytest.mark.parametrize("collision", ["cell_id", "floor_fraction", "n", "rep", "actual_cell_id",
                                       "seed_anchor_cell_id", "seed_anchor_floor_fraction", "source_coupling_policy"])
def test_upstream_diagnostic_collisions_cannot_overwrite_actual_identity(monkeypatch, collision):
    cell = floor_design().cells[0]
    monkeypatch.setattr(OuSampleSource, "sample", lambda *args: (object(), {collision: "unexpected"}))
    with pytest.raises(ValueError, match="collide"):
        PairedFloorSource(VALIDATION_SEED).sample(cell.to_dict(), 0)


def test_three_regular_variants_have_identical_full_fits_but_distinct_actual_rows():
    # This is the exact negative control at one separate-namespace replication,
    # not three independent datasets or a coverage study.
    cells = [cell for cell in floor_design().cells if cell.cv == 0]
    source = PairedFloorSource(VALIDATION_SEED)
    records = [source_records(cell, 0, source, {name: registered_wald(name) for name in ESTIMATORS})
               for cell in cells]
    assert len({cell.identity for cell in cells}) == 3
    for cell, rows in zip(cells, records):
        assert len(rows) == 3
        for row in rows:
            assert row["cell_id"] == row["actual_cell_id"] == cell.identity
            assert row["floor_fraction"] == cell.floor_fraction
            assert row["seed_anchor_cell_id"] == anchor_cell(cell).identity
            assert row["source_coupling_policy"] == COUPLING_POLICY
            assert row["source_valid"] is True
            assert not row["pipeline_error"]
            assert row["clip_count"] == 0 and row["floor_applied"] is False
    stripped = [[{key: value for key, value in row.items()
                  if key not in {"cell_id", "actual_cell_id", "floor_fraction"}}
                 for row in rows] for rows in records]
    assert canonical(stripped[0]) == canonical(stripped[1]) == canonical(stripped[2])
    assert all(json.loads(row["source_spawn_key"])[0] == 4 for rows in records for row in rows)
