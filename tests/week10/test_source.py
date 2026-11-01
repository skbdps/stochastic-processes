from dataclasses import replace

import numpy as np
import pytest

from ou_irregular.week8.source import generate_gaps
from ou_irregular.week10.design import (HorizonCell, full_design, parse_cell,
                                        parse_design, validation_design)
from ou_irregular.week10.source import HorizonSource, generate_schedule


def test_exact_horizon_design_and_seed_separation():
    full, validation = full_design(), validation_design()
    assert len(full.cells) == 7 and full.replications == 2000
    assert validation.seed == 21261101 and validation.replications == 4
    assert full.seed == 20261101
    assert parse_design(full.to_dict()) == full
    assert len({c.identity for c in full.cells}) == 7
    for cell in full.cells:
        assert cell.horizon == 499.5
        assert parse_cell(cell.to_dict()) == cell
        assert cell.identity != cell.parent_cell().identity
        assert replace(cell, horizon=500.).identity != cell.identity
        broken = cell.to_dict()
        broken["nominal_observations"] += 1
        with pytest.raises(ValueError):
            parse_cell(broken)


@pytest.mark.parametrize("cell", full_design().cells)
def test_uncensored_parent_sequence_is_exact_unchanged_e4_prefix(cell):
    horizon = generate_schedule(cell, np.random.default_rng(123))
    parent = generate_gaps(cell.parent_cell(), np.random.default_rng(123),
                           n_gaps=len(horizon.raw_parent_gaps))
    assert np.array_equal(horizon.raw_parent_gaps, parent.raw_gaps)
    assert np.array_equal(horizon.transformed_parent_gaps, parent.transformed_gaps)
    assert np.array_equal(horizon.start_phases, parent.start_phases)
    assert np.array_equal(horizon.end_phases, parent.end_phases)
    assert np.array_equal(horizon.switch_counts, parent.switch_counts)
    assert np.array_equal(horizon.times[:-1], parent.times[:-1])
    assert horizon.times[-1] == 499.5
    assert parent.times[-1] >= 499.5
    assert np.all(np.diff(horizon.times) > 0)
    with pytest.raises(ValueError):
        horizon.times.setflags(write=True)


def test_regular_endpoint_is_not_duplicated():
    cell = full_design().cells[-1]
    result = generate_schedule(cell, np.random.default_rng(0))
    assert result.terminal_is_arrival
    assert len(result.times) == 1000
    assert np.array_equal(result.times, np.arange(1000) * .5)


def test_source_replays_independent_ou_and_schedule_children():
    source = HorizonSource(21261101)
    for cell in full_design().cells:
        sample, diagnostic = source.sample(cell.to_dict(), 0)
        again, second = source.sample(cell.to_dict(), 0)
        assert diagnostic == second
        assert np.array_equal(sample.values, again.values)
        assert np.array_equal(source.schedule(cell, 0).times, sample.times)
        assert diagnostic["realized_span"] == cell.horizon
        assert diagnostic["actual_n"] == len(sample.values)
        assert diagnostic["actual_transitions"] == len(sample.values) - 1
        assert diagnostic["gap_spawn_key"] != diagnostic["path_spawn_key"]
    assert len({source.address(c, 0)["source_spawn_key"] for c in full_design().cells}) == 7
    assert source.address(full_design().cells[0], 0) != source.address(full_design().cells[0], 1)


class ConstructedRace:
    def __init__(self):
        self.waits = iter([1e-8, 1.])

    def random(self):
        return 0.

    def exponential(self, scale):
        return next(self.waits)


def test_terminal_gap_is_not_floored_or_clock_renormalized():
    cell = HorizonCell(n=1000, q=.5, mechanism="mmpp", ratio=20., horizon=7e-7)
    schedule = generate_schedule(cell, ConstructedRace())
    assert np.array_equal(schedule.raw_parent_gaps, [1e-8, 1.])
    assert np.array_equal(schedule.transformed_parent_gaps, [5e-7, 1.])
    assert np.array_equal(schedule.times, [0., 5e-7, 7e-7])
    assert schedule.times[-1] - schedule.times[-2] < schedule.floor
    assert not schedule.terminal_is_arrival
    assert np.array_equal(schedule.switch_counts, [0, 0])


def test_r1_count_has_poisson_mean_and_variance():
    # Independent development sequences, not a fitted research result.
    # r=1 has a Poisson arrival process even though the hidden phase switches.
    source = HorizonSource(22261101)
    for cell in (full_design().cells[0], full_design().cells[3]):
        counts = np.array([len(source.schedule(cell, rep).times) - 2 for rep in range(256)])
        expected = cell.horizon / (cell.q / cell.theta)
        assert abs(counts.mean() - expected) <= 6 * np.sqrt(expected / len(counts))
        variance_se = np.sqrt((expected + 2 * expected**2) / (len(counts) - 1))
        assert abs(counts.var(ddof=1) - expected) <= 6 * variance_se


def test_source_rejects_offgrid_or_invalid_addresses():
    source = HorizonSource(20261101)
    cell = full_design().cells[0]
    with pytest.raises(ValueError):
        source.sample(replace(cell, horizon=1.).to_dict(), 0)
    with pytest.raises(ValueError):
        source.sample(cell.to_dict(), True)
    with pytest.raises(ValueError):
        HorizonSource(-1)
