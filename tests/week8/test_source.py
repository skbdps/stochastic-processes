"""Constructed source checks and a few fixed development addresses; no full run."""
from dataclasses import replace
from fractions import Fraction
import json
import math

import numpy as np
import pytest

from ou_irregular.ou import simulate_ou
from ou_irregular.week4.simulation import array_digest
from ou_irregular.week8.design import E4Cell, VALIDATION_SEED, full_design
from ou_irregular.week8 import source as e4


class ScriptedRng:
    def __init__(self, uniforms, waits):
        self.uniforms, self.waits = list(uniforms), list(waits)
        self.scales = []
        self.uniform_calls = 0

    def random(self):
        self.uniform_calls += 1
        return self.uniforms.pop(0)

    def exponential(self, scale):
        self.scales.append(scale)
        return self.waits.pop(0)


def test_joint_race_retains_actual_arrival_phase_for_mmpp():
    cell = full_design().cells[1]
    rng = ScriptedRng([.1, .95, .1, .2], [.02, .03, .04])
    result = e4.generate_gaps(cell, rng, n_gaps=2)
    np.testing.assert_allclose(result.raw_gaps, [.05, .04], rtol=0, atol=1e-17)
    np.testing.assert_array_equal(result.start_phases, [0, 1])
    np.testing.assert_array_equal(result.end_phases, [1, 1])
    np.testing.assert_array_equal(result.switch_counts, [1, 0])
    assert result.initial_phase == 0
    assert rng.uniform_calls == 4 and not rng.uniforms and not rng.waits
    assert rng.scales == pytest.approx([1 / 1.4, 1 / 7, 1 / 7])


def test_iid_resets_palm_phase_at_every_gap_but_keeps_internal_switches():
    cell = full_design().cells[4]
    rng = ScriptedRng([.1, .95, .1, .1, .2], [.02, .03, .04])
    result = e4.generate_gaps(cell, rng, n_gaps=2)
    np.testing.assert_allclose(result.raw_gaps, [.05, .04], rtol=0, atol=1e-17)
    np.testing.assert_array_equal(result.start_phases, [0, 0])
    np.testing.assert_array_equal(result.end_phases, [1, 0])
    np.testing.assert_array_equal(result.switch_counts, [1, 0])
    assert rng.uniform_calls == 5
    assert rng.scales == pytest.approx([1 / 1.4, 1 / 7, 1 / 1.4])
    assert e4.gap_statistics(cell, result)["source_phase_reset_count"] == 2


def test_floor_is_posthoc_and_never_evolves_the_latent_state_through_padding():
    cell = full_design().cells[2]
    uniforms, waits = [.05, .1, .9, .2, .2], [1e-8, 1e-8, 1e-8, 3e-6]
    rng = ScriptedRng(uniforms, waits)
    result = e4.generate_gaps(cell, rng, n_gaps=3)
    other_rng = ScriptedRng(uniforms, waits)
    tiny_floor = e4.generate_gaps(replace(cell, floor_fraction=1e-9), other_rng, n_gaps=3)
    assert result.floor == 5e-7
    np.testing.assert_array_equal(result.raw_gaps, [1e-8, 2e-8, 3e-6])
    np.testing.assert_array_equal(result.transformed_gaps, [5e-7, 5e-7, 3e-6])
    for name in ("raw_gaps", "start_phases", "end_phases", "switch_counts"):
        np.testing.assert_array_equal(getattr(result, name), getattr(tiny_floor, name))
    assert rng.scales == other_rng.scales and rng.uniform_calls == other_rng.uniform_calls
    stats = e4.gap_statistics(cell, result)
    assert stats["clip_count"] == stats["raw_below_floor_count"] == stats["transformed_at_floor_count"] == 2
    assert stats["transformed_below_floor_count"] == 0
    assert stats["source_switch_count"] == 1


def test_regular_source_uses_no_schedule_draws_and_handles_analytics_separately():
    cell = full_design().cells[-1]
    rng = ScriptedRng([], [])
    result = e4.generate_gaps(cell, rng)
    np.testing.assert_array_equal(result.raw_gaps, np.full(999, .5))
    np.testing.assert_array_equal(result.times, np.arange(1000) * .5)
    assert result.initial_phase == -1 and not result.switch_counts.any()
    assert rng.uniform_calls == 0 and rng.scales == []
    assert e4.moments(cell) == {"mean": .5, "second_moment": .25, "variance": 0., "cv": 0.}
    assert e4.adjacent_product(cell) == .25
    assert (e4.cdf(cell, .49), e4.cdf(cell, .5), e4.survival(cell, .5)) == (0., 1., 0.)
    with pytest.raises(ValueError, match="no two-state"):
        e4.phase_type(cell)


@pytest.mark.parametrize("index,second,adjacent", [
    (0, Fraction(1, 2), Fraction(1, 4)),
    (1, Fraction(2, 3), Fraction(5, 16)),
    (2, Fraction(290, 219), Fraction(1, 4) + Fraction(1444, 5329)),
])
def test_independent_rational_raw_moment_and_resolvent_identities(index, second, adjacent):
    cell = full_design().cells[index]
    ph = e4.phase_type(cell)
    np.testing.assert_allclose(ph.pi @ ph.Q, np.zeros(2), atol=1e-15, rtol=0)
    np.testing.assert_allclose(ph.P @ np.ones(2), np.ones(2), atol=1e-15, rtol=0)
    np.testing.assert_allclose(ph.alpha @ ph.P, ph.alpha, atol=1e-15, rtol=0)
    np.testing.assert_allclose(ph.alpha @ ph.inverse_A, .5 * ph.pi, atol=1e-15, rtol=0)
    expected = e4.moments(cell)
    assert expected["mean"] == pytest.approx(.5, abs=2e-15)
    assert expected["second_moment"] == pytest.approx(float(second), rel=2e-14)
    assert e4.adjacent_product(cell) == pytest.approx(float(adjacent), rel=2e-14)
    assert e4.adjacent_product(replace(cell, mechanism="iid")) == pytest.approx(.25, abs=2e-15)
    beta = Fraction(60 * int(cell.ratio), int(cell.ratio)**2 + 70 * int(cell.ratio) + 25)
    covariance = float(beta**2) * .25 * (cell.cv**2 - 1) / 2
    assert e4.adjacent_product(cell, lag=2) == pytest.approx(.25 + covariance, rel=2e-14)
    with pytest.raises(ValueError):
        ph.Q.setflags(write=True)
    expected["mean"] = 999.
    assert e4.moments(cell)["mean"] == pytest.approx(.5, abs=2e-15)


def independent_two_exponential_law(cell):
    # Derive the two decay rates and weights from scalar trace/determinant and
    # initial hazard, independently of the producer's matrix exponentials.
    m, r = cell.q / cell.theta, cell.ratio
    switch_c, switch_b = 1 / (10 * m), 1 / (2 * m)
    arrival_c, arrival_b = 6 / ((r + 5) * m), 6 * r / ((r + 5) * m)
    trace = switch_c + arrival_c + switch_b + arrival_b
    determinant = arrival_c * arrival_b + switch_c * arrival_b + switch_b * arrival_c
    delta = math.sqrt(trace * trace - 4 * determinant)
    slow, fast = (trace - delta) / 2, (trace + delta) / 2
    hazard = (5 * arrival_c + r * arrival_b) / (r + 5)
    weight = (fast - hazard) / (fast - slow)
    return ((weight, slow), (1 - weight, fast))


@pytest.mark.parametrize("index", [1, 2, 4, 5])
def test_scalar_spectral_oracle_for_raw_and_floored_distribution(index):
    cell = full_design().cells[index]
    for fraction in (1e-6, .25):
        trial = replace(cell, floor_fraction=fraction, floor_cv_threshold=.5)
        f = fraction * .5
        terms = independent_two_exponential_law(trial)
        mean = f + sum(w * math.exp(-rate * f) / rate for w, rate in terms)
        second = f*f + 2 * sum(w * math.exp(-rate*f) * (f/rate + 1/rate**2) for w, rate in terms)
        result = e4.moments(trial, transformed=True)
        assert result["mean"] == pytest.approx(mean, rel=2e-14)
        assert result["second_moment"] == pytest.approx(second, rel=2e-14)
        for t in (0., 1e-15, f, .05, .25, .5, 1., 2.5):
            expected = sum(w * -math.expm1(-rate*t) for w, rate in terms)
            assert e4.cdf(trial, t) == pytest.approx(expected, abs=3e-15)
            transformed = 0. if t < f else expected
            assert e4.cdf(trial, t, transformed=True) == pytest.approx(transformed, abs=3e-15)
        assert e4.cdf(trial, np.nextafter(f, 0.), transformed=True) == 0.
        assert e4.cdf(trial, f, transformed=True) > 0.


def test_ratio_one_exponential_limit_and_cdf_endpoints():
    cell = full_design().cells[0]
    for t in (0., .05, .25, .5, 1., 2.5):
        assert e4.cdf(cell, t) == pytest.approx(-math.expm1(-2*t), abs=2e-15)
        assert e4.survival(cell, t) == pytest.approx(math.exp(-2*t), abs=2e-15)
    assert e4.cdf(cell, -math.inf) == 0. and e4.cdf(cell, math.inf) == 1.
    assert e4.survival(cell, -1.) == 1.


def test_stored_time_diagnostics_do_not_substitute_the_pre_cumsum_vector():
    cell = full_design().cells[0]
    result = e4.generate_gaps(cell, ScriptedRng([.1] * 6, [.1] * 5), n_gaps=5)
    stats = e4.gap_statistics(cell, result)
    observed = np.diff(result.times).astype(np.longdouble)
    assert stats["raw_cv"] == stats["transformed_cv"] == 0.
    assert stats["realized_cv"] > 0.
    assert stats["realized_mean_gap"] == float(observed.mean())
    assert stats["min_gap"] == float(observed.min()) < stats["transformed_min_gap"]
    assert stats["realized_span"] == result.times[-1]
    for array in (result.raw_gaps, result.transformed_gaps, result.times, result.start_phases):
        with pytest.raises(ValueError):
            array.setflags(write=True)


def test_nonpositive_stored_interval_is_not_repaired_or_redrawn():
    cell = full_design().cells[0]
    rng = ScriptedRng([.1, .1, .1], [1., 2**-54])
    with pytest.raises(ValueError, match="strictly positive"):
        e4.generate_gaps(cell, rng, n_gaps=2)
    assert rng.uniform_calls == 3 and len(rng.scales) == 2


@pytest.mark.parametrize("bad", [0., -1., math.nan, math.inf])
def test_invalid_wait_is_retained_as_failure_without_redraw(bad):
    rng = ScriptedRng([.1], [bad])
    with pytest.raises(ValueError, match="no redraw"):
        e4.generate_gaps(full_design().cells[2], rng, n_gaps=1)
    assert len(rng.scales) == 1 and rng.uniform_calls == 1


@pytest.mark.parametrize("bad", [-.1, 1., math.nan, math.inf])
def test_invalid_phase_uniform_rejected(bad):
    with pytest.raises(ValueError, match="uniform"):
        e4.generate_gaps(full_design().cells[0], ScriptedRng([bad], []), n_gaps=1)


@pytest.mark.parametrize("bad", [True, np.bool_(False), -1, 1., "1", None])
def test_invalid_seed_and_replication_rejected(bad):
    with pytest.raises(ValueError):
        e4.E4Source(bad)
    with pytest.raises(ValueError):
        e4.E4Source(VALIDATION_SEED).gaps(full_design().cells[0].to_dict(), bad)


def test_strict_cell_mapping_and_fixed_source_scope():
    source = e4.E4Source(VALIDATION_SEED)
    mapping = full_design().cells[0].to_dict()
    for bad in ({**mapping, "source_namespace": 4}, {**mapping, "extra": 1},
                replace(full_design().cells[0], n=4).to_dict()):
        with pytest.raises(ValueError):
            source.gaps(bad, 0)
    with pytest.raises(ValueError):
        source.gaps(full_design().cells[0], 0)


def test_development_address_replay_and_independent_ou_child_agree_exactly():
    cell = full_design().cells[2]
    source = e4.E4Source(VALIDATION_SEED)
    sample, diagnostics = source.sample(cell.to_dict(), 0)
    gaps = source.gaps(cell.to_dict(), 0)
    words = tuple(int(cell.identity[i:i+8], 16) for i in range(0, 64, 8))
    sequence = np.random.SeedSequence(21261011, spawn_key=(8, *words, 0))
    children = sequence.spawn(2)
    expected_values = simulate_ou(1., 0., .5, gaps.times, np.random.default_rng(children[1]))
    np.testing.assert_array_equal(sample.values, expected_values)
    np.testing.assert_array_equal(sample.times, gaps.times)
    assert set(diagnostics) == e4.SOURCE_FIELDS
    assert json.loads(diagnostics["source_spawn_key"]) == [8, *words, 0]
    assert json.loads(diagnostics["gap_spawn_key"]) == [8, *words, 0, 0]
    assert json.loads(diagnostics["path_spawn_key"]) == [8, *words, 0, 1]
    assert diagnostics["times_sha256"] == array_digest(sample.times)
    assert diagnostics["observations_sha256"] == array_digest(sample.values)
    repeat, repeated = source.sample(cell.to_dict(), 0)
    assert diagnostics == repeated
    np.testing.assert_array_equal(repeat.values, sample.values)
    assert not sample.times.flags.writeable and not sample.values.flags.writeable
    assert not set(diagnostics) & (set(cell.to_dict()) | {"rep", "point_valid", "pipeline_error"})


def test_mechanism_and_replication_streams_are_distinct_at_development_seed():
    source = e4.E4Source(VALIDATION_SEED)
    first, iid = full_design().cells[0], full_design().cells[3]
    addresses = [source.address(cell.to_dict(), rep) for cell, rep in ((first, 0), (first, 1), (iid, 0))]
    assert len({row["gap_spawn_key"] for row in addresses}) == 3
    hashes = {array_digest(source.gaps(cell.to_dict(), rep).raw_gaps)
              for cell, rep in ((first, 0), (first, 1), (iid, 0))}
    assert len(hashes) == 3


def test_schedule_only_api_never_calls_ou_simulator(monkeypatch):
    monkeypatch.setattr(e4, "simulate_ou", lambda *a, **k: pytest.fail("schedule-only call drew OU"))
    e4.E4Source(VALIDATION_SEED).gaps(full_design().cells[-1].to_dict(), 0)


def test_scientific_or_fit_diagnostic_collision_is_rejected(monkeypatch):
    original = e4.gap_statistics
    monkeypatch.setattr(e4, "gap_statistics", lambda cell, result: {**original(cell, result), "rep": 0})
    monkeypatch.setattr(e4, "simulate_ou", lambda theta, mu, sigma, times, rng: np.zeros(len(times)))
    with pytest.raises(ValueError, match="collide"):
        e4.E4Source(VALIDATION_SEED).sample(full_design().cells[-1].to_dict(), 0)


@pytest.mark.parametrize("call", [
    lambda c: e4.cdf(c, math.nan), lambda c: e4.cdf(c, True),
    lambda c: e4.moments(c, transformed=1), lambda c: e4.adjacent_product(c, lag=0),
    lambda c: e4.adjacent_product(c, lag=True),
    lambda c: e4.generate_gaps(c, ScriptedRng([], []), n_gaps=True),
])
def test_analytical_and_constructed_api_domains(call):
    with pytest.raises(ValueError):
        call(full_design().cells[0])
