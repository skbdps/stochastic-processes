"""Port-level bootstrap checks independent of Week 3 optimizers and truths."""
from dataclasses import FrozenInstanceError, replace

import numpy as np
import pytest

from ou_irregular.inference.bootstrap import ParametricBootstrap
from ou_irregular.inference.contracts import PointFit, Sample


def fit(parameters=(1., 0., .5), *, name="exact", valid=True,
        status="valid_interior", success=True):
    return PointFit(name, *parameters, 4., valid, success, status,
                    "test fit", {"mutable_source": [1, 2]})


@pytest.fixture
def sample():
    return Sample([.4, -.1, .2, .05], [0., .3, .8, 1.2])


class RecordingSimulator:
    """The second observation is a known deterministic or random marker."""
    def __init__(self, markers=None):
        self.markers = markers
        self.calls = []

    def simulate(self, parameters, times, rng, *, x0):
        self.calls.append((parameters.copy(), times.copy(), x0))
        marker = rng.uniform(1, 2) if self.markers is None else self.markers[len(self.calls) - 1]
        return Sample([x0, marker, .2, .05], times)


class MarkerEstimator:
    def __init__(self, name="exact", invalid=(), errors=()):
        self.name = name
        self.invalid = set(invalid)
        self.errors = set(errors)
        self.calls = 0

    def fit(self, sample):
        self.calls += 1
        if self.calls in self.errors:
            raise RuntimeError("controlled failed fit")
        z = sample.values[1]
        if self.calls in self.invalid:
            return fit((z, -z, 2*z), name=self.name, valid=False,
                       status="boundary", success=True)
        return fit((z, -z, 2*z), name=self.name)


def test_percentile_endpoints_and_sample_sd_have_hand_computable_values(sample):
    engine = ParametricBootstrap(MarkerEstimator(), RecordingSimulator([1, 2, 3, 4]),
                                 draws=4, level=.5)
    result = engine.run(sample, fit(), seed=21)
    # Linear .25/.75 quantiles of 1,2,3,4 are 1.75 and 3.25;
    # sign reversal changes which tail applies to mu.
    np.testing.assert_allclose(result.lower, [1.75, -3.25, 3.5])
    np.testing.assert_allclose(result.upper, [3.25, -1.75, 6.5])
    np.testing.assert_allclose(result.standard_errors,
                               np.sqrt(5 / 3) * np.array([1, 1, 2]))
    assert result.valid and result.status == "valid"
    assert result.planned_count == result.attempted_count == result.valid_count == 4
    assert result.invalid_count == 0
    assert result.diagnostics["quantile_method"] == "linear"


def test_fixes_observed_timestamps_and_initial_value_at_plugin_parameters(sample):
    simulator = RecordingSimulator([1, 2, 3])
    parameters = (2., 1., 3.)
    result = ParametricBootstrap(MarkerEstimator(), simulator, draws=3).run(sample, fit(parameters))
    assert result.valid
    for simulated_parameters, simulated_times, x0 in simulator.calls:
        np.testing.assert_array_equal(simulated_parameters, parameters)
        np.testing.assert_array_equal(simulated_times, sample.times)
        assert x0 == sample.values[0]
    assert result.diagnostics["condition_on_times"] is True
    assert result.diagnostics["condition_on_x0"] is True


def test_seed_reproduces_draws_and_extension_preserves_prefix(sample):
    def run(draws, seed):
        return ParametricBootstrap(MarkerEstimator(), RecordingSimulator(), draws=draws).run(sample, fit(), seed=seed)
    first, again, longer = run(3, 581), run(3, 581), run(5, 581)
    assert first.draw_records == again.draw_records == longer.draw_records[:3]
    assert first.draw_records != run(3, 582).draw_records
    assert [row["seed_spawn_key"] for row in longer.draw_records] == [(0,), (1,), (2,), (3,), (4,)]
    assert all(row["seed_entropy"] == 581 for row in longer.draw_records)
    assert len({row["theta"] for row in longer.draw_records}) == 5


def test_invalid_draws_are_retained_and_gate_withholds_intervals_without_replacement(sample):
    estimator = MarkerEstimator(invalid={2}, errors={4})
    simulator = RecordingSimulator([1, 2, 3, 4, 5])
    result = ParametricBootstrap(estimator, simulator, draws=5, min_valid_fraction=.8).run(sample, fit())
    assert result.status == "insufficient_valid_draws" and not result.valid
    assert result.planned_count == result.attempted_count == estimator.calls == len(simulator.calls) == 5
    assert result.valid_count == 3 and result.invalid_count == 2
    assert result.diagnostics["required_valid_count"] == 4
    assert result.diagnostics["valid_fraction"] == .6
    assert np.isnan(result.lower).all() and np.isnan(result.upper).all()
    assert np.isfinite(result.standard_errors).all()  # descriptive valid-draw SD
    assert result.draw_records[1]["status"] == "boundary"
    assert result.draw_records[1]["optimizer_success"] is True
    assert result.draw_records[1]["point_valid"] is False
    assert result.draw_records[3]["reason"] == "RuntimeError: controlled failed fit"


def test_threshold_is_inclusive_and_uses_ceiling(sample):
    result = ParametricBootstrap(MarkerEstimator(invalid={20}), RecordingSimulator(), draws=20).run(sample, fit())
    assert result.valid and result.valid_count == result.diagnostics["required_valid_count"] == 19
    result = ParametricBootstrap(MarkerEstimator(invalid={3}), RecordingSimulator(),
                                 draws=3, min_valid_fraction=.7).run(sample, fit())
    assert result.valid_count == 2 and result.diagnostics["required_valid_count"] == 3
    assert not result.valid


@pytest.mark.parametrize("draws,invalid", [(1, ()), (3, (1, 2)), (3, (1, 2, 3))])
def test_at_least_two_valid_draws_needed_even_with_low_fraction(sample, draws, invalid):
    result = ParametricBootstrap(MarkerEstimator(invalid=invalid), RecordingSimulator(),
                                 draws=draws, min_valid_fraction=.1).run(sample, fit())
    assert not result.valid and result.valid_count < 2
    assert np.isnan(result.lower).all() and np.isnan(result.standard_errors).all()


@pytest.mark.parametrize("point", [
    fit(valid=False, status="boundary"),
    fit(valid=True, status="boundary"),
    fit(valid=True, status="boundary_suspected"),
    fit(valid=False, status="degenerate"),
    fit(success=False),
    fit((-1., 0., .5)),
    fit((1., 0., -1.)),
    fit((0., 0., .5)),
    fit((1., np.nan, .5)),
    replace(fit(), nll=np.inf),
])
def test_ineligible_original_fit_is_skipped_with_no_simulation_or_refit(sample, point):
    estimator, simulator = MarkerEstimator(), RecordingSimulator()
    result = ParametricBootstrap(estimator, simulator, draws=500).run(sample, point)
    assert result.status == "skipped_invalid_point" and not result.valid
    assert result.planned_count == 500 and result.attempted_count == result.valid_count == result.invalid_count == 0
    assert not simulator.calls and estimator.calls == 0
    assert np.isnan(result.lower).all()


@pytest.mark.parametrize("kwargs", [
    {"draws": True}, {"draws": np.bool_(False)}, {"draws": 0}, {"draws": -1}, {"draws": 2.0},
    {"level": True}, {"level": 0}, {"level": 1}, {"level": np.nan}, {"level": np.inf},
    {"min_valid_fraction": False}, {"min_valid_fraction": 0}, {"min_valid_fraction": 1.1},
    {"min_valid_fraction": "0.95"},
])
def test_invalid_configuration_is_rejected(kwargs):
    with pytest.raises(ValueError):
        ParametricBootstrap(MarkerEstimator(), RecordingSimulator(), **kwargs)


@pytest.mark.parametrize("seed", [True, np.bool_(True), -1, 1., "1", None])
def test_invalid_seed_is_rejected_before_any_work(sample, seed):
    simulator = RecordingSimulator()
    with pytest.raises(ValueError, match="seed"):
        ParametricBootstrap(MarkerEstimator(), simulator).run(sample, fit(), seed=seed)
    assert simulator.calls == []


def test_numpy_integer_counts_and_seed_are_accepted(sample):
    result = ParametricBootstrap(MarkerEstimator(), RecordingSimulator(), draws=np.int64(2)).run(sample, fit(), seed=np.uint64(5))
    assert result.valid


def test_invalid_dependency_sample_or_point_type_rejected(sample):
    with pytest.raises(ValueError, match="estimator"):
        ParametricBootstrap(object(), RecordingSimulator())
    with pytest.raises(ValueError, match="simulator"):
        ParametricBootstrap(MarkerEstimator(), object())
    engine = ParametricBootstrap(MarkerEstimator(), RecordingSimulator(), draws=2)
    with pytest.raises(ValueError, match="Sample"):
        engine.run(object(), fit())
    with pytest.raises(ValueError, match="PointFit"):
        engine.run(sample, {})
    with pytest.raises(ValueError, match="match"):
        engine.run(sample, fit(name="euler"))


@pytest.mark.parametrize("changed", ["times", "x0"])
def test_simulator_violating_conditioning_is_recorded_as_failure(sample, changed):
    class BrokenSimulator:
        def simulate(self, parameters, times, rng, *, x0):
            if changed == "times":
                times = times * 2
            else:
                x0 += 1
            return Sample([x0, .1, .2, .3], times)
    estimator = MarkerEstimator()
    result = ParametricBootstrap(estimator, BrokenSimulator(), draws=2).run(sample, fit())
    assert not result.valid and result.invalid_count == 2 and estimator.calls == 0
    assert all("preserve" in row["reason"] for row in result.draw_records)


@pytest.mark.parametrize("name", ["pfml", "euler"])
def test_approximate_estimators_are_explicitly_diagnostic(sample, name):
    result = ParametricBootstrap(MarkerEstimator(name), RecordingSimulator(), draws=2).run(sample, fit(name=name))
    assert result.valid
    assert result.diagnostics["interpretation"] == "diagnostic_model_based_plugin_bootstrap"
    assert "not a robust interval" in result.reason
    assert result.diagnostics["coverage_guarantee"] is False


def test_result_is_an_immutable_snapshot(sample):
    point = fit()
    result = ParametricBootstrap(MarkerEstimator(), RecordingSimulator(), draws=2).run(sample, point)
    point.diagnostics["mutable_source"].append(3)
    assert result.point_fit.diagnostics["mutable_source"] == (1, 2)
    with pytest.raises(FrozenInstanceError):
        result.valid = False
    with pytest.raises(ValueError):
        result.lower[0] = 1
    with pytest.raises(TypeError):
        result.draw_records[0]["theta"] = 1
    with pytest.raises(TypeError):
        result.diagnostics["level"] = .8


def test_optional_original_fit_is_called_once_before_bootstrap(sample):
    class Estimator(MarkerEstimator):
        def fit(self, value):
            if self.calls == 0:
                self.calls += 1
                return fit()
            return super().fit(value)
    estimator = Estimator()
    result = ParametricBootstrap(estimator, RecordingSimulator(), draws=2).run(sample)
    assert result.valid and estimator.calls == 3


def test_original_fit_exception_is_visible_and_skips_draws(sample):
    estimator, simulator = MarkerEstimator(errors={1}), RecordingSimulator()
    result = ParametricBootstrap(estimator, simulator, draws=2).run(sample)
    assert result.status == "skipped_invalid_point" and result.attempted_count == 0
    assert "controlled failed fit" in result.reason and not simulator.calls
