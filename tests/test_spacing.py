"""Spacing checks against distributional identities and literal floor semantics."""

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.stats import gamma

from ou_irregular.spacing import SpacingSample, generate_times


def test_regular_schedule_reports_exact_design_moments():
    sample = generate_times(100, 0.1, 0, np.random.default_rng(1))
    assert isinstance(sample, SpacingSample)
    np.testing.assert_array_equal(sample.times, np.arange(100) * 0.1)
    stats = sample.diagnostics
    assert stats["family"] == "equidistant"
    assert stats["raw_mean_gap"] == stats["realized_mean_gap"] == 0.1
    assert stats["raw_cv"] == stats["realized_cv"] == 0.0
    assert stats["clip_count"] == stats["clip_fraction"] == 0
    assert not stats["floor_applied"]
    assert stats["expected_mean_gap_after_floor"] == 0.1
    assert stats["expected_clip_fraction"] == 0.0


@pytest.mark.parametrize("cv,family", [(0.5, "gamma"), (1.0, "exponential"), (1.4, "gamma")])
def test_unfloored_gamma_moments(cv, family):
    sample = generate_times(200_001, 0.7, cv, np.random.default_rng(9))
    stats = sample.diagnostics
    assert stats["family"] == family
    assert not stats["floor_applied"]
    assert abs(stats["realized_mean_gap"] / 0.7 - 1) < 0.015
    assert abs(stats["realized_cv"] / cv - 1) < 0.025
    assert stats["expected_mean_gap_after_floor"] == 0.7
    assert stats["expected_clip_fraction"] == 0.0


def test_high_cv_floor_is_literal_and_never_rescales():
    n, mean, cv, epsilon, seed = 1001, 0.3, 2.0, 0.1, 24
    raw = np.random.default_rng(seed).gamma(cv ** -2, mean * cv ** 2, size=n - 1)
    floor = mean * epsilon
    expected = np.maximum(raw, floor)
    sample = generate_times(n, mean, cv, np.random.default_rng(seed), floor_fraction=epsilon)
    np.testing.assert_array_equal(sample.times, np.r_[0, np.cumsum(expected)])
    stats = sample.diagnostics
    assert stats["floor_applied"]
    assert stats["floor_value"] == floor
    assert stats["clip_count"] == np.count_nonzero(raw < floor)
    assert stats["clip_fraction"] == stats["clip_count"] / (n - 1)
    assert stats["raw_mean_gap"] == pytest.approx(raw.mean())
    assert stats["raw_cv"] == pytest.approx(raw.std(ddof=0) / raw.mean())
    assert stats["realized_mean_gap"] == pytest.approx(expected.mean())
    assert stats["realized_cv"] == pytest.approx(expected.std(ddof=0) / expected.mean())
    assert stats["realized_span"] == sample.times[-1]
    assert stats["min_gap"] >= floor * (1 - 1e-10)
    assert stats["realized_mean_gap"] > stats["raw_mean_gap"]
    assert not np.isclose(stats["realized_mean_gap"], mean)


def test_floor_analytic_mean_matches_independent_tail_integral_and_sample():
    # For any nonnegative D, E[max(D,a)] = a + integral_a^infinity P(D>t) dt.
    # This independently checks the incomplete-Gamma identity in production.
    mean, cv, epsilon = 2.0, 2.0, 0.3
    sample = generate_times(200_001, mean, cv, np.random.default_rng(123), epsilon)
    distribution = gamma(a=cv ** -2, scale=mean * cv ** 2)
    floor = epsilon * mean
    integral, _ = quad(distribution.sf, floor, np.inf, epsabs=1e-10)
    expected = floor + integral
    stats = sample.diagnostics
    assert stats["expected_mean_gap_after_floor"] == pytest.approx(expected, rel=1e-10)
    assert stats["expected_clip_fraction"] == pytest.approx(distribution.cdf(floor))
    assert abs(stats["realized_mean_gap"] / expected - 1) < 0.015
    assert abs(stats["clip_fraction"] - distribution.cdf(floor)) < 0.005
    assert expected > mean


def test_floor_is_only_applied_at_or_above_threshold():
    lower = generate_times(100, 1, 1, np.random.default_rng(4), floor_fraction=100)
    threshold = generate_times(100, 1, 1.5, np.random.default_rng(4), floor_fraction=100)
    assert not lower.diagnostics["floor_applied"]
    assert lower.diagnostics["clip_count"] == 0
    assert threshold.diagnostics["floor_applied"]
    assert threshold.diagnostics["clip_count"] == 99
    np.testing.assert_array_equal(np.diff(threshold.times), np.full(99, 100.0))


def test_none_floor_is_allowed_only_when_inactive():
    generate_times(10, 1, 1, np.random.default_rng(4), floor_fraction=None)
    with pytest.raises(ValueError, match="positive floor_fraction"):
        generate_times(10, 1, 2, np.random.default_rng(4), floor_fraction=None)


@pytest.mark.parametrize("n", [0, 1, 2, 3.5, "3", True, np.bool_(False)])
def test_invalid_n(n):
    with pytest.raises(ValueError, match="n must"):
        generate_times(n, 1, 1, np.random.default_rng(1))


@pytest.mark.parametrize("name,value", [
    ("mean_gap", 0), ("mean_gap", -1), ("mean_gap", np.inf), ("mean_gap", np.nan),
    ("cv", -1), ("cv", np.inf), ("cv", np.nan), ("cv", True),
    ("floor_fraction", 0), ("floor_fraction", -1), ("floor_fraction", np.nan),
    ("floor_cv_threshold", 0), ("floor_cv_threshold", np.inf),
])
def test_invalid_distribution_inputs(name, value):
    kwargs = dict(n=10, mean_gap=1, cv=2, rng=np.random.default_rng(1))
    kwargs[name] = value
    with pytest.raises(ValueError, match=name):
        generate_times(**kwargs)


@pytest.mark.parametrize("cv", [1e-200, 1e200])
def test_unrepresentable_gamma_parameters_raise(cv):
    with pytest.raises(ValueError, match="Gamma shape and scale"):
        generate_times(10, 1, cv, np.random.default_rng(1))


def test_cumulative_roundoff_must_not_silently_merge_timestamps():
    class FixedGenerator:
        def gamma(self, *args, **kwargs):
            return np.array([1e16, 1.0])

    with pytest.raises(ValueError, match="floating-point accumulation"):
        generate_times(3, 1, 1, FixedGenerator())


def test_time_overflow_raises():
    with pytest.raises(ValueError, match="finite and strictly increasing"):
        generate_times(5, 1e308, 0, np.random.default_rng(1))
