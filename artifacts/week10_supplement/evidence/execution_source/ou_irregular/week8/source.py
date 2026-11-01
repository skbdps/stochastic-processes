"""Event-stationary E4 schedules and unchanged exact OU simulation.

E4Source(seed).sample(mapping, rep) returns (Sample, scalar diagnostics).
E4Source.gaps(mapping, rep) returns an inspectable immutable GapResult without
OU draws; .address(mapping, rep) returns its deterministic address fields.
These methods accept only the seven published scientific cells. Low-level
generate_gaps(cell, rng, n_gaps=None) also supports strictly parsed off-grid
cells for constructed tests. It implements the joint arrival/switch race.

phase_type, survival, cdf, moments, adjacent_product and floor_value expose
the analytical law. gap_statistics separates raw gaps, transformed gaps before
cumsum, and the actual differences of stored times. SOURCE_FIELDS is the
complete diagnostic schema; it excludes scientific and fit-record fields.
"""
from dataclasses import dataclass
from functools import lru_cache
import hashlib
import math

import numpy as np
from scipy.linalg import expm

from ..inference.contracts import Sample
from ..ou import simulate_ou
from ..week4.config import canonical_json
from ..week4.simulation import array_digest
from .design import E4Cell, FILES as DESIGN_FILES, SOURCE_NAMESPACE, full_design, parse_cell


FILES = tuple(sorted(set(DESIGN_FILES) | {
    "ou_irregular/week8/source.py", "docs/week8_source.md",
}))
CDF_GRID = ((.1, "0p1"), (.5, "0p5"), (1., "1"), (2., "2"), (5., "5"))
_PLANNED_IDS = frozenset(c.identity for c in full_design().cells)
ADDRESS_FIELDS = frozenset({"source_seed", "source_spawn_key", "gap_spawn_key", "path_spawn_key"})
GAP_FIELDS = frozenset({
    "family", "target_mean_gap", "requested_cv", "floor_applied", "floor_value",
    "expected_mean_gap_after_floor", "expected_clip_fraction", "clip_count", "clip_fraction",
    "raw_below_floor_count", "transformed_below_floor_count", "transformed_at_floor_count",
    "realized_mean_gap", "realized_cv", "realized_span", "min_gap",
    "raw_mean_gap", "raw_cv", "raw_second_moment", "raw_adjacent_product", "raw_min_gap", "raw_span",
    "transformed_mean_gap", "transformed_cv", "transformed_second_moment", "transformed_adjacent_product",
    "transformed_min_gap", "transformed_span", "raw_gaps_sha256", "transformed_gaps_sha256",
    "phase_trace_sha256", "source_initial_phase", "source_switch_count", "source_max_switches_per_gap",
    "source_phase_reset_count", "source_calm_start_count", "source_burst_start_count",
    "source_calm_end_count", "source_burst_end_count",
    "expected_raw_mean_gap", "expected_raw_second_moment", "expected_raw_cv", "expected_raw_adjacent_product",
    "expected_transformed_mean_gap", "expected_transformed_second_moment", "expected_transformed_cv",
    "expected_transformed_adjacent_product",
    *(f"{kind}_cdf_t_over_m_{label}" for kind in ("raw", "transformed") for _, label in CDF_GRID),
})
SOURCE_FIELDS = ADDRESS_FIELDS | GAP_FIELDS | {"times_sha256", "observations_sha256"}
_PARAMETERS = ("theta", "mu", "sigma")
_FIT_FIELDS = frozenset({
    "rep", "estimator", "source_valid", "source_error", "pipeline_error",
    "point_valid", "optimizer_success", "point_status", "point_reason", "nll", "point_diagnostics_json",
    "interval_valid", "interval_status", "interval_reason", "interval_diagnostics_json", "covariance_type",
    "hessian_condition_number", "hessian_relative_difference", "hessian_whitened_difference", "standardized_score",
    *(_PARAMETERS), *(f"{p}_{suffix}" for p in _PARAMETERS for suffix in ("lower", "upper", "se")),
    *(f"gradient_{i}" for i in range(3)),
    *(f"cov_{a}_{b}" for i, a in enumerate(_PARAMETERS) for b in _PARAMETERS[i:]),
    *(f"working_{kind}_{i}_{j}" for kind in ("cov", "hessian") for i in range(3) for j in range(i, 3)),
})


def _integer(value, label, minimum=0):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) or value < minimum:
        raise ValueError(f"{label} must be an integer >= {minimum}")
    return int(value)


def _cell(value):
    return parse_cell(value.to_dict() if type(value) is E4Cell else value)


def _boolean(value, label):
    if type(value) is not bool:
        raise ValueError(f"{label} must be Boolean")
    return value


def _immutable(value, dtype):
    array = np.asarray(value, dtype=dtype)
    # A bytes-backed view cannot have writeability re-enabled by a caller.
    return np.frombuffer(array.tobytes(), dtype=dtype).reshape(array.shape)


@dataclass(frozen=True, eq=False)
class PhaseType:
    Q: np.ndarray
    Lambda: np.ndarray
    D0: np.ndarray
    A: np.ndarray
    inverse_A: np.ndarray
    P: np.ndarray
    alpha: np.ndarray
    pi: np.ndarray


@lru_cache(maxsize=64)
def _phase_type(cell):
    if cell.mechanism == "regular":
        raise ValueError("regular gaps have no two-state phase-type representation")
    data = cell.to_dict()
    calm, burst = data["calm_switch_rate"], data["burst_switch_rate"]
    q = np.array([[-calm, calm], [burst, -burst]])
    arrivals = np.diag([data["calm_arrival_rate"], data["burst_arrival_rate"]])
    a = arrivals - q
    inverse = np.linalg.solve(a, np.eye(2))
    values = (q, arrivals, -a, a, inverse, inverse @ arrivals,
              [data["calm_event_probability"], data["burst_event_probability"]],
              [data["calm_time_probability"], data["burst_time_probability"]])
    return PhaseType(*(_immutable(value, np.float64) for value in values))


def phase_type(cell):
    """Immutable row-vector PH matrices for a strictly parsed stochastic cell."""
    return _phase_type(_cell(cell))


def floor_value(cell):
    cell = _cell(cell)
    return cell.floor_fraction * (cell.q / cell.theta) if cell.cv >= cell.floor_cv_threshold else 0.0


def _time(value):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, float, np.integer, np.floating)):
        raise ValueError("distribution time must be real")
    result = float(value)
    if math.isnan(result):
        raise ValueError("distribution time must not be NaN")
    return result


def _probability(value):
    # Bound only matrix-exponential probability roundoff, never generated data.
    tolerance = 128 * np.finfo(float).eps
    if not math.isfinite(value) or not -tolerance <= value <= 1 + tolerance:
        raise ValueError("analytical probability is outside its numerical domain")
    return float(min(1.0, max(0.0, value)))


@lru_cache(maxsize=512)
def _distribution(cell, t, transformed):
    f = floor_value(cell) if transformed else 0.0
    if t < f or t < 0:
        return 1.0, 0.0
    if cell.mechanism == "regular":
        hit = t >= max(cell.q / cell.theta, f)
        return float(not hit), float(hit)
    if math.isinf(t):
        return 0.0, 1.0
    ph = _phase_type(cell)
    s = float(ph.alpha @ expm(ph.D0 * t) @ np.ones(2))
    # The augmented absorbing generator computes tiny CDFs without 1-S loss.
    generator = np.zeros((3, 3))
    generator[:2, :2] = ph.D0
    generator[:2, 2] = np.diag(ph.Lambda)
    probability = float(ph.alpha @ expm(generator * t)[:2, 2])
    return _probability(s), _probability(probability)


def survival(cell, t, *, transformed=False):
    """P(gap>t), including the correct transformed jump at the floor."""
    return _distribution(_cell(cell), _time(t), _boolean(transformed, "transformed"))[0]


def cdf(cell, t, *, transformed=False):
    """P(gap<=t); transformed CDF equals the raw CDF at and above its floor."""
    return _distribution(_cell(cell), _time(t), _boolean(transformed, "transformed"))[1]


@lru_cache(maxsize=128)
def _moments(cell, transformed):
    f = floor_value(cell) if transformed else 0.0
    if cell.mechanism == "regular":
        mean = max(cell.q / cell.theta, f)
        return mean, mean * mean, 0.0, 0.0
    ph = _phase_type(cell)
    transition = expm(ph.D0 * f)
    first = ph.inverse_A @ np.ones(2)
    second = ph.inverse_A @ first
    mean = float(f + ph.alpha @ transition @ first)
    second_moment = float(f * f + 2 * ph.alpha @ transition @ (f * first + second))
    variance = second_moment - mean * mean
    if not np.isfinite([mean, second_moment, variance]).all() or mean <= 0 or variance <= 0:
        raise ValueError("analytical moments are outside finite admissible arithmetic")
    return mean, second_moment, variance, math.sqrt(variance) / mean


def moments(cell, *, transformed=False):
    """First/second moments, variance and CV; returned mappings are fresh."""
    return dict(zip(("mean", "second_moment", "variance", "cv"),
                    _moments(_cell(cell), _boolean(transformed, "transformed"))))


@lru_cache(maxsize=256)
def _adjacent_product(cell, transformed, lag):
    mean = _moments(cell, transformed)[0]
    if cell.mechanism in {"iid", "regular"}:
        return mean * mean
    ph = _phase_type(cell)
    f = floor_value(cell) if transformed else 0.0
    r = f * ph.P + expm(ph.D0 * f) @ ph.inverse_A @ ph.inverse_A @ ph.Lambda
    value = float(ph.alpha @ r @ np.linalg.matrix_power(ph.P, lag - 1) @ r @ np.ones(2))
    if not math.isfinite(value) or value <= 0:
        raise ValueError("analytical cross moment is outside finite positive arithmetic")
    return value


def adjacent_product(cell, *, transformed=False, lag=1):
    """E[gap_j gap_(j+lag)]; the default is the adjacent-product target."""
    return _adjacent_product(_cell(cell), _boolean(transformed, "transformed"),
                             _integer(lag, "lag", 1))


@dataclass(frozen=True, eq=False)
class GapResult:
    cell_id: str
    raw_gaps: np.ndarray
    transformed_gaps: np.ndarray
    times: np.ndarray
    start_phases: np.ndarray
    end_phases: np.ndarray
    switch_counts: np.ndarray
    floor: float

    def __post_init__(self):
        raw, transformed, times = (np.asarray(getattr(self, name), dtype=float)
                                   for name in ("raw_gaps", "transformed_gaps", "times"))
        if (raw.ndim != 1 or not len(raw) or transformed.shape != raw.shape
                or times.shape != (len(raw) + 1,) or not np.isfinite(raw).all()
                or not np.isfinite(transformed).all() or not np.isfinite(times).all()
                or np.any(raw <= 0) or np.any(transformed <= 0)):
            raise ValueError("gap arrays must be finite, positive and conformable")
        if (type(self.floor) is not float or not math.isfinite(self.floor) or self.floor < 0
                or not np.array_equal(transformed, np.maximum(raw, self.floor))):
            raise ValueError("transformed gaps must equal the declared posthoc floor")
        expected_times = np.concatenate(([0.0], np.cumsum(transformed)))
        if not np.array_equal(times, expected_times) or np.any(np.diff(times) <= 0):
            raise ValueError("stored times must have strictly positive unmodified accumulated gaps")
        for name in ("start_phases", "end_phases", "switch_counts"):
            array = np.asarray(getattr(self, name))
            if (array.shape != raw.shape or not np.issubdtype(array.dtype, np.integer)
                    or (np.any(array < 0) if name == "switch_counts" else np.any((array < -1) | (array > 1)))):
                raise ValueError("invalid phase or switch-count trace")
            object.__setattr__(self, name, _immutable(array, np.int64))
        for name, array in (("raw_gaps", raw), ("transformed_gaps", transformed), ("times", times)):
            object.__setattr__(self, name, _immutable(array, np.float64))

    @property
    def initial_phase(self):
        return int(self.start_phases[0])


def _uniform(rng):
    value = float(rng.random())
    if not math.isfinite(value) or not 0 <= value < 1:
        raise ValueError("race uniform draw must lie in [0,1)")
    return value


def generate_gaps(cell, rng, *, n_gaps=None):
    """Generate a complete raw phase trace, then floor without CTMC padding.

    The RNG port needs random() and exponential(scale). There are no retries,
    burn-in, conditional resampling, transition caps or clock repairs.
    """
    cell = _cell(cell)
    count = cell.n - 1 if n_gaps is None else _integer(n_gaps, "n_gaps", 1)
    raw = np.empty(count)
    starts, ends, switches = (np.zeros(count, dtype=np.int64) for _ in range(3))
    if cell.mechanism == "regular":
        raw.fill(cell.q / cell.theta)
        starts.fill(-1)
        ends.fill(-1)
    else:
        data = cell.to_dict()
        arrivals = (data["calm_arrival_rate"], data["burst_arrival_rate"])
        exits = (data["calm_switch_rate"], data["burst_switch_rate"])
        phase = None
        for index in range(count):
            if index == 0 or cell.mechanism == "iid":
                phase = 0 if _uniform(rng) < data["calm_event_probability"] else 1
            starts[index] = phase
            elapsed = 0.0
            while True:
                rate = arrivals[phase] + exits[phase]
                wait = float(rng.exponential(1.0 / rate))
                if not math.isfinite(wait) or wait <= 0:
                    raise ValueError("race wait must be finite and strictly positive; no redraw")
                elapsed += wait
                if not math.isfinite(elapsed):
                    raise ValueError("raw race duration overflow; no redraw")
                if _uniform(rng) < arrivals[phase] / rate:
                    raw[index], ends[index] = elapsed, phase
                    break
                phase = 1 - phase
                switches[index] += 1
    f = floor_value(cell)
    transformed = np.maximum(raw, f)
    times = np.concatenate(([0.0], np.cumsum(transformed)))
    return GapResult(cell.identity, raw, transformed, times, starts, ends, switches, float(f))


def _vector_statistics(values):
    extended = np.asarray(values, dtype=np.longdouble)
    mean = extended.mean()
    result = {
        "mean_gap": float(mean), "cv": float(np.sqrt(np.mean((extended - mean) ** 2)) / mean),
        "second_moment": float(np.mean(extended ** 2)),
        "adjacent_product": float(np.mean(extended[:-1] * extended[1:])) if len(extended) > 1 else None,
        "min_gap": float(np.min(extended)), "span": float(np.sum(extended)),
    }
    if any(value is not None and not math.isfinite(value) for value in result.values()):
        raise ValueError("gap statistics overflow")
    return result


@lru_cache(maxsize=64)
def _theory(cell):
    raw, transformed = moments(cell), moments(cell, transformed=True)
    result = {}
    for name, item in (("raw", raw), ("transformed", transformed)):
        result.update({f"expected_{name}_mean_gap": item["mean"],
                       f"expected_{name}_second_moment": item["second_moment"],
                       f"expected_{name}_cv": item["cv"],
                       f"expected_{name}_adjacent_product": adjacent_product(cell, transformed=name == "transformed")})
    result.update(expected_mean_gap_after_floor=transformed["mean"],
                  expected_clip_fraction=cdf(cell, floor_value(cell)))
    return result


def gap_statistics(cell, result):
    """Scalar diagnostics; realized_* uses differences of the exact stored times."""
    cell = _cell(cell)
    if not isinstance(result, GapResult) or result.cell_id != cell.identity or result.floor != floor_value(cell):
        raise ValueError("gap result does not belong to this scientific cell")
    raw, transformed = result.raw_gaps, result.transformed_gaps
    observed = np.diff(result.times)
    count = len(raw)
    clipped = int(np.count_nonzero(raw < result.floor))
    diagnostics = {
        "family": {"mmpp": "event_stationary_mmpp", "iid": "iid_full_palm_marginal", "regular": "equidistant"}[cell.mechanism],
        "target_mean_gap": cell.q / cell.theta, "requested_cv": cell.cv,
        "floor_applied": cell.cv >= cell.floor_cv_threshold, "floor_value": result.floor,
        "clip_count": clipped, "clip_fraction": clipped / count,
        "raw_below_floor_count": clipped,
        "transformed_below_floor_count": int(np.count_nonzero(transformed < result.floor)),
        "transformed_at_floor_count": int(np.count_nonzero(transformed == result.floor)),
        "raw_gaps_sha256": array_digest(raw), "transformed_gaps_sha256": array_digest(transformed),
        "phase_trace_sha256": hashlib.sha256(b"".join(np.asarray(array, dtype="<i8").tobytes()
            for array in (result.start_phases, result.end_phases, result.switch_counts))).hexdigest(),
        "source_initial_phase": result.initial_phase,
        "source_switch_count": int(result.switch_counts.sum()),
        "source_max_switches_per_gap": int(result.switch_counts.max()),
        "source_phase_reset_count": count if cell.mechanism == "iid" else 1 if cell.mechanism == "mmpp" else 0,
        "source_calm_start_count": int(np.count_nonzero(result.start_phases == 0)),
        "source_burst_start_count": int(np.count_nonzero(result.start_phases == 1)),
        "source_calm_end_count": int(np.count_nonzero(result.end_phases == 0)),
        "source_burst_end_count": int(np.count_nonzero(result.end_phases == 1)),
        **_theory(cell),
    }
    for name, values in (("raw", raw), ("transformed", transformed)):
        diagnostics.update({f"{name}_{key}": value for key, value in _vector_statistics(values).items()})
        for t_over_m, label in CDF_GRID:
            diagnostics[f"{name}_cdf_t_over_m_{label}"] = float(np.mean(values <= t_over_m * cell.q / cell.theta))
    actual = _vector_statistics(observed)
    diagnostics.update(realized_mean_gap=actual["mean_gap"], realized_cv=actual["cv"],
                       realized_span=float(result.times[-1] - result.times[0]), min_gap=actual["min_gap"])
    if set(diagnostics) != GAP_FIELDS:
        raise ValueError("gap diagnostic schema differs from declared GAP_FIELDS")
    return diagnostics


@dataclass(frozen=True)
class E4Source:
    seed: int

    def __post_init__(self):
        object.__setattr__(self, "seed", _integer(self.seed, "seed"))

    def _prepare(self, mapping, rep):
        cell = parse_cell(mapping)
        if cell.identity not in _PLANNED_IDS:
            raise ValueError("E4 source requires one of the seven published scientific cells")
        rep = _integer(rep, "rep")
        words = tuple(int(cell.identity[i:i + 8], 16) for i in range(0, 64, 8))
        sequence = np.random.SeedSequence(self.seed, spawn_key=(SOURCE_NAMESPACE, *words, rep))
        return cell, sequence, sequence.spawn(2)

    def address(self, mapping, rep):
        """Pure address allocation; neither child produces a random observation."""
        _, sequence, children = self._prepare(mapping, rep)
        return {"source_seed": self.seed, "source_spawn_key": canonical_json(sequence.spawn_key),
                "gap_spawn_key": canonical_json(children[0].spawn_key),
                "path_spawn_key": canonical_json(children[1].spawn_key)}

    def gaps(self, mapping, rep):
        """Schedule-only generation; the OU child is never drawn."""
        cell, _, children = self._prepare(mapping, rep)
        return generate_gaps(cell, np.random.default_rng(children[0]))

    def sample(self, mapping, rep):
        cell, sequence, children = self._prepare(mapping, rep)
        gaps = generate_gaps(cell, np.random.default_rng(children[0]))
        values = simulate_ou(cell.theta, cell.mu, cell.sigma, gaps.times, np.random.default_rng(children[1]))
        sample = Sample(values, gaps.times)
        diagnostics = {
            **gap_statistics(cell, gaps), "source_seed": self.seed,
            "source_spawn_key": canonical_json(sequence.spawn_key),
            "gap_spawn_key": canonical_json(children[0].spawn_key),
            "path_spawn_key": canonical_json(children[1].spawn_key),
            "times_sha256": array_digest(sample.times), "observations_sha256": array_digest(sample.values),
        }
        if set(diagnostics) & (set(cell.to_dict()) | _FIT_FIELDS):
            raise ValueError("source diagnostics collide with scientific or fit-record fields")
        if set(diagnostics) != SOURCE_FIELDS:
            raise ValueError("source diagnostic schema differs from declared SOURCE_FIELDS")
        return sample, diagnostics
