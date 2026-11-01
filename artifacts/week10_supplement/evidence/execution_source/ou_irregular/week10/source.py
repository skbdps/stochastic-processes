"""Palm-origin schedules with an exact common right endpoint and random count."""
from dataclasses import dataclass
import hashlib
import math

import numpy as np

from ..inference.contracts import Sample
from ..ou import simulate_ou
from ..week4.config import canonical_json
from ..week4.simulation import array_digest
from ..week8.source import floor_value, generate_gaps, moments
from .design import HorizonCell, SOURCE_NAMESPACE, full_design, parse_cell


def _integer(value, name):
    if type(value) is not int or value < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return value


def _immutable(values, dtype):
    array = np.asarray(values, dtype=dtype)
    return np.frombuffer(array.tobytes(), dtype=dtype).reshape(array.shape)


def _uniform(rng):
    value = float(rng.random())
    if not math.isfinite(value) or not 0 <= value < 1:
        raise ValueError("race uniform outside [0,1)")
    return value


@dataclass(frozen=True, eq=False)
class HorizonSchedule:
    cell_id: str
    horizon: float
    times: np.ndarray
    raw_parent_gaps: np.ndarray
    transformed_parent_gaps: np.ndarray
    start_phases: np.ndarray
    end_phases: np.ndarray
    switch_counts: np.ndarray
    floor: float
    terminal_is_arrival: bool

    def __post_init__(self):
        raw = np.asarray(self.raw_parent_gaps, dtype=float)
        transformed = np.asarray(self.transformed_parent_gaps, dtype=float)
        times = np.asarray(self.times, dtype=float)
        if (raw.ndim != 1 or not len(raw) or transformed.shape != raw.shape
                or times.shape != (len(raw) + 1,) or not np.isfinite(raw).all()
                or np.any(raw <= 0) or not np.isfinite(times).all()):
            raise ValueError("invalid fixed-horizon schedule arrays")
        if not np.array_equal(transformed, np.maximum(raw, self.floor)):
            raise ValueError("parent transform differs from declared floor")
        if times[0] != 0. or times[-1] != self.horizon or np.any(np.diff(times) <= 0):
            raise ValueError("schedule must have strictly increasing exact endpoints")
        parent_times = np.cumsum(transformed)
        if (not np.array_equal(times[1:-1], parent_times[:-1])
                or parent_times[-1] < self.horizon
                or self.terminal_is_arrival != bool(parent_times[-1] == self.horizon)):
            raise ValueError("terminal censoring or retained parent arrival times disagree")
        for name in ("start_phases", "end_phases", "switch_counts"):
            values = np.asarray(getattr(self, name))
            if (values.shape != raw.shape or not np.issubdtype(values.dtype, np.integer)
                    or (np.any(values < 0) if name == "switch_counts"
                        else np.any((values < -1) | (values > 1)))):
                raise ValueError("invalid schedule phase trace")
            object.__setattr__(self, name, _immutable(values, np.int64))
        for name, values in (("times", times), ("raw_parent_gaps", raw),
                             ("transformed_parent_gaps", transformed)):
            object.__setattr__(self, name, _immutable(values, np.float64))


def generate_schedule(cell, rng):
    """Stop the parent sequence on its first arrival at/after T; never redraw.

    Store the complete crossing parent gap for audit, but observe only X(T).
    The deterministic terminal observation is not an extra parent arrival.
    A schedule with fewer than three observations remains a possible source
    outcome; inference will retain its unavailable result rather than resample.
    """
    if type(cell) is not HorizonCell:
        cell = parse_cell(cell)
    parent = cell.parent_cell()
    data = parent.to_dict()
    floor = floor_value(parent)
    raw, transformed, starts, ends, switches = [], [], [], [], []
    times = [0.]
    elapsed_clock = 0.
    phase = None
    arrivals = (data["calm_arrival_rate"], data["burst_arrival_rate"])
    exits = (data["calm_switch_rate"], data["burst_switch_rate"])
    while elapsed_clock < cell.horizon:
        if cell.mechanism == "regular":
            duration, start, ending, count = cell.q / cell.theta, -1, -1, 0
        else:
            if phase is None or cell.mechanism == "iid":
                phase = 0 if _uniform(rng) < data["calm_event_probability"] else 1
            start, duration, count = phase, 0., 0
            while True:
                rate = arrivals[phase] + exits[phase]
                wait = float(rng.exponential(1. / rate))
                if not math.isfinite(wait) or wait <= 0:
                    raise ValueError("race wait must be positive and finite; no redraw")
                duration += wait
                if not math.isfinite(duration):
                    raise ValueError("nonfinite raw parent gap")
                if _uniform(rng) < arrivals[phase] / rate:
                    ending = phase
                    break
                phase = 1 - phase
                count += 1
        used = max(duration, floor)
        next_clock = elapsed_clock + used
        if not math.isfinite(next_clock) or next_clock <= elapsed_clock:
            raise ValueError("parent arrival time cannot be represented; no clock repair")
        raw.append(duration)
        transformed.append(used)
        starts.append(start)
        ends.append(ending)
        switches.append(count)
        if next_clock >= cell.horizon:
            times.append(cell.horizon)
            break
        times.append(next_clock)
        elapsed_clock = next_clock
    return HorizonSchedule(cell.identity, cell.horizon, np.asarray(times), np.asarray(raw),
                           np.asarray(transformed), np.asarray(starts), np.asarray(ends),
                           np.asarray(switches), float(floor), bool(next_clock == cell.horizon))


def schedule_diagnostics(cell, result):
    gaps = np.diff(result.times)
    raw, transformed = result.raw_parent_gaps, result.transformed_parent_gaps
    count = len(raw)
    clip_count = int(np.count_nonzero(raw < result.floor))
    # The final parent draw is crossing/censored; all parent-draw summaries
    # include that draw and are descriptive stopped-sample diagnostics only.
    parent = cell.parent_cell()
    return {
        "actual_n": len(result.times), "actual_transitions": len(gaps),
        "interior_arrival_count": len(result.times) - 2,
        "terminal_is_arrival": result.terminal_is_arrival,
        "terminal_gap": float(gaps[-1]), "terminal_parent_raw_gap": float(raw[-1]),
        "terminal_parent_transformed_gap": float(transformed[-1]),
        "terminal_parent_overshoot": float(result.times[-2] + transformed[-1] - cell.horizon),
        "terminal_gap_below_floor": bool(gaps[-1] < result.floor),
        "family": "fixed_horizon_" + cell.mechanism,
        "target_mean_gap": cell.q / cell.theta, "requested_cv": cell.cv,
        "floor_applied": result.floor > 0., "floor_value": result.floor,
        "clip_count": clip_count, "clip_fraction": clip_count / count,
        "parent_draw_count_including_crossing": count,
        "realized_mean_gap": float(np.mean(gaps)), "realized_cv": float(np.std(gaps) / np.mean(gaps)),
        "realized_span": float(result.times[-1]), "min_gap": float(np.min(gaps)),
        "raw_mean_gap": float(np.mean(raw)), "raw_cv": float(np.std(raw) / np.mean(raw)),
        "raw_gaps_sha256": array_digest(raw), "transformed_gaps_sha256": array_digest(transformed),
        "phase_trace_sha256": hashlib.sha256(b"".join(np.asarray(a, dtype="<i8").tobytes()
            for a in (result.start_phases, result.end_phases, result.switch_counts))).hexdigest(),
        "source_initial_phase": int(result.start_phases[0]),
        "source_switch_count": int(np.sum(result.switch_counts)),
        "source_phase_reset_count": count if cell.mechanism == "iid" else int(cell.mechanism == "mmpp"),
        "expected_mean_gap_after_floor": moments(parent, transformed=True)["mean"],
    }


@dataclass(frozen=True)
class HorizonSource:
    seed: int

    def __post_init__(self):
        _integer(self.seed, "seed")

    def _prepare(self, mapping, rep):
        cell = mapping if type(mapping) is HorizonCell else parse_cell(mapping)
        if cell.identity not in {c.identity for c in full_design().cells}:
            raise ValueError("HorizonSource accepts only the seven specified cells")
        _integer(rep, "rep")
        words = tuple(int(cell.identity[i:i+8], 16) for i in range(0, 64, 8))
        sequence = np.random.SeedSequence(self.seed, spawn_key=(SOURCE_NAMESPACE, *words, rep))
        return cell, sequence, sequence.spawn(2)

    def address(self, mapping, rep):
        _, sequence, children = self._prepare(mapping, rep)
        return {"source_seed": self.seed, "source_spawn_key": canonical_json(sequence.spawn_key),
                "gap_spawn_key": canonical_json(children[0].spawn_key),
                "path_spawn_key": canonical_json(children[1].spawn_key)}

    def schedule(self, mapping, rep):
        cell, _, children = self._prepare(mapping, rep)
        return generate_schedule(cell, np.random.default_rng(children[0]))

    def sample(self, mapping, rep):
        cell, sequence, children = self._prepare(mapping, rep)
        schedule = generate_schedule(cell, np.random.default_rng(children[0]))
        values = simulate_ou(cell.theta, cell.mu, cell.sigma, schedule.times,
                             np.random.default_rng(children[1]))
        sample = Sample(values, schedule.times)
        diagnostic = {**schedule_diagnostics(cell, schedule),
                      "source_seed": self.seed, "source_spawn_key": canonical_json(sequence.spawn_key),
                      "gap_spawn_key": canonical_json(children[0].spawn_key),
                      "path_spawn_key": canonical_json(children[1].spawn_key),
                      "times_sha256": array_digest(sample.times),
                      "observations_sha256": array_digest(sample.values)}
        if set(diagnostic) & set(cell.to_dict()):
            raise ValueError("source diagnostics collide with design fields")
        return sample, diagnostic


def validate_sources():
    """Fixed independent source-only diagnostics; no OU samples or fits.

    Six-MCSE limits are prespecified implementation checks, not research
    confidence statements. Non-Poisson finite-window count means are reported
    without incorrectly assuming calendar-stationary initialization.
    """
    seed, replications = 22261101, 256
    source = HorizonSource(seed)
    rows, checks = [], []
    for cell in full_design().cells:
        if cell.mechanism == "regular":
            continue
        first_phases, counts = [], []
        for rep in range(replications):
            schedule = source.schedule(cell, rep)
            first_phases.append(int(schedule.start_phases[0]))
            counts.append(len(schedule.times) - 2)
            rows.append({"cell_id": cell.identity, "mechanism": cell.mechanism, "ratio": cell.ratio,
                         "rep": rep, "first_phase": first_phases[-1], "interior_arrivals": counts[-1],
                         "actual_n": len(schedule.times), "span": float(schedule.times[-1]),
                         "terminal_gap": float(np.diff(schedule.times)[-1]),
                         "times_sha256": array_digest(schedule.times)})
            if rep in {0, replications - 1}:
                _, _, children = source._prepare(cell, rep)
                parent = generate_gaps(cell.parent_cell(), np.random.default_rng(children[0]),
                                       n_gaps=len(schedule.raw_parent_gaps))
                equivalent = all(np.array_equal(a, b) for a, b in (
                    (schedule.raw_parent_gaps, parent.raw_gaps),
                    (schedule.transformed_parent_gaps, parent.transformed_gaps),
                    (schedule.start_phases, parent.start_phases),
                    (schedule.end_phases, parent.end_phases),
                    (schedule.switch_counts, parent.switch_counts)))
                checks.append({"cell_id": cell.identity, "check": "unchanged_uncensored_E4_prefix",
                               "rep": rep, "passed": bool(equivalent)})
        probability = cell.parent_cell().to_dict()["burst_event_probability"]
        actual = float(np.mean(first_phases))
        tolerance = 6 * math.sqrt(probability * (1 - probability) / replications)
        checks.append({"cell_id": cell.identity, "check": "Palm_initial_burst_phase",
                       "observed": actual, "expected": probability, "tolerance": tolerance,
                       "passed": bool(abs(actual - probability) <= tolerance)})
        if cell.ratio == 1.:
            counts = np.asarray(counts, dtype=float)
            expected = cell.horizon / (cell.q / cell.theta)
            for label, observed, tolerance in (
                ("Poisson_count_mean", float(counts.mean()), 6 * math.sqrt(expected / replications)),
                ("Poisson_count_variance", float(counts.var(ddof=1)),
                 6 * math.sqrt((expected + 2 * expected**2) / (replications - 1))),
            ):
                checks.append({"cell_id": cell.identity, "check": label, "observed": observed,
                               "expected": expected, "tolerance": tolerance,
                               "passed": bool(abs(observed - expected) <= tolerance)})
    return {"seed": seed, "replications_per_stochastic_cell": replications,
            "source_sequences": len(rows), "new_ou_paths": 0, "new_fits": 0,
            "passed": all(c["passed"] for c in checks), "checks": checks, "sequences": rows}
