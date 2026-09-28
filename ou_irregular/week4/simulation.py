"""Sampling adapter reusing stable model primitives, never the Week 3 runner."""
from dataclasses import dataclass
import hashlib

import numpy as np

from ..inference.contracts import Sample
from ..ou import simulate_ou
from ..spacing import generate_times
from .config import canonical_json


def array_digest(values):
    return hashlib.sha256(np.asarray(values, dtype="<f8").tobytes()).hexdigest()


@dataclass(frozen=True)
class OuSampleSource:
    seed: int

    def __post_init__(self):
        _seed(self.seed)

    def sample(self, cell, rep):
        if isinstance(rep, bool) or not isinstance(rep, (int, np.integer)) or rep < 0:
            raise ValueError("rep must be a nonnegative integer")
        words = tuple(int(cell["cell_id"][i:i + 8], 16) for i in range(0, 64, 8))
        sequence = np.random.SeedSequence(self.seed, spawn_key=(4, *words, int(rep)))
        gap_seed, path_seed = sequence.spawn(2)
        schedule = generate_times(cell["n"], cell["mean_gap"], cell["cv"],
                                  np.random.default_rng(gap_seed), floor_fraction=cell["floor_fraction"],
                                  floor_cv_threshold=cell["floor_cv_threshold"])
        values = simulate_ou(cell["true_theta"], cell["true_mu"], cell["true_sigma"],
                             schedule.times, np.random.default_rng(path_seed))
        sample = Sample(values, schedule.times)
        return sample, {**schedule.diagnostics, "source_seed": self.seed,
                        "source_spawn_key": canonical_json(sequence.spawn_key),
                        "gap_spawn_key": canonical_json(gap_seed.spawn_key),
                        "path_spawn_key": canonical_json(path_seed.spawn_key),
                        "times_sha256": array_digest(sample.times),
                        "observations_sha256": array_digest(sample.values)}


def bootstrap_seed(root_seed, cell_id, estimator):
    """Stable child identity; estimator ordering does not change a bootstrap."""
    _seed(root_seed)
    words = [int(cell_id[i:i + 8], 16) for i in range(0, 64, 8)]
    index = {"exact": 0, "pfml": 1, "euler": 2}[estimator]
    sequence = np.random.SeedSequence(root_seed, spawn_key=(4, *words, index))
    return int(sequence.generate_state(1, dtype=np.uint64)[0])


def _seed(value):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) or value < 0:
        raise ValueError("seed must be a nonnegative integer")
