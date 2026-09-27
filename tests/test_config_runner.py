"""Integration properties that protect the simulation design, not just syntax."""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from ou_irregular.config import ExperimentConfig, load_config

ROOT = Path(__file__).resolve().parents[1]


def minimal_config():
    data = load_config(ROOT / "configs/week3_smoke.yaml").to_dict()
    data["grid"] = {"theta_mean_gap": [0.5], "cv": [0.0, 1.0], "n": [100]}
    data["replications"] = 2
    return ExperimentConfig.from_dict(data)


def test_config_rejects_silent_typos_and_headline_phase():
    data = minimal_config().to_dict()
    data["optimizer"]["n_start"] = 3
    with pytest.raises(ValueError):
        ExperimentConfig.from_dict(data)
    data = minimal_config().to_dict()
    data["grid"]["theta_mean_gap"] = [0.1, 0.5, 1, 2]
    with pytest.raises(ValueError, match="at most 3"):
        ExperimentConfig.from_dict(data)
    data = minimal_config().to_dict()
    data["phase"] = "headline"
    with pytest.raises(ValueError, match="Week 4"):
        ExperimentConfig.from_dict(data)


def test_duplicate_yaml_keys_are_rejected(tmp_path):
    path = tmp_path / "bad.yaml"
    path.write_text("seed: 1\nseed: 2\n")
    with pytest.raises(ValueError, match="unique"):
        load_config(path)


def test_existing_cell_and_replication_streams_survive_grid_extension():
    from ou_irregular.runner import replication_streams
    config = minimal_config()
    first = next(config.cells())
    data = config.to_dict()
    data["grid"]["cv"] = [2.0, 1.0, 0.0]
    data["replications"] = 20
    other = next(cell for cell in ExperimentConfig.from_dict(data).cells() if cell["cv"] == 0)
    assert first == other
    gap, path, diagnostics = replication_streams(config.seed, first["cell_id"], 1)
    gap2, path2, diagnostics2 = replication_streams(config.seed, other["cell_id"], 1)
    np.testing.assert_array_equal(gap.normal(size=16), gap2.normal(size=16))
    np.testing.assert_array_equal(path.normal(size=16), path2.normal(size=16))
    assert diagnostics == diagnostics2
    assert diagnostics["gap_spawn_key"] != diagnostics["path_spawn_key"]


def test_end_to_end_replay_and_replot_from_raw_rows(tmp_path):
    from ou_irregular.runner import replay, replot, run_experiment
    config = minimal_config()
    directory = tmp_path / "run"
    manifest = run_experiment(config, directory)
    assert manifest["planned_fits"] == manifest["recorded_fits"] == 12
    raw = pd.read_csv(directory / "replications.csv", float_precision="round_trip")
    assert not raw.duplicated(["cell_id", "rep", "estimator"]).any()
    cell_id = next(config.cells())["cell_id"]
    regenerated = replay(directory, cell_id, 0).sort_values("estimator")
    saved = raw[(raw.cell_id == cell_id) & (raw.rep == 0)].sort_values("estimator")
    np.testing.assert_allclose(regenerated[["theta", "mu", "sigma"]],
                               saved[["theta", "mu", "sigma"]], rtol=0, atol=0)
    summary, figures = replot(directory)
    assert summary.coverage.isna().all()
    assert set(summary.ci_status) == {"not_computed_week3"}
    assert len(figures) == 4 and all(path.is_file() for path in figures)
    with pytest.raises(FileExistsError):
        run_experiment(config, directory)
    # Truncated data cannot silently change simulation denominators on replot.
    raw.iloc[:-1].to_csv(directory / "replications.csv", index=False)
    with pytest.raises(ValueError, match="hash"):
        replot(directory)


@pytest.mark.parametrize("mode", ["nonfinite", "overflow"])
def test_bad_simulator_output_keeps_failed_rows(monkeypatch, mode):
    import ou_irregular.runner as runner
    config = minimal_config()
    def bad_simulator(theta, mu, sigma, times, rng):
        if mode == "overflow":
            raise OverflowError("arithmetic overflow")
        return np.full(times.size, np.nan)
    monkeypatch.setattr(runner, "simulate_ou", bad_simulator)
    rows, spacing = runner.run_replication(config, next(config.cells()), 0)
    assert len(rows) == 3 and not any(row["success"] for row in rows)
    assert not spacing["simulation_success"]


def test_simulation_failure_preserves_every_estimator_denominator(monkeypatch):
    import ou_irregular.runner as runner
    config = minimal_config()
    def failing_generator(*args, **kwargs):
        raise FloatingPointError("unrepresentable gap")
    monkeypatch.setattr(runner, "generate_times", failing_generator)
    rows, spacing = runner.run_replication(config, next(config.cells()), 0)
    assert len(rows) == 3 and not any(row["success"] for row in rows)
    assert all("unrepresentable gap" in row["error"] for row in rows)
    assert not spacing["simulation_success"]
