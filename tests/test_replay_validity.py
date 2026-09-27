"""End-to-end safeguards for changed fit semantics and historical reproduction."""
import copy
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from ou_irregular.config import ExperimentConfig, load_config
from ou_irregular import runner

ROOT = Path(__file__).resolve().parents[1]


def small_config():
    data = load_config(ROOT / "configs/week3_smoke.yaml").to_dict()
    data["grid"] = {"theta_mean_gap": [0.5], "cv": [0], "n": [100]}
    data["replications"] = 2
    return ExperimentConfig.from_dict(data)


@pytest.fixture(scope="module")
def saved_run(tmp_path_factory):
    path = tmp_path_factory.mktemp("classified_run")
    runner.run_experiment(small_config(), path)
    return path


def test_every_outcome_reconciles_across_exports_summary_and_manifest(tmp_path, monkeypatch):
    statuses = iter(["valid_interior", "boundary", "boundary_suspected", "degenerate",
                     "numerical_failure", "simulation_failure"])

    def fake_fit(estimator, x, times, **kwargs):
        status = next(statuses)
        converged = status in {"valid_interior", "boundary", "boundary_suspected"}
        return {"estimator": estimator, "theta": 1., "mu": 0., "sigma": .5,
                "success": converged, "optimizer_success": converged,
                "valid_for_point_summary": status == "valid_interior",
                "fit_status": status, "reason": "injected outcome", "solver": "test",
                "nll": 1.}

    monkeypatch.setattr(runner, "fit_estimator", fake_fit)
    manifest = runner.run_experiment(small_config(), tmp_path / "run")
    assert manifest["status"] == "completed_with_excluded_fits"
    assert manifest["planned_fits"] == manifest["recorded_fits"] == 6
    counts = manifest["validity_counts"]
    assert counts["n_valid"] == 1 and counts["n_excluded"] == 5
    assert counts["n_optimizer_success"] == 3
    for status in ("valid_interior", "boundary", "boundary_suspected", "degenerate",
                   "numerical_failure", "simulation_failure"):
        assert counts[f"n_{status}"] == 1
    summary = pd.read_csv(tmp_path / "run/summary.csv")
    theta = summary[summary.parameter == "theta"]
    assert theta.n_valid.sum() == 1 and theta.n_total.sum() == 6
    assert len(pd.read_csv(tmp_path / "run/failures.csv")) == 5
    assert len(pd.read_csv(tmp_path / "run/boundary_fits.csv")) == 2
    assert summary.coverage.isna().all()


@pytest.mark.parametrize("changed", ["versions", "source_sha256", "python"])
def test_replay_rejects_provenance_mismatch_before_fitting(saved_run, monkeypatch, changed):
    original = runner._provenance

    def altered(config):
        value = copy.deepcopy(original(config))
        if changed == "python":
            value[changed] = "different"
        else:
            value[changed][next(iter(value[changed]))] = "different"
        return value

    monkeypatch.setattr(runner, "_provenance", altered)
    monkeypatch.setattr(runner, "run_replication", lambda *args: pytest.fail("must check before fitting"))
    cell = next(small_config().cells())["cell_id"]
    with pytest.raises(ValueError, match="provenance mismatch"):
        runner.replay(saved_run, cell, 0)


def test_deliberate_patched_replay_records_both_provenances(saved_run, tmp_path, monkeypatch):
    original = runner._provenance

    def altered(config):
        result = original(config)
        result["source_sha256"]["estimators.py"] = "test-patched-source"
        return result

    monkeypatch.setattr(runner, "_provenance", altered)
    cell = next(small_config().cells())["cell_id"]
    with pytest.raises(ValueError, match="separate --output"):
        runner.replay(saved_run, cell, 0, allow_provenance_mismatch=True)
    with pytest.warns(RuntimeWarning, match="PATCHED-CODE"):
        results = runner.replay(saved_run, cell, 0, allow_provenance_mismatch=True,
                                output_dir=tmp_path / "replay")
    evidence = json.loads((tmp_path / "replay/replay_metadata.json").read_text())
    assert evidence["kind"] == "patched_code_replay" and evidence["input_identity"] is True
    assert evidence["saved_provenance"]["source_sha256"] != evidence["executing_provenance"]["source_sha256"]
    assert len(results) == 3


def test_replay_checks_arrays_not_just_seed(saved_run, monkeypatch):
    original = runner.simulate_ou
    monkeypatch.setattr(runner, "simulate_ou", lambda *a, **kw: original(*a, **kw) + 1e-3)
    with pytest.raises(ValueError, match="identical input"):
        runner.replay(saved_run, next(small_config().cells())["cell_id"], 0)


def test_changed_replot_requires_new_derived_directory(saved_run, tmp_path, monkeypatch):
    original = runner._provenance
    before = (saved_run / "run_metadata.json").read_bytes()

    def altered(config):
        result = original(config)
        result["versions"]["numpy"] = "changed"
        return result

    monkeypatch.setattr(runner, "_provenance", altered)
    with pytest.raises(ValueError, match="separate --output"):
        runner.replot(saved_run)
    summary, paths = runner.replot(saved_run, tmp_path / "derived")
    assert (saved_run / "run_metadata.json").read_bytes() == before
    evidence = json.loads((tmp_path / "derived/run_metadata.json").read_text())
    assert evidence["derived_from"]["manifest_sha256"] == hashlib.sha256(before).hexdigest()
    assert len(summary) == 9 and len(paths) == 4


def test_failed_simulation_has_explicit_classification(monkeypatch):
    monkeypatch.setattr(runner, "simulate_ou", lambda *a, **kw: np.full(100, np.nan))
    rows, _ = runner.run_replication(small_config(), next(small_config().cells()), 0)
    assert len(rows) == 3
    assert all(row["fit_status"] == "simulation_failure" and not row["valid_for_point_summary"] for row in rows)
