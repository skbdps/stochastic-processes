"""Small workflow executions, including failures that must retain denominators."""
import hashlib
import json

import numpy as np
import pandas as pd
import pytest
import yaml

from ou_irregular.inference.bootstrap import BootstrapResult, ParametricBootstrap
from ou_irregular.inference.adapters import ExactOuSimulator
from ou_irregular.inference.contracts import PointFit
from ou_irregular.inference.wald import WaldResult
from ou_irregular.week4.config import ValidationConfig, load_config
from ou_irregular.week4.runner import PROJECT_ROOT, ValidationExperiment, main
from ou_irregular.week4.simulation import array_digest, bootstrap_seed


@pytest.fixture
def tiny_config():
    data = load_config(PROJECT_ROOT / "configs/week4_validation.yaml").to_dict()
    data["name"] = "tiny_test_configuration"
    data["certification"].update(n=60, replications=2)
    data["bootstrap"]["draws"] = 3
    for cell in data["bootstrap"]["cells"]:
        cell["n"] = 60
    return ValidationConfig.from_dict(data)


class FigureStub:
    def write(self, summary, bootstrap_summary, output_dir):
        output_dir.mkdir(parents=True, exist_ok=True)
        artifact = output_dir / "figure_stub.txt"
        artifact.write_text("The test injected a figure writer.\n")
        return [artifact]


class EstimatorStub:
    def __init__(self, name, *, valid=True):
        self.name, self.valid = name, valid

    def fit(self, sample):
        return PointFit(self.name, 1., 0., .5, 5., self.valid, True,
                        "valid_interior" if self.valid else "boundary", "test point")


class WaldStub:
    def __init__(self, name, *, valid_point=True, raise_inference=False):
        self.estimator = EstimatorStub(name, valid=valid_point)
        self.raise_inference = raise_inference

    def infer(self, sample, point_fit=None):
        if self.raise_inference:
            raise RuntimeError("controlled inference exception")
        point = point_fit if point_fit is not None else self.estimator.fit(sample)
        if point.valid:
            return WaldResult(point, True, "valid", "test interval",
                              point.parameters - .2, point.parameters + .2,
                              .01 * np.eye(3), .01 * np.eye(3), {})
        missing = np.full(3, np.nan)
        return WaldResult(point, False, "invalid_point_fit", "test boundary",
                          missing, missing, np.full((3, 3), np.nan), np.full((3, 3), np.nan), {})


class BootstrapStub:
    def __init__(self, name, draws, calls=None):
        self.name, self.draws, self.calls = name, draws, calls

    def run(self, sample, point_fit=None, *, seed):
        if self.calls is not None:
            self.calls.append({"name": self.name, "sample_id": id(sample), "seed": seed,
                               "times_hash": array_digest(sample.times),
                               "values_hash": array_digest(sample.values),
                               "x0": float(sample.values[0])})
        rows = tuple({"draw": i, "seed_entropy": seed, "seed_spawn_key": (i,),
                      "theta": point_fit.theta, "mu": point_fit.mu, "sigma": point_fit.sigma,
                      "nll": point_fit.nll, "point_valid": True, "optimizer_success": True,
                      "status": "valid_interior", "reason": "test draw"} for i in range(self.draws))
        return BootstrapResult(point_fit, True, "valid", "test bootstrap",
                               point_fit.parameters - .1, point_fit.parameters + .1, np.zeros(3),
                               rows, self.draws, self.draws, {"seed_entropy": seed})


def stubbed_experiment(config, **kwargs):
    draws = config.to_dict()["bootstrap"]["draws"]
    options = {"wald_factory": WaldStub,
               "bootstrap_factory": lambda name: BootstrapStub(name, draws),
               "figure_writer": FigureStub()}
    options.update(kwargs)
    return ValidationExperiment(config, **options)


def read_table(directory, filename):
    return pd.read_csv(directory / filename)


def test_default_services_complete_tiny_end_to_end_run_and_manifest(tiny_config, tmp_path):
    output = tmp_path / "default_run"
    metadata = ValidationExperiment(tiny_config).run(output)
    certification = read_table(output, "certification.csv")
    summary = read_table(output, "coverage_summary.csv")
    checks = read_table(output, "bootstrap_summary.csv")
    draws = read_table(output, "bootstrap_draws.csv")
    assert metadata["status"] == "completed"
    assert metadata["unexpected_pipeline_exceptions"] == 0
    assert len(certification) == 2 and set(certification.rep) == {0, 1}
    assert len(summary) == 3 and (summary.n_planned == 2).all()
    assert len(checks) == 6 and set(checks.estimator) == {"exact", "pfml", "euler"}
    assert (checks.planned_draws == 3).all()
    assert len(draws) == checks.attempted_draws.sum()
    assert metadata["bootstrap"]["planned_draws"] == 18
    assert metadata["bootstrap"]["attempted_draws"] == len(draws)
    assert metadata["registered_configuration_matches"] is False
    assert metadata["registered_procedure_matches"] is False
    assert metadata["study_role"] == "diagnostic_unregistered"
    assert metadata["overridden_scientific_dependencies"] == []
    assert metadata["headline_experiments_run"] is False
    assert metadata["resolved_config_sha256"] == tiny_config.digest
    assert json.loads((output / "run_metadata.json").read_text()) == metadata
    for relative, expected_digest in metadata["artifact_sha256"].items():
        assert hashlib.sha256((output / relative).read_bytes()).hexdigest() == expected_digest
    assert "run_metadata.json" not in metadata["artifact_sha256"]
    assert "ou_irregular/week4/runner.py" in metadata["source_sha256"]
    assert "ou_irregular/inference/bootstrap.py" in metadata["source_sha256"]
    assert any(name.startswith("figures/") and name.endswith(".png") for name in metadata["artifact_sha256"])
    assert not read_table(output, "inference_failures.csv").interval_valid.any()


def test_simulation_failures_retain_every_planned_row_and_no_fabricated_bootstrap_draws(tiny_config, tmp_path):
    class BrokenSource:
        def sample(self, cell, rep):
            raise RuntimeError("controlled simulation exception")
    output = tmp_path / "simulation_failures"
    metadata = stubbed_experiment(tiny_config, sample_source=BrokenSource()).run(output)
    points = read_table(output, "certification.csv")
    summary = read_table(output, "coverage_summary.csv")
    checks = read_table(output, "bootstrap_summary.csv")
    draws = read_table(output, "bootstrap_draws.csv")
    assert len(points) == 2 and set(points.rep) == {0, 1}
    assert (points.interval_status == "simulation_failure").all()
    assert not points.point_valid.any() and not points.interval_valid.any()
    assert points.pipeline_error.str.contains("controlled simulation exception").all()
    assert (summary.n_planned == 2).all() and (summary.coverage_operational == 0).all()
    assert summary.coverage_conditional.isna().all()
    assert len(checks) == 6 and (checks.planned_draws == 3).all()
    assert (checks.attempted_draws == 0).all() and (checks.invalid_attempted_draws == 0).all()
    assert draws.empty and {"draw", "status", "seed_spawn_key"}.issubset(draws.columns)
    assert metadata["status"] == "implementation_review_required"
    assert metadata["unexpected_pipeline_exceptions"] == 8
    assert metadata["bootstrap"]["attempted_draws"] == 0


def test_inference_exception_preserves_valid_point_estimates_and_planned_coverage(tiny_config, tmp_path):
    output = tmp_path / "inference_failure"
    factory = lambda name: WaldStub(name, raise_inference=True)
    metadata = stubbed_experiment(tiny_config, wald_factory=factory).run(output)
    points = read_table(output, "certification.csv")
    summary = read_table(output, "coverage_summary.csv")
    checks = read_table(output, "bootstrap_summary.csv")
    assert points.point_valid.all() and not points.interval_valid.any()
    np.testing.assert_allclose(points.theta, 1.)
    assert (points.point_status == "valid_interior").all()
    assert (points.interval_status == "pipeline_exception").all()
    assert points.interval_reason.str.contains("controlled inference exception").all()
    assert (summary.n_point_valid == 2).all() and (summary.n_interval_valid == 0).all()
    assert np.isfinite(summary.bias).all() and (summary.coverage_operational == 0).all()
    assert checks.point_valid.all() and not checks.wald_valid.any()
    assert checks.bootstrap_valid.all() and (checks.attempted_draws == 3).all()
    assert metadata["status"] == "implementation_review_required"
    assert metadata["unexpected_pipeline_exceptions"] == 8


def test_boundary_original_skips_actual_bootstrap_without_marking_undrawn_samples_failed(tiny_config, tmp_path):
    class MustNotSimulate:
        def simulate(self, *args, **kwargs):
            pytest.fail("bootstrap must not simulate from a boundary original fit")
    factory = lambda name: ParametricBootstrap(EstimatorStub(name), MustNotSimulate(), draws=3)
    output = tmp_path / "boundary_original"
    metadata = stubbed_experiment(
        tiny_config, wald_factory=lambda name: WaldStub(name, valid_point=False),
        bootstrap_factory=factory).run(output)
    checks = read_table(output, "bootstrap_summary.csv")
    assert len(checks) == 6 and (checks.bootstrap_status == "skipped_invalid_point").all()
    assert (checks.planned_draws == 3).all() and (checks.attempted_draws == 0).all()
    assert (checks.valid_draws == 0).all() and (checks.invalid_attempted_draws == 0).all()
    assert read_table(output, "bootstrap_draws.csv").empty
    # A diagnosed boundary is an ordinary retained statistical outcome, not
    # an unexpected code/dependency exception.
    assert metadata["status"] == "completed" and metadata["unexpected_pipeline_exceptions"] == 0


def test_bootstrap_all_estimators_share_source_data_and_use_documented_seeds(tiny_config, tmp_path):
    calls = []
    output = tmp_path / "shared_sources"
    factory = lambda name: BootstrapStub(name, 3, calls)
    metadata = stubbed_experiment(tiny_config, bootstrap_factory=factory).run(output)
    checks = read_table(output, "bootstrap_summary.csv")
    draws = read_table(output, "bootstrap_draws.csv")
    cert = read_table(output, "certification.csv")
    assert len(calls) == len(checks) == 6
    for cell_id, group in checks.groupby("cell_id"):
        assert group.times_sha256.nunique() == 1 and group.observations_sha256.nunique() == 1
        expected = {name: bootstrap_seed(tiny_config.to_dict()["bootstrap"]["seed"], cell_id, name)
                    for name in ("exact", "pfml", "euler")}
        for row in group.itertuples():
            call = next(call for call in calls if call["name"] == row.estimator
                        and call["values_hash"] == row.observations_sha256)
            assert call["seed"] == row.bootstrap_seed == expected[row.estimator]
            assert call["times_hash"] == row.times_sha256
            subset = draws[(draws.cell_id == cell_id) & (draws.estimator == row.estimator)]
            assert set(subset.draw) == {0, 1, 2}
            assert [json.loads(key) for key in subset.seed_spawn_key] == [[0], [1], [2]]
        matching_calls = [call for call in calls if call["values_hash"] == group.observations_sha256.iloc[0]]
        assert len({call["sample_id"] for call in matching_calls}) == 1
    regular_checks = checks.loc[checks.cv == 0]
    original = cert.loc[cert.rep == 0].iloc[0]
    assert (regular_checks.times_sha256 == original.times_sha256).all()
    assert (regular_checks.observations_sha256 == original.observations_sha256).all()
    assert metadata["registered_configuration_matches"] is False
    assert metadata["registered_procedure_matches"] is False
    assert set(metadata["overridden_scientific_dependencies"]) == {"wald_factory", "bootstrap_factory"}


def test_bootstrap_draw_dependency_exceptions_require_implementation_review(tiny_config, tmp_path):
    class BrokenEstimator(EstimatorStub):
        def fit(self, sample):
            raise RuntimeError("controlled bootstrap refit exception")
    factory = lambda name: ParametricBootstrap(BrokenEstimator(name), ExactOuSimulator(), draws=3)
    output = tmp_path / "draw_exceptions"
    metadata = stubbed_experiment(tiny_config, bootstrap_factory=factory).run(output)
    checks = read_table(output, "bootstrap_summary.csv")
    draws = read_table(output, "bootstrap_draws.csv")
    assert len(draws) == 18 and draws.dependency_exception.all()
    assert (draws.exception_type == "RuntimeError").all()
    assert draws.reason.str.contains("controlled bootstrap refit exception").all()
    assert (checks.attempted_draws == 3).all() and (checks.invalid_attempted_draws == 3).all()
    assert (checks.valid_draws == 0).all() and not checks.bootstrap_valid.any()
    assert metadata["bootstrap"]["dependency_exceptions"] == 18
    assert metadata["unexpected_pipeline_exceptions"] == 18
    assert metadata["status"] == "implementation_review_required"


def test_interrupted_figure_writer_leaves_incomplete_manifest_and_completed_tables(tiny_config, tmp_path):
    class InterruptedFigures:
        def write(self, *args):
            raise KeyboardInterrupt("controlled figure interruption")
    output = tmp_path / "interrupted"
    with pytest.raises(KeyboardInterrupt, match="controlled figure interruption"):
        stubbed_experiment(tiny_config, figure_writer=InterruptedFigures()).run(output)
    metadata = json.loads((output / "run_metadata.json").read_text())
    assert metadata["status"] == "incomplete"
    assert metadata["termination_type"] == "KeyboardInterrupt"
    assert metadata["termination_reason"] == "controlled figure interruption"
    assert metadata["resolved_config_sha256"] == tiny_config.digest
    assert len(read_table(output, "certification.csv")) == 2
    assert len(read_table(output, "bootstrap_summary.csv")) == 6
    assert len(read_table(output, "bootstrap_draws.csv")) == 18


@pytest.mark.parametrize("existing_file", [False, True])
def test_existing_evidence_is_never_overwritten(tiny_config, tmp_path, existing_file):
    output = tmp_path / "evidence"
    if existing_file:
        output.write_text("keep me")
        evidence = output
    else:
        output.mkdir()
        evidence = output / "original.txt"
        evidence.write_text("keep me")
    with pytest.raises(FileExistsError, match="never overwritten"):
        stubbed_experiment(tiny_config).run(output)
    assert evidence.read_text() == "keep me"


def test_empty_output_directory_is_allowed_but_second_run_is_rejected(tiny_config, tmp_path):
    output = tmp_path / "empty"
    output.mkdir()
    experiment = stubbed_experiment(tiny_config)
    metadata = experiment.run(output)
    assert metadata["status"] == "completed"
    with pytest.raises(FileExistsError):
        experiment.run(output)


@pytest.mark.parametrize("prereg_commit", ["main", "abc123", "x" * 40, "A" * 40])
def test_prereg_reference_requires_full_lowercase_git_sha(tiny_config, tmp_path, prereg_commit):
    output = tmp_path / "bad_reference"
    with pytest.raises(ValueError, match="full Git commit SHA"):
        stubbed_experiment(tiny_config).run(output, prereg_commit=prereg_commit)
    assert not output.exists()


def test_cli_loads_config_passes_arguments_and_reports_success(tiny_config, tmp_path, monkeypatch, capsys):
    from ou_irregular.week4 import runner
    config_path, output = tmp_path / "tiny.yaml", tmp_path / "cli_output"
    config_path.write_text(yaml.safe_dump(tiny_config.to_dict()))
    seen = {}
    class ExperimentStub:
        def __init__(self, config):
            seen["digest"] = config.digest
        def run(self, path, *, prereg_commit, progress):
            seen.update(path=path, prereg_commit=prereg_commit)
            progress("retained planned rows")
            return {"status": "completed", "unexpected_pipeline_exceptions": 0}
    monkeypatch.setattr(runner, "ValidationExperiment", ExperimentStub)
    reference = "a" * 40
    assert main(["run", "--config", str(config_path), "--output", str(output),
                 "--prereg-commit", reference]) == 0
    assert seen == {"digest": tiny_config.digest, "path": output, "prereg_commit": reference}
    assert "retained planned rows" in capsys.readouterr().out


def test_cli_returns_failure_code_for_unexpected_pipeline_exception(tiny_config, tmp_path, monkeypatch):
    from ou_irregular.week4 import runner
    config_path = tmp_path / "tiny.yaml"
    config_path.write_text(yaml.safe_dump(tiny_config.to_dict()))
    class ExperimentStub:
        def __init__(self, config):
            pass
        def run(self, *args, **kwargs):
            return {"status": "implementation_review_required", "unexpected_pipeline_exceptions": 1}
    monkeypatch.setattr(runner, "ValidationExperiment", ExperimentStub)
    assert main(["run", "--config", str(config_path), "--output", str(tmp_path / "output")]) == 1
