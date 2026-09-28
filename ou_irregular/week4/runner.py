"""Isolated Week 4 validation orchestration and command-line composition root.

The inference layer sees observations, never simulation truth. This layer owns
designs, denominators, files and plots. No Week 3 runner/config/metric schema is
imported or extended. Factories make numerical services independently replaceable
in tests; the default composition uses the registered adapters and procedures.
"""
import argparse
from collections.abc import Mapping
from datetime import datetime, timezone
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import platform
import re
import time

import numpy as np
import pandas as pd

from ..inference.adapters import ExactOuSimulator, OuLikelihood, Week3Estimator
from ..inference.bootstrap import BootstrapResult, ParametricBootstrap
from ..inference.contracts import PointFit
from ..inference.hessian import FiniteDifferenceHessian
from ..inference.wald import WaldInference, WaldResult
from .config import ValidationConfig, canonical_json, load_config
from .plots import DiagnosticFigureWriter
from .simulation import OuSampleSource, bootstrap_seed
from .summaries import CoverageSummarizer, PARAMETERS


PROJECT_ROOT = Path(__file__).resolve().parents[2]
PREREG_COMMIT = "1b1255b0650ad7efcf700316bd1cccf704af4d2d"
# These digests belong to that published, pre-execution freeze. Comparing
# against today's editable YAML alone would falsely bless later amendments.
PREREG_SHA256 = "d7b1a7b6a6afa6351d3115ade91100402f846eb931952b394ecec4b68dd2813a"
FROZEN_CONFIG_SHA256 = "a02f99aa82898e22dd4dde76f9b1eb12b4acc940e7f695cab0d170e1ee34ac25"
REGISTERED_CONFIG_DIGEST = "66844e989d1e9dfddcb13fe8eef5c3d4fad70c7a01215e34e560a2592196a24b"


def _json_safe(value):
    """Serialize immutable diagnostics; JSON null is an explicit missing value."""
    if isinstance(value, Mapping):
        return {str(k): _json_safe(v) for k, v in value.items()}
    if isinstance(value, (tuple, list, np.ndarray)):
        return [_json_safe(v) for v in value]
    if isinstance(value, np.generic):
        return _json_safe(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def _json(value):
    return json.dumps(_json_safe(value), sort_keys=True, allow_nan=False)


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _failed_wald(name, status, reason, point_fit=None):
    fit = point_fit if point_fit is not None else PointFit(
        name, np.nan, np.nan, np.nan, np.nan, False, False, status, reason)
    return WaldResult(fit, False, status, reason, np.full(3, np.nan),
                      np.full(3, np.nan), np.full((3, 3), np.nan),
                      np.full((3, 3), np.nan), {})


def _wald_record(result):
    """Keep the point candidate even if the full-information interval fails."""
    fit, diagnostics = result.point_fit, result.diagnostics
    hessian = diagnostics.get("hessian_diagnostics", {})
    row = {"estimator": fit.estimator, "point_valid": bool(fit.valid),
           "optimizer_success": bool(fit.optimizer_success),
           "point_status": fit.status, "point_reason": fit.reason, "nll": fit.nll,
           "point_diagnostics_json": _json(fit.diagnostics),
           "interval_valid": bool(result.valid), "interval_status": result.status,
           "interval_reason": result.reason,
           "interval_diagnostics_json": _json(diagnostics),
           "covariance_type": "inverse_observed_hessian_model_based",
           "hessian_condition_number": diagnostics.get("hessian_condition_number", np.nan),
           "hessian_relative_difference": diagnostics.get("hessian_relative_difference", np.nan),
           "hessian_whitened_difference": hessian.get("whitened_relative_difference", np.nan),
           "standardized_score": diagnostics.get("standardized_score", np.nan)}
    for i, parameter in enumerate(PARAMETERS):
        row[parameter] = float(fit.parameters[i])
        row[f"{parameter}_lower"], row[f"{parameter}_upper"] = result.lower[i], result.upper[i]
        row[f"{parameter}_se"] = np.sqrt(result.covariance[i, i])
        for j in range(i, 3):
            row[f"cov_{parameter}_{PARAMETERS[j]}"] = result.covariance[i, j]
            row[f"working_cov_{i}_{j}"] = result.working_covariance[i, j]
    return row


class ValidationExperiment:
    """Coordinate a fixed plan using small, injected service interfaces.

    A failed simulation, fit or interval retains its planned row. Unexpected
    dependency exceptions are marked separately and block an implementation
    validation claim, rather than being hidden as ordinary sampling outcomes.
    """

    def __init__(self, config, *, sample_source=None, wald_factory=None,
                 bootstrap_factory=None, summarizer=None, figure_writer=None):
        self.config = ValidationConfig.from_dict(config.to_dict())
        self.settings = config.to_dict()
        self.overridden_scientific_dependencies = [name for name, dependency in (
            ("sample_source", sample_source), ("wald_factory", wald_factory),
            ("bootstrap_factory", bootstrap_factory), ("summarizer", summarizer)) if dependency is not None]
        self.sample_source = sample_source or OuSampleSource(self.settings["seed"])
        self.wald_factory = wald_factory or self._wald
        self.bootstrap_factory = bootstrap_factory or self._bootstrap
        self.summarizer = summarizer or CoverageSummarizer()
        self.figure_writer = figure_writer or DiagnosticFigureWriter()

    def _wald(self, name):
        settings = self.settings["wald"]
        return WaldInference(
            Week3Estimator(name), OuLikelihood(name),
            hessian=FiniteDifferenceHessian(
                relative_step=settings["relative_step"],
                max_relative_difference=settings["hessian_relative_tolerance"],
                max_whitened_difference=settings["max_whitened_difference"]),
            level=settings["level"], max_condition_number=settings["max_condition"],
            max_standardized_score=settings["max_standardized_score"])

    def _bootstrap(self, name):
        settings = self.settings["bootstrap"]
        return ParametricBootstrap(Week3Estimator(name), ExactOuSimulator(),
                                   draws=settings["draws"], level=settings["level"],
                                   min_valid_fraction=settings["min_valid_fraction"])

    def _sample(self, cell, rep):
        try:
            sample, diagnostics = self.sample_source.sample(cell, rep)
            return sample, diagnostics, ""
        except Exception as exc:
            return None, {}, f"{type(exc).__name__}: {exc}"

    @staticmethod
    def _infer(service, sample, source_error):
        fit = None
        if sample is None:
            return _failed_wald(service.estimator.name, "simulation_failure", source_error), source_error
        try:
            # Separating fitting preserves an available point if inference raises.
            fit = service.estimator.fit(sample)
            return service.infer(sample, point_fit=fit), ""
        except Exception as exc:
            reason = f"{type(exc).__name__}: {exc}"
            return _failed_wald(service.estimator.name, "pipeline_exception", reason, fit), reason

    def run(self, output_dir, *, prereg_commit=PREREG_COMMIT, progress=None):
        if not re.fullmatch(r"[0-9a-f]{40}", prereg_commit):
            raise ValueError("prereg_commit must be a full Git commit SHA")
        output = Path(output_dir)
        if output.exists() and (not output.is_dir() or any(output.iterdir())):
            raise FileExistsError("Use a fresh output directory; existing evidence is never overwritten")
        output.mkdir(parents=True, exist_ok=True)
        started, clock = datetime.now(timezone.utc).isoformat(), time.perf_counter()
        if canonical_json(self.settings) != canonical_json(self.config.to_dict()):
            raise ValueError("configuration changed after experiment construction")
        manifest_path = output / "run_metadata.json"
        manifest_path.write_text(_json({"project": "ou-irregular-sampling-week4", "schema_version": 1,
                                      "status": "running_incomplete", "started_utc": started,
                                      "resolved_config_sha256": self.config.digest,
                                      "preregistration_commit": prereg_commit}) + "\n")
        (output / "resolved_config.json").write_text(_json(self.settings) + "\n")
        try:
            return self._execute(output, prereg_commit, progress, started, clock)
        except BaseException as exc:
            # Interrupted runs retain completed batches but cannot masquerade as
            # completed evidence. No automatic retries change the planned draw set.
            metadata = json.loads(manifest_path.read_text())
            metadata.update(status="incomplete", termination_type=type(exc).__name__,
                            termination_reason=str(exc), finished_utc=datetime.now(timezone.utc).isoformat())
            manifest_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
            raise

    def _execute(self, output, prereg_commit, progress, started, clock):
        source_hashes = {str(p.relative_to(PROJECT_ROOT)): _sha(p)
                         for p in sorted((PROJECT_ROOT / "ou_irregular").rglob("*.py"))}
        cert, boot = self.settings["certification"], self.settings["bootstrap"]
        cell = self.config.scientific_cell(cert)
        service = self.wald_factory(cert["estimator"])
        records = []
        for rep in range(cert["replications"]):
            sample, diagnostics, error = self._sample(cell, rep)
            result, error = self._infer(service, sample, error)
            records.append({**cell, "rep": rep, **diagnostics, **_wald_record(result),
                            "pipeline_error": error})
            if (rep + 1) % 250 == 0 or rep + 1 == cert["replications"]:
                pd.DataFrame(records).to_csv(output / "certification.csv", index=False)
                if progress:
                    progress(f"Checkpoint: {rep + 1}/{cert['replications']} planned paths retained")
        certification = pd.DataFrame(records)
        coverage = self.summarizer.summarize(certification, planned_replications=cert["replications"])
        checks, draws, cells = [], [], {cell["cell_id"]: cell}
        for specification in boot["cells"]:
            source_cell = self.config.scientific_cell(specification)
            cells[source_cell["cell_id"]] = source_cell
            sample, diagnostics, source_error = self._sample(source_cell, boot["sample_rep"])
            for name in boot["estimators"]:
                wald, error = self._infer(self.wald_factory(name), sample, source_error)
                seed = bootstrap_seed(boot["seed"], source_cell["cell_id"], name)
                try:
                    if sample is None:
                        raise ValueError(f"source path unavailable: {source_error}")
                    result = self.bootstrap_factory(name).run(sample, point_fit=wald.point_fit, seed=seed)
                except Exception as exc:
                    reason = f"{type(exc).__name__}: {exc}"
                    error = "; ".join(x for x in (error, reason) if x)
                    result = BootstrapResult(wald.point_fit, False, "pipeline_exception", reason,
                                             np.full(3, np.nan), np.full(3, np.nan), np.full(3, np.nan),
                                             (), boot["draws"], 0, {})
                common = {**source_cell, "source_rep": boot["sample_rep"],
                          **diagnostics, "estimator": name, "bootstrap_seed": seed}
                check = {**common, "point_valid": bool(wald.point_fit.valid),
                         "point_status": wald.point_fit.status, "point_reason": wald.point_fit.reason,
                         "wald_valid": bool(wald.valid), "wald_status": wald.status,
                         "wald_reason": wald.reason, "wald_diagnostics_json": _json(wald.diagnostics),
                         "bootstrap_valid": bool(result.valid), "bootstrap_status": result.status,
                         "bootstrap_reason": result.reason, "planned_draws": result.planned_count,
                         "attempted_draws": result.attempted_count, "valid_draws": result.valid_count,
                         "invalid_attempted_draws": result.invalid_count,
                         "bootstrap_diagnostics_json": _json(result.diagnostics), "pipeline_error": error}
                for i, parameter in enumerate(PARAMETERS):
                    check[parameter] = wald.point_fit.parameters[i]
                    check[f"wald_{parameter}_lower"] = wald.lower[i]
                    check[f"wald_{parameter}_upper"] = wald.upper[i]
                    check[f"wald_{parameter}_se"] = np.sqrt(wald.covariance[i, i])
                    check[f"bootstrap_{parameter}_lower"] = result.lower[i]
                    check[f"bootstrap_{parameter}_upper"] = result.upper[i]
                    check[f"bootstrap_{parameter}_se"] = result.standard_errors[i]
                checks.append(check)
                for record in result.draw_records:
                    # Complete per-draw records; no resampling to replace failures.
                    row = {**common, **dict(record)}
                    row["seed_spawn_key"] = _json(row["seed_spawn_key"])
                    draws.append(row)
                pd.DataFrame(checks).to_csv(output / "bootstrap_summary.csv", index=False)
                if draws:
                    pd.DataFrame(draws).to_csv(output / "bootstrap_draws.csv", index=False)
                if progress:
                    progress(f"Bootstrap CV={source_cell['cv']:g} {name}: "
                             f"{result.valid_count}/{result.planned_count} valid; {result.status}")
        bootstrap_summary = pd.DataFrame(checks)
        if draws:
            bootstrap_draws = pd.DataFrame(draws)
        else:
            bootstrap_draws = pd.DataFrame(columns=["cell_id", "source_rep", "estimator", "draw",
                                                   "seed_entropy", "seed_spawn_key", *PARAMETERS,
                                                   "nll", "point_valid", "optimizer_success", "status", "reason",
                                                   "dependency_exception", "exception_type"])
        tables = {"certification.csv": certification, "coverage_summary.csv": coverage,
                  "bootstrap_summary.csv": bootstrap_summary, "bootstrap_draws.csv": bootstrap_draws,
                  "inference_failures.csv": certification.loc[~certification.interval_valid]}
        for filename, data in tables.items():
            data.to_csv(output / filename, index=False)
        (output / "resolved_config.json").write_text(_json(self.settings) + "\n")
        (output / "cells.json").write_text(_json(list(cells.values())) + "\n")
        self.figure_writer.write(coverage, bootstrap_summary, output / "figures")
        draw_exceptions = int(bootstrap_draws.get("dependency_exception", pd.Series(dtype=bool)).sum())
        unexpected = int(certification.pipeline_error.ne("").sum() + bootstrap_summary.pipeline_error.ne("").sum()) + draw_exceptions
        frozen_path = PROJECT_ROOT / "configs/week4_validation.yaml"
        registered = (prereg_commit == PREREG_COMMIT
                      and _sha(PROJECT_ROOT / "preregistration.md") == PREREG_SHA256
                      and _sha(frozen_path) == FROZEN_CONFIG_SHA256
                      and self.config.digest == REGISTERED_CONFIG_DIGEST)
        metadata = {
            "project": "ou-irregular-sampling-week4", "schema_version": 1,
            "started_utc": started, "finished_utc": datetime.now(timezone.utc).isoformat(),
            "elapsed_seconds": time.perf_counter() - clock,
            "status": "implementation_review_required" if unexpected else "completed",
            "unexpected_pipeline_exceptions": unexpected,
            "preregistration_commit": prereg_commit,
            "preregistration_sha256": _sha(PROJECT_ROOT / "preregistration.md"),
            "frozen_config_file_sha256": _sha(frozen_path),
            "resolved_config_sha256": self.config.digest,
            "registered_configuration_matches": registered,
            "registered_procedure_matches": registered and not self.overridden_scientific_dependencies,
            "overridden_scientific_dependencies": self.overridden_scientific_dependencies,
            "study_role": "registered_week4_validation" if registered and not self.overridden_scientific_dependencies else "diagnostic_unregistered",
            "source_sha256": source_hashes,
            "python_version": platform.python_version(), "platform": platform.platform(),
            "dependencies": {name: version(name) for name in ("numpy", "scipy", "pandas", "matplotlib", "PyYAML")},
            "certification": {"planned_paths": cert["replications"], "retained_rows": len(certification),
                              "point_valid": int(certification.point_valid.sum()),
                              "interval_valid": int(certification.interval_valid.sum()),
                              "point_status_counts": certification.point_status.value_counts().to_dict(),
                              "interval_status_counts": certification.interval_status.value_counts().to_dict()},
            "bootstrap": {"planned_checks": len(checks), "planned_draws": int(bootstrap_summary.planned_draws.sum()),
                          "attempted_draws": len(bootstrap_draws), "valid_draws": int(bootstrap_summary.valid_draws.sum()),
                          "dependency_exceptions": draw_exceptions,
                          "status_counts": bootstrap_summary.bootstrap_status.value_counts().to_dict(),
                          "draw_status_counts": bootstrap_draws.status.value_counts().to_dict()},
            "checkpoint_reference_band": [0.94, 0.96],
            "checkpoint_in_reference_band": dict(zip(coverage.parameter, coverage.checkpoint_in_reference_band)),
            "interpretation": "Reference-band misses trigger investigation, never seed or interval retuning; see Week 4 progress report.",
            "headline_experiments_run": False,
            "artifact_sha256": {str(p.relative_to(output)): _sha(p) for p in sorted(output.rglob("*"))
                                if p.is_file() and p.name != "run_metadata.json"},
        }
        # Source mutation during execution invalidates a single-version provenance claim.
        if source_hashes != {str(p.relative_to(PROJECT_ROOT)): _sha(p)
                             for p in sorted((PROJECT_ROOT / "ou_irregular").rglob("*.py"))}:
            raise RuntimeError("source changed during execution; retain outputs as incomplete evidence")
        (output / "run_metadata.json").write_text(json.dumps(_json_safe(metadata), indent=2, sort_keys=True) + "\n")
        return metadata


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    run = commands.add_parser("run", help="Run checkpoint, bootstrap diagnostics, raw exports and figures")
    run.add_argument("--config", required=True, type=Path)
    run.add_argument("--output", required=True, type=Path)
    run.add_argument("--prereg-commit", default=PREREG_COMMIT)
    args = parser.parse_args(argv)
    metadata = ValidationExperiment(load_config(args.config)).run(
        args.output, prereg_commit=args.prereg_commit, progress=lambda message: print(message, flush=True))
    print(f"Saved {args.output}: {metadata['status']}", flush=True)
    return 1 if metadata["unexpected_pipeline_exceptions"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
