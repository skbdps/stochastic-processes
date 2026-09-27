"""Run and replay the Week 3 OU experiment from one YAML configuration.

Each cell/replication gets its own SeedSequence, then independent gap/path
children. Estimators share the same generated observations for a fair comparison;
optimizer order or retries consume no simulation randomness. Cell identities use
SHA-256 rather than Python's process-randomized hash(). See NumPy's primary docs:
https://numpy.org/doc/stable/reference/random/parallel.html
"""
import argparse
from datetime import datetime, timezone
import hashlib
from importlib import metadata
import json
import os
from pathlib import Path
import platform
import shutil
import time
import warnings

# The likelihoods operate on small vectors. Avoid BLAS thread oversubscription
# in the command-line process, consistently with the Week 2 reproduction runner.
for _variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_variable] = "1"

import numpy as np
import pandas as pd
import yaml

from .config import ExperimentConfig, canonical_json, load_config
from .estimators import fit_estimator
from .metrics import (add_twin_comparisons, summarize_results, point_validity,
                      validity_counts, paired_estimator_comparisons)
from .ou import simulate_ou
from .plots import plot_summary
from .spacing import generate_times


PROJECT = "ou-irregular-sampling"


def _json_write(path, value):
    # Atomic replacement prevents a truncated JSON file being mistaken for a run.
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
    temporary.replace(path)


def _csv_write(path, frame):
    temporary = path.with_suffix(path.suffix + ".tmp")
    frame.to_csv(temporary, index=False)
    temporary.replace(path)


def replication_streams(seed, cell_id, rep):
    """Reconstruct experiment -> cell -> replication -> gap/path seed branches.

    The fixed-length SHA-256 words identify the complete data-generating cell.
    Replication IDs occupy a separate final branch coordinate. This explicit
    branch address makes existing draws independent of grid traversal, additional
    cells, estimator order and total replication count. The two spawned children
    separate observation-time generation from Brownian innovations.
    """
    if isinstance(rep, bool) or not isinstance(rep, (int, np.integer)) or rep < 0:
        raise ValueError("rep must be a nonnegative integer")
    words = tuple(int(cell_id[i:i + 8], 16) for i in range(0, 64, 8))
    sequence = np.random.SeedSequence(seed, spawn_key=(*words, int(rep)))
    gap_seed, path_seed = sequence.spawn(2)
    diagnostics = {"seed_entropy": str(seed), "seed_spawn_key": canonical_json(sequence.spawn_key),
                   "gap_spawn_key": canonical_json(gap_seed.spawn_key),
                   "path_spawn_key": canonical_json(path_seed.spawn_key), "bit_generator": "PCG64"}
    return np.random.default_rng(gap_seed), np.random.default_rng(path_seed), diagnostics


def _failed_record(estimator, message):
    return {"estimator": estimator, "theta": np.nan, "mu": np.nan, "sigma": np.nan,
            "success": False, "initial_success": False, "retried": False,
            "optimizer_success": False, "valid_for_point_summary": False,
            "fit_status": "simulation_failure", "reason": message, "solver": "not_run",
            "attempted_starts": 0, "nll": np.nan, "status": -1, "message": message,
            "error": message}


def array_digest(values):
    """Identify the actual input array, not just the random seed that made it.

    Canonical little-endian float64 bytes make hashes independent of native
    byte order. This catches changed simulation code or lost input resolution.
    """
    return hashlib.sha256(np.asarray(values, dtype="<f8").tobytes()).hexdigest()


def run_replication(config, cell, rep):
    """Simulate once and fit all estimators, preserving every planned result row.

    Estimation uses the conditional likelihood given x0 although the simulator
    draws x0 from the stationary law. Data always use the exact OU recursion;
    the Euler approximation is applied only to its estimator's likelihood.
    """
    if rep >= config.replications:
        raise ValueError("rep is outside the configured replication range")
    gap_rng, path_rng, seeds = replication_streams(config.seed, cell["cell_id"], rep)
    common = {**cell, "rep": int(rep), **seeds, "phase": config.phase}
    try:
        spacing = generate_times(cell["n"], cell["mean_gap"], cell["cv"], gap_rng,
                                 floor_fraction=cell["floor_fraction"],
                                 floor_cv_threshold=cell["floor_cv_threshold"])
        observations = simulate_ou(cell["true_theta"], cell["true_mu"], cell["true_sigma"],
                                   spacing.times, path_rng)
        if observations.shape != (cell["n"],) or not np.isfinite(observations).all():
            raise FloatingPointError("simulator returned nonfinite or incorrectly shaped observations")
        diagnostics = {**spacing.diagnostics, "simulation_success": True, "simulation_error": "",
                       "times_sha256": array_digest(spacing.times),
                       "observations_sha256": array_digest(observations)}
    except (ValueError, FloatingPointError, RuntimeError, OverflowError) as error:
        # A failed simulation is visible in every estimator's denominator. It
        # is not silently replaced with a fresh seed or omitted from the cell.
        message = f"simulation failed: {type(error).__name__}: {error}"
        spacing_row = {**common, "simulation_success": False, "simulation_error": message}
        return [{**common, **_failed_record(e, message)} for e in config.estimators], spacing_row
    fits = []
    for estimator in config.estimators:
        result = fit_estimator(estimator, observations, spacing.times, **config.optimizer)
        fits.append({**common, **diagnostics, **result})
    return fits, {**common, **diagnostics}


def _provenance(config):
    package = Path(__file__).resolve().parent
    versions = {name: metadata.version(name) for name in ("numpy", "scipy", "pandas", "matplotlib", "PyYAML")}
    return {"project": PROJECT, "schema_version": 1, "result_schema_version": 2,
            "phase": config.phase,
            "config_sha256": hashlib.sha256(canonical_json(config.to_dict()).encode()).hexdigest(),
            "python": platform.python_version(), "platform": platform.platform(),
            "versions": versions,
            "dependency_versions_sha256": hashlib.sha256(canonical_json(versions).encode()).hexdigest(),
            "source_sha256": {
                path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in sorted(package.glob("*.py"))}}


def provenance_differences(saved, current):
    """Compare scientific code and runtime versions, not machine-specific paths.

    Dependency fingerprints cover version identities, not installed binary-file
    contents. Platform is recorded separately; bitwise cross-platform equality
    is not promised. Missing historical evidence also counts as a mismatch.
    """
    differences = {}
    for key in ("source_sha256", "versions", "python"):
        if saved.get(key) != current.get(key):
            differences[key] = {"saved": saved.get(key), "executing": current.get(key)}
    return differences


def _validated_saved_rows(run_dir, run_metadata):
    """Check the planned design before using a CSV's row count as denominator.

    A truncated CSV must never masquerade as a smaller, successful experiment.
    Check both the full Cartesian key set and every cell's generating settings.
    If a completed manifest records a raw-data digest, verify that as well.
    """
    config = ExperimentConfig.from_dict(json.loads((run_dir / "resolved_config.json").read_text()))
    digest = hashlib.sha256(canonical_json(config.to_dict()).encode()).hexdigest()
    if digest != run_metadata.get("config_sha256"):
        raise ValueError("Resolved configuration does not match the run manifest")
    raw_path = run_dir / "replications.csv"
    recorded = run_metadata.get("output_sha256", {}).get("replications.csv")
    if recorded and hashlib.sha256(raw_path.read_bytes()).hexdigest() != recorded:
        raise ValueError("Raw replication CSV hash does not match the completed run")
    frame = pd.read_csv(raw_path, float_precision="round_trip")
    keys = ["cell_id", "rep", "estimator"]
    if not set(keys).issubset(frame) or frame.duplicated(keys).any():
        raise ValueError("Replication keys are missing or duplicated")
    cells = list(config.cells())
    expected = {(cell["cell_id"], rep, estimator) for cell in cells
                for rep in range(config.replications) for estimator in config.estimators}
    if set(frame[keys].itertuples(index=False, name=None)) != expected:
        raise ValueError("Replication CSV does not contain every planned cell/rep/estimator")
    for cell in cells:
        selected = frame.loc[frame.cell_id == cell["cell_id"]]
        for key, value in cell.items():
            if key not in selected or not selected[key].eq(value).all():
                raise ValueError(f"Replication metadata differs from config: {key}")
    return config, frame


def replot(run_dir, output_dir=None):
    """Recompute aggregates and figures solely from saved replication rows.

    A changed implementation must write a separate derived run. Old rows keep
    their original success/finite policy; no new boundary diagnosis is invented
    without a refit. Classified rows use the shared point-validity policy.
    """
    run_dir = Path(run_dir)
    run_metadata = json.loads((run_dir / "run_metadata.json").read_text())
    if run_metadata.get("project") != PROJECT:
        raise ValueError("Not an OU experiment output directory")
    config, frame = _validated_saved_rows(run_dir, run_metadata)
    current = _provenance(config)
    differences = provenance_differences(run_metadata, current)
    destination = run_dir if output_dir is None else Path(output_dir)
    if differences and destination.resolve() == run_dir.resolve():
        raise ValueError("Replot provenance mismatch: use a separate --output directory to preserve historical outputs")
    if destination.resolve() != run_dir.resolve():
        if destination.exists() and any(destination.iterdir()):
            raise FileExistsError("Derived output must be a new or empty directory")
        destination.mkdir(parents=True, exist_ok=True)
        for name in ("replications.csv", "resolved_config.json", "config.yaml", "cells.json",
                     "spacing.csv", "spacing_summary.csv"):
            if (run_dir / name).exists():
                shutil.copyfile(run_dir / name, destination / name)
        run_metadata["derived_from"] = {
            "run_directory": str(run_dir),
            "manifest_sha256": hashlib.sha256((run_dir / "run_metadata.json").read_bytes()).hexdigest(),
            "provenance_differences": differences}
    run_dir = destination
    summary = add_twin_comparisons(summarize_results(frame))
    _csv_write(run_dir / "summary.csv", summary)
    _csv_write(run_dir / "failures.csv", frame.loc[~point_validity(frame)])
    # Boundaries may be numerical successes, so give them their own export as
    # well as retaining them among all point-summary exclusions.
    boundary = frame["fit_status"].isin(["boundary", "boundary_suspected"]) if "fit_status" in frame else pd.Series(False, index=frame.index)
    _csv_write(run_dir / "boundary_fits.csv", frame.loc[boundary])
    _csv_write(run_dir / "estimator_comparisons.csv", paired_estimator_comparisons(frame))
    paths = plot_summary(summary, run_dir / "figures")
    # Raw-estimate provenance stays unchanged when only aggregation/plots rerun.
    # Record current derived-code hashes and refresh digests of regenerated files.
    package = Path(__file__).resolve().parent
    run_metadata["derived_source_sha256"] = {
        name: hashlib.sha256((package / name).read_bytes()).hexdigest()
        for name in ("runner.py", "metrics.py", "plots.py")}
    run_metadata["derived_environment"] = {key: current[key] for key in ("python", "versions", "dependency_versions_sha256")}
    run_metadata["validity_counts"] = validity_counts(frame)
    hashes = run_metadata.setdefault("output_sha256", {})
    for path in [*run_dir.glob("*.csv"), *paths]:
        hashes[str(path.relative_to(run_dir))] = hashlib.sha256(path.read_bytes()).hexdigest()
    _json_write(run_dir / "run_metadata.json", run_metadata)
    return summary, paths


def run_experiment(config, output_dir, *, overwrite=False):
    """Serialize a complete smoke run; failed fits are results, not lost rows.

    Existing unrelated/nonempty directories are never overwritten. Repeating an
    identified OU run requires the explicit overwrite flag. No directory tree
    is deleted; only the known generated files are replaced.
    """
    output_dir = Path(output_dir)
    if output_dir.exists() and any(output_dir.iterdir()):
        known = output_dir / "run_metadata.json"
        if not overwrite or not known.exists() or json.loads(known.read_text()).get("project") != PROJECT:
            raise FileExistsError("Output is nonempty; choose a new directory or --overwrite for an existing OU run")
    output_dir.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    run_metadata = _provenance(config)
    run_metadata.update(status="running", started_utc=datetime.now(timezone.utc).isoformat())
    _json_write(output_dir / "run_metadata.json", run_metadata)
    _json_write(output_dir / "resolved_config.json", config.to_dict())
    (output_dir / "config.yaml").write_text(yaml.safe_dump(config.to_dict(), sort_keys=False))
    cells = list(config.cells())
    _json_write(output_dir / "cells.json", cells)
    rows, spacing_rows = [], []
    for cell in cells:
        print(f"cell {cell['cell_id'][:12]}: n={cell['n']}, theta*gap={cell['theta_mean_gap']}, CV={cell['cv']}", flush=True)
        for rep in range(config.replications):
            fits, spacing = run_replication(config, cell, rep)
            rows.extend(fits)
            spacing_rows.append(spacing)
    frame = pd.DataFrame(rows)
    _csv_write(output_dir / "replications.csv", frame)
    _csv_write(output_dir / "spacing.csv", pd.DataFrame(spacing_rows))
    failures = frame.loc[~point_validity(frame)]
    _csv_write(output_dir / "failures.csv", failures)
    summary, paths = replot(output_dir)
    spacing_frame = pd.DataFrame(spacing_rows)
    spacing_summary = spacing_frame.groupby("cell_id", sort=True).agg(
        replications=("rep", "size"), successful_simulations=("simulation_success", "sum"))
    for column in ("clip_fraction", "expected_clip_fraction", "realized_mean_gap", "expected_mean_gap_after_floor",
                   "raw_mean_gap", "realized_cv", "realized_span", "min_gap"):
        if column in spacing_frame:
            spacing_summary[f"mean_{column}"] = spacing_frame.groupby("cell_id")[column].mean()
    _csv_write(output_dir / "spacing_summary.csv", spacing_summary.reset_index())
    run_metadata = json.loads((output_dir / "run_metadata.json").read_text())
    counts = validity_counts(frame)
    run_metadata.update(status="completed" if failures.empty else "completed_with_excluded_fits",
                        cells=len(cells), replications_per_cell=config.replications,
                        planned_fits=len(cells) * config.replications * len(config.estimators),
                        recorded_fits=len(frame), failed_fits=len(failures),
                        excluded_fits=len(failures), validity_counts=counts,
                        ci_status="not_computed_week3", elapsed_seconds=round(time.perf_counter() - started, 3),
                        figures=[str(path.relative_to(output_dir)) for path in paths])
    # File hashes permit checking saved raw data before interpreting or replaying it.
    run_metadata.setdefault("output_sha256", {}).update({
        path.name: hashlib.sha256(path.read_bytes()).hexdigest()
        for path in sorted(output_dir.glob("*.csv"))})
    _json_write(output_dir / "run_metadata.json", run_metadata)
    return run_metadata


def replay(run_dir, cell_id, rep, *, allow_provenance_mismatch=False, output_dir=None):
    """Strict historical replay, or an explicitly recorded patched-code replay.

    Raw-file integrity does not prove historical code identity. By default a
    source/dependency mismatch stops before any fitting. An intentional replay
    with changed code requires a separate output directory and records both
    provenance records, the original raw hash and exact input hashes.
    """
    run_dir = Path(run_dir)
    run_metadata = json.loads((run_dir / "run_metadata.json").read_text())
    config, saved = _validated_saved_rows(run_dir, run_metadata)
    current = _provenance(config)
    differences = provenance_differences(run_metadata, current)
    if differences:
        if not allow_provenance_mismatch or output_dir is None:
            raise ValueError("Replay provenance mismatch: deliberate patched-code replay requires --allow-provenance-mismatch and a separate --output directory")
        warnings.warn("PATCHED-CODE REPLAY: executing source/runtime differs from the historical run", RuntimeWarning)
    if output_dir is not None:
        output_dir = Path(output_dir)
        if output_dir.resolve() == run_dir.resolve():
            raise ValueError("Replay output must be separate from the source run")
        if output_dir.exists() and any(output_dir.iterdir()):
            raise FileExistsError("Replay output must be a new or empty directory")
    matches = [cell for cell in config.cells() if cell["cell_id"] == cell_id]
    if len(matches) != 1:
        raise ValueError("cell_id must be an exact ID from cells.json")
    result = pd.DataFrame(run_replication(config, matches[0], rep)[0])
    historical = saved[(saved.cell_id == cell_id) & (saved.rep == rep)]
    input_identity = "not_recorded_in_historical_run"
    if {"times_sha256", "observations_sha256"}.issubset(historical):
        hash_columns = ["times_sha256", "observations_sha256"]
        if historical[hash_columns].notna().all().all():
            input_identity = all(key in result and historical[key].iloc[0] == result[key].iloc[0]
                                 for key in hash_columns)
            if not input_identity:
                raise ValueError("Replay did not reconstruct identical input arrays")
        else:
            input_identity = "unavailable_for_historical_simulation_failure"
    if output_dir is not None:
        output_dir.mkdir(parents=True, exist_ok=True)
        _csv_write(output_dir / "replications.csv", result)
        _json_write(output_dir / "replay_metadata.json", {
            "project": PROJECT, "kind": "patched_code_replay" if differences else "strict_replay",
            "source_run": str(run_dir), "cell_id": cell_id, "rep": rep,
            "source_manifest_sha256": hashlib.sha256((run_dir / "run_metadata.json").read_bytes()).hexdigest(),
            "source_raw_sha256": hashlib.sha256((run_dir / "replications.csv").read_bytes()).hexdigest(),
            "saved_provenance": run_metadata, "executing_provenance": current,
            "provenance_differences": differences, "input_identity": input_identity,
            "output_sha256": hashlib.sha256((output_dir / "replications.csv").read_bytes()).hexdigest()})
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    run = commands.add_parser("run", help="one YAML -> raw estimates -> summaries -> figures")
    run.add_argument("--config", required=True, type=Path)
    run.add_argument("--output", required=True, type=Path)
    run.add_argument("--overwrite", action="store_true")
    plotting = commands.add_parser("replot", help="rebuild summaries and figures from saved rows")
    plotting.add_argument("--run-dir", required=True, type=Path)
    plotting.add_argument("--output", type=Path, help="separate derived run; required if provenance differs")
    repeating = commands.add_parser("replay", help="regenerate one replication as CSV on stdout")
    repeating.add_argument("--run-dir", required=True, type=Path)
    repeating.add_argument("--cell-id", required=True)
    repeating.add_argument("--rep", required=True, type=int)
    repeating.add_argument("--allow-provenance-mismatch", action="store_true")
    repeating.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    if arguments.command == "run":
        print(json.dumps(run_experiment(load_config(arguments.config), arguments.output,
                                        overwrite=arguments.overwrite), indent=2))
    elif arguments.command == "replot":
        summary, paths = replot(arguments.run_dir, arguments.output)
        print(f"Regenerated {len(summary)} summary rows and {len(paths)} figures")
    else:
        print(replay(arguments.run_dir, arguments.cell_id, arguments.rep,
                     allow_provenance_mismatch=arguments.allow_provenance_mismatch,
                     output_dir=arguments.output).to_csv(index=False), end="")


if __name__ == "__main__":
    main()
