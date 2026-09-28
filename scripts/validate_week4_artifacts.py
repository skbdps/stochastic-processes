"""Read-only audit of the registered Week 4 evidence, with a separate report.

Checks provenance, planned keys, seeded input replay and independently recomputed
summary arithmetic. It does not refit models, recompute Hessians, certify nominal
coverage, or establish a theorem about an estimator. Historical files are never
rewritten. Run from any directory with the repository's pinned environment.
"""
import argparse
from datetime import datetime, timezone
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import platform
import sys

import numpy as np
import pandas as pd
from scipy.stats import norm

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from ou_irregular.week4.simulation import OuSampleSource  # noqa: E402

PARAMETERS = ("theta", "mu", "sigma")
FREEZE_COMMIT = "1b1255b0650ad7efcf700316bd1cccf704af4d2d"
PREREG_SHA = "d7b1a7b6a6afa6351d3115ade91100402f846eb931952b394ecec4b68dd2813a"
CONFIG_SHA = "a02f99aa82898e22dd4dde76f9b1eb12b4acc940e7f695cab0d170e1ee34ac25"
RESOLVED_SHA = "66844e989d1e9dfddcb13fe8eef5c3d4fad70c7a01215e34e560a2592196a24b"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def table(path):
    return pd.read_csv(path, float_precision="round_trip")


def numeric(actual, expected, message):
    require(bool(np.allclose(actual, expected, rtol=2e-12, atol=2e-14, equal_nan=True)), message)


def booleans(frame, names):
    for name in names:
        require(frame[name].map(lambda x: isinstance(x, (bool, np.bool_))).all(),
                f"{name} must contain actual Boolean values")


def expected_cell(settings, specification):
    truth, floor = settings["truth"], settings["gap_floor"]
    cell = {"n": int(specification["n"]), "theta_mean_gap": float(specification["q"]),
            "mean_gap": float(specification["q"] / truth["theta"]), "cv": float(specification["cv"]),
            **{f"true_{key}": float(value) for key, value in truth.items()},
            "floor_fraction": float(floor["fraction"]), "floor_cv_threshold": float(floor["cv_threshold"])}
    return {"cell_id": hashlib.sha256(canonical(cell).encode()).hexdigest(), **cell}


def wilson(successes, total):
    if not total:
        return np.nan, np.nan
    p, z = successes / total, norm.ppf(.975)
    center = (p + z * z / (2 * total)) / (1 + z * z / total)
    radius = z * np.sqrt(p * (1 - p) / total + z * z / (4 * total ** 2)) / (1 + z * z / total)
    return max(0., min(p, center - radius)), min(1., max(p, center + radius))


def audit(run_dir):
    run_dir = Path(run_dir).resolve()
    manifest_path = run_dir / "run_metadata.json"
    manifest_before = sha(manifest_path)
    metadata = json.loads(manifest_path.read_text())
    require(metadata["status"] == "completed", "run must have completed")
    require(metadata["registered_configuration_matches"] and metadata["registered_procedure_matches"],
            "run does not claim the frozen procedure/configuration")
    require(metadata["preregistration_commit"] == FREEZE_COMMIT, "unexpected preregistration commit")
    require(sha(ROOT / "preregistration.md") == metadata["preregistration_sha256"] == PREREG_SHA,
            "preregistration bytes differ from the pre-execution freeze")
    require(sha(ROOT / "configs/week4_validation.yaml") == metadata["frozen_config_file_sha256"] == CONFIG_SHA,
            "frozen configuration bytes changed")
    settings = json.loads((run_dir / "resolved_config.json").read_text())
    require(hashlib.sha256(canonical(settings).encode()).hexdigest()
            == metadata["resolved_config_sha256"] == RESOLVED_SHA, "resolved configuration differs from freeze")
    actual_sources = {str(p.relative_to(ROOT)): sha(p) for p in sorted((ROOT / "ou_irregular").rglob("*.py"))}
    require(actual_sources == metadata["source_sha256"], "executing source tree differs from recorded run")
    actual_artifacts = {str(p.relative_to(run_dir)): sha(p) for p in sorted(run_dir.rglob("*"))
                        if p.is_file() and p.name != "run_metadata.json"}
    require(actual_artifacts == metadata["artifact_sha256"], "artifact set or bytes differ from run manifest")
    require(platform.python_version() == metadata["python_version"], "input replay requires historical Python")
    require({name: version(name) for name in metadata["dependencies"]} == metadata["dependencies"],
            "input replay requires historical dependency versions")

    cert = table(run_dir / "certification.csv")
    coverage = table(run_dir / "coverage_summary.csv")
    checks = table(run_dir / "bootstrap_summary.csv")
    draws = table(run_dir / "bootstrap_draws.csv")
    cert_settings, boot = settings["certification"], settings["bootstrap"]
    certification_cell = expected_cell(settings, cert_settings)
    boot_cells = [expected_cell(settings, s) for s in boot["cells"]]
    expected_cells = {c["cell_id"]: c for c in [certification_cell, *boot_cells]}
    saved_cells = json.loads((run_dir / "cells.json").read_text())
    require(len(saved_cells) == len(expected_cells)
            and {c["cell_id"]: c for c in saved_cells} == expected_cells, "cell identities differ from configuration")
    cert_keys = ["cell_id", "rep", "estimator"]
    require(not cert.duplicated(cert_keys).any(), "duplicate certification keys")
    require(set(cert[cert_keys].itertuples(index=False, name=None)) == {
        (certification_cell["cell_id"], rep, cert_settings["estimator"])
        for rep in range(cert_settings["replications"])}, "certification did not retain all planned keys")
    check_keys = ["cell_id", "source_rep", "estimator"]
    expected_checks = {(c["cell_id"], boot["sample_rep"], name) for c in boot_cells for name in boot["estimators"]}
    require(not checks.duplicated(check_keys).any()
            and set(checks[check_keys].itertuples(index=False, name=None)) == expected_checks, "bootstrap check keys differ")
    draw_keys = [*check_keys, "draw"]
    require(not draws.duplicated(draw_keys).any()
            and set(draws[draw_keys].itertuples(index=False, name=None)) == {
                (*key, draw) for key in expected_checks for draw in range(boot["draws"])},
            "bootstrap draw keys are incomplete or duplicated")
    booleans(cert, ["point_valid", "interval_valid", "optimizer_success"])
    booleans(checks, ["point_valid", "wald_valid", "bootstrap_valid"])
    booleans(draws, ["point_valid", "optimizer_success", "dependency_exception"])
    require(not (cert.interval_valid & ~cert.point_valid).any(), "valid intervals have invalid points")
    require(not draws.dependency_exception.any() and metadata["unexpected_pipeline_exceptions"] == 0,
            "unexpected scientific dependency exception recorded")
    require(metadata["certification"]["planned_paths"] == metadata["certification"]["retained_rows"] == len(cert),
            "manifest certification denominator mismatch")
    require(metadata["certification"]["point_valid"] == int(cert.point_valid.sum())
            and metadata["certification"]["interval_valid"] == int(cert.interval_valid.sum()), "manifest validity mismatch")
    require(metadata["bootstrap"]["planned_checks"] == len(checks)
            and metadata["bootstrap"]["planned_draws"] == metadata["bootstrap"]["attempted_draws"] == len(draws),
            "manifest bootstrap denominator mismatch")

    # Replay ALL registered checkpoint inputs plus both fixed bootstrap paths.
    # This reconstructs historical inputs only; it does not refit or create new
    # experiment outcomes, and hashes must agree exactly rather than within tol.
    source = OuSampleSource(settings["seed"])
    for row in cert.itertuples(index=False):
        _, reproduced = source.sample(certification_cell, int(row.rep))
        for key in ("times_sha256", "observations_sha256", "source_spawn_key", "gap_spawn_key", "path_spawn_key"):
            require(getattr(row, key) == reproduced[key], f"checkpoint input replay mismatch at rep {row.rep}: {key}")
    for cell in boot_cells:
        _, reproduced = source.sample(cell, boot["sample_rep"])
        for frame in (checks, draws):
            subset = frame[frame.cell_id == cell["cell_id"]]
            for key in ("times_sha256", "observations_sha256", "source_spawn_key", "gap_spawn_key", "path_spawn_key"):
                require(subset[key].eq(reproduced[key]).all(), f"bootstrap source-path identity mismatch: {key}")

    require(not coverage.duplicated(["cell_id", "estimator", "parameter"]).any()
            and set(coverage.parameter) == set(PARAMETERS) and len(coverage) == 3,
            "coverage summary keys differ")
    coverage_counts = {}
    for row in coverage.itertuples(index=False):
        parameter = row.parameter
        truth = settings["truth"][parameter]
        valid = cert.interval_valid.to_numpy()
        point_valid = cert.point_valid.to_numpy()
        lower, upper = cert[f"{parameter}_lower"].to_numpy(), cert[f"{parameter}_upper"].to_numpy()
        require(np.isfinite(lower[valid]).all() and np.isfinite(upper[valid]).all()
                and (lower[valid] <= upper[valid]).all(), "valid interval has invalid endpoints")
        count = int((valid & (lower <= truth) & (truth <= upper)).sum())
        total, available = len(cert), int(valid.sum())
        require((row.n_planned, row.n_point_valid, row.n_interval_valid, row.n_covered)
                == (total, int(point_valid.sum()), available, count), f"coverage counts differ for {parameter}")
        numeric([row.coverage_operational, row.coverage_conditional],
                [count / total, count / available if available else np.nan], f"coverage proportions differ: {parameter}")
        numeric([row.operational_wilson_lower, row.operational_wilson_upper], wilson(count, total), "operational Wilson interval differs")
        numeric([row.conditional_wilson_lower, row.conditional_wilson_upper], wilson(count, available), "conditional Wilson interval differs")
        errors = cert[parameter].to_numpy()[point_valid] - truth
        numeric([row.bias, row.rmse, row.mcse_bias],
                [errors.mean(), np.sqrt(np.mean(errors ** 2)), errors.std(ddof=1) / np.sqrt(len(errors))],
                f"point summary arithmetic differs: {parameter}")
        coverage_counts[parameter] = {"covered": count, "planned": total, "valid_intervals": available}

    for check in checks.itertuples(index=False):
        group = draws[(draws.cell_id == check.cell_id) & (draws.estimator == check.estimator)
                      & (draws.source_rep == check.source_rep)]
        eligible = (group.point_valid & group.optimizer_success & group.status.eq("valid_interior")
                    & np.isfinite(group[list(PARAMETERS)]).all(axis=1) & group.theta.gt(0) & group.sigma.gt(0))
        require(group.point_valid.equals(eligible), "bootstrap point validity inconsistent with classification")
        require((check.planned_draws, check.attempted_draws, check.valid_draws, check.invalid_attempted_draws)
                == (boot["draws"], len(group), int(eligible.sum()), int((~eligible).sum())), "bootstrap counts differ")
        minimum = max(2, int(np.ceil(boot["min_valid_fraction"] * boot["draws"])))
        require(check.bootstrap_valid == (int(eligible.sum()) >= minimum), "bootstrap availability gate differs")
        require(group.seed_entropy.map(int).eq(int(check.bootstrap_seed)).all(), "bootstrap seed entropy mismatch")
        require(all(json.loads(row.seed_spawn_key) == [int(row.draw)] for row in group.itertuples(index=False)),
                "bootstrap draw stream labels differ")
        for p in PARAMETERS:
            values = group.loc[eligible, p].to_numpy()
            expected_bounds = np.quantile(values, [(1 - boot["level"]) / 2, (1 + boot["level"]) / 2], method="linear")
            numeric([getattr(check, f"bootstrap_{p}_lower"), getattr(check, f"bootstrap_{p}_upper")],
                    expected_bounds, f"bootstrap percentile endpoints differ: {check.cell_id}/{check.estimator}/{p}")
            numeric(getattr(check, f"bootstrap_{p}_se"), values.std(ddof=1), f"bootstrap sample SD differs: {p}")
    require(metadata["bootstrap"]["valid_draws"] == int(draws.point_valid.sum()), "manifest bootstrap-valid count differs")
    failures = table(run_dir / "inference_failures.csv")
    require(set(failures[cert_keys].itertuples(index=False, name=None))
            == set(cert.loc[~cert.interval_valid, cert_keys].itertuples(index=False, name=None)), "inference-failure export differs")

    week3 = {}
    for name in ("week3_smoke", "week3_validation_rerun"):
        path = ROOT / "artifacts" / name
        archived = json.loads((path / "run_metadata.json").read_text())
        for relative, expected in archived["output_sha256"].items():
            require(sha(path / relative) == expected, f"Week 3 artifact changed: {name}/{relative}")
        week3[name] = {"artifact_hashes_verified": len(archived["output_sha256"])}
        if name == "week3_validation_rerun":
            for relative, expected in archived["source_sha256"].items():
                require(sha(ROOT / "ou_irregular" / relative) == expected, f"Week 3 source changed: {relative}")
            week3[name]["source_hashes_verified"] = len(archived["source_sha256"])
    require(sha(manifest_path) == manifest_before, "audit changed the source manifest")
    require(actual_artifacts == {relative: sha(run_dir / relative) for relative in actual_artifacts},
            "audit changed a source artifact")
    return {"audit_passed": True, "audited_utc": datetime.now(timezone.utc).isoformat(),
            "script_sha256": sha(__file__), "source_manifest_sha256": manifest_before,
            "preregistration_commit": FREEZE_COMMIT, "freeze_hashes_verified": True,
            "source_hashes_verified": len(actual_sources), "artifact_hashes_verified": len(actual_artifacts),
            "certification_keys_verified": len(cert), "bootstrap_check_keys_verified": len(checks),
            "bootstrap_draw_keys_verified": len(draws), "checkpoint_input_replays": len(cert),
            "bootstrap_source_path_replays": len(boot_cells), "input_replay_hashes_exact": True,
            "coverage_counts": coverage_counts, "bootstrap_parameter_quantile_pairs_verified": len(checks) * 3,
            "bootstrap_sample_standard_deviations_verified": len(checks) * 3,
            "arithmetic_rtol": 2e-12, "arithmetic_atol": 2e-14, "week3_preservation": week3,
            "original_artifacts_unchanged": True, "new_model_fits_performed": 0,
            "scope": "Provenance, full planned keys, exact seeded input replay, raw-to-summary arithmetic, and Week 3 preservation. Does not independently refit models, recompute Hessians, or certify nominal coverage."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, default=ROOT / "artifacts/week4_validation")
    parser.add_argument("--output", type=Path, default=ROOT / "artifacts/week4_validation_evidence/audit.json")
    args = parser.parse_args()
    output, run_dir = args.output.resolve(), args.run_dir.resolve()
    require(not output.is_relative_to(run_dir), "audit output must be outside the original run directory")
    for historical_manifest in (ROOT / "artifacts").rglob("run_metadata.json"):
        require(not output.is_relative_to(historical_manifest.parent.resolve()),
                "audit output cannot overwrite any historical run directory")
    require(not output.is_relative_to(ROOT / "ou_irregular"), "audit output cannot overwrite source code")
    require(output.suffix == ".json", "audit output must be a JSON report")
    report = audit(run_dir)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print(f"Audit passed: {report['certification_keys_verified']} checkpoint rows, "
          f"{report['bootstrap_draw_keys_verified']} bootstrap draws; report {output}")


if __name__ == "__main__":
    main()
