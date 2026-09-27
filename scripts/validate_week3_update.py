"""Reconcile a Week 3 solver update on the original 180 paths and 540 fit keys.

This is an audit utility, not an experiment selector. It never changes seeds,
truth, grids, or archived results. Independent references live in tests and do
not call production fitting code. Gates concern numerical agreement with the
same objective, never closeness of estimates to the simulation truth.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tests"))
from ou_irregular.config import load_config
from ou_irregular.runner import array_digest, replication_streams
from ou_irregular.ou import simulate_ou
from ou_irregular.spacing import generate_times
from reference_estimators import (reference_fit, conditional_nll, ReferenceInvalid,
                                  PARAMETER_RTOL, PARAMETER_ATOL, NLL_ATOL, NLL_RTOL)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-run", type=Path, required=True)
    parser.add_argument("--updated-run", type=Path, required=True)
    parser.add_argument("--evidence", type=Path, required=True)
    args = parser.parse_args()
    frozen = json.loads((args.evidence / "baseline_provenance.json").read_text())
    for name, expected in frozen["baseline_artifacts_sha256"].items():
        actual = hashlib.sha256((args.baseline_run / name).read_bytes()).hexdigest()
        if actual != expected:
            raise ValueError(f"Historical artifact changed: {name}")
    before = pd.read_csv(args.baseline_run / "replications.csv", float_precision="round_trip")
    after = pd.read_csv(args.updated_run / "replications.csv", float_precision="round_trip")
    keys = ["cell_id", "rep", "estimator"]
    if before.duplicated(keys).any() or after.duplicated(keys).any():
        raise ValueError("Duplicate replication keys")
    if set(map(tuple, before[keys].values)) != set(map(tuple, after[keys].values)):
        raise ValueError("Changed planned replication keys")
    old_config = json.loads((args.baseline_run / "resolved_config.json").read_text())
    new_config = json.loads((args.updated_run / "resolved_config.json").read_text())
    if old_config != new_config:
        raise ValueError("The original smoke configuration must be preserved exactly")
    inputs = pd.read_csv(args.evidence / "baseline_inputs.csv")
    spacing = pd.read_csv(args.updated_run / "spacing.csv")
    input_check = inputs.merge(spacing, on=["cell_id", "rep"], suffixes=("_baseline", "_updated"),
                               how="outer", validate="one_to_one", indicator=True)
    identical = input_check["_merge"].eq("both")
    for column in ("times_sha256", "observations_sha256"):
        identical &= input_check[f"{column}_baseline"].eq(input_check[f"{column}_updated"])
    if len(input_check) != len(inputs) or not identical.all():
        raise ValueError("Changed input arrays; seeds alone are not sufficient evidence")

    config = load_config(args.updated_run / "config.yaml")
    cells = {cell["cell_id"]: cell for cell in config.cells()}
    new_index = after.set_index(keys)
    comparisons = []
    for (cell_id, rep), old_group in before.groupby(["cell_id", "rep"], sort=True):
        cell = cells[cell_id]
        gap_rng, path_rng, _ = replication_streams(config.seed, cell_id, int(rep))
        times = generate_times(cell["n"], cell["mean_gap"], cell["cv"], gap_rng,
                               floor_fraction=cell["floor_fraction"],
                               floor_cv_threshold=cell["floor_cv_threshold"]).times
        x = simulate_ou(cell["true_theta"], cell["true_mu"], cell["true_sigma"], times, path_rng)
        state_scale, time_scale = np.std(x), np.mean(np.diff(times))

        def coordinates(parameters):
            theta, mu, sigma = parameters
            return np.array([theta * time_scale, (mu - np.mean(x)) / state_scale,
                             sigma * np.sqrt(time_scale) / state_scale])

        for _, old in old_group.iterrows():
            estimator = old.estimator
            new = new_index.loc[(cell_id, rep, estimator)]
            if new.times_sha256 != array_digest(times) or new.observations_sha256 != array_digest(x):
                raise ValueError("Reference audit regenerated different arrays")
            old_parameters = np.array([old[k] for k in ("theta", "mu", "sigma")])
            new_parameters = np.array([new[k] for k in ("theta", "mu", "sigma")])
            row = {"cell_id": cell_id, "rep": rep, "estimator": estimator,
                   "theta_mean_gap": cell["theta_mean_gap"], "cv": cell["cv"],
                   "before_success": bool(old.success), "before_policy": "legacy_success_finite_v1",
                   "after_optimizer_success": bool(new.optimizer_success),
                   "after_fit_status": new.fit_status, "after_valid": bool(new.valid_for_point_summary),
                   "after_reason": new.reason, "after_solver": new.solver,
                   "before_nll": old.nll, "after_nll": new.nll, "nll_change": new.nll - old.nll,
                   "before_status": old.status, "before_message": old.message}
            for k in ("theta", "mu", "sigma"):
                row.update({f"before_{k}": old[k], f"after_{k}": new[k], f"change_{k}": new[k] - old[k]})
            try:
                ref = reference_fit(estimator, x, times)
            except ReferenceInvalid as exc:
                row.update(reference_status="no_interior_reference", reference_reason=str(exc))
                if new.valid_for_point_summary:
                    raise AssertionError(f"Valid production fit lacks independent interior reference: {row}")
            else:
                reference_coordinates = coordinates(ref.parameters)
                difference = np.abs(coordinates(new_parameters) - reference_coordinates)
                tolerance = PARAMETER_ATOL + PARAMETER_RTOL * np.abs(reference_coordinates)
                nll_tolerance = NLL_ATOL + NLL_RTOL * abs(ref.nll)
                candidate_nll = conditional_nll(estimator, new_parameters, x, times)
                old_nll = conditional_nll(estimator, old_parameters, x, times)
                row.update(reference_status="valid_interior", reference_nll=ref.nll,
                           before_nll_recheck=old_nll, after_nll_recheck=candidate_nll,
                           before_nll_gap=old.nll - ref.nll, after_nll_gap=new.nll - ref.nll,
                           nll_tolerance=nll_tolerance,
                           max_parameter_tolerance_fraction=float(np.max(difference / tolerance)),
                           max_before_parameter_tolerance_fraction=float(np.max(np.abs(coordinates(old_parameters) - reference_coordinates) / tolerance)))
                assert new.valid_for_point_summary, row
                assert np.all(difference <= tolerance), row
                assert abs(new.nll - ref.nll) <= nll_tolerance, row
                assert abs(candidate_nll - ref.nll) <= nll_tolerance, row
                row["material_baseline_discrepancy"] = (abs(old.nll - ref.nll) > nll_tolerance
                                                         or row["max_before_parameter_tolerance_fraction"] > 1)
            comparisons.append(row)

    frame = pd.DataFrame(comparisons)
    frame.to_csv(args.evidence / "before_after.csv", index=False)
    report = {"baseline_commit": frozen["baseline_commit"],
              "configuration_unchanged": True, "historical_artifacts_unchanged": True,
              "identical_input_pairs": len(input_check), "matched_fit_rows": len(frame),
              "after_outcomes": frame.after_fit_status.value_counts().to_dict(),
              "after_max_absolute_nll_reference_difference": float(frame.after_nll_gap.abs().max()),
              "after_max_parameter_tolerance_fraction": float(frame.max_parameter_tolerance_fraction.max()),
              "material_baseline_discrepancies": int(frame.material_baseline_discrepancy.fillna(False).sum()),
              "nll_change_min": float(frame.nll_change.min()), "nll_change_max": float(frame.nll_change.max()),
              "parameter_rtol": PARAMETER_RTOL, "parameter_atol": PARAMETER_ATOL,
              "nll_atol": NLL_ATOL, "nll_rtol": NLL_RTOL}
    (args.evidence / "reconciliation.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
