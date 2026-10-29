"""Complete E2 classification and numerical-failure views of audited tables.

These companion reports live under figures/, the frozen workflow's derived
product namespace. They do not change the six core tables or their manifest.
The figure manifest binds their bytes and the core source-table hashes.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

from ..experiments.design import ESTIMATORS
from ..experiments.provenance import canonical, file_digest, verify
from ..experiments.storage import atomic_bytes, atomic_json


PARAMETERS = ("theta", "mu", "sigma")
CLASSES = ("within_criteria", "supported_breach", "unresolved", "not_applicable_zero_truth")
INPUTS = ("parameter_summary.csv", "outcome_counts.csv", "numerical_diagnostics.csv")


def load_tables(run_dir, design):
    """Require the full planned grid, even if a missing group leaves others valid."""
    root = Path(run_dir)
    metadata = json.loads((root / "run_metadata.json").read_text())
    verify(metadata["provenance"], expected_design=design)
    manifest = json.loads((root / "table_manifest.json").read_text())
    tables = {}
    for name in INPUTS:
        if file_digest(root / name) != manifest[name]["sha256"]:
            raise ValueError(f"grid report input digest mismatch: {name}")
        tables[name] = pd.read_csv(root / name, float_precision="round_trip")
    summary, outcomes, numerical = (tables[n] for n in INPUTS)
    expected = {(c.identity, m, p) for c in design.cells for m in ESTIMATORS for p in PARAMETERS}
    keys = list(summary[["cell_id", "estimator", "parameter"]].itertuples(index=False, name=None))
    if len(keys) != len(expected) or set(keys) != expected:
        raise ValueError("parameter report does not contain the full planned key set")
    expected_pairs = {(c.identity, m) for c in design.cells for m in ESTIMATORS}
    numeric_keys = list(numerical[["cell_id", "estimator"]].itertuples(index=False, name=None))
    if len(numeric_keys) != len(expected_pairs) or set(numeric_keys) != expected_pairs:
        raise ValueError("numerical report does not contain the full planned cell/method set")
    outcome_keys = list(outcomes[["cell_id", "estimator", "point_status", "interval_status"]]
                        .itertuples(index=False, name=None))
    if (len(outcome_keys) != len(set(outcome_keys))
            or set(outcomes[["cell_id", "estimator"]].itertuples(index=False, name=None)) != expected_pairs):
        raise ValueError("outcome report does not contain unique statuses for all planned cell/methods")
    cells = {c.identity: c for c in design.cells}
    for row in summary.itertuples():
        cell = cells[row.cell_id]
        if (row.n, row.theta_mean_gap, row.cv, row.truth) != (
                cell.n, cell.q, cell.cv, getattr(cell, row.parameter)):
            raise ValueError("summary design metadata mismatch")
        counts = (row.n_planned, row.n_point_valid, row.n_interval_valid)
        if (any(not np.isfinite(x) or x != int(x) for x in counts)
                or counts[0] != design.replications or not 0 <= counts[2] <= counts[1] <= counts[0]):
            raise ValueError("summary validity counts are inconsistent")
        if (row.practical_classification not in CLASSES
                or (row.parameter == "mu") != (row.practical_classification == CLASSES[-1])):
            raise ValueError("classification scope is inconsistent with zero truth")
    return tables, {name: manifest[name]["sha256"] for name in INPUTS}


def write_grid_reports(run_dir, design, *, output_dir=None):
    tables, inputs = load_tables(run_dir, design)
    summary, outcomes, numerical = (tables[n] for n in INPUTS)
    output = Path(output_dir) if output_dir is not None else Path(run_dir) / "figures"
    output.mkdir(parents=True, exist_ok=True)
    classifications = summary.sort_values(["n", "parameter", "estimator", "theta_mean_gap", "cv"])
    counts = []
    for (n, parameter, method), group in classifications.groupby(["n", "parameter", "estimator"], sort=True):
        for label in CLASSES:
            counts.append({"n": int(n), "parameter": parameter, "estimator": method,
                           "practical_classification": label, "n_cells": int((group.practical_classification == label).sum()),
                           "n_planned_cells": len(group)})
    failures = []
    for cell in design.cells:
        for method in ESTIMATORS:
            group = summary[(summary.cell_id == cell.identity) & (summary.estimator == method)]
            if any(group[column].nunique() != 1 for column in ("n_planned", "n_point_valid", "n_interval_valid")):
                raise ValueError("joint validity counts disagree between parameters")
            planned, point, interval = (int(group.iloc[0][c]) for c in ("n_planned", "n_point_valid", "n_interval_valid"))
            statuses = outcomes[(outcomes.cell_id == cell.identity) & (outcomes.estimator == method)]
            values = statuses["count"].to_numpy()
            if (not np.all(np.isfinite(values)) or np.any(values != values.astype(int))
                    or np.any(values <= 0) or int(values.sum()) != planned
                    or np.any(statuses.n_planned != planned)):
                raise ValueError("outcome counts do not equal all planned paths")
            valid_points = int(statuses.loc[statuses.point_status == "valid_interior", "count"].sum())
            valid_intervals = int(statuses.loc[statuses.interval_status == "valid", "count"].sum())
            if (valid_points, valid_intervals) != (point, interval):
                raise ValueError("outcome status accounting disagrees with validity summaries")
            composition = statuses.sort_values(["point_status", "interval_status"])[
                ["point_status", "interval_status", "count"]].to_dict(orient="records")
            diagnostics = numerical[(numerical.cell_id == cell.identity) & (numerical.estimator == method)].iloc[0]
            if diagnostics.n_planned != planned or diagnostics.n_interval_valid != interval:
                raise ValueError("numerical validity count mismatch")
            # Point-invalid records are already among interval-unavailable
            # records. Summing those two columns would double-count failures.
            # Statistical noncoverage is not a numerical failure.
            failures.append({"cell_id": cell.identity, "n": cell.n, "theta_mean_gap": cell.q,
                "cv": cell.cv, "estimator": method, "n_planned": planned,
                "n_point_valid": point, "n_interval_valid": interval,
                "n_point_invalid": planned-point, "n_interval_unavailable": planned-interval,
                "n_point_valid_without_interval": point-interval,
                "has_failure": interval < planned, "failure_fraction": (planned-interval)/planned,
                "status_counts_json": canonical(composition),
                **{name: diagnostics[name] for name in numerical.columns if name.endswith("_max")}})
    reports = {"classification_table.csv": classifications,
               "classification_counts.csv": pd.DataFrame(counts), "failed_cells.csv": pd.DataFrame(failures)}
    saved = {}
    for name, table in reports.items():
        atomic_bytes(output / name, table.to_csv(index=False, float_format="%.17g").encode(), immutable=True)
        # Read-back checks the complete exported values and labels. Scientific
        # classifications themselves come from the independently audited core.
        reread = pd.read_csv(output / name, float_precision="round_trip")
        pd.testing.assert_frame_equal(reread, table.reset_index(drop=True), check_dtype=False, check_exact=True)
        saved[name] = {"rows": len(table), "sha256": file_digest(output / name)}
    report = {"design_sha256": design.digest, "input_table_sha256": inputs, "reports": saved,
              "failed_cell_estimator_rows": sum(r["has_failure"] for r in failures),
              "failed_scientific_cells": len({r["cell_id"] for r in failures if r["has_failure"]}),
              "failure_log_unit": "all planned scientific cells times estimators, including zero-failure rows",
              "classification_scope": "pointwise; copied from audited parameter summaries without rounding or rethresholding"}
    atomic_json(output / "report_manifest.json", report, immutable=True)
    return report
