"""Reproducible complete-grid tables from immutable saved records.

Read one cell at a time; the 60-cell experiment need not fit all raw diagnostics
in memory. Every cell still undergoes complete Cartesian-key validation before
contributing to any table. References are matched by full scientific identity.
"""
from dataclasses import replace
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .comparisons import paired_estimator_comparisons, regular_reference_comparisons
from .provenance import file_digest, verify
from .storage import CellStore, atomic_bytes, atomic_json
from .summaries import ExperimentSummarizer, outcome_counts, spacing_summary, validate_records


TABLE_NAMES = ("parameter_summary.csv", "regular_comparisons.csv", "paired_comparisons.csv",
               "spacing_summary.csv", "outcome_counts.csv", "numerical_diagnostics.csv")


def read_cell(output, design, cell, identity):
    store = CellStore(Path(output) / "cells" / cell.identity, cell, design, identity, resume=True)
    if store.next_replication != design.replications:
        raise ValueError("all planned source paths are required before analysis")
    return validate_records(list(store.records()), replace(design, cells=(cell,)))


def csv_bytes(table):
    return table.to_csv(index=False, float_format="%.17g", lineterminator="\n").encode()


def numerical_diagnostics(data, design):
    data = validate_records(data, design)
    rows = []
    fields = ("hessian_condition_number", "hessian_relative_difference",
              "hessian_whitened_difference", "standardized_score")
    for (cid, name), group in data.groupby(["cell_id", "estimator"], sort=True):
        row = {"cell_id": cid, "estimator": name, "n_planned": len(group),
               "n_interval_valid": int(group.interval_valid.sum()),
               **{key: group[key].iloc[0] for key in ("n", "theta_mean_gap", "cv")}}
        for field in fields:
            values = pd.to_numeric(group[field], errors="raise").to_numpy(dtype=float)
            finite = values[np.isfinite(values)]
            row.update({f"{field}_n_finite": len(finite), f"{field}_n_nonfinite": len(values)-len(finite),
                        f"{field}_max": float(finite.max()) if len(finite) else np.nan,
                        f"{field}_median": float(np.median(finite)) if len(finite) else np.nan,
                        f"{field}_p95": float(np.quantile(finite, .95)) if len(finite) else np.nan})
        rows.append(row)
    return pd.DataFrame(rows)


def write_tables(output_dir, design):
    output = Path(output_dir)
    metadata = json.loads((output / "run_metadata.json").read_text())
    verify(metadata["provenance"], expected_design=design)
    expected_cells = {c.identity for c in design.cells}
    actual_cells = {p.name for p in (output / "cells").iterdir() if p.is_dir()}
    if actual_cells != expected_cells:
        raise ValueError("cell directories do not equal the full planned cell set")
    identity = metadata["provenance"]["identity_sha256"]
    collections = {name: [] for name in TABLE_NAMES}
    failures = []
    for cell in design.cells:
        cell_design = replace(design, cells=(cell,))
        data = read_cell(output, design, cell, identity)
        collections["parameter_summary.csv"].append(ExperimentSummarizer().summarize(data, cell_design))
        collections["paired_comparisons.csv"].append(paired_estimator_comparisons(data, cell_design))
        collections["spacing_summary.csv"].append(spacing_summary(data, cell_design))
        collections["outcome_counts.csv"].append(outcome_counts(data, cell_design))
        collections["numerical_diagnostics.csv"].append(numerical_diagnostics(data, cell_design))
        reference = replace(cell, cv=0)
        if reference.identity not in expected_cells:
            raise ValueError("missing planned regular reference")
        if cell.cv == 0:
            reference_data, comparison_design = data, cell_design
        else:
            reference_data = pd.concat([data, read_cell(output, design, reference, identity)], ignore_index=True)
            comparison_design = replace(design, cells=(reference, cell))
        comparisons = regular_reference_comparisons(reference_data, comparison_design)
        collections["regular_comparisons.csv"].append(comparisons[comparisons.cell_id == cell.identity])
        mask = ~data.point_valid | ~data.interval_valid | data.pipeline_error.fillna("").ne("")
        failures.append(data.loc[mask])
    result = {}
    for name, pieces in collections.items():
        table = pd.concat(pieces, ignore_index=True)
        atomic_bytes(output / name, csv_bytes(table), immutable=True)
        result[name] = {"rows": len(table), "sha256": file_digest(output / name)}
    failed = pd.concat(failures, ignore_index=True)
    name = "failed_fits.csv.gz"
    atomic_bytes(output / name, gzip.compress(csv_bytes(failed), compresslevel=6, mtime=0), immutable=True)
    result[name] = {"rows": len(failed), "sha256": file_digest(output / name)}
    verify(metadata["provenance"], expected_design=design)
    atomic_json(output / "table_manifest.json", result, immutable=True)
    return result
