"""Run the fixed E4 source-only validation after its implementation is frozen.

Exactly six cells x 256 independent 999-gap sequences; no OU draws or fits.
There are no seed, replication, tolerance, resume or overwrite overrides.
The validation identity below describes this source-only experiment; it is
not a fitted RunDesign. Saved per-sequence values and all checks are retained
for both passing and failing validations in a new immutable output directory.
"""
import argparse
import csv
from dataclasses import dataclass
import gzip
import hashlib
import io
import math
import os
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from ou_irregular.experiments.provenance import CORE_FILES, canonical, capture, file_digest, utc_now, verify
from ou_irregular.experiments.storage import atomic_bytes, atomic_json
from ou_irregular.week8.design import (SOURCE_ADDRESS, SOURCE_NAMESPACE, SOURCE_POLICY_VERSION,
                                      SOURCE_VALIDATION_REPLICATIONS, SOURCE_VALIDATION_SEED, full_design)
from ou_irregular.week8.source import (ADDRESS_FIELDS, CDF_GRID, FILES as SOURCE_FILES, GAP_FIELDS,
                                      E4Source, adjacent_product, cdf, gap_statistics, moments, phase_type)


FILES = tuple(sorted(set(SOURCE_FILES) | {
    "scripts/week8_validate_source.py", "ou_irregular/experiments/storage.py",
}))
REPLICATIONS = SOURCE_VALIDATION_REPLICATIONS
GAPS = 999
ADJACENT_PAIRS = 998
MULTIPLIER = 6.0
ABSOLUTE_TOLERANCE = 1e-10
CHECK_COLUMNS = (
    "cell_id", "mechanism", "ratio", "metric", "n_sequences", "terms_per_sequence",
    "observed", "expected", "standard_error", "absolute_error", "tolerance", "passed",
    "standard_error_method", "status",
)


def validation_cells():
    return full_design().cells[:6]


def validation_plan():
    """Fresh source-only protocol record, without fitted-run or estimator fields."""
    return {
        "schema_version": 1, "name": "week8_e4_source_only_validation",
        "scope": "source_only_validation", "source_policy_version": SOURCE_POLICY_VERSION,
        "seed": SOURCE_VALIDATION_SEED, "source_namespace": SOURCE_NAMESPACE,
        "source_address": SOURCE_ADDRESS, "schedule_child": 0,
        "replications_per_cell": REPLICATIONS, "gaps_per_sequence": GAPS,
        "adjacent_pairs_per_sequence": ADJACENT_PAIRS,
        "planned_sequences": 6 * REPLICATIONS, "planned_raw_gaps": 6 * REPLICATIONS * GAPS,
        "planned_adjacent_pairs": 6 * REPLICATIONS * ADJACENT_PAIRS,
        "cells": [cell.to_dict() for cell in validation_cells()],
        "cdf_t_over_m": [value for value, _ in CDF_GRID], "cdf_indicator": "gap<=t",
        "moment_cdf_gate": {"mcse_multiplier": MULTIPLIER, "absolute_tolerance": ABSOLUTE_TOLERANCE,
                            "mcse": "across_sequence_sample_sd_ddof1_div_sqrt256"},
        "phase_gate": {"observations_per_cell": REPLICATIONS, "phase_observation": "first_phase_only",
                       "standard_error": "theoretical_binomial", "multiplier": MULTIPLIER},
        "clipping_gate": "none_report_only", "new_ou_observations": 0, "new_fits": 0,
    }


@dataclass(frozen=True)
class _SourceValidationIdentity:
    """Adapter for provenance.capture, not an executable fitted-run design."""

    def to_dict(self):
        return validation_plan()

    @property
    def digest(self):
        return hashlib.sha256(canonical(self.to_dict()).encode()).hexdigest()


def expected_statistics(cell):
    targets = {}
    m = cell.q / cell.theta
    for kind, transformed in (("raw", False), ("transformed", True)):
        values = moments(cell, transformed=transformed)
        targets[f"{kind}_mean_gap"] = (values["mean"], GAPS)
        targets[f"{kind}_second_moment"] = (values["second_moment"], GAPS)
        targets[f"{kind}_adjacent_product"] = (adjacent_product(cell, transformed=transformed), ADJACENT_PAIRS)
        for ratio, label in CDF_GRID:
            targets[f"{kind}_cdf_t_over_m_{label}"] = (cdf(cell, ratio * m, transformed=transformed), GAPS)
    return targets


def comparison(values, expected, *, theoretical_probability=None):
    """The published sequence-level or initial-phase gate, on exactly N=256."""
    values = np.asarray(values, dtype=float)
    if values.shape != (REPLICATIONS,) or not np.isfinite(values).all() or not math.isfinite(expected):
        raise ValueError("comparison requires 256 finite independent sequence statistics")
    observed = float(values.mean())
    if theoretical_probability is None:
        se = float(values.std(ddof=1) / math.sqrt(REPLICATIONS))
        tolerance = MULTIPLIER * se + ABSOLUTE_TOLERANCE
        method = "across_sequence_sample_sd_ddof1_over_sqrt256"
    else:
        p = float(theoretical_probability)
        if not 0 <= p <= 1 or not np.isin(values, (0., 1.)).all():
            raise ValueError("phase check requires first-phase indicators and a probability")
        se = math.sqrt(p * (1 - p) / REPLICATIONS)
        tolerance = MULTIPLIER * se
        method = "theoretical_binomial_first_phase_N256"
    error = abs(observed - expected)
    return {"observed": observed, "expected": expected, "standard_error": se,
            "absolute_error": error, "tolerance": tolerance, "passed": error <= tolerance,
            "standard_error_method": method, "status": "evaluated"}


def evaluate_rows(rows):
    """All 108 fixed check rows, including unavailable checks after a source failure."""
    cells = validation_cells()
    expected_keys = {(cell.identity, rep) for cell in cells for rep in range(REPLICATIONS)}
    keys = [(row.get("cell_id"), row.get("rep")) for row in rows]
    if len(set(keys)) != len(keys) or set(keys) != expected_keys:
        raise ValueError("sequence rows must equal the complete fixed 1536-key plan")
    by_key = dict(zip(keys, rows))
    checks = []
    for cell in cells:
        group = [by_key[cell.identity, rep] for rep in range(REPLICATIONS)]
        targets = expected_statistics(cell)
        alpha = phase_type(cell).alpha
        for phase, label in enumerate(("calm", "burst")):
            targets[f"initial_phase_{label}"] = (float(alpha[phase]), 1)
        complete = all(row.get("status") == "evaluated" for row in group)
        for metric, (expected, terms) in targets.items():
            row = {"cell_id": cell.identity, "mechanism": cell.mechanism, "ratio": cell.ratio,
                   "metric": metric, "n_sequences": REPLICATIONS, "terms_per_sequence": terms}
            if not complete:
                row.update(observed=None, expected=expected, standard_error=None, absolute_error=None,
                           tolerance=None, passed=False, standard_error_method="unavailable",
                           status="source_failure_or_unfinished_sequence")
            else:
                if metric.startswith("initial_phase_"):
                    phase = 0 if metric.endswith("calm") else 1
                    values = [float(item["source_initial_phase"] == phase) for item in group]
                    result = comparison(values, expected, theoretical_probability=expected)
                else:
                    result = comparison([item[metric] for item in group], expected)
                row.update(result)
            checks.append(row)
    return checks


def _csv_bytes(rows, columns):
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=columns, lineterminator="\n", extrasaction="raise")
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue().encode()


def _output_path(output):
    supplied = Path(output).absolute()
    if any(path.is_symlink() for path in (supplied, *supplied.parents)):
        raise ValueError("validation output must not use symbolic links")
    markers = {"run_metadata.json", "cell_manifest.json", "resolved_design.json", "diagnostic_design.json",
               "table_manifest.json", "figure_manifest.json", "validation.json", "validation_plan.json"}
    if any((ancestor / name).is_file() for ancestor in supplied.parents for name in markers):
        raise ValueError("validation output must not be nested in existing scientific evidence")
    if supplied.exists():
        raise FileExistsError("source validation requires a fresh output directory; no resume or replacement")
    return supplied


def run_validation(output):
    """Execute the exact source-only protocol once, preserving all outcomes."""
    output = _output_path(output)
    identity = _SourceValidationIdentity()
    provenance = capture(identity, extra_files=FILES, require_frozen=True)
    if set(provenance["source_files"]) != set(CORE_FILES) | set(FILES):
        raise ValueError("source-validation dependency membership changed")
    started = utc_now()
    output.mkdir(parents=True, exist_ok=False)
    atomic_json(output / "validation_plan.json", identity.to_dict(), immutable=True)
    atomic_json(output / "provenance.json", provenance, immutable=True)
    atomic_json(output / "started.json", {"started_utc": started, "status": "started",
                                         "new_ou_observations": 0, "new_fits": 0}, immutable=True)
    source = E4Source(SOURCE_VALIDATION_SEED)
    columns = tuple(sorted(set(validation_cells()[0].to_dict()) | GAP_FIELDS | ADDRESS_FIELDS | {
        "rep", "status", "error", "n_gaps", "adjacent_pairs", "times_sha256",
    }))
    rows = []
    for cell in validation_cells():
        for rep in range(REPLICATIONS):
            row = dict.fromkeys(columns)
            row.update(cell.to_dict())
            row.update(source.address(cell.to_dict(), rep))
            row.update(rep=rep, status="not_evaluated", error="", n_gaps=GAPS, adjacent_pairs=ADJACENT_PAIRS)
            rows.append(row)
    fatal = None
    try:
        for index, row in enumerate(rows):
            cell = validation_cells()[index // REPLICATIONS]
            try:
                gaps = source.gaps(cell.to_dict(), row["rep"])
                if len(gaps.raw_gaps) != GAPS or len(gaps.raw_gaps[:-1]) != ADJACENT_PAIRS:
                    raise ValueError("generated sequence dimensions differ from the protocol")
                row.update(gap_statistics(cell, gaps))
                from ou_irregular.week4.simulation import array_digest
                row.update(times_sha256=array_digest(gaps.times), status="evaluated")
            except Exception as exc:
                row.update(status="source_failure", error=f"{type(exc).__name__}: {exc}")
            if (index + 1) % REPLICATIONS == 0:
                group = rows[index + 1 - REPLICATIONS:index + 1]
                payload = ("\n".join(canonical(item) for item in group) + "\n").encode()
                atomic_bytes(output / f"sequences-{cell.identity}.jsonl.gz",
                             gzip.compress(payload, compresslevel=6, mtime=0), immutable=True)
    except BaseException as exc:
        fatal = f"{type(exc).__name__}: {exc}"
    checks = evaluate_rows(rows)
    try:
        verify(provenance, expected_design=identity)
    except Exception as exc:
        fatal = (fatal + "; " if fatal else "") + f"provenance: {type(exc).__name__}: {exc}"
    atomic_bytes(output / "sequence_summary.csv", _csv_bytes(rows, columns), immutable=True)
    atomic_bytes(output / "checks.csv", _csv_bytes(checks, CHECK_COLUMNS), immutable=True)
    atomic_json(output / "checks.json", checks, immutable=True)
    failures = [{"cell_id": row["cell_id"], "rep": row["rep"], "status": row["status"], "error": row["error"]}
                for row in rows if row["status"] != "evaluated"]
    failed_checks = [{"cell_id": row["cell_id"], "metric": row["metric"], "status": row["status"]}
                     for row in checks if not row["passed"]]
    report = {
        "schema_version": 1, "validation_contract": "week8_source_only_validation_v1",
        "status": "passed" if not failures and not failed_checks and fatal is None else "failed",
        "started_utc": started, "finished_utc": utc_now(), "plan_sha256": identity.digest,
        "provenance_identity_sha256": provenance["identity_sha256"], "code_commit": provenance["code_commit"],
        "planned_sequences": 6 * REPLICATIONS, "retained_sequence_rows": len(rows),
        "evaluated_sequences": sum(row["status"] == "evaluated" for row in rows), "check_rows": len(checks),
        "source_failures": failures, "failed_checks": failed_checks, "fatal_error": fatal,
        "clipping_test": "reported_only_no_empirical_normal_gate",
        "new_ou_observations": 0, "new_fits": 0,
        "artifacts": {path.name: file_digest(path) for path in sorted(output.iterdir()) if path.is_file()},
    }
    atomic_json(output / "validation.json", report, immutable=True)
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    report = run_validation(args.output)
    print(canonical({key: report[key] for key in ("status", "planned_sequences", "evaluated_sequences",
                                                 "check_rows", "new_ou_observations", "new_fits")}), flush=True)
    return 0 if report["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
