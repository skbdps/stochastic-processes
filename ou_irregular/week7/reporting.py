"""Immutable Week 7 tables from complete saved raw evidence, never new paths.

The historical six-table/failed-fit writer is composed unchanged.  Floor
contrasts are additionally computed one complete CV triple at a time.  Root
lifecycle bytes are neither advanced nor repaired here; independent audit,
figures, review and sealing remain separate stages.
"""
import csv
from dataclasses import dataclass, replace
import gzip
from pathlib import Path

import pandas as pd

from ..experiments.analysis import (TABLE_NAMES, csv_bytes, read_cell,
                                     write_tables as standard_write_tables)
from ..experiments.design import ESTIMATORS, integer
from ..experiments.provenance import canonical, digest_bytes, file_digest
from ..experiments.storage import atomic_bytes, atomic_json, run_lock
from ..week4.summaries import PARAMETERS
from . import execution
from .floor_comparisons import COLUMNS, epsilon_comparisons
from .floor_source import ANCHOR_FLOOR_FRACTION, FLOOR_FRACTIONS
from .path_diagnostics import _ManifestGuard, _check_frozen


FILES = tuple(sorted(set(execution.FILES) | {
    "ou_irregular/experiments/analysis.py", "ou_irregular/experiments/comparisons.py",
    "ou_irregular/experiments/summaries.py", "ou_irregular/week7/floor_comparisons.py",
    "ou_irregular/week7/reporting.py", "docs/week7_reporting.md"}))
REPORTING_CONTRACT = "week7_saved_raw_tables_v1"
STANDARD_FILES = (*TABLE_NAMES, "failed_fits.csv.gz")
EPSILON_FILE = "epsilon_comparisons.csv"
EPSILON_MANIFEST = "epsilon_table_manifest.json"
_WRITABLE_STATES = {"raw_complete_pending_audit", "implementation_review_required"}


def _digest(value):
    return digest_bytes(canonical(value).encode())


def raw_input_binding(output, design, metadata):
    """Bind scientific raw inputs without root lifecycle bytes that change later."""
    output = Path(output)
    kind = metadata["execution_kind"]
    execution.require_design(kind, design)
    if canonical(metadata["source_policy"]) != canonical(execution.source_policy(kind, design)):
        raise ValueError("raw source policy differs from the fixed design")
    resolved, resolved_bytes = execution._json_file(output / "resolved_design.json")
    if canonical(resolved) != canonical(design.to_dict()):
        raise ValueError("raw resolved design changed")
    expected_cells = {cell.identity for cell in design.cells}
    directory = output / "cells"
    if {path.name for path in directory.iterdir()} != expected_cells:
        raise ValueError("raw cell directories do not equal the complete design")
    cells = {}
    identity = metadata["provenance"]["identity_sha256"]
    for cell in design.cells:
        root = directory / cell.identity
        state, payload = execution._json_file(root / "cell_manifest.json")
        if (state["run_identity"] != identity or canonical(state["cell"]) != canonical(cell.to_dict())
                or type(state["replications"]) is not int or state["replications"] != design.replications
                or type(state["completed_replications"]) is not int
                or state["completed_replications"] != design.replications or state["status"] != "complete"):
            raise ValueError("reporting requires a complete unchanged raw cell")
        partitions, start = [], 0
        for part in state["partitions"]:
            for key in ("start", "stop", "rows"):
                integer(part[key], "raw partition " + key, 0)
            if (part["filename"] != f"part-{part['start']:06d}-{part['stop']:06d}.jsonl.gz"
                    or part["start"] != start or not start < part["stop"] <= design.replications
                    or part["rows"] != 3 * (part["stop"] - start)):
                raise ValueError("raw partition map is noncontiguous or incomplete")
            path = root / part["filename"]
            if path.is_symlink() or not path.is_file() or file_digest(path) != part["sha256"]:
                raise ValueError("raw partition bytes changed")
            partitions.append(dict(part))
            start = part["stop"]
        if start != design.replications or {path.name for path in root.iterdir()} != {
                "cell_manifest.json", *(part["filename"] for part in partitions)}:
            raise ValueError("raw cell layout has missing or surplus evidence")
        cells[cell.identity] = {"manifest_sha256": digest_bytes(payload), "partitions": partitions}
    return {"run_identity": identity, "design_sha256": design.digest,
            "resolved_design_sha256": digest_bytes(resolved_bytes),
            "source_policy": metadata["source_policy"], "cells": cells}


def _csv_shape(path):
    opener = gzip.open if Path(path).suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as stream:
        reader = csv.reader(stream)
        header = next(reader, None)
        if not header or len(set(header)) != len(header) or any(not name for name in header):
            raise ValueError("table has a missing or duplicate column header")
        count = 0
        for row in reader:
            if len(row) != len(header):
                raise ValueError("table row width differs from its header")
            count += 1
    return header, count


def _verify_standard_manifest(output, design, *, required=False):
    path = output / "table_manifest.json"
    if not path.exists():
        if required:
            raise ValueError("standard table manifest is missing")
        return None
    manifest, _ = execution._json_file(path)
    if set(manifest) != set(STANDARD_FILES):
        raise ValueError("standard table manifest must contain exactly the seven declared artifacts")
    fixed_counts = {"parameter_summary.csv": 9 * len(design.cells),
                    "regular_comparisons.csv": 9 * len(design.cells),
                    "paired_comparisons.csv": 9 * len(design.cells),
                    "spacing_summary.csv": len(design.cells),
                    "numerical_diagnostics.csv": 3 * len(design.cells)}
    for name, entry in manifest.items():
        if not isinstance(entry, dict) or set(entry) != {"rows", "sha256"}:
            raise ValueError("standard artifact manifest entry changed")
        integer(entry["rows"], "standard table row count", 0)
        if file_digest(output / name) != entry["sha256"]:
            raise ValueError("declared standard table is missing or changed")
        _, rows = _csv_shape(output / name)
        if rows != entry["rows"] or (name in fixed_counts and rows != fixed_counts[name]):
            raise ValueError("declared standard table row count is inconsistent")
        if name == "failed_fits.csv.gz" and rows > 3 * len(design.cells) * design.replications:
            raise ValueError("failed-fit table exceeds the complete planned fit count")
        if name == "outcome_counts.csv" and not 3 * len(design.cells) <= rows <= 3 * len(design.cells) * design.replications:
            raise ValueError("outcome table count is outside the planned grid")
    return manifest


def _expected_epsilon_keys(design):
    return {(cell.identity, method, parameter) for cell in design.cells
            if cell.floor_fraction != ANCHOR_FLOOR_FRACTION
            for method in ESTIMATORS for parameter in PARAMETERS}


def _epsilon_csv(output, design):
    path = output / EPSILON_FILE
    header, count = _csv_shape(path)
    expected = _expected_epsilon_keys(design)
    if tuple(header) != COLUMNS or count != len(expected):
        raise ValueError("epsilon CSV must retain all 71 columns and planned comparison rows")
    with path.open(newline="", encoding="utf-8") as stream:
        keys = [(row["cell_id"], row["estimator"], row["parameter"]) for row in csv.DictReader(stream)]
    if len(set(keys)) != len(keys) or set(keys) != expected:
        raise ValueError("epsilon CSV does not equal the complete planned comparison key set")
    return {"filename": EPSILON_FILE, "rows": count, "columns": list(COLUMNS), "sha256": file_digest(path)}


def _epsilon_manifest(design, binding, table):
    return {"schema_version": 1, "project": "ou-week7-paired-floor-reporting",
            "reporting_contract": REPORTING_CONTRACT, "run_identity": binding["run_identity"],
            "design_sha256": design.digest, "raw_inputs": binding, "raw_input_sha256": _digest(binding),
            "table": table, "difference_direction": "epsilon_minus_registered_floor",
            "ratio_direction": "epsilon_over_registered_floor"}


def _verify_epsilon_manifest(output, kind, design, binding, *, required=False):
    csv_path, path = output / EPSILON_FILE, output / EPSILON_MANIFEST
    if kind != "floor":
        if path.exists() or csv_path.exists():
            raise ValueError("follow-up output must not contain paired-floor artifacts")
        return None
    if not path.exists():
        if required:
            raise ValueError("epsilon table manifest is missing")
        return None
    manifest, _ = execution._json_file(path)
    # Check declared existence/hash first.  A declared missing file must never
    # be silently regenerated and described as uninterrupted existing evidence.
    table = manifest.get("table")
    if not isinstance(table, dict) or table.get("filename") != EPSILON_FILE:
        raise ValueError("epsilon manifest has an unknown table descriptor")
    if not csv_path.is_file() or file_digest(csv_path) != table.get("sha256"):
        raise ValueError("declared epsilon table is missing or changed")
    expected = _epsilon_manifest(design, binding, _epsilon_csv(output, design))
    if canonical(manifest) != canonical(expected):
        raise ValueError("epsilon manifest source binding, schema or table identity changed")
    return manifest


@dataclass(frozen=True)
class ReportingPreflight:
    raw: execution.Preflight
    binding: dict
    standard: dict | None
    epsilon: dict | None


def preflight_reporting(kind, design, output_dir, *, extra_files=()):
    """Validate a proposed table write without output or lifecycle mutation."""
    raw = execution.preflight(kind, design, output_dir, extra_files=(*FILES, *tuple(extra_files)), resume=True)
    if (not raw.cached_completion or raw.metadata is None
            or raw.metadata["status"] not in _WRITABLE_STATES):
        raise ValueError("reporting requires complete zero-pipeline-error raw evidence")
    binding = raw_input_binding(raw.output, design, raw.metadata)
    standard = _verify_standard_manifest(raw.output, design)
    epsilon = _verify_epsilon_manifest(raw.output, kind, design, binding)
    if (raw.output / "run_metadata.json").read_bytes() != raw.expected_metadata_bytes:
        raise ValueError("root metadata changed during read-only reporting preflight")
    return ReportingPreflight(raw, binding, standard, epsilon)


def iter_epsilon_tables(output, design, identity):
    """Load exactly three cells per CV; keep only the 18-row result between groups."""
    for cv in sorted({cell.cv for cell in design.cells}):
        cells = tuple(cell for cell in design.cells if cell.cv == cv)
        if len(cells) != 3 or {cell.floor_fraction for cell in cells} != set(FLOOR_FRACTIONS):
            raise ValueError("bounded epsilon reporting requires complete three-floor CV groups")
        records = pd.concat([read_cell(output, design, cell, identity) for cell in cells], ignore_index=True)
        result = epsilon_comparisons(records, replace(design, cells=cells))
        del records
        yield result


def _verify_stage(prepared, kind, design, guard, stage):
    guard.check(stage + ":root")
    execution._root_layout(prepared.raw.output, design)
    _check_frozen(prepared.raw.metadata["provenance"], design, prepared.raw.provenance,
                  prepared.raw.require_frozen)
    current = raw_input_binding(prepared.raw.output, design, prepared.raw.metadata)
    if canonical(current) != canonical(prepared.binding):
        raise ValueError("raw input binding changed during reporting")
    guard.check(stage + ":after_raw_verification")


def write_tables(kind, design, output_dir, *, extra_files=()):
    """Write/adopt immutable derived tables; root lifecycle bytes never change."""
    extra_files = tuple(extra_files)
    prepared = preflight_reporting(kind, design, output_dir, extra_files=extra_files)
    output = prepared.raw.output
    with run_lock(output):
        prepared = preflight_reporting(kind, design, output, extra_files=extra_files)
        guard = _ManifestGuard(output / "run_metadata.json", prepared.raw.expected_metadata_bytes)
        try:
            _verify_stage(prepared, kind, design, guard, "before_calculations")
            epsilon_bytes = None
            if kind == "floor":
                pieces = []
                for table in iter_epsilon_tables(output, design, prepared.binding["run_identity"]):
                    guard.check("after_cv_comparison")
                    pieces.append(table)
                epsilon = pd.concat(pieces, ignore_index=True)
                expected = _expected_epsilon_keys(design)
                keys = list(epsilon[["cell_id", "estimator", "parameter"]].itertuples(index=False, name=None))
                if (tuple(epsilon.columns) != COLUMNS or len(keys) != len(expected)
                        or len(set(keys)) != len(keys) or set(keys) != expected):
                    raise ValueError("bounded epsilon calculations lost planned rows or columns")
                epsilon_bytes = csv_bytes(epsilon)
            _verify_stage(prepared, kind, design, guard, "before_standard_tables")
            standard_write_tables(output, design)
            _verify_stage(prepared, kind, design, guard, "after_standard_tables")
            standard = _verify_standard_manifest(output, design, required=True)
            epsilon_manifest = None
            if kind == "floor":
                guard.check("before_epsilon_table")
                atomic_bytes(output / EPSILON_FILE, epsilon_bytes, immutable=True)
                _verify_stage(prepared, kind, design, guard, "before_epsilon_manifest")
                epsilon_manifest = _epsilon_manifest(design, prepared.binding, _epsilon_csv(output, design))
                atomic_json(output / EPSILON_MANIFEST, epsilon_manifest, immutable=True)
                _verify_stage(prepared, kind, design, guard, "after_epsilon_manifest")
                _verify_epsilon_manifest(output, kind, design, prepared.binding, required=True)
            guard.check("reporting_complete")
            return {"standard": standard, "epsilon": epsilon_manifest}
        except BaseException:
            # No root-state write or stale finally block is appropriate here.
            # If a foreign writer changed the root, preserve that incident and
            # its expected/observed bytes before returning control to a runner.
            if not guard.drifted:
                guard.check("reporting_failure")
            raise
