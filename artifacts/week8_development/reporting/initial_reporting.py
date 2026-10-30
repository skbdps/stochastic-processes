"""Immutable E4 tables from complete saved records, without new sources or fits.

The public boundary comprises FILES, TABLE_KEYS, TABLE_COLUMNS, SOURCE_METRICS,
RAW_FIELDS, FAILURE_COLUMNS, raw_input_binding, independent_comparisons,
source_summary, preflight_reporting, write_tables and verify_tables. Independent
contrasts use the explicit regular/matched_iid maps and each cell's own eligible
summary. The old CV-replacement reference calculator/writer is never called.

Pure source_summary permits canonical internal cell views; filesystem stages
and independent_comparisons require the complete fixed full or validation root.
verify_tables checks saved schemas, planned keys/counts and raw bindings, not
independent statistical arithmetic. That remains a separate audit stage.
"""
import csv
from dataclasses import dataclass, replace
import gzip
import io
from itertools import combinations
import math
from pathlib import Path

import numpy as np
import pandas as pd

from ..experiments.analysis import csv_bytes, numerical_diagnostics, read_cell
from ..experiments.comparisons import _interval, _ratio, _relative, paired_estimator_comparisons
from ..experiments.design import ESTIMATORS, integer
from ..experiments.provenance import canonical, digest_bytes, file_digest
from ..experiments.storage import atomic_bytes, atomic_json, run_lock
from ..experiments.summaries import (ExperimentSummarizer, mean_mcse, outcome_counts,
                                      spacing_summary, validate_records)
from ..week4.summaries import PARAMETERS
from . import execution
from .design import E4Cell, SOURCE_NAMESPACE, design_for, full_design, reference_map, source_policy, validate_design
from .source import CDF_GRID, GAP_FIELDS, SOURCE_FIELDS, _FIT_FIELDS


FILES = tuple(sorted(set(execution.FILES) | {
    "ou_irregular/experiments/analysis.py", "ou_irregular/experiments/comparisons.py",
    "ou_irregular/experiments/summaries.py", "ou_irregular/week8/reporting.py", "docs/week8_reporting.md",
}))
REPORTING_CONTRACT = "week8_saved_raw_tables_v1"
_CELL_FIELDS = frozenset(full_design().cells[0].to_dict())
RAW_FIELDS = _CELL_FIELDS | SOURCE_FIELDS | _FIT_FIELDS
FAILURE_COLUMNS = tuple(sorted(RAW_FIELDS))
TABLE_KEYS = {
    "parameter_summary.csv": ("cell_id", "estimator", "parameter"),
    "paired_comparisons.csv": ("cell_id", "estimator_a", "estimator_b", "parameter"),
    "regular_comparisons.csv": ("cell_id", "estimator", "parameter"),
    "matched_iid_comparisons.csv": ("cell_id", "estimator", "parameter"),
    "spacing_summary.csv": ("cell_id",), "source_summary.csv": ("cell_id",),
    "outcome_counts.csv": ("cell_id", "estimator", "point_status", "interval_status"),
    "numerical_diagnostics.csv": ("cell_id", "estimator"),
}
_INTERVAL_SUFFIXES = ("", "_mcse", "_mc_lower", "_mc_upper")
_RATIO_COLUMNS = {"rmse_ratio", "rmse_ratio_mcse", "rmse_ratio_log_mcse", "rmse_ratio_mc_lower",
                  "rmse_ratio_mc_upper", "rmse_ratio_uncertainty"}
_PARAMETER_COLUMNS = {
    "estimator", "parameter", "truth", "n_planned", "n_point_valid", "n_interval_valid", "n_covered",
    "point_availability", "interval_availability", "bias", "rmse", "mcse_bias", "relative_bias",
    "relative_mcse_bias", "relative_bias_mc_lower", "relative_bias_mc_upper", "coverage_conditional",
    "coverage_operational", "conditional_mcse", "operational_mcse", "conditional_wilson_lower",
    "conditional_wilson_upper", "operational_wilson_lower", "operational_wilson_upper",
    "negative_lower_count", "practical_classification", "point_metric_scope", "primary_coverage_scope",
    "mcse_rmse", "relative_rmse", "relative_mcse_rmse", "empirical_sd_point_valid", "bias_mc_lower",
    "bias_mc_upper", "empirical_sd_interval_valid", "mean_reported_se", "rms_reported_se",
    "nominal_dimensionless_span",
    *(f"{kind}_availability_{suffix}" for kind in ("point", "interval")
      for suffix in ("mcse", "wilson_lower", "wilson_upper")),
}
_PAIRED_COLUMNS = {
    "estimator_a", "estimator_b", "parameter", "truth", "comparison_type", "direction", "n_planned",
    "n_point_valid_a", "n_point_valid_b", "n_common_point_valid", "n_common_interval_valid",
    "point_subset", "coverage_subset", "bias_a_common", "bias_b_common", "rmse_a_common", "rmse_b_common",
    "squared_error_covariance_common", "n_covered_a", "n_covered_b", "n_both_cover", "n_a_only_cover",
    "n_b_only_cover", "n_neither_cover", "operational_coverage_a", "operational_coverage_b",
    *(name + suffix for name in ("bias_difference", "mse_difference", "relative_bias_difference",
                                "rmse_difference", "operational_coverage_difference")
      for suffix in _INTERVAL_SUFFIXES), *_RATIO_COLUMNS,
}
COMPONENT_FIELDS = (
    "n_planned", "n_point_valid", "n_interval_valid", "n_covered", "point_availability", "interval_availability",
    "bias", "rmse", "mcse_bias", "mcse_rmse", "relative_bias", "relative_mcse_bias",
    "relative_rmse", "relative_mcse_rmse", "coverage_operational", "coverage_conditional",
    "operational_mcse", "conditional_mcse",
)
_COMPARISON_COLUMNS = {
    "estimator", "parameter", "truth", "comparison_kind", "comparison_type", "point_subset", "coverage_subset",
    "matching", "difference_direction", "ratio_direction",
    *(f"{side}_{field}" for side in ("cell", "reference") for field in COMPONENT_FIELDS),
    *("reference_" + field for field in _CELL_FIELDS),
    *(name + suffix for name in ("excess_bias", "relative_excess_bias", "operational_coverage_difference",
                                "conditional_coverage_difference") for suffix in _INTERVAL_SUFFIXES),
    *_RATIO_COLUMNS,
}
SOURCE_METRICS = (
    *(f"{kind}_{suffix}" for kind in ("raw", "transformed")
      for suffix in ("mean_gap", "cv", "second_moment", "adjacent_product", "min_gap", "span",
                     *(f"cdf_t_over_m_{label}" for _, label in CDF_GRID))),
    "realized_mean_gap", "realized_cv", "realized_span", "min_gap",
    "clip_count", "clip_fraction", "raw_below_floor_count", "transformed_below_floor_count", "transformed_at_floor_count",
    "source_switch_count", "source_max_switches_per_gap", "source_phase_reset_count",
    "source_calm_start_count", "source_burst_start_count", "source_calm_end_count", "source_burst_end_count",
)
SOURCE_COUNT_METRICS = frozenset(name for name in SOURCE_METRICS if name.endswith("_count")) | {"source_max_switches_per_gap"}
SOURCE_EXPECTATIONS = tuple(sorted(name for name in GAP_FIELDS if name.startswith("expected_")))
_SOURCE_CONSTANTS = ("family", "target_mean_gap", "requested_cv", "floor_applied", "floor_value", *SOURCE_EXPECTATIONS)
_SOURCE_PHASES = (("calm", 0), ("burst", 1), ("regular", -1))
_SOURCE_COLUMNS = {
    "n_planned", "n_unique_sources", "n_source_valid", "n_fit_records", "source_metric_mcse",
    "raw_metric_scope", "transformed_metric_scope", "realized_metric_scope", *_SOURCE_CONSTANTS,
    *(f"{name}_{suffix}" for name in SOURCE_METRICS for suffix in ("n_available", "mean", "mcse", "min", "max")),
    *(name + "_total" for name in SOURCE_COUNT_METRICS),
    *(f"initial_phase_{name}_{suffix}" for name, _ in _SOURCE_PHASES for suffix in ("count", "fraction", "mcse")),
    *("expected_initial_phase_" + name for name, _ in _SOURCE_PHASES),
}
_SPACING_METRICS = ("realized_mean_gap", "realized_cv", "realized_span", "clip_fraction", "min_gap")
_NUMERICAL_METRICS = ("hessian_condition_number", "hessian_relative_difference", "hessian_whitened_difference", "standardized_score")
TABLE_COLUMNS = {name: tuple(sorted(_CELL_FIELDS | columns)) for name, columns in {
    "parameter_summary.csv": _PARAMETER_COLUMNS, "paired_comparisons.csv": _PAIRED_COLUMNS,
    "regular_comparisons.csv": _COMPARISON_COLUMNS, "matched_iid_comparisons.csv": _COMPARISON_COLUMNS,
    "spacing_summary.csv": {"n_planned", *(f"{field}_{suffix}" for field in _SPACING_METRICS
                                            for suffix in ("n_available", "mean", "mcse", "min", "max"))},
    "source_summary.csv": _SOURCE_COLUMNS,
    "outcome_counts.csv": {"estimator", "point_status", "interval_status", "count", "n_planned", "fraction", "wilson_lower", "wilson_upper"},
    "numerical_diagnostics.csv": {"estimator", "n_planned", "n_interval_valid",
                                   *(f"{field}_{suffix}" for field in _NUMERICAL_METRICS
                                     for suffix in ("n_finite", "n_nonfinite", "max", "median", "p95"))},
}.items()}
TABLE_COLUMNS["failed_fits.csv.gz"] = FAILURE_COLUMNS
_READABLE_STATES = {"raw_complete_pending_audit", "audited_pending_visual_review", "completed"}


def _digest(value):
    return digest_bytes(canonical(value).encode())


def _view_design(design):
    root = design_for(validation=design.scope == "validation")
    validate_design(replace(design, cells=root.cells))
    planned = {cell.identity: cell.to_dict() for cell in root.cells}
    if (not design.cells or len({cell.identity for cell in design.cells}) != len(design.cells)
            or any(type(cell) is not E4Cell or cell.identity not in planned
                   or canonical(cell.to_dict()) != canonical(planned[cell.identity]) for cell in design.cells)):
        raise ValueError("reporting cell view must retain fixed E4 scientific identities")


def _enrich(table, design):
    data = pd.DataFrame(table).copy()
    mapping = {cell.identity: cell.to_dict() for cell in design.cells}
    if "cell_id" not in data or not set(data.cell_id) <= set(mapping):
        raise ValueError("table contains an unplanned cell")
    for field in _CELL_FIELDS - {"cell_id"}:
        expected = data.cell_id.map(lambda cid: mapping[cid][field])
        if field in data and not data[field].eq(expected).all():
            raise ValueError("table scientific metadata differs: " + field)
        data[field] = expected
    return data


def _table(name, data, design):
    data = _enrich(data, design)
    if not data.columns.is_unique or set(data) != set(TABLE_COLUMNS[name]):
        raise ValueError(f"unexpected complete table schema: {name}; difference={sorted(set(data) ^ set(TABLE_COLUMNS[name]))}")
    for field in data.select_dtypes(include=[np.number]):
        if np.isinf(data[field].to_numpy(dtype=float)).any():
            raise ValueError("table arithmetic overflow: " + field)
    return data.loc[:, TABLE_COLUMNS[name]].sort_values(list(TABLE_KEYS[name]), kind="stable").reset_index(drop=True)


def _summary_input(summary, design):
    validate_design(design)
    data = pd.DataFrame(summary).copy()
    if tuple(data.columns) != TABLE_COLUMNS["parameter_summary.csv"]:
        raise ValueError("independent contrasts require the exact parameter-summary schema")
    data = _enrich(data, design)
    keys = list(data[list(TABLE_KEYS["parameter_summary.csv"])].itertuples(index=False, name=None))
    expected = {(cell.identity, estimator, parameter) for cell in design.cells for estimator in ESTIMATORS for parameter in PARAMETERS}
    if len(set(keys)) != len(keys) or set(keys) != expected:
        raise ValueError("independent contrasts require the complete planned parameter-summary keys")
    cells = {cell.identity: cell for cell in design.cells}
    for row in data.to_dict("records"):
        for name in ("n_planned", "n_point_valid", "n_interval_valid", "n_covered", "negative_lower_count"):
            integer(row[name], "summary " + name, 0)
        r, p, i, c = (row[name] for name in ("n_planned", "n_point_valid", "n_interval_valid", "n_covered"))
        if r != design.replications or not 0 <= c <= i <= p <= r or row["negative_lower_count"] > i:
            raise ValueError("summary eligibility/count semantics differ from the fixed plan")
        if row["truth"] != cells[row["cell_id"]].to_dict()["true_" + row["parameter"]]:
            raise ValueError("summary truth differs from the fixed scientific cell")
        for name, value in {"point_availability": p/r, "interval_availability": i/r,
                            "coverage_operational": c/r, "coverage_conditional": c/i if i else np.nan,
                            "operational_mcse": math.sqrt((c/r)*(1-c/r)/r),
                            "conditional_mcse": math.sqrt((c/i)*(1-c/i)/i) if i else np.nan}.items():
            if isinstance(row[name], (bool, np.bool_)) or not np.isclose(row[name], value, rtol=2e-13, atol=1e-15, equal_nan=True):
                raise ValueError("summary counts and coverage components disagree: " + name)
        for name in COMPONENT_FIELDS:
            value = row[name]
            if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, float, np.integer, np.floating)) or np.isinf(value):
                raise ValueError("summary component must be numeric or genuinely undefined: " + name)
        if ((p and (not np.isfinite([row["bias"], row["rmse"]]).all() or row["rmse"] < 0))
                or (not p and not np.isnan([row["bias"], row["rmse"], row["mcse_bias"], row["mcse_rmse"]]).all())):
            raise ValueError("summary point components disagree with point eligibility")
        for field, available in (("mcse_bias", p > 1), ("mcse_rmse", p > 1 and row["rmse"] > 0)):
            value = row[field]
            if (available and (not np.isfinite(value) or value < 0)) or (not available and not np.isnan(value)):
                raise ValueError("summary point uncertainty has invalid eligibility: " + field)
        if row["truth"] == 0 and any(not np.isnan(row[name]) for name in COMPONENT_FIELDS if name.startswith("relative_")):
            raise ValueError("relative zero-truth summaries must remain undefined")
    return data


def independent_comparisons(summary, design, kind):
    """Cell-minus-reference contrasts on each independent cell's own eligible fits."""
    if kind not in {"regular", "matched_iid"}:
        raise ValueError("comparison kind must be regular or matched_iid")
    data = _summary_input(summary, design)
    index = data.set_index(["cell_id", "estimator", "parameter"])
    cells = {cell.identity: cell for cell in design.cells}
    rows = []
    with np.errstate(over="raise", invalid="ignore", divide="ignore"):
        for cid, reference_id in reference_map(design)[kind].items():
            for estimator in ESTIMATORS:
                for parameter in PARAMETERS:
                    a, b = index.loc[(cid, estimator, parameter)], index.loc[(reference_id, estimator, parameter)]
                    row = {**cells[cid].to_dict(), **{"reference_" + key: value for key, value in cells[reference_id].to_dict().items()},
                           "estimator": estimator, "parameter": parameter, "truth": a.truth,
                           "comparison_kind": kind, "comparison_type": "independent_cells",
                           "point_subset": "each_cell_own_valid_point_fits",
                           "coverage_subset": "all_planned_paths_unavailable_counts_as_not_covering",
                           "difference_direction": "cell_minus_reference", "ratio_direction": "cell_over_reference",
                           "matching": "same_full_event_gap_marginal_and_floor_transform_independent_paths" if kind == "matched_iid"
                           else "same_n_truth_nominal_mean_gap_and_floor_settings_not_realized_span"}
                    for field in COMPONENT_FIELDS:
                        row["cell_" + field], row["reference_" + field] = a[field], b[field]
                    for name, field, se in (("excess_bias", "bias", "mcse_bias"),
                                            ("operational_coverage_difference", "coverage_operational", "operational_mcse"),
                                            ("conditional_coverage_difference", "coverage_conditional", "conditional_mcse")):
                        _interval(row, name, a[field] - b[field], np.hypot(a[se], b[se]))
                    _relative(row, "excess_bias", a.truth)
                    log_se = np.hypot(a.mcse_rmse/a.rmse, b.mcse_rmse/b.rmse) if a.rmse > 0 and b.rmse > 0 else np.nan
                    _ratio(row, a.rmse, b.rmse, log_se)
                    rows.append(row)
    return _table("regular_comparisons.csv" if kind == "regular" else "matched_iid_comparisons.csv", rows, design)


def _source_rows(records, design):
    _view_design(design)
    data = pd.DataFrame(records)
    if not data.columns.is_unique or set(data) != RAW_FIELDS:
        raise ValueError("reporting requires the complete exact raw record schema")
    data = validate_records(data, design)
    if (not data.source_valid.map(lambda value: isinstance(value, (bool, np.bool_)) and bool(value)).all()
            or not data.source_error.map(lambda value: type(value) is str and value == "").all()
            or not data.pipeline_error.map(lambda value: type(value) is str and value == "").all()):
        raise ValueError("source/pipeline failures require review before scientific reporting")
    columns = tuple(sorted(SOURCE_FIELDS))
    for _, group in data.groupby(["cell_id", "rep"], sort=False):
        values = [canonical(item) for item in group.loc[:, columns].to_dict("records")]
        if len(values) != len(ESTIMATORS) or len(set(values)) != 1:
            raise ValueError("all SOURCE_FIELDS must agree across the complete estimator triplet")
    sources = data[data.estimator == "exact"].copy()
    cells = {cell.identity: cell for cell in design.cells}
    hashes = {name for name in SOURCE_FIELDS if name.endswith("_sha256")}
    for row in sources.to_dict("records"):
        cell, rep = cells[row["cell_id"]], row["rep"]
        m, k = cell.q/cell.theta, cell.n-1
        f = cell.floor_fraction*m if cell.cv >= cell.floor_cv_threshold else 0.
        parent = [SOURCE_NAMESPACE, *(int(cell.identity[i:i+8], 16) for i in range(0, 64, 8)), rep]
        expected = {"source_seed": design.seed, "source_spawn_key": canonical(parent),
                    "gap_spawn_key": canonical(parent+[0]), "path_spawn_key": canonical(parent+[1]),
                    "family": {"mmpp":"event_stationary_mmpp", "iid":"iid_full_palm_marginal", "regular":"equidistant"}[cell.mechanism],
                    "target_mean_gap": m, "requested_cv": cell.cv, "floor_value": f,
                    "floor_applied": cell.cv >= cell.floor_cv_threshold,
                    "source_phase_reset_count": k if cell.mechanism == "iid" else 1 if cell.mechanism == "mmpp" else 0}
        if any(canonical(row[name]) != canonical(value) for name, value in expected.items()):
            raise ValueError("saved source identity or policy diagnostic changed")
        if not isinstance(row["floor_applied"], (bool, np.bool_)):
            raise ValueError("source floor flag must be Boolean")
        for name in hashes:
            value = row[name]
            if type(value) is not str or len(value) != 64 or any(c not in "0123456789abcdef" for c in value):
                raise ValueError("invalid saved source hash: " + name)
        for name in (*SOURCE_COUNT_METRICS, "source_seed", "source_initial_phase"):
            integer(row[name], "source " + name, -1 if name == "source_initial_phase" else 0)
        if row["source_initial_phase"] not in ((-1,) if cell.mechanism == "regular" else (0, 1)):
            raise ValueError("saved initial phase is invalid for this mechanism")
        total = 0 if cell.mechanism == "regular" else k
        if any(row[f"source_calm_{end}_count"]+row[f"source_burst_{end}_count"] != total for end in ("start", "end")):
            raise ValueError("saved phase counts do not equal the sequence dimensions")
        if (not 0 <= row["clip_count"] <= row["transformed_at_floor_count"] <= k
                or row["clip_count"] != row["raw_below_floor_count"] or row["transformed_below_floor_count"] != 0
                or row["clip_fraction"] != row["clip_count"]/k):
            raise ValueError("saved clipping counts/fraction are inconsistent")
        for name in (*SOURCE_METRICS, *SOURCE_EXPECTATIONS):
            value = row[name]
            if (isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, float, np.integer, np.floating))
                    or not math.isfinite(value) or value < 0):
                raise ValueError("source summary requires finite nonnegative diagnostics: " + name)
            if ("_cdf_" in name or "fraction" in name) and value > 1:
                raise ValueError("source probability/fraction exceeds one: " + name)
    for _, group in sources.groupby("cell_id", sort=False):
        for name in _SOURCE_CONSTANTS:
            if len({canonical(value) for value in group[name]}) != 1:
                raise ValueError("saved source expectation/policy must be constant within a cell: " + name)
    return data, sources


def _summarize_sources(sources, design):
    rows = []
    for cell in design.cells:
        group = sources[sources.cell_id == cell.identity].sort_values("rep")
        row = {**cell.to_dict(), **{name: group[name].iloc[0] for name in _SOURCE_CONSTANTS},
               "n_planned": design.replications, "n_unique_sources": len(group), "n_source_valid": len(group),
               "n_fit_records": len(group)*len(ESTIMATORS),
               "source_metric_mcse": "across_independent_paths_sample_sd_ddof1_over_sqrtR",
               "raw_metric_scope": "generated_raw_gap_vector_before_cumsum",
               "transformed_metric_scope": "posthoc_floored_gap_vector_before_cumsum",
               "realized_metric_scope": "differences_of_stored_cumulative_times"}
        for name in SOURCE_METRICS:
            values = group[name].to_numpy(dtype=float)
            row.update({name+"_n_available": len(values), name+"_mean": float(values.mean()),
                        name+"_mcse": mean_mcse(values), name+"_min": float(values.min()), name+"_max": float(values.max())})
            if name in SOURCE_COUNT_METRICS:
                row[name+"_total"] = sum(int(value) for value in group[name])
        for name, value in _SOURCE_PHASES:
            indicators = (group.source_initial_phase.to_numpy() == value).astype(float)
            row.update({f"initial_phase_{name}_count": int(indicators.sum()),
                        f"initial_phase_{name}_fraction": float(indicators.mean()),
                        f"initial_phase_{name}_mcse": mean_mcse(indicators)})
            row["expected_initial_phase_"+name] = (float(value == -1) if cell.mechanism == "regular" else
                0. if value == -1 else cell.to_dict()[f"{name}_event_probability"])
        rows.append(row)
    return _table("source_summary.csv", rows, design)


def source_summary(records, design):
    """One copy per path, with complete triplet agreement and across-path MCSEs."""
    _, sources = _source_rows(records, design)
    return _summarize_sources(sources, design)


def raw_input_binding(output, design, metadata=None):
    """Bind complete immutable raw inputs, excluding later root lifecycle bytes."""
    validate_design(design)
    output = execution._ordinary_path(output)
    execution._root_layout(output, design)
    if metadata is None:
        metadata, _ = execution._json_file(output / "run_metadata.json")
    if canonical(metadata.get("source_policy")) != canonical(source_policy(design)):
        raise ValueError("raw source policy differs from the fixed E4 design")
    resolved, payload = execution._json_file(output / "resolved_design.json")
    if canonical(resolved) != canonical(design.to_dict()):
        raise ValueError("raw resolved design changed")
    identity = metadata["provenance"]["identity_sha256"]
    directory, cells = output / "cells", {}
    if {path.name for path in directory.iterdir()} != {cell.identity for cell in design.cells}:
        raise ValueError("raw cell directories do not equal the complete E4 design")
    for cell in design.cells:
        root = directory / cell.identity
        state, saved = execution._json_file(root / "cell_manifest.json")
        if (set(state) != {"schema_version", "run_identity", "cell", "replications", "status",
                           "completed_replications", "partitions", "recovered_orphans"}
                or type(state["schema_version"]) is not int or state["schema_version"] != 1
                or state["run_identity"] != identity or canonical(state["cell"]) != canonical(cell.to_dict())
                or type(state["replications"]) is not int or state["replications"] != design.replications
                or type(state["completed_replications"]) is not int or state["completed_replications"] != design.replications
                or state["status"] != "complete" or not isinstance(state["partitions"], list)
                or not isinstance(state["recovered_orphans"], list)):
            raise ValueError("reporting requires a complete unchanged raw cell")
        partitions, start = [], 0
        for part in state["partitions"]:
            if not isinstance(part, dict) or set(part) != {"filename", "start", "stop", "rows", "sha256"}:
                raise ValueError("raw partition map has unexpected fields")
            for key in ("start", "stop", "rows"):
                integer(part[key], "raw partition " + key, 0)
            if (part["filename"] != f"part-{part['start']:06d}-{part['stop']:06d}.jsonl.gz"
                    or part["start"] != start or not start < part["stop"] <= design.replications
                    or part["rows"] != len(ESTIMATORS)*(part["stop"]-start)):
                raise ValueError("raw partition map is noncontiguous or incomplete")
            path = root / part["filename"]
            if path.is_symlink() or not path.is_file() or file_digest(path) != part["sha256"]:
                raise ValueError("raw partition bytes changed")
            partitions.append(dict(part))
            start = part["stop"]
        if start != design.replications or {path.name for path in root.iterdir()} != {
                "cell_manifest.json", *(part["filename"] for part in partitions)}:
            raise ValueError("raw cell layout has missing or surplus evidence")
        cells[cell.identity] = {"manifest_sha256": digest_bytes(saved), "partitions": partitions}
    return {"run_identity": identity, "design_sha256": design.digest,
            "resolved_design_sha256": digest_bytes(payload), "source_policy": source_policy(design), "cells": cells}


def _binding_record(design, binding):
    return {"schema_version": 1, "project": "ou-week8-e4-reporting", "reporting_contract": REPORTING_CONTRACT,
            "run_identity": binding["run_identity"], "design_sha256": design.digest,
            "raw_inputs": binding, "raw_input_sha256": _digest(binding)}


def _raw_cells(output, design, identity):
    for cell in sorted(design.cells, key=lambda item: item.identity):
        view = replace(design, cells=(cell,))
        data, sources = _source_rows(read_cell(output, design, cell, identity), view)
        yield cell, view, data, sources


def _failure_part(data):
    failed = data.loc[~data.point_valid | ~data.interval_valid | data.pipeline_error.ne(""), FAILURE_COLUMNS]
    failed = failed.sort_values(["cell_id", "rep", "estimator"], kind="stable")
    return len(failed), csv_bytes(failed).split(b"\n", 1)[1]


def _failure_payload(parts):
    header = csv_bytes(pd.DataFrame(columns=FAILURE_COLUMNS))
    return gzip.compress(header+b"".join(parts), compresslevel=6, mtime=0)


def _raw_expectations(output, design, identity):
    outcomes, points, constants, parts, failed_count = {}, {}, {}, [], 0
    for cell, _, data, sources in _raw_cells(output, design, identity):
        for key, group in data.groupby(["cell_id", "estimator", "point_status", "interval_status"], dropna=False):
            outcomes[key] = len(group)
        for key, group in data.groupby(["cell_id", "estimator"]):
            points[key] = (int(group.point_valid.sum()), int(group.interval_valid.sum()))
        constants[cell.identity] = {name: sources[name].iloc[0] for name in _SOURCE_CONSTANTS}
        count, part = _failure_part(data)
        failed_count += count
        parts.append(part)
    return {"outcomes": outcomes, "points": points, "source_constants": constants,
            "failure_count": failed_count, "failure_payload": _failure_payload(parts)}


def _build_tables(output, design, identity):
    pieces = {name: [] for name in TABLE_KEYS if name not in {"regular_comparisons.csv", "matched_iid_comparisons.csv"}}
    parts, failed_count = [], 0
    for _, view, data, sources in _raw_cells(output, design, identity):
        for name, table in (
            ("parameter_summary.csv", ExperimentSummarizer().summarize(data, view)),
            ("paired_comparisons.csv", paired_estimator_comparisons(data, view)),
            ("spacing_summary.csv", spacing_summary(data, view)),
            ("source_summary.csv", _summarize_sources(sources, view)),
            ("outcome_counts.csv", outcome_counts(data, view)),
            ("numerical_diagnostics.csv", numerical_diagnostics(data, view)),
        ):
            pieces[name].append(_table(name, table, view))
        count, part = _failure_part(data)
        failed_count += count
        parts.append(part)
    tables = {name: _table(name, pd.concat(items, ignore_index=True), design) for name, items in pieces.items()}
    for kind, name in (("regular", "regular_comparisons.csv"), ("matched_iid", "matched_iid_comparisons.csv")):
        tables[name] = independent_comparisons(tables["parameter_summary.csv"], design, kind)
    payloads = {name: csv_bytes(table) for name, table in tables.items()}
    payloads["failed_fits.csv.gz"] = _failure_payload(parts)
    counts = {name: len(table) for name, table in tables.items()}
    counts["failed_fits.csv.gz"] = failed_count
    manifest = {name: {"rows": counts[name], "sha256": digest_bytes(payloads[name])} for name in sorted(payloads)}
    return payloads, manifest


def _csv_rows(path, name):
    raw = path.read_bytes()
    if name.endswith(".gz"):
        raw = gzip.decompress(raw)
    reader = csv.reader(io.StringIO(raw.decode("utf-8"), newline=""))
    header = next(reader, None)
    if header != list(TABLE_COLUMNS[name]):
        raise ValueError("saved table has a missing, duplicate, reordered or unknown column: " + name)
    rows = []
    for values in reader:
        if len(values) != len(header):
            raise ValueError("saved table row width differs from its exact schema")
        rows.append(dict(zip(header, values)))
    return rows


def _csv_integer(value, label):
    if type(value) is not str or not value or any(character not in "0123456789" for character in value):
        raise ValueError("noncanonical nonnegative table count: " + label)
    number = int(value)
    if str(number) != value:
        raise ValueError("noncanonical integer table representation: " + label)
    return number


def _csv_value(value, expected, label):
    if isinstance(expected, (bool, np.bool_)):
        matches = value == str(bool(expected))
    elif isinstance(expected, str):
        matches = value == expected
    elif isinstance(expected, (int, np.integer)):
        matches = _csv_integer(value, label) == expected
    else:
        try:
            number = float(value)
            matches = math.isfinite(number) and number == float(expected)
        except (ValueError, TypeError):
            matches = False
    if not matches:
        raise ValueError("table scientific/count metadata changed: " + label)


def _planned_keys(name, design, expected):
    if name == "outcome_counts.csv":
        return set(expected["outcomes"])
    if name == "paired_comparisons.csv":
        return {(cell.identity, a, b, parameter) for cell in design.cells
                for a, b in combinations(ESTIMATORS, 2) for parameter in PARAMETERS}
    if name in {"regular_comparisons.csv", "matched_iid_comparisons.csv"}:
        kind = "regular" if name.startswith("regular") else "matched_iid"
        ids = reference_map(design)[kind]
    else:
        ids = [cell.identity for cell in design.cells]
    if name in {"spacing_summary.csv", "source_summary.csv"}:
        return {(cid,) for cid in ids}
    if name == "numerical_diagnostics.csv":
        return {(cid, method) for cid in ids for method in ESTIMATORS}
    return {(cid, method, parameter) for cid in ids for method in ESTIMATORS for parameter in PARAMETERS}


def _verify_csv(output, name, design, expected):
    rows = _csv_rows(output/name, name)
    if name == "failed_fits.csv.gz":
        if len(rows) != expected["failure_count"] or (output/name).read_bytes() != expected["failure_payload"]:
            raise ValueError("failed-fit export must preserve every complete raw failure record exactly")
        return len(rows)
    keys = [tuple(row[field] for field in TABLE_KEYS[name]) for row in rows]
    if len(set(keys)) != len(keys) or set(keys) != _planned_keys(name, design, expected):
        raise ValueError("saved table does not equal its complete planned key set: " + name)
    cells = {cell.identity: cell.to_dict() for cell in design.cells}
    for key, row in zip(keys, rows):
        cid = row["cell_id"]
        for field, value in cells[cid].items():
            _csv_value(row[field], value, name+":"+field)
        if "parameter" in row:
            _csv_value(row["truth"], cells[cid]["true_"+row["parameter"]], "truth")
        if "n_planned" in row:
            _csv_value(row["n_planned"], design.replications, "n_planned")
        if name in {"regular_comparisons.csv", "matched_iid_comparisons.csv"}:
            kind = "regular" if name.startswith("regular") else "matched_iid"
            ref = reference_map(design)[kind][cid]
            for field, value in cells[ref].items():
                _csv_value(row["reference_"+field], value, "reference_"+field)
            _csv_value(row["comparison_kind"], kind, "comparison_kind")
            _csv_value(row["comparison_type"], "independent_cells", "comparison_type")
            _csv_value(row["point_subset"], "each_cell_own_valid_point_fits", "point_subset")
            _csv_value(row["difference_direction"], "cell_minus_reference", "difference_direction")
            _csv_value(row["ratio_direction"], "cell_over_reference", "ratio_direction")
            for side, other in (("cell", cid), ("reference", ref)):
                _csv_value(row[side+"_n_planned"], design.replications, side+"_n_planned")
                for field, value in zip(("n_point_valid", "n_interval_valid"), expected["points"][(other, row["estimator"])]):
                    _csv_value(row[side+"_"+field], value, side+"_"+field)
        elif name == "outcome_counts.csv":
            _csv_value(row["count"], expected["outcomes"][key], "outcome count")
        elif name in {"parameter_summary.csv", "numerical_diagnostics.csv"}:
            point, interval = expected["points"][(cid, row["estimator"])]
            _csv_value(row["n_interval_valid"], interval, "n_interval_valid")
            if name == "parameter_summary.csv":
                _csv_value(row["n_point_valid"], point, "n_point_valid")
        elif name == "source_summary.csv":
            for field, value in {**expected["source_constants"][cid], "n_unique_sources": design.replications,
                                 "n_source_valid": design.replications, "n_fit_records": design.replications*len(ESTIMATORS)}.items():
                _csv_value(row[field], value, field)
        for field in row:
            if field.endswith("_n_available"):
                _csv_value(row[field], design.replications, field)
    return len(rows)


def _saved_tables(output, design, binding, expected, *, required):
    binding_path, manifest_path = output/"table_binding.json", output/"table_manifest.json"
    saved_binding = manifest = None
    if binding_path.exists():
        saved_binding, _ = execution._json_file(binding_path)
        if canonical(saved_binding) != canonical(_binding_record(design, binding)):
            raise ValueError("immutable table binding differs from current raw inputs")
    elif required or manifest_path.exists():
        raise ValueError("saved tables require their immutable raw-input binding")
    if manifest_path.exists():
        manifest, _ = execution._json_file(manifest_path)
        if set(manifest) != set(TABLE_COLUMNS):
            raise ValueError("table manifest must declare exactly the nine E4 artifacts")
        for name, entry in manifest.items():
            if not isinstance(entry, dict) or set(entry) != {"rows", "sha256"}:
                raise ValueError("table manifest entry has unexpected fields")
            integer(entry["rows"], "table row count", 0)
            if file_digest(output/name) != entry["sha256"] or _verify_csv(output, name, design, expected) != entry["rows"]:
                raise ValueError("saved table bytes/rows differ from the declared manifest")
    elif required:
        raise ValueError("complete table manifest is missing")
    else:
        for name in TABLE_COLUMNS:
            if (output/name).exists():
                _verify_csv(output, name, design, expected)
    return saved_binding, manifest


@dataclass(frozen=True)
class ReportingPreflight:
    raw: execution.Preflight
    binding: dict
    table_binding: dict | None
    manifest: dict | None


def _extra_files(extra_files):
    if isinstance(extra_files, (str, bytes)):
        raise ValueError("extra_files must be an iterable of repository-relative file names")
    return execution._dependencies((*FILES, *tuple(extra_files)))


def _prepare_reporting(design, output, *, extra_files, mutable, required):
    validate_design(design)
    output = execution._output_path(output) if mutable else execution._ordinary_path(output)
    raw = execution.readback(design, output, extra_files=_extra_files(extra_files))
    if (raw.metadata["status"] not in _READABLE_STATES
            or (mutable and raw.metadata["status"] != "raw_complete_pending_audit")
            or raw.metadata.get("pipeline_errors") != 0 or not raw.cached_completion):
        raise ValueError("reporting requires complete zero-pipeline-error raw evidence at the permitted stage")
    binding = raw_input_binding(output, design, raw.metadata)
    expected = _raw_expectations(output, design, binding["run_identity"])
    saved_binding, manifest = _saved_tables(output, design, binding, expected, required=required)
    if canonical(raw_input_binding(output, design, raw.metadata)) != canonical(binding):
        raise ValueError("raw inputs changed during read-only reporting preflight")
    execution._check_frozen(raw.metadata["provenance"], design, raw.provenance, raw.require_frozen)
    if (output/"run_metadata.json").read_bytes() != raw.expected_metadata_bytes:
        raise ValueError("root metadata changed during read-only reporting preflight")
    return ReportingPreflight(raw, binding, saved_binding, manifest)


def preflight_reporting(design, output, *, extra_files=()):
    """Pure proposed-write preflight, including all source triplets and partial files."""
    return _prepare_reporting(design, output, extra_files=extra_files, mutable=True, required=False)


def verify_tables(design, output, *, extra_files=()):
    """Read-only schemas/keys/raw-binding verification, also valid for closed runs."""
    prepared = _prepare_reporting(design, output, extra_files=extra_files, mutable=False, required=True)
    return {"tables": prepared.manifest, "binding": prepared.table_binding}


def _verify_stage(prepared, design, guard, stage):
    guard.check(stage+":root")
    execution._root_layout(prepared.raw.output, design)
    execution._check_frozen(prepared.raw.metadata["provenance"], design, prepared.raw.provenance, prepared.raw.require_frozen)
    if canonical(raw_input_binding(prepared.raw.output, design, prepared.raw.metadata)) != canonical(prepared.binding):
        raise ValueError("raw input binding changed during table publication")
    guard.check(stage+":after_input_verification")


def write_tables(design, output, *, extra_files=()):
    """Publish/adopt complete immutable tables, without changing root lifecycle."""
    extra_files = _extra_files(extra_files)
    prepared = preflight_reporting(design, output, extra_files=extra_files)
    output = prepared.raw.output
    with run_lock(output):
        prepared = preflight_reporting(design, output, extra_files=extra_files)
        guard = execution._ManifestGuard(output/"run_metadata.json", prepared.raw.expected_metadata_bytes)
        try:
            _verify_stage(prepared, design, guard, "before_calculations")
            payloads, manifest = _build_tables(output, design, prepared.binding["run_identity"])
            record = _binding_record(design, prepared.binding)
            candidates = {**payloads, "table_binding.json": (canonical(record)+"\n").encode(),
                          "table_manifest.json": (canonical(manifest)+"\n").encode()}
            # Validate every existing candidate before publishing any missing one.
            for name, payload in candidates.items():
                if (output/name).exists() and (output/name).read_bytes() != payload:
                    raise FileExistsError("refusing to replace different immutable table evidence: " + name)
            _verify_stage(prepared, design, guard, "before_first_publication")
            guard.check("before_table_binding")
            atomic_json(output/"table_binding.json", record, immutable=True)
            guard.check("after_table_binding")
            for name in sorted(payloads):
                guard.check("before_table:"+name)
                atomic_bytes(output/name, payloads[name], immutable=True)
                guard.check("after_table:"+name)
            _verify_stage(prepared, design, guard, "before_final_manifest")
            for name, payload in candidates.items():
                if name != "table_manifest.json" and (output/name).read_bytes() != payload:
                    raise ValueError("table bytes changed before manifest publication: " + name)
            guard.check("before_table_manifest")
            atomic_json(output/"table_manifest.json", manifest, immutable=True)
            _verify_stage(prepared, design, guard, "after_table_manifest")
            verified = verify_tables(design, output, extra_files=extra_files)
            guard.check("reporting_complete")
            return verified
        except BaseException:
            if not guard.drifted:
                guard.check("reporting_failure")
            raise
