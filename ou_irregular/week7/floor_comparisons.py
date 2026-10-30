"""Paired epsilon contrasts on verified shared innovation addresses.

All point metrics are recomputed on the same common-valid replication subset.
Operational coverage uses every planned pair, with unavailable intervals coded
as not covering.  Shared normal innovations do not imply identical irregular
sample arrays; this is an explicitly exploratory floor coupling, not the
independent-cell or within-source-estimator comparison in earlier weeks.
"""
import json
import math

import numpy as np
import pandas as pd

from ..experiments.comparisons import _interval, _ratio, _relative
from ..experiments.design import ESTIMATORS, RunDesign
from ..experiments.provenance import canonical
from ..experiments.summaries import ExperimentSummarizer, covered, mean_mcse, validate_records
from ..week4.summaries import PARAMETERS
from .floor_source import (ANCHOR_FLOOR_FRACTION, COUPLING_POLICY, FLOOR_FRACTIONS,
                           anchor_cell, validate_cell)


_SOURCE_FIELDS = {"actual_cell_id", "seed_anchor_cell_id", "seed_anchor_floor_fraction",
                  "source_coupling_policy", "source_seed", "source_spawn_key", "gap_spawn_key",
                  "path_spawn_key", "times_sha256", "observations_sha256"}
_FIT_FIELDS = {"point_valid", "point_status", "interval_valid", "interval_status",
               "optimizer_success", "nll", *PARAMETERS,
               *(f"{p}_{suffix}" for p in PARAMETERS for suffix in ("se", "lower", "upper"))}
COLUMNS = (
    "cell_id", "reference_cell_id", "seed_anchor_cell_id", "floor_fraction", "reference_floor_fraction",
    "estimator", "parameter", "truth", "n", "theta_mean_gap", "mean_gap", "cv", "source_seed",
    "source_coupling_policy", "comparison_type", "difference_direction", "ratio_direction", "point_subset",
    "coverage_subset", "n_planned", "n_point_valid_a", "n_point_valid_b", "n_interval_valid_a", "n_interval_valid_b",
    "n_common_point_valid", "n_common_interval_valid", "bias_a_common", "bias_b_common", "mse_a_common",
    "mse_b_common", "rmse_a_common", "rmse_b_common", "relative_bias_a_common", "relative_bias_b_common",
    "relative_rmse_a_common", "relative_rmse_b_common", "squared_error_covariance_common",
    "bias_difference", "bias_difference_mcse", "bias_difference_mc_lower", "bias_difference_mc_upper",
    "mse_difference", "mse_difference_mcse", "mse_difference_mc_lower", "mse_difference_mc_upper",
    "relative_bias_difference", "relative_bias_difference_mcse", "relative_bias_difference_mc_lower",
    "relative_bias_difference_mc_upper", "rmse_difference", "rmse_difference_mcse", "rmse_difference_mc_lower",
    "rmse_difference_mc_upper", "rmse_ratio", "rmse_ratio_mcse", "rmse_ratio_log_mcse", "rmse_ratio_mc_lower",
    "rmse_ratio_mc_upper", "rmse_ratio_uncertainty", "n_covered_a", "n_covered_b", "n_both_cover",
    "n_a_only_cover", "n_b_only_cover", "n_neither_cover", "operational_coverage_a", "operational_coverage_b",
    "operational_coverage_difference", "operational_coverage_difference_mcse",
    "operational_coverage_difference_mc_lower", "operational_coverage_difference_mc_upper",
)


def _design_cells(design):
    if not isinstance(design, RunDesign) or design.scope not in {"validation", "exploratory"}:
        raise ValueError("paired-floor comparisons require an exploratory or validation RunDesign")
    groups = {}
    for cell in design.cells:
        validate_cell(cell.to_dict())
        groups.setdefault(cell.cv, {})[cell.floor_fraction] = cell
    if not groups or any(set(group) != set(FLOOR_FRACTIONS) for group in groups.values()):
        raise ValueError("each included CV requires all three fixed floor variants")
    return groups


def _spawn_key(value, expected, field):
    if not isinstance(value, str):
        raise ValueError(f"paired source {field} must retain its JSON spawn-key string")
    try:
        parsed = json.loads(value)
    except (ValueError, TypeError) as exc:
        raise ValueError(f"invalid paired source {field}") from exc
    if (not isinstance(parsed, list) or any(type(v) is not int for v in parsed)
            or tuple(parsed) != expected):
        raise ValueError(f"paired source {field} differs from the canonical anchor address")


def _hash(value, field):
    if (type(value) is not str or len(value) != 64
            or any(character not in "0123456789abcdef" for character in value)):
        raise ValueError(f"paired source {field} must be a SHA-256 hexadecimal digest")


def _coupling_records(records, design, groups):
    """Validate actual cells and expected addresses before claiming any pairing."""
    raw = pd.DataFrame(records)
    required = _SOURCE_FIELDS | _FIT_FIELDS | set(design.cells[0].to_dict()) | {"rep", "estimator"}
    if not required <= set(raw):
        raise ValueError(f"paired records lack required fields: {sorted(required - set(raw))}")
    # Check canonical scientific types before the historical numerical-column
    # conversion can hide an invalid representation.  Only validation columns
    # are iterated, avoiding a second full dictionary copy of the raw archive.
    cell_fields = tuple(design.cells[0].to_dict())
    planned = {cell.identity: cell for cell in design.cells}
    checked_fields = sorted(required | ({"source_valid", "source_error", "pipeline_error"} & set(raw)))
    for values in raw[checked_fields].itertuples(index=False, name=None):
        row = dict(zip(checked_fields, values))
        actual = validate_cell({name: row[name] for name in cell_fields})
        if actual.identity not in planned:
            raise ValueError("paired record contains an unplanned actual cell")
        if type(row["rep"]) is not int or not 0 <= row["rep"] < design.replications:
            raise ValueError("paired replication is outside the complete planned range")
        anchor = anchor_cell(actual)
        if (type(row["actual_cell_id"]) is not str or row["actual_cell_id"] != actual.identity
                or type(row["seed_anchor_cell_id"]) is not str or row["seed_anchor_cell_id"] != anchor.identity
                or type(row["seed_anchor_floor_fraction"]) is not float
                or row["seed_anchor_floor_fraction"] != ANCHOR_FLOOR_FRACTION
                or type(row["source_coupling_policy"]) is not str or row["source_coupling_policy"] != COUPLING_POLICY
                or type(row["source_seed"]) is not int or row["source_seed"] != design.seed):
            raise ValueError("paired source actual/anchor identity, seed or coupling policy changed")
        words = tuple(int(anchor.identity[i:i + 8], 16) for i in range(0, 64, 8))
        source_key = (4, *words, row["rep"])
        for field, expected in (("source_spawn_key", source_key),
                                ("gap_spawn_key", (*source_key, 0)),
                                ("path_spawn_key", (*source_key, 1))):
            _spawn_key(row[field], expected, field)
        for field in ("times_sha256", "observations_sha256"):
            _hash(row[field], field)
        for field in ("point_status", "interval_status"):
            if type(row[field]) is not str or not row[field]:
                raise ValueError("point and interval statuses must retain nonempty diagnostic strings")
        numeric_fields = {"nll", *PARAMETERS,
                          *(f"{p}_{suffix}" for p in PARAMETERS for suffix in ("se", "lower", "upper"))}
        for field in numeric_fields:
            value = row[field]
            if value is not None and (isinstance(value, (bool, np.bool_))
                    or not isinstance(value, (int, float, np.integer, np.floating))):
                raise ValueError(f"paired numerical field {field} must be numeric or missing")
        if "source_valid" in row and (type(row["source_valid"]) is not bool or not row["source_valid"]):
            raise ValueError("unavailable or invalid source evidence cannot certify paired paths")
        if "source_error" in row and (type(row["source_error"]) is not str or row["source_error"]):
            raise ValueError("source-error evidence cannot certify paired paths")
        if "pipeline_error" in row and type(row["pipeline_error"]) is not str:
            raise ValueError("pipeline errors must retain their diagnostic strings")
        if row.get("pipeline_error") and row["interval_valid"]:
            raise ValueError("a pipeline error cannot accompany a valid interval")
    data = validate_records(raw, design)
    # This historical component enforces valid-point, valid-interval, endpoint,
    # SE, NLL and status semantics.  Its own-valid summaries are NOT substituted
    # for the common-valid point metrics calculated below.
    ExperimentSummarizer().summarize(data, design)
    if 0. in groups:
        columns = [name for name in data.columns if name not in {"cell_id", "actual_cell_id", "floor_fraction"}]
        reference = groups[0.][ANCHOR_FLOOR_FRACTION].identity
        baseline = data[data.cell_id == reference].sort_values(["rep", "estimator"])
        reference_bytes = canonical(baseline[columns].to_dict(orient="records"))
        for epsilon, cell in groups[0.].items():
            current = data[data.cell_id == cell.identity].sort_values(["rep", "estimator"])
            if canonical(current[columns].to_dict(orient="records")) != reference_bytes:
                raise ValueError("CV=0 paired negative-control arrays or numerical evidence differ across floors")
    return data


def _point_contrasts(row, errors_a, errors_b, truth):
    count = len(errors_a)
    squares_a, squares_b = errors_a ** 2, errors_b ** 2
    mean_a = float(np.mean(squares_a)) if count else np.nan
    mean_b = float(np.mean(squares_b)) if count else np.nan
    rmse_a, rmse_b = float(np.sqrt(mean_a)), float(np.sqrt(mean_b))
    bias_a = float(np.mean(errors_a)) if count else np.nan
    bias_b = float(np.mean(errors_b)) if count else np.nan
    row.update(bias_a_common=bias_a, bias_b_common=bias_b,
               mse_a_common=mean_a, mse_b_common=mean_b,
               rmse_a_common=rmse_a, rmse_b_common=rmse_b,
               relative_bias_a_common=bias_a / abs(truth) if truth else np.nan,
               relative_bias_b_common=bias_b / abs(truth) if truth else np.nan,
               relative_rmse_a_common=rmse_a / abs(truth) if truth else np.nan,
               relative_rmse_b_common=rmse_b / abs(truth) if truth else np.nan,
               squared_error_covariance_common=float(np.cov(squares_a, squares_b, ddof=1)[0, 1]) if count > 1 else np.nan)
    for name, differences in (("bias_difference", errors_a - errors_b),
                              ("mse_difference", squares_a - squares_b)):
        _interval(row, name, float(np.mean(differences)) if count else np.nan, mean_mcse(differences))
    _relative(row, "bias_difference", truth)
    # Each expression is one paired influence vector.  Its sample variance
    # includes the cross-variant covariance; independent variances cannot be
    # added for these common-random-number samples.
    rmse_se = mean_mcse(squares_a / (2 * rmse_a) - squares_b / (2 * rmse_b)) if rmse_a > 0 and rmse_b > 0 else np.nan
    _interval(row, "rmse_difference", rmse_a - rmse_b, rmse_se)
    log_se = mean_mcse((squares_a / mean_a - squares_b / mean_b) / 2) if mean_a > 0 and mean_b > 0 else np.nan
    # _ratio preserves numerator-zero ratio=0 when the denominator is positive,
    # while withholding its singular square-root/log delta uncertainty.
    _ratio(row, rmse_a, rmse_b, log_se)


def epsilon_comparisons(records, design):
    """Epsilon A versus registered-floor B on verified paired replications.

    The full twelve-cell design produces 72 rows: 4 CVs * 2 nonreference floors
    * 3 estimators * 3 parameters.  Complete three-floor CV subsets are supported
    for bounded reporting, giving 18 rows per included CV.  No reference-to-self
    row is generated.  CV=0 is an exact shared negative control, not an
    independent sample or three independent control datasets.

    The returned MC intervals concern simulation-summary uncertainty.  They do
    not replace the registered parameter Wald intervals.  Statistical failures
    survive the planned denominators; invalid pairing evidence or overflowing
    arithmetic blocks the table instead of silently discarding any record.
    """
    groups = _design_cells(design)
    try:
        with np.errstate(over="raise", invalid="raise", divide="raise"):
            data = _coupling_records(records, design, groups)
            rows = []
            for cv, variants in sorted(groups.items()):
                reference = variants[ANCHOR_FLOOR_FRACTION]
                for epsilon in FLOOR_FRACTIONS:
                    if epsilon == ANCHOR_FLOOR_FRACTION:
                        continue
                    cell = variants[epsilon]
                    for method in ESTIMATORS:
                        a = data[(data.cell_id == cell.identity) & (data.estimator == method)].set_index("rep").sort_index()
                        b = data[(data.cell_id == reference.identity) & (data.estimator == method)].set_index("rep").sort_index()
                        common = a.point_valid.to_numpy(dtype=bool) & b.point_valid.to_numpy(dtype=bool)
                        common_intervals = a.interval_valid.to_numpy(dtype=bool) & b.interval_valid.to_numpy(dtype=bool)
                        for parameter in PARAMETERS:
                            truth = cell.to_dict()[f"true_{parameter}"]
                            row = {"cell_id": cell.identity, "reference_cell_id": reference.identity,
                                   "seed_anchor_cell_id": reference.identity, "floor_fraction": epsilon,
                                   "reference_floor_fraction": ANCHOR_FLOOR_FRACTION,
                                   "estimator": method, "parameter": parameter, "truth": truth,
                                   "n": cell.n, "theta_mean_gap": cell.q, "mean_gap": cell.q / cell.theta, "cv": cv,
                                   "source_seed": design.seed, "source_coupling_policy": COUPLING_POLICY,
                                   "comparison_type": "paired_floor_shared_standardized_innovations",
                                   "difference_direction": "epsilon_minus_registered_floor",
                                   "ratio_direction": "epsilon_over_registered_floor",
                                   "point_subset": "intersection_of_valid_point_fits",
                                   "coverage_subset": "all_planned_paths_unavailable_counts_as_not_covering",
                                   "n_planned": design.replications,
                                   "n_point_valid_a": int(a.point_valid.sum()), "n_point_valid_b": int(b.point_valid.sum()),
                                   "n_interval_valid_a": int(a.interval_valid.sum()), "n_interval_valid_b": int(b.interval_valid.sum()),
                                   "n_common_point_valid": int(common.sum()),
                                   "n_common_interval_valid": int(common_intervals.sum())}
                            errors_a = a[parameter].to_numpy(dtype=float)[common] - truth
                            errors_b = b[parameter].to_numpy(dtype=float)[common] - truth
                            _point_contrasts(row, errors_a, errors_b, truth)
                            coverage_a, coverage_b = covered(a, parameter), covered(b, parameter)
                            differences = coverage_a.astype(float) - coverage_b.astype(float)
                            row.update(n_covered_a=int(coverage_a.sum()), n_covered_b=int(coverage_b.sum()),
                                       n_both_cover=int((coverage_a & coverage_b).sum()),
                                       n_a_only_cover=int((coverage_a & ~coverage_b).sum()),
                                       n_b_only_cover=int((~coverage_a & coverage_b).sum()),
                                       n_neither_cover=int((~coverage_a & ~coverage_b).sum()),
                                       operational_coverage_a=float(coverage_a.mean()),
                                       operational_coverage_b=float(coverage_b.mean()))
                            _interval(row, "operational_coverage_difference", float(differences.mean()), mean_mcse(differences))
                            if any(isinstance(v, (float, np.floating)) and math.isinf(v) for v in row.values()):
                                raise ValueError("nonfinite paired-floor contrast arithmetic")
                            if set(row) != set(COLUMNS):
                                raise ValueError("paired-floor comparison column schema changed")
                            rows.append(row)
    except (FloatingPointError, OverflowError) as exc:
        raise ValueError("paired-floor metric or Monte Carlo arithmetic overflow") from exc
    result = pd.DataFrame(rows, columns=COLUMNS)
    expected = len(groups) * 2 * len(ESTIMATORS) * len(PARAMETERS)
    if len(result) != expected or result.duplicated(["cell_id", "estimator", "parameter"]).any():
        raise ValueError("paired-floor comparison row accounting failed")
    return result
