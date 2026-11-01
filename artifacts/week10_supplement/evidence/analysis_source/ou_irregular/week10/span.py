"""Post-result span sensitivity using the immutable original E4 records.

Strata and support rules are fixed independently of estimator outcomes. The
standardized comparison is descriptive and conditional on observed stratum
counts and shared weights. It is not a fixed-horizon experiment or a causal
adjustment, and it does not replace the original unconditional E4 estimand.
"""
from __future__ import annotations

import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

from ..experiments.analysis import csv_bytes, read_cell
from ..experiments.provenance import file_digest
from ..experiments.storage import atomic_bytes, atomic_json
from ..week4.summaries import PARAMETERS, wilson_interval
from ..week8.design import full_design, parse_design


SPAN_EDGES = (0.0, 400.0, 450.0, 475.0, 500.0, 525.0, 550.0, 600.0, np.inf)
MIN_STRATUM_COUNT = 30
ESTIMATORS = ("exact", "pfml", "euler")
Z = float(norm.ppf(0.975))
CONTRACT = "week10_post_result_span_standardization_v1"
_BASE_COLUMNS = (
    "cell_id", "rep", "estimator", "mechanism", "ratio", "n", "realized_span",
    "point_valid", "interval_valid", "point_status", "interval_status", "pipeline_error",
    "source_valid", "source_seed", "source_spawn_key", "times_sha256", "observations_sha256",
)
COMPACT_COLUMNS = (*_BASE_COLUMNS,
                   *(name for p in PARAMETERS for name in
                     (p, "true_" + p, p + "_se", p + "_lower", p + "_upper")))


def read_e4_records(root):
    """Read all 42,000 planned fits after partition, key and source-pair checks.

    Existing E4 CellStore verification is read-only for an existing archive.
    The source binding contains every raw partition and its cell manifest,
    plus the design, lifecycle, numerical audit and original visual review.
    This readback does not claim a fresh numerical replay of old fits.
    """
    root = Path(root).resolve()
    design = parse_design(json.loads((root / "resolved_design.json").read_text()))
    if design.to_dict() != full_design().to_dict():
        raise ValueError("span sensitivity requires the complete original full E4 design")
    metadata = json.loads((root / "run_metadata.json").read_text())
    if metadata.get("status") != "completed":
        raise ValueError("span sensitivity requires the completed original E4 archive")
    names = ["resolved_design.json", "run_metadata.json", "audit.json", "visual_review.json",
             "parameter_summary.csv", "matched_iid_comparisons.csv"]
    frames = []
    for cell in design.cells:
        directory = root / "cells" / cell.identity
        manifest = json.loads((directory / "cell_manifest.json").read_text())
        names.append(str((directory / "cell_manifest.json").relative_to(root)))
        names.extend(str((directory / item["filename"]).relative_to(root))
                     for item in manifest["partitions"])
        data = read_cell(root, design, cell, metadata["provenance"]["identity_sha256"])
        frames.append(data.loc[:, COMPACT_COLUMNS])
    data = pd.concat(frames, ignore_index=True).sort_values(
        ["cell_id", "rep", "estimator"], kind="stable").reset_index(drop=True)
    _validate(data)
    binding = {"contract": CONTRACT, "source_design": design.to_dict(),
               "fit_records": len(data), "source_paths": len(data.drop_duplicates(["cell_id", "rep"])),
               "source_files": {name: {"sha256": file_digest(root / name),
                                        "bytes": (root / name).stat().st_size}
                                for name in sorted(names)}}
    return data, binding


def _validate(data):
    if data.empty or not set(COMPACT_COLUMNS).issubset(data):
        raise ValueError("span analysis requires complete compact fit records")
    if data.duplicated(["cell_id", "rep", "estimator"]).any():
        raise ValueError("duplicate fit key")
    if not data.estimator.isin(ESTIMATORS).all() or not data.mechanism.isin(("mmpp", "iid", "regular")).all():
        raise ValueError("unrecognized estimator or mechanism")
    if not np.isfinite(data.realized_span).all() or (data.realized_span <= 0).any():
        raise ValueError("all planned sources require a finite positive realized span")
    for field in ("point_valid", "interval_valid", "source_valid"):
        if not data[field].map(lambda x: isinstance(x, (bool, np.bool_))).all():
            raise ValueError("validity flags must be Boolean")
    if not data.source_valid.all() or data.pipeline_error.fillna("").ne("").any():
        raise ValueError("pipeline/source errors must be resolved before span reporting")
    if (data.interval_valid & ~data.point_valid).any():
        raise ValueError("interval-valid fits must be point-valid")
    if not data.point_valid.eq(data.point_status.eq("valid_interior")).all():
        raise ValueError("point status mismatch")
    if not data.interval_valid.eq(data.interval_status.eq("valid")).all():
        raise ValueError("interval status mismatch")
    for parameter in PARAMETERS:
        for field in (parameter, "true_" + parameter, parameter + "_se",
                      parameter + "_lower", parameter + "_upper"):
            if data[field].map(lambda x: isinstance(x, (bool, np.bool_))).any():
                raise ValueError("Boolean numeric fit field")
        if not np.isfinite(data["true_" + parameter]).all():
            raise ValueError("nonfinite truth")
        valid = data.point_valid
        if not np.isfinite(data.loc[valid, parameter]).all():
            raise ValueError("point-valid estimates must be finite")
        valid = data.interval_valid
        values = data.loc[valid, [parameter + "_se", parameter + "_lower", parameter + "_upper"]]
        if (not np.isfinite(values.to_numpy(dtype=float)).all()
                or (values[parameter + "_se"] <= 0).any()
                or (values[parameter + "_lower"] > values[parameter + "_upper"]).any()):
            raise ValueError("interval-valid fits require positive SEs and ordered finite limits")
    source_columns = ["mechanism", "ratio", "n", "realized_span", "source_seed",
                      "source_spawn_key", "times_sha256", "observations_sha256"]
    groups = data.groupby(["cell_id", "rep"], sort=False)
    if groups[source_columns].nunique(dropna=False).gt(1).any().any():
        raise ValueError("estimator fits do not share the same source")
    if not groups.estimator.agg(lambda values: set(values) == set(ESTIMATORS)).all():
        raise ValueError("all three estimators must survive for every planned source")
    for (_, _), group in data.groupby(["mechanism", "ratio"]):
        if group.cell_id.nunique() != 1:
            raise ValueError("multiple source cells for one mechanism/ratio")
        for p in PARAMETERS:
            if group["true_" + p].nunique() != 1:
                raise ValueError("truth must be constant within each cell")


def stratum_indices(spans, edges=SPAN_EDGES):
    """Left-closed, right-open positive-span strata, including the upper tail."""
    edges = np.asarray(edges, dtype=float)
    spans = np.asarray(spans, dtype=float)
    if (edges.ndim != 1 or len(edges) < 2 or edges[0] != 0 or not np.isposinf(edges[-1])
            or not np.isfinite(edges[:-1]).all() or not np.all(np.diff(edges) > 0)):
        raise ValueError("span edges must strictly increase from zero to positive infinity")
    if not np.isfinite(spans).all() or (spans <= 0).any():
        raise ValueError("spans must be finite and positive")
    return np.searchsorted(edges, spans, side="right") - 1


def _mean_variance(values):
    values = np.asarray(values, dtype=float)
    return float(values.var(ddof=1) / len(values)) if len(values) > 1 else np.nan


def _bounds(value, variance, *, degenerate_missing=False):
    if not np.isfinite([value, variance]).all() or variance < 0 or (degenerate_missing and variance == 0):
        return np.nan, np.nan, np.nan if not np.isfinite(variance) else float(np.sqrt(variance))
    se = float(np.sqrt(variance))
    return value - Z * se, value + Z * se, se


def _summary(group, parameter):
    n = len(group)
    points = group.point_valid.to_numpy(dtype=bool)
    intervals = group.interval_valid.to_numpy(dtype=bool)
    truth = group["true_" + parameter].to_numpy(dtype=float)
    errors = group[parameter].to_numpy(dtype=float)[points] - truth[points]
    covers = (intervals & (group[parameter + "_lower"].to_numpy(dtype=float) <= truth)
              & (truth <= group[parameter + "_upper"].to_numpy(dtype=float)))
    n_point, n_interval, n_covered = int(points.sum()), int(intervals.sum()), int(covers.sum())
    bias = float(errors.mean()) if n_point else np.nan
    mse = float(np.mean(errors ** 2)) if n_point else np.nan
    bias_variance, mse_variance = _mean_variance(errors), _mean_variance(errors ** 2)
    coverage = n_covered / n if n else np.nan
    coverage_variance = coverage * (1 - coverage) / n if n else np.nan
    low, high = wilson_interval(n_covered, n)
    result = {"n_planned": n, "n_point_valid": n_point, "n_interval_valid": n_interval,
              "n_covered": n_covered, "point_availability": n_point / n if n else np.nan,
              "interval_availability": n_interval / n if n else np.nan,
              "bias": bias, "bias_mean_variance": bias_variance,
              "mse": mse, "mse_mean_variance": mse_variance,
              "rmse": float(np.sqrt(mse)), "coverage_operational": coverage,
              "coverage_mean_variance": coverage_variance,
              "coverage_conditional": n_covered / n_interval if n_interval else np.nan,
              "coverage_wilson_lower": low, "coverage_wilson_upper": high,
              "span_mean": float(group.realized_span.mean()) if n else np.nan,
              "span_min": float(group.realized_span.min()) if n else np.nan,
              "span_max": float(group.realized_span.max()) if n else np.nan}
    for metric, value, variance in (("bias", bias, bias_variance),
                                     ("mse", mse, mse_variance),
                                     ("rmse", result["rmse"], mse_variance / (4 * mse) if mse > 0 else np.nan)):
        low, high, se = _bounds(value, variance)
        result.update({metric + "_mc_lower": low, metric + "_mc_upper": high, metric + "_mcse": se})
    result["coverage_mcse"] = float(np.sqrt(coverage_variance))
    return result


def _contrast(mmpp, iid):
    result = {}
    for metric, field in (("bias", "bias"), ("mse", "mse"), ("coverage", "coverage_operational")):
        difference = mmpp[field] - iid[field]
        variance = mmpp[metric + "_mean_variance"] + iid[metric + "_mean_variance"]
        low, high, se = _bounds(difference, variance, degenerate_missing=metric == "coverage")
        result.update({metric + "_difference": difference, metric + "_difference_mcse": se,
                       metric + "_difference_mc_lower": low, metric + "_difference_mc_upper": high})
    mse_m, mse_i = mmpp["mse"], iid["mse"]
    if np.isfinite([mse_m, mse_i]).all() and mse_m > 0 and mse_i > 0:
        ratio = float(np.sqrt(mse_m / mse_i))
        log_variance = 0.25 * (mmpp["mse_mean_variance"] / mse_m ** 2
                              + iid["mse_mean_variance"] / mse_i ** 2)
        log_se = float(np.sqrt(log_variance))
        low, high = ratio * np.exp(-Z * log_se), ratio * np.exp(Z * log_se)
    else:
        ratio, log_se, low, high = (np.nan,) * 4
    result.update({"rmse_ratio": ratio, "rmse_ratio_log_mcse": log_se,
                   "rmse_ratio_mc_lower": low, "rmse_ratio_mc_upper": high})
    for side, summary in (("mmpp", mmpp), ("iid", iid)):
        result.update({side + "_" + key: value for key, value in summary.items()})
    return result


def _standardize(summaries, weights):
    """Fixed-weight mixture of stratum-conditional valid-point error means."""
    weights = np.asarray(weights, dtype=float)
    if not len(weights):
        result = {key: np.nan for key in ("bias", "mse", "rmse", "bias_mean_variance", "mse_mean_variance",
                                          "coverage_operational", "coverage_mean_variance", "point_availability",
                                          "interval_availability", "span_mean")}
        result.update({key: 0 for key in ("n_planned", "n_point_valid", "n_interval_valid", "n_covered")})
        return result
    if not np.isclose(weights.sum(), 1) or (weights <= 0).any():
        raise ValueError("standardization weights must be positive and sum to one")
    result = {key: sum(item[key] for item in summaries)
              for key in ("n_planned", "n_point_valid", "n_interval_valid", "n_covered")}
    for field in ("bias", "mse", "coverage_operational", "point_availability", "interval_availability", "span_mean"):
        result[field] = float(np.dot(weights, [item[field] for item in summaries]))
    for metric in ("bias", "mse", "coverage"):
        field = metric + "_mean_variance"
        result[field] = float(np.dot(weights ** 2, [item[field] for item in summaries]))
    result["rmse"] = float(np.sqrt(result["mse"]))
    return result


def analyze_spans(records, *, edges=SPAN_EDGES, minimum_count=MIN_STRATUM_COUNT):
    """Return complete strata and raw/standardized independent-arm contrasts.

    Optional arguments exist for synthetic unit checks; production uses the
    frozen module defaults. No fitted-outcome selection changes the weights.
    """
    data = pd.DataFrame(records).copy()
    _validate(data)
    if isinstance(minimum_count, (bool, np.bool_)) or not isinstance(minimum_count, (int, np.integer)) or minimum_count < 2:
        raise ValueError("minimum stratum count must be an integer at least two")
    data["stratum"] = stratum_indices(data.realized_span, edges)
    data = data[data.mechanism.isin(("mmpp", "iid"))].copy()
    if data.empty:
        raise ValueError("no matched stochastic cells")
    edges = tuple(float(value) for value in edges)
    strata_rows, contrast_rows, standardized_rows, support_rows, raw_rows = [], [], [], [], []
    for ratio, ratio_data in data.groupby("ratio", sort=True):
        if set(ratio_data.mechanism) != {"mmpp", "iid"}:
            raise ValueError("every ratio requires both independent arms")
        sources = ratio_data.drop_duplicates(["cell_id", "rep"])
        counts = {(arm, h): len(sources[(sources.mechanism == arm) & (sources.stratum == h)])
                  for arm in ("mmpp", "iid") for h in range(len(edges) - 1)}
        totals = {arm: int((sources.mechanism == arm).sum()) for arm in ("mmpp", "iid")}
        eligible = [h for h in range(len(edges) - 1)
                    if min(counts["mmpp", h], counts["iid", h]) >= minimum_count]
        pooled = sum(counts[arm, h] for arm in ("mmpp", "iid") for h in eligible)
        weight = {h: (counts["mmpp", h] + counts["iid", h]) / pooled if h in eligible else 0.0
                  for h in range(len(edges) - 1)}
        for h in range(len(edges) - 1):
            support_rows.append({"ratio": ratio, "stratum": h, "span_lower_inclusive": edges[h],
                                 "span_upper_exclusive": edges[h + 1], "mmpp_count": counts["mmpp", h],
                                 "iid_count": counts["iid", h], "common_support": h in eligible,
                                 "pooled_standardization_weight": weight[h]})
        excluded = {arm + "_excluded_fraction": 1 - sum(counts[arm, h] for h in eligible) / totals[arm]
                    for arm in ("mmpp", "iid")}
        for estimator in ESTIMATORS:
            method_data = ratio_data[ratio_data.estimator == estimator]
            for parameter in PARAMETERS:
                base = {"ratio": ratio, "estimator": estimator, "parameter": parameter}
                arm_summaries = {}
                for arm in ("mmpp", "iid"):
                    arm_data = method_data[method_data.mechanism == arm]
                    arm_summaries[arm] = {}
                    for h in range(len(edges) - 1):
                        summary = _summary(arm_data[arm_data.stratum == h], parameter)
                        arm_summaries[arm][h] = summary
                        strata_rows.append({**base, "mechanism": arm, "stratum": h,
                                            "span_lower_inclusive": edges[h], "span_upper_exclusive": edges[h + 1],
                                            "common_support": h in eligible, "weight": weight[h], **summary})
                for h in range(len(edges) - 1):
                    contrast_rows.append({**base, "stratum": h, "common_support": h in eligible,
                                          **_contrast(arm_summaries["mmpp"][h], arm_summaries["iid"][h])})
                standardized = {arm: _standardize([arm_summaries[arm][h] for h in eligible],
                                                    [weight[h] for h in eligible])
                                for arm in ("mmpp", "iid")}
                standardized_rows.append({**base, "support_strata": len(eligible), **excluded,
                                          "uncertainty_scope": "conditional_on_observed_stratum_counts_and_weights",
                                          **_contrast(standardized["mmpp"], standardized["iid"])})
                raw_rows.append({**base, **_contrast(*(_summary(method_data[method_data.mechanism == arm], parameter)
                                                       for arm in ("mmpp", "iid")))})
    return {"span_strata.csv": pd.DataFrame(strata_rows),
            "span_stratum_contrasts.csv": pd.DataFrame(contrast_rows),
            "span_standardized.csv": pd.DataFrame(standardized_rows),
            "span_support.csv": pd.DataFrame(support_rows),
            "span_unadjusted.csv": pd.DataFrame(raw_rows)}


def write_span_analysis(root, output):
    """Create a separate, immutable saved-data supplement; never replace files."""
    root, output = Path(root).resolve(), Path(output).resolve()
    if output == root or output.is_relative_to(root) or root.is_relative_to(output):
        raise ValueError("span output must be separate from the original E4 archive")
    if output.exists():
        raise FileExistsError("span output must be a new directory")
    records, binding = read_e4_records(root)
    tables = analyze_spans(records)
    original = pd.read_csv(root / "matched_iid_comparisons.csv", float_precision="round_trip")
    checks = 0
    compared = {"mmpp_bias": "cell_bias", "iid_bias": "reference_bias",
                "mmpp_rmse": "cell_rmse", "iid_rmse": "reference_rmse",
                "mmpp_coverage_operational": "cell_coverage_operational",
                "iid_coverage_operational": "reference_coverage_operational",
                "rmse_ratio": "rmse_ratio", "coverage_difference": "operational_coverage_difference"}
    for row in tables["span_unadjusted.csv"].to_dict("records"):
        matching = original[(original.ratio == row["ratio"])
                            & (original.estimator == row["estimator"])
                            & (original.parameter == row["parameter"])]
        if len(matching) != 1:
            raise ValueError("original matched contrast is absent or ambiguous")
        prior = matching.iloc[0]
        for new_field, prior_field in compared.items():
            if not np.isclose(row[new_field], prior[prior_field], rtol=2e-13, atol=1e-15, equal_nan=True):
                raise ValueError("saved-record arithmetic differs from the original contrast: " + new_field)
            checks += 1
    tables["span_fit_records.csv.gz"] = records
    output.mkdir(parents=True)
    outputs = {}
    for name, frame in tables.items():
        payload = csv_bytes(frame)
        if name.endswith(".gz"):
            payload = gzip.compress(payload, compresslevel=6, mtime=0)
        atomic_bytes(output / name, payload, immutable=True)
        outputs[name] = {"sha256": file_digest(output / name), "rows": len(frame), "bytes": len(payload)}
    for name, entry in binding["source_files"].items():
        if file_digest(root / name) != entry["sha256"]:
            raise ValueError("original source changed during span analysis")
    manifest = {**binding, "span_edges": ["infinity" if np.isposinf(x) else x for x in SPAN_EDGES],
                "stratum_convention": "left_closed_right_open", "minimum_count_each_arm": MIN_STRATUM_COUNT,
                "weighting": "pooled_counts_normalized_on_common_support",
                "uncertainty": "independent_arm_delta_method_conditional_on_observed_stratum_counts_and_weights",
                "zero_coverage_variance": "MCSE_zero_but_normal_confidence_bounds_undefined;_stratum_Wilson_bounds_retained",
                "scope": "post_result_descriptive_sensitivity_not_causal_adjustment",
                "original_unadjusted_arithmetic_checks": checks,
                "new_simulated_sources": 0, "new_fits": 0, "outputs": outputs}
    atomic_json(output / "span_manifest.json", manifest, immutable=True)
    return manifest


def verify_span_analysis(root, output):
    """Recheck complete settings, source hashes and every exported table.

    This is a read-only arithmetic replay of saved E4 fits, not a new simulation
    or a refit. An incomplete manifest cannot shrink the required inventory.
    """
    root, output = Path(root).resolve(), Path(output).resolve()
    manifest = json.loads((output / "span_manifest.json").read_text())
    settings = {
        "contract": CONTRACT,
        "span_edges": ["infinity" if np.isposinf(x) else x for x in SPAN_EDGES],
        "stratum_convention": "left_closed_right_open",
        "minimum_count_each_arm": MIN_STRATUM_COUNT,
        "weighting": "pooled_counts_normalized_on_common_support",
        "uncertainty": "independent_arm_delta_method_conditional_on_observed_stratum_counts_and_weights",
        "zero_coverage_variance": "MCSE_zero_but_normal_confidence_bounds_undefined;_stratum_Wilson_bounds_retained",
        "scope": "post_result_descriptive_sensitivity_not_causal_adjustment",
        "original_unadjusted_arithmetic_checks": 216,
        "new_simulated_sources": 0, "new_fits": 0,
    }
    for field, value in settings.items():
        if manifest.get(field) != value:
            raise ValueError("span analysis setting changed: " + field)
    records, binding = read_e4_records(root)
    for field, value in binding.items():
        if manifest.get(field) != value:
            raise ValueError("original span input binding changed: " + field)
    tables = analyze_spans(records)
    tables["span_fit_records.csv.gz"] = records
    expected = set(tables) | {"span_manifest.json"}
    if (any(p.is_symlink() or not p.is_file() for p in output.iterdir())
            or {p.name for p in output.iterdir()} != expected
            or set(manifest.get("outputs", {})) != set(tables)):
        raise ValueError("span output inventory differs from the complete planned exports")
    for name, table in tables.items():
        path = output / name
        entry = manifest["outputs"][name]
        if (set(entry) != {"sha256", "rows", "bytes"} or file_digest(path) != entry["sha256"]
                or entry["rows"] != len(table) or entry["bytes"] != path.stat().st_size):
            raise ValueError("span table digest or accounting changed: " + name)
        payload = path.read_bytes()
        if name.endswith(".gz"):
            payload = gzip.decompress(payload)
        if payload != csv_bytes(table):
            raise ValueError("span raw-to-table arithmetic changed: " + name)
    for name, entry in binding["source_files"].items():
        if file_digest(root / name) != entry["sha256"]:
            raise ValueError("original source changed during span verification")
    return {"status": "passed", "fit_records": len(records), "source_paths": binding["source_paths"],
            "source_files": len(binding["source_files"]), "rederived_exports": len(tables),
            "new_simulated_sources": 0, "new_fits": 0}
