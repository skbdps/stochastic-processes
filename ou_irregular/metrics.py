"""Monte Carlo summaries and equidistant controls for the Week 3 harness.

The runner must emit a row for *every* planned replication/estimator, including
failed fits. Bias and RMSE below are conditional on the recorded, versioned
point-validity policy; exclusion rates always use all emitted/planned rows. Neither
discarding a failed row nor replacing its estimate with the truth is valid.

Week 3 does not construct confidence intervals. Its coverage columns therefore
remain missing unless callers explicitly provide actual interval endpoints.
"""

from __future__ import annotations

from itertools import combinations

import numpy as np
import pandas as pd


PARAMETERS = ("theta", "mu", "sigma")
_META = ("theta_mean_gap", "cv", "n", "mean_gap", "true_theta", "true_mu", "true_sigma")
FIT_STATUSES = (
    "valid_interior", "boundary", "boundary_suspected", "degenerate",
    "numerical_failure", "simulation_failure",
)
VALIDITY_POLICY = "classified_interior_v2"
LEGACY_VALIDITY_POLICY = "legacy_success_finite_v1"
_CLASSIFICATION = ("optimizer_success", "fit_status", "valid_for_point_summary")
_COUNT_COLUMNS = [
    "n_planned", "n_valid", "n_optimizer_success", "n_excluded", "validity_policy",
    *(f"n_{status}" for status in FIT_STATUSES), "n_legacy_usable", "n_legacy_excluded",
]
_SUMMARY_COLUMNS = [
    "cell_id", *_META, "estimator", "parameter", "truth", "n_total", "n_success",
    "n_failed", "failure_rate", *_COUNT_COLUMNS, "bias", "relative_bias", "rmse", "relative_rmse",
    "mcse_bias", "mcse_rmse", "relative_mcse_bias", "relative_mcse_rmse",
    "coverage", "coverage_mcse", "n_ci_valid", "ci_status",
    "metric_scope", "coverage_scope",
]


def _boolean_column(results: pd.DataFrame, column: str) -> np.ndarray:
    if not results[column].map(lambda value: isinstance(value, (bool, np.bool_))).all():
        raise ValueError(f"{column} must contain non-missing Boolean values")
    return results[column].to_numpy(dtype=bool)


def point_validity(results: pd.DataFrame) -> pd.Series:
    """Return the one fit-level mask used by exports, manifests and summaries.

    New records must contain all three classification columns. A usable row is
    explicitly ``valid_interior``, has optimizer success, finite theta/mu/sigma
    (and finite NLL when recorded), and positive theta/sigma. The recorded
    ``valid_for_point_summary`` must agree. Contradictory or incomplete records
    raise instead of silently changing denominators. ``success`` remains an
    alias for optimizer termination, never a replacement for validity.

    A table with none of these columns is a historical v1 table: retain exactly
    its original success-and-all-three-finite rule, without applying new domain
    checks or retrospectively inferring boundaries. Mixed schemas are rejected.
    No estimate is overwritten and invalid candidate values are preserved.
    """
    required = {"success", *PARAMETERS}
    missing = required.difference(results.columns)
    if missing:
        raise ValueError(f"results missing validity columns: {sorted(missing)}")
    success = _boolean_column(results, "success")
    finite = np.isfinite(results[list(PARAMETERS)].to_numpy(dtype=float)).all(axis=1)
    present = set(_CLASSIFICATION).intersection(results.columns)
    if not present:
        return pd.Series(success & finite, index=results.index, name="point_valid")
    if present != set(_CLASSIFICATION):
        raise ValueError("all classification columns are required; do not mix legacy and classified records")
    optimizer_success = _boolean_column(results, "optimizer_success")
    declared = _boolean_column(results, "valid_for_point_summary")
    if not np.array_equal(success, optimizer_success):
        raise ValueError("success must remain an alias of optimizer_success")
    if not results["fit_status"].isin(FIT_STATUSES).all():
        raise ValueError(f"fit_status must be one of {FIT_STATUSES}")
    interior = results["fit_status"].eq("valid_interior").to_numpy()
    admissible = finite & results["theta"].gt(0).to_numpy() & results["sigma"].gt(0).to_numpy()
    if "nll" in results:
        admissible &= np.isfinite(results["nll"].to_numpy(dtype=float))
    if np.any(interior & ~(optimizer_success & admissible)):
        raise ValueError("valid_interior requires optimizer success and finite admissible parameters/NLL")
    if not np.array_equal(declared, interior):
        raise ValueError("valid_for_point_summary must agree with fit_status=valid_interior")
    return pd.Series(declared, index=results.index, name="point_valid")


def validity_counts(results: pd.DataFrame) -> dict:
    """Count every emitted fit once; caller verifies the planned key set.

    ``n_planned`` is the emitted-row denominator, justified only after the
    runner's complete-key-set check. Legacy outcome buckets are explicitly
    separate: an old usable record is not certified as a new interior fit.
    """
    valid = point_validity(results)
    classified = all(column in results for column in _CLASSIFICATION)
    count = int(valid.sum())
    counts = {"n_planned": len(results), "n_valid": count,
              "n_optimizer_success": int(results["success"].sum()),
              "n_excluded": len(results) - count,
              "validity_policy": VALIDITY_POLICY if classified else LEGACY_VALIDITY_POLICY}
    counts.update({f"n_{status}": int(results["fit_status"].eq(status).sum()) if classified else 0
                   for status in FIT_STATUSES})
    counts.update(n_legacy_usable=0 if classified else count,
                  n_legacy_excluded=0 if classified else len(results) - count)
    return counts


def _validate_results(results: pd.DataFrame) -> pd.Series:
    required = {"cell_id", *_META, "rep", "estimator", "success", *PARAMETERS}
    missing = required.difference(results.columns)
    if missing:
        raise ValueError(f"results missing required columns: {sorted(missing)}")
    usable = point_validity(results)
    if results[["cell_id", "estimator", "rep", *_META]].isna().any().any():
        raise ValueError("cell metadata, estimator, and replication IDs cannot be missing")
    if results.duplicated(["cell_id", "estimator", "rep"]).any():
        raise ValueError("duplicate cell_id/estimator/rep rows would double-count replications")
    for column in _META:
        if not np.isfinite(results[column].to_numpy(dtype=float)).all():
            raise ValueError(f"metadata {column} must be finite")
        if results.groupby("cell_id")[column].nunique().gt(1).any():
            raise ValueError(f"metadata {column} is inconsistent within a cell_id")
    return usable


def summarize_results(results: pd.DataFrame) -> pd.DataFrame:
    """Return one row per cell, estimator, and parameter.

    For S point-valid fits, bias = sum(estimate - truth)/S and
    RMSE = sqrt(sum((estimate - truth)**2)/S). The Monte Carlo standard error
    of estimated bias is the sample SD of estimates divided by sqrt(S), using
    S-1 in the SD, and is undefined when S < 2. Applying the delta method to
    sqrt(mean(squared errors)) gives MCSE(RMSE) = sample_SD(squared errors) /
    (2*sqrt(S)*RMSE), undefined for S < 2 or RMSE = 0. This approximation need
    not be accurate in a tiny smoke sample. These are simulation-error measures,
    not confidence intervals for individual OU parameter estimates.

    ``point_validity`` supplies a shared versioned eligibility rule. Every
    parameter row has identical fit counts. ``n_success`` is a compatibility
    alias for ``n_valid``; optimizer termination has its own explicit count.
    ``n_failed``/``failure_rate`` include every excluded boundary or failed fit.
    Relative quantities divide by abs(truth) and are undefined at zero;
    in particular mu=0 must be assessed through absolute bias and RMSE.

    Optional ``{parameter}_lower`` / ``{parameter}_upper`` endpoints are used
    only for point-valid fits with finite, correctly ordered endpoints.
    Coverage uses n_ci_valid as its denominator (not n_total or n_success),
    includes the interval endpoints, and reports binomial plug-in MCSE. Missing
    pairs leave coverage explicitly uncomputed; a half-supplied pair is an
    input error. No interval estimator is implemented here.
    """
    _validate_results(results)
    interval_parameters = set()
    for parameter in PARAMETERS:
        lower, upper = f"{parameter}_lower", f"{parameter}_upper"
        if (lower in results) != (upper in results):
            raise ValueError(f"supply both {lower} and {upper}, or neither")
        if lower in results:
            interval_parameters.add(parameter)
    if results.empty:
        return pd.DataFrame(columns=_SUMMARY_COLUMNS)
    rows = []
    for (cell_id, estimator), group in results.groupby(["cell_id", "estimator"], sort=True):
        usable = point_validity(group).to_numpy()
        counts = validity_counts(group)
        legacy = counts["validity_policy"] == LEGACY_VALIDITY_POLICY
        count = counts["n_valid"]
        n_total = len(group)
        metadata = {column: group[column].iloc[0] for column in _META}
        for parameter in PARAMETERS:
            truth = float(metadata[f"true_{parameter}"])
            estimates = group[parameter].to_numpy(dtype=float)[usable]
            errors = estimates - truth
            squared_errors = errors**2
            bias = float(errors.mean()) if count else np.nan
            rmse = float(np.sqrt(np.mean(squared_errors))) if count else np.nan
            # Repeated paths are independent Monte Carlo draws within a cell.
            # Sharing paths across estimators does not change this marginal SE.
            mcse_bias = float(estimates.std(ddof=1) / np.sqrt(count)) if count >= 2 else np.nan
            # Delta method: for g(MSE)=sqrt(MSE), g'(MSE)=1/(2*RMSE).
            # At RMSE=0 this derivative is undefined, so do not report zero SE.
            mcse_rmse = (float(squared_errors.std(ddof=1) / (2.0 * np.sqrt(count) * rmse))
                         if count >= 2 and rmse > 0 else np.nan)
            relative_bias = bias / abs(truth) if truth != 0 else np.nan
            relative_rmse = rmse / abs(truth) if truth != 0 else np.nan
            relative_mcse_bias = mcse_bias / abs(truth) if truth != 0 else np.nan
            relative_mcse_rmse = mcse_rmse / abs(truth) if truth != 0 else np.nan
            coverage, coverage_mcse, n_ci_valid = np.nan, np.nan, 0
            ci_status = "not_computed_week3"
            coverage_scope = "not_computed"
            if parameter in interval_parameters:
                lower = group[f"{parameter}_lower"].to_numpy(dtype=float)
                upper = group[f"{parameter}_upper"].to_numpy(dtype=float)
                valid = usable & np.isfinite(lower) & np.isfinite(upper) & (lower <= upper)
                n_ci_valid = int(valid.sum())
                ci_status = "computed" if n_ci_valid else "no_valid_intervals"
                coverage_scope = ("conditional_on_success_and_valid_interval" if legacy
                                  else "conditional_on_point_valid_and_valid_interval")
                if n_ci_valid:
                    covered = (lower[valid] <= truth) & (truth <= upper[valid])
                    coverage = float(covered.mean())
                    coverage_mcse = float(np.sqrt(coverage * (1.0 - coverage) / n_ci_valid))
            rows.append({
                "cell_id": cell_id, **metadata, "estimator": estimator,
                "parameter": parameter, "truth": truth,
                "n_total": n_total, "n_success": count, "n_failed": n_total - count,
                "failure_rate": (n_total - count) / n_total, **counts,
                "bias": bias, "relative_bias": relative_bias, "rmse": rmse,
                "relative_rmse": relative_rmse, "mcse_bias": mcse_bias,
                "mcse_rmse": mcse_rmse, "relative_mcse_bias": relative_mcse_bias,
                "relative_mcse_rmse": relative_mcse_rmse,
                "coverage": coverage, "coverage_mcse": coverage_mcse,
                "n_ci_valid": n_ci_valid, "ci_status": ci_status,
                "metric_scope": ("conditional_on_success" if legacy else "conditional_on_point_valid"),
                "coverage_scope": coverage_scope,
            })
    return pd.DataFrame(rows, columns=_SUMMARY_COLUMNS)


def add_twin_comparisons(summary: pd.DataFrame) -> pd.DataFrame:
    """Attach the CV=0 control with the same configured design and truth.

    The equidistant twin matches n, nominal mean_gap, theta*mean_gap, estimator,
    parameter, and available truths. It does not match realized total span.
    Flooring very small Gamma gaps raises their expected mean above the nominal
    mean, so this control is a *nominal-design* comparison, not a guaranteed
    equal-expected-span experiment. Realized-gap diagnostics belong in results.
    It is also not the 2003 paper's formal FIML-versus-IOML information-loss
    decomposition: these twins change the sampling design at fixed nominal
    mean gap rather than integrate over unobserved observation times.

    ``excess_bias`` subtracts the equidistant bias; ``rmse_ratio`` divides by
    the equidistant RMSE. The ratio is undefined when the twin RMSE is zero.
    Coverage gaps stay missing when intervals were not computed. These are
    descriptive comparisons; no uncertainty estimate for a paired difference
    is inferred from marginal Monte Carlo SEs.
    """
    keys = ["n", "theta_mean_gap", "mean_gap", "estimator", "parameter", "truth"]
    keys.extend(column for column in ("true_theta", "true_mu", "true_sigma") if column in summary)
    values = ["cell_id", "bias", "relative_bias", "rmse", "coverage"]
    missing = set(keys + values + ["cv"]).difference(summary.columns)
    if missing:
        raise ValueError(f"summary missing required columns: {sorted(missing)}")
    if summary.duplicated(["cell_id", "estimator", "parameter"]).any():
        raise ValueError("summary must contain one row per cell/estimator/parameter")
    baseline = summary.loc[summary["cv"].eq(0), keys + values].copy()
    if baseline.duplicated(keys).any():
        raise ValueError("multiple equidistant twins match the same nominal design")
    names = {column: f"twin_{column}" for column in values}
    baseline = baseline.rename(columns=names)
    # validate='many_to_one' prevents accidental row multiplication on merge.
    result = summary.merge(baseline, on=keys, how="left", validate="many_to_one", sort=False)
    result["excess_bias"] = result["bias"] - result["twin_bias"]
    result["excess_relative_bias"] = result["relative_bias"] - result["twin_relative_bias"]
    positive_twin_rmse = result["twin_rmse"].where(result["twin_rmse"] > 0)
    result["rmse_ratio"] = result["rmse"] / positive_twin_rmse
    result["coverage_gap"] = result["coverage"] - result["twin_coverage"]
    result["comparison_status"] = np.where(result["twin_cell_id"].notna(), "matched", "missing_twin")
    return result


_PAIR_COLUMNS = [
    "cell_id", *_META, "estimator_a", "estimator_b", "parameter", "truth",
    "n_planned_pairs", "n_recorded_pairs", "n_valid_a", "n_valid_b", "n_common_valid",
    "n_excluded_pairs", "bias_a_common", "bias_b_common", "bias_difference",
    "mcse_bias_difference", "mse_a_common", "mse_b_common", "mse_difference",
    "mcse_mse_difference", "difference_direction", "pairing_scope", "validity_policy",
]


def paired_estimator_comparisons(results: pd.DataFrame) -> pd.DataFrame:
    """Compare estimators on shared paths within each scientific cell.

    For every pair A,B and parameter, restrict to paths valid for BOTH fits.
    Report mean(e_A-e_B) and mean(e_A**2-e_B**2), with sample-SD/sqrt(S)
    Monte Carlo SEs of these *paired* quantities. Pairing uses (cell_id, rep),
    never rep alone: CV cells have independent random streams. The two marginal
    biases/MSEs are recomputed on the same common-valid subset, so subtracting
    means from the main summary (which can have different denominators) is not
    substituted for this calculation. No RMSE-ratio uncertainty is invented.

    ``n_planned_pairs`` is the union of recorded replication IDs for A and B;
    ``n_recorded_pairs`` exposes unmatched rows. The runner checks complete keys
    before invoking this helper. All-invalid pairs remain as rows with missing
    metrics; one common valid path yields a difference but no estimable MCSE.
    """
    valid = _validate_results(results)
    if results.empty:
        return pd.DataFrame(columns=_PAIR_COLUMNS)
    data = results.copy()
    # Position-based assignment also handles a caller's repeated dataframe index.
    data["_point_valid"] = valid.to_numpy()
    rows = []
    for cell_id, cell in data.groupby("cell_id", sort=True):
        metadata = {column: cell[column].iloc[0] for column in _META}
        policy = validity_counts(cell)["validity_policy"]
        for a, b in combinations(sorted(cell["estimator"].unique()), 2):
            first = cell.loc[cell.estimator.eq(a)].set_index("rep")
            second = cell.loc[cell.estimator.eq(b)].set_index("rep")
            union = first.index.union(second.index)
            shared = first.index.intersection(second.index)
            common = shared[first.loc[shared, "_point_valid"].to_numpy()
                            & second.loc[shared, "_point_valid"].to_numpy()]
            count = len(common)
            for parameter in PARAMETERS:
                truth = float(metadata[f"true_{parameter}"])
                e_a = first.loc[common, parameter].to_numpy(dtype=float) - truth
                e_b = second.loc[common, parameter].to_numpy(dtype=float) - truth
                error_diff = e_a - e_b
                squared_diff = e_a**2 - e_b**2
                rows.append({
                    "cell_id": cell_id, **metadata, "estimator_a": a, "estimator_b": b,
                    "parameter": parameter, "truth": truth,
                    "n_planned_pairs": len(union), "n_recorded_pairs": len(shared),
                    "n_valid_a": int(first["_point_valid"].sum()),
                    "n_valid_b": int(second["_point_valid"].sum()),
                    "n_common_valid": count, "n_excluded_pairs": len(union) - count,
                    "bias_a_common": float(e_a.mean()) if count else np.nan,
                    "bias_b_common": float(e_b.mean()) if count else np.nan,
                    "bias_difference": float(error_diff.mean()) if count else np.nan,
                    "mcse_bias_difference": (float(error_diff.std(ddof=1) / np.sqrt(count))
                                             if count >= 2 else np.nan),
                    "mse_a_common": float(np.mean(e_a**2)) if count else np.nan,
                    "mse_b_common": float(np.mean(e_b**2)) if count else np.nan,
                    "mse_difference": float(squared_diff.mean()) if count else np.nan,
                    "mcse_mse_difference": (float(squared_diff.std(ddof=1) / np.sqrt(count))
                                            if count >= 2 else np.nan),
                    "difference_direction": "estimator_a_minus_estimator_b",
                    "pairing_scope": "common_valid_shared_paths_within_cell",
                    "validity_policy": policy,
                })
    return pd.DataFrame(rows, columns=_PAIR_COLUMNS)
