"""Monte Carlo summaries and equidistant controls for the Week 3 harness.

The runner must emit a row for *every* planned replication/estimator, including
failed fits. Bias and RMSE below are conditional on a successful finite fit;
the accompanying failure rate always uses all emitted/planned rows. Neither
discarding a failed row nor replacing its estimate with the truth is valid.

Week 3 does not construct confidence intervals. Its coverage columns therefore
remain missing unless callers explicitly provide actual interval endpoints.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


PARAMETERS = ("theta", "mu", "sigma")
_META = ("theta_mean_gap", "cv", "n", "mean_gap", "true_theta", "true_mu", "true_sigma")
_SUMMARY_COLUMNS = [
    "cell_id", *_META, "estimator", "parameter", "truth", "n_total", "n_success",
    "n_failed", "failure_rate", "bias", "relative_bias", "rmse", "relative_rmse",
    "mcse_bias", "mcse_rmse", "relative_mcse_bias", "relative_mcse_rmse",
    "coverage", "coverage_mcse", "n_ci_valid", "ci_status",
    "metric_scope", "coverage_scope",
]


def summarize_results(results: pd.DataFrame) -> pd.DataFrame:
    """Return one row per cell, estimator, and parameter.

    For S successful fits, bias = sum(estimate - truth)/S and
    RMSE = sqrt(sum((estimate - truth)**2)/S). The Monte Carlo standard error
    of estimated bias is the sample SD of estimates divided by sqrt(S), using
    S-1 in the SD, and is undefined when S < 2. Applying the delta method to
    sqrt(mean(squared errors)) gives MCSE(RMSE) = sample_SD(squared errors) /
    (2*sqrt(S)*RMSE), undefined for S < 2 or RMSE = 0. This approximation need
    not be accurate in a tiny smoke sample. These are simulation-error measures,
    not confidence intervals for individual OU parameter estimates.

    A fit is usable only if ``success`` is True and *all three* parameter
    estimates are finite. This keeps failure counts identical across parameter
    rows. Relative quantities divide by abs(truth) and are undefined at zero;
    in particular mu=0 must be assessed through absolute bias and RMSE.

    Optional ``{parameter}_lower`` / ``{parameter}_upper`` endpoints are used
    only for successful finite fits with finite, correctly ordered endpoints.
    Coverage uses n_ci_valid as its denominator (not n_total or n_success),
    includes the interval endpoints, and reports binomial plug-in MCSE. Missing
    pairs leave coverage explicitly uncomputed; a half-supplied pair is an
    input error. No interval estimator is implemented here.
    """
    required = {"cell_id", *_META, "rep", "estimator", "success", *PARAMETERS}
    missing = required.difference(results.columns)
    if missing:
        raise ValueError(f"results missing required columns: {sorted(missing)}")
    interval_parameters = set()
    for parameter in PARAMETERS:
        lower, upper = f"{parameter}_lower", f"{parameter}_upper"
        if (lower in results) != (upper in results):
            raise ValueError(f"supply both {lower} and {upper}, or neither")
        if lower in results:
            interval_parameters.add(parameter)
    if results.empty:
        return pd.DataFrame(columns=_SUMMARY_COLUMNS)
    if results[["cell_id", "estimator", "rep", *_META]].isna().any().any():
        raise ValueError("cell metadata, estimator, and replication IDs cannot be missing")
    if not results["success"].map(lambda value: isinstance(value, (bool, np.bool_))).all():
        raise ValueError("success must contain non-missing Boolean values")
    if results.duplicated(["cell_id", "estimator", "rep"]).any():
        raise ValueError("duplicate cell_id/estimator/rep rows would double-count replications")
    for column in _META:
        if not np.isfinite(results[column].to_numpy(dtype=float)).all():
            raise ValueError(f"metadata {column} must be finite")
        if results.groupby("cell_id")[column].nunique().gt(1).any():
            raise ValueError(f"metadata {column} is inconsistent within a cell_id")

    rows = []
    for (cell_id, estimator), group in results.groupby(["cell_id", "estimator"], sort=True):
        finite_fit = np.isfinite(group[list(PARAMETERS)].to_numpy(dtype=float)).all(axis=1)
        usable = group["success"].to_numpy(dtype=bool) & finite_fit
        count = int(usable.sum())
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
                coverage_scope = "conditional_on_success_and_valid_interval"
                if n_ci_valid:
                    covered = (lower[valid] <= truth) & (truth <= upper[valid])
                    coverage = float(covered.mean())
                    coverage_mcse = float(np.sqrt(coverage * (1.0 - coverage) / n_ci_valid))
            rows.append({
                "cell_id": cell_id, **metadata, "estimator": estimator,
                "parameter": parameter, "truth": truth,
                "n_total": n_total, "n_success": count, "n_failed": n_total - count,
                "failure_rate": (n_total - count) / n_total,
                "bias": bias, "relative_bias": relative_bias, "rmse": rmse,
                "relative_rmse": relative_rmse, "mcse_bias": mcse_bias,
                "mcse_rmse": mcse_rmse, "relative_mcse_bias": relative_mcse_bias,
                "relative_mcse_rmse": relative_mcse_rmse,
                "coverage": coverage, "coverage_mcse": coverage_mcse,
                "n_ci_valid": n_ci_valid, "ci_status": ci_status,
                "metric_scope": "conditional_on_success", "coverage_scope": coverage_scope,
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
