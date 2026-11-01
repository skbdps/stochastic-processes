"""Additive calibration checks using only the sealed Week 10 saved records.

For fixed weights and a known exogenous schedule t, the linear diagnostic is
exactly N(mu, v(t)). Its coverage with the known width z*sqrt(w(t)) is exactly
2*Phi(z*sqrt(w(t)/v(t)))-1. Averaging these probabilities integrates a Gaussian
mixture; it does not require a marginal normal approximation. Transfer to the
fitted mean or its fitted Hessian width remains a separate approximation.

The schedules are the same ones used in the observed errors. These response-
free conditional benchmarks are not independent validation samples. Paired
coverage differences retain this dependence. No sources or fits are rerun.
"""
from pathlib import Path
import gzip
import json
import math

import numpy as np
import pandas as pd
from scipy.special import ndtr
from scipy.stats import beta, norm

from ..experiments.analysis import csv_bytes
from ..experiments.provenance import file_digest
from ..experiments.storage import atomic_bytes, atomic_json


Z95 = float(norm.ppf(.975))
INPUTS = ("analysis/schedule_predictions.csv.gz", "analysis/parameter_summary.csv",
          "analysis/original_e4_mechanism.csv")


def known_width_coverage(variance, working_variance):
    """Exact conditional coverage for the fixed-coefficient Gaussian diagnostic.

Both variances must be positive and fixed given the schedule. Using a fitted
response-dependent width here would not give its conditional coverage.
    """
    variance, working = np.broadcast_arrays(np.asarray(variance, dtype=float),
                                            np.asarray(working_variance, dtype=float))
    if (not np.isfinite(variance).all() or not np.isfinite(working).all()
            or np.any(variance <= 0) or np.any(working <= 0)):
        raise ValueError("diagnostic and working variances must be positive and finite")
    return 2 * ndtr(Z95 * np.sqrt(working / variance)) - 1


def paired_difference(first, second):
    """Paired coverage difference with a conservative finite-sample interval.

The discordances have probabilities p10 and p01, whose difference is the
coverage difference. Each marginal count is binomial. Two-sided 97.5%
Clopper-Pearson bounds on each probability hold simultaneously with at least
95% confidence by the union bound, even though the counts are dependent.
Subtracting their endpoints therefore gives a conservative pointwise 95%
interval. The plug-in MCSE and normal bounds are descriptive approximations;
zero observed discordances do not imply zero population uncertainty.
    """
    first, second = np.asarray(first, dtype=float), np.asarray(second, dtype=float)
    if (first.ndim != 1 or first.shape != second.shape or len(first) < 2
            or not np.isin(first, [0., 1.]).all() or not np.isin(second, [0., 1.]).all()):
        raise ValueError("paired coverage needs at least two matching binary indicators")
    difference = first - second
    mean = float(difference.mean())
    se = float(difference.std(ddof=1) / math.sqrt(len(difference)))
    count10, count01 = int(np.sum(difference == 1)), int(np.sum(difference == -1))
    n = len(difference)
    tail = .0125
    def bounds(count):
        lower = 0. if count == 0 else float(beta.ppf(tail, count, n-count+1))
        upper = 1. if count == n else float(beta.ppf(1-tail, count+1, n-count))
        return lower, upper
    lower10, upper10 = bounds(count10)
    lower01, upper01 = bounds(count01)
    return {"coverage_difference": mean, "paired_mcse": se,
            "pointwise_mc_lower": lower10-upper01,
            "pointwise_mc_upper": upper10-lower01,
            "normal_approx_mc_lower": mean-Z95*se,
            "normal_approx_mc_upper": mean+Z95*se,
            "first_only_covered": count10, "second_only_covered": count01,
            "p10_cp_lower": lower10, "p10_cp_upper": upper10,
            "p01_cp_lower": lower01, "p01_cp_upper": upper01}


def _mcse(values):
    return float(np.std(values, ddof=1) / math.sqrt(len(values)))


def diagnostic_tables(paths, summary):
    """Calculate conditional calibration and paired center/width diagnostics.

This review deliberately requires the complete point/interval availability
of the sealed study. It fails on unavailable fits instead of silently changing
denominators. A uses fixed-rho centers and population widths; B fitted centers
and population widths; C fitted centers and actual reported Hessian widths.
    """
    required = {"cell_id", "rep", "mechanism", "ratio", "fixed_rho_mu", "pfml_mu",
                "pfml_mu_se", "fixed_rho_conditional_variance", "population_working_variance",
                "point_valid", "interval_valid"}
    if not required.issubset(paths.columns) or paths.empty:
        raise ValueError("missing saved path diagnostics")
    if paths.duplicated(["cell_id", "rep"]).any():
        raise ValueError("duplicate source diagnostic")
    for field in ("point_valid", "interval_valid"):
        if not pd.api.types.is_bool_dtype(paths[field]) or not paths[field].all():
            raise ValueError("review requires complete point and interval availability")
    fitted_summary = summary[(summary.estimator == "pfml") & (summary.parameter == "mu")]
    if (fitted_summary.cell_id.duplicated().any()
            or set(fitted_summary.cell_id) != set(paths.cell_id)):
        raise ValueError("missing or duplicate fitted mean summary")
    cells, paired = [], []
    for cid, group in paths.groupby("cell_id", sort=True):
        fitted = fitted_summary[fitted_summary.cell_id == cid].iloc[0]
        n = len(group)
        if n < 2 or any(int(fitted[x]) != n for x in
                        ("n_planned", "n_point_valid", "n_interval_valid")):
            raise ValueError("source count or availability differs from summary")
        if group.mechanism.nunique() != 1 or group.ratio.nunique() != 1:
            raise ValueError("source cell mixes mechanisms or ratios")
        truth = float(fitted.truth)
        v = group.fixed_rho_conditional_variance.to_numpy()
        w = group.population_working_variance.to_numpy()
        probability = known_width_coverage(v, w)
        fixed_error = group.fixed_rho_mu.to_numpy() - truth
        fitted_error = group.pfml_mu.to_numpy() - truth
        fitted_se = group.pfml_mu_se.to_numpy()
        if (not np.isfinite(truth) or not np.isfinite(fixed_error).all()
                or not np.isfinite(fitted_error).all() or not np.isfinite(fitted_se).all()
                or np.any(fitted_se <= 0)):
            raise ValueError("nonfinite centers or invalid fitted widths")
        if ("normal_working_coverage_approximation" in group and not np.allclose(
                group.normal_working_coverage_approximation, probability, rtol=1e-13, atol=1e-15)):
            raise ValueError("stored conditional coverage calculation changed")
        coverage = {"A": np.abs(fixed_error) <= Z95 * np.sqrt(w),
                    "B": np.abs(fitted_error) <= Z95 * np.sqrt(w),
                    "C": np.abs(fitted_error) <= Z95 * fitted_se}
        if (int(coverage["C"].sum()) != int(fitted.n_covered)
                or not math.isclose(float(coverage["C"].mean()), fitted.coverage_operational,
                                    rel_tol=0, abs_tol=1e-14)):
            raise ValueError("actual fitted interval coverage differs from summary")
        observed = {"observed_pfml_bias": float(fitted_error.mean()),
                    "observed_pfml_sd": float(fitted_error.std(ddof=1)),
                    "observed_pfml_rmse": float(np.sqrt(np.mean(fitted_error**2)))}
        for name, field in (("observed_pfml_bias", "bias"),
                            ("observed_pfml_sd", "empirical_sd_point_valid"),
                            ("observed_pfml_rmse", "rmse")):
            if not math.isclose(observed[name], float(fitted[field]), rel_tol=1e-12, abs_tol=1e-14):
                raise ValueError("fitted mean metric differs from summary: " + field)
        mean_v, mean_w, mean_c = float(v.mean()), float(w.mean()), float(probability.mean())
        variance_se, coverage_se = _mcse(v), _mcse(probability)
        sd = math.sqrt(mean_v)
        log_sd_se = variance_se / (2 * mean_v)
        identity = {"cell_id": cid, "mechanism": group.mechanism.iloc[0],
                    "ratio": float(group.ratio.iloc[0]), "n_paths": n, "truth": truth}
        row = {**identity, "predicted_variance": mean_v, "predicted_variance_mcse": variance_se,
               "predicted_sd": sd, "predicted_sd_mcse": variance_se / (2 * sd),
               "predicted_sd_mc_lower": sd * math.exp(-Z95 * log_sd_se),
               "predicted_sd_mc_upper": sd * math.exp(Z95 * log_sd_se),
               "exact_mixture_known_width_coverage_estimate": mean_c,
               "mixture_coverage_schedule_mcse": coverage_se,
               "mixture_coverage_mc_lower": max(0., mean_c - Z95 * coverage_se),
               "mixture_coverage_mc_upper": min(1., mean_c + Z95 * coverage_se),
               "marginal_normal_rms_width_coverage_approximation": float(known_width_coverage(mean_v, mean_w)),
               "population_rms_working_se": math.sqrt(mean_w),
               "fixed_rho_rmse": float(np.sqrt(np.mean(fixed_error**2))),
               "fixed_rho_bias": float(fixed_error.mean()),
               "fixed_rho_sd": float(fixed_error.std(ddof=1)),
               "fitted_minus_fixed_mu_rms": float(np.sqrt(np.mean((fitted_error-fixed_error)**2))),
               "A_coverage_minus_exact_conditional_expectation": float(coverage["A"].mean()-mean_c),
               "A_coverage_residual_conditional_mcse": float(np.sqrt(np.sum(probability*(1-probability))) / n),
               **observed}
        for label, indicators in coverage.items():
            row[label + "_covered"] = int(indicators.sum())
            row[label + "_coverage"] = float(indicators.mean())
        cells.append(row)
        for first, second in (("B", "A"), ("C", "B"), ("C", "A")):
            paired.append({**identity, "contrast": first + "_minus_" + second,
                           **paired_difference(coverage[first], coverage[second])})
    cell_table = pd.DataFrame(cells)
    ratios = []
    for a in cell_table[cell_table.mechanism == "mmpp"].itertuples(index=False):
        control = cell_table[(cell_table.mechanism == "iid") & (cell_table.ratio == a.ratio)]
        if len(control) != 1:
            raise ValueError("missing unique IID reference")
        b = control.iloc[0]
        ratio = a.predicted_sd / b.predicted_sd
        logse = .5 * math.hypot(a.predicted_variance_mcse / a.predicted_variance,
                                b.predicted_variance_mcse / b.predicted_variance)
        ratios.append({"ratio": a.ratio, "mmpp_cell_id": a.cell_id, "iid_cell_id": b.cell_id,
                       "predicted_sd_ratio": ratio, "predicted_sd_ratio_log_mcse": logse,
                       "predicted_sd_ratio_mc_lower": ratio * math.exp(-Z95 * logse),
                       "predicted_sd_ratio_mc_upper": ratio * math.exp(Z95 * logse),
                       "observed_pfml_sd_ratio": a.observed_pfml_sd / b.observed_pfml_sd,
                       "observed_pfml_rmse_ratio": a.observed_pfml_rmse / b.observed_pfml_rmse,
                       "mixture_coverage_difference": a.exact_mixture_known_width_coverage_estimate - b.exact_mixture_known_width_coverage_estimate,
                       "mixture_coverage_difference_mcse": math.hypot(a.mixture_coverage_schedule_mcse, b.mixture_coverage_schedule_mcse)})
    return {"cell_calibration.csv": cell_table,
            "matched_prediction_comparisons.csv": pd.DataFrame(ratios).sort_values("ratio").reset_index(drop=True),
            "paired_coverage.csv": pd.DataFrame(paired)}


def original_comparisons(original):
    """Keep fixed-count exact diagnostic SD distinct from fitted SD and RMSE."""
    rows = []
    for a in original[original.mechanism == "mmpp"].itertuples(index=False):
        selected = original[(original.mechanism == "iid") & (original.ratio == a.ratio)]
        if len(selected) != 1:
            raise ValueError("missing original IID reference")
        b = selected.iloc[0]
        rows.append({"ratio": a.ratio, "predicted_sd_ratio": a.fixed_rho_mean_sd / b.fixed_rho_mean_sd,
                     "observed_pfml_sd_ratio": a.observed_sd / b.observed_sd,
                     "observed_pfml_rmse_ratio": a.observed_rmse / b.observed_rmse,
                     "mmpp_observed_bias": a.observed_bias, "iid_observed_bias": b.observed_bias,
                     "prediction_is_analytic_not_schedule_monte_carlo": True})
    return pd.DataFrame(rows).sort_values("ratio").reset_index(drop=True)


def _sealed_inputs(package):
    """Verify the entire existing seal without changing its files or inventory."""
    package = Path(package)
    seal_file = package / "completion_seal.json"
    seal = json.loads(seal_file.read_text())
    files = seal.get("files", {})
    if seal.get("status") != "completed" or not files:
        raise ValueError("original supplement is not sealed complete")
    # Match the original seal policy: these two advisory locks were never
    # research payload. Other additions, including nested seal names, fail.
    excluded = {"completion_seal.json", "horizon_full/.execution.lock",
                "horizon_validation/.execution.lock"}
    actual = set()
    for path in package.rglob("*"):
        if path.is_symlink():
            raise ValueError("original sealed package contains a symbolic link")
        relative = str(path.relative_to(package))
        if path.is_file() and relative not in excluded:
            actual.add(relative)
    if actual != set(files):
        raise ValueError("original sealed inventory changed")
    for name, digest in files.items():
        if Path(name).is_absolute() or ".." in Path(name).parts or file_digest(package/name) != digest:
            raise ValueError("original sealed bytes changed: " + name)
    return {"seal_sha256": file_digest(seal_file), "verified_files": len(files)}


def _raw_coverage(package, paths):
    rows, partitions, total = [], [], 0
    for path in sorted((package/"horizon_full/cells").glob("*/part-*.jsonl.gz")):
        partitions.append(str(path.relative_to(package)))
        with gzip.open(path, "rt", encoding="utf-8") as stream:
            for line in stream:
                if line.strip():
                    record = json.loads(line)
                    total += 1
                    if record["estimator"] == "pfml":
                        rows.append(record)
    raw = pd.DataFrame(rows)
    if raw.empty or raw.duplicated(["cell_id", "rep"]).any() or len(raw) != len(paths):
        raise ValueError("raw PFML sources missing or duplicated")
    merged = paths.merge(raw, on=["cell_id", "rep"], suffixes=("_saved", "_raw"),
                         validate="one_to_one", how="outer", indicator=True)
    if not merged._merge.eq("both").all():
        raise ValueError("raw and diagnostic source identities differ")
    for saved, direct in (("pfml_mu", "mu"), ("pfml_mu_se", "mu_se"),
                          ("point_valid_saved", "point_valid_raw"),
                          ("interval_valid_saved", "interval_valid_raw"),
                          ("times_sha256_saved", "times_sha256_raw"),
                          ("observations_sha256_saved", "observations_sha256_raw")):
        if not np.array_equal(merged[saved], merged[direct]):
            raise ValueError("raw fitted diagnostic mismatch: " + saved)
    actual = (merged.mu_lower <= merged.true_mu) & (merged.true_mu <= merged.mu_upper)
    derived = np.abs(merged.pfml_mu-merged.true_mu) <= Z95*merged.pfml_mu_se
    if not np.array_equal(actual, derived):
        raise ValueError("reconstructed C indicators differ from stored full-Hessian intervals")
    counts = actual.groupby(merged.cell_id).sum()
    return {"raw_records_read": total, "pfml_sources_checked": len(raw),
            "every_C_indicator_matches_raw_interval": True,
            "raw_C_covered_by_cell": {str(cid): int(value) for cid, value in counts.items()}}, partitions


def _calculate(package):
    package = Path(package)
    seal = _sealed_inputs(package)
    paths, summary, original = [pd.read_csv(package/name, float_precision="round_trip") for name in INPUTS]
    tables = diagnostic_tables(paths, summary)
    tables["original_sd_comparisons.csv"] = original_comparisons(original)
    raw_check, partitions = _raw_coverage(package, paths)
    for row in tables["cell_calibration.csv"].itertuples(index=False):
        if row.C_covered != raw_check["raw_C_covered_by_cell"][row.cell_id]:
            raise ValueError("review C count differs from raw interval count")
    inputs = {name: file_digest(package/name) for name in (*INPUTS, *partitions)}
    return tables, {"original_seal": seal, "input_files": inputs, "raw_readback": raw_check}


def build(package, output):
    """Write a new review directory; existing evidence is never overwritten."""
    package, output = Path(package), Path(output)
    if output.exists() or output.resolve().is_relative_to(package.resolve()):
        raise FileExistsError("review output must be new and outside the sealed supplement")
    tables, evidence = _calculate(package)
    output.mkdir(parents=True)
    for name, table in tables.items():
        atomic_bytes(output/name, csv_bytes(table), immutable=True)
    manifest = {"schema_version": 1, "scope": "post_result_analytical_critique_saved_data_only",
                "uncertainty": "pointwise Monte Carlo; source replication is the independent unit",
                "paired_interval_method": "conservative >=95%: subtract simultaneous 97.5% two-sided Clopper-Pearson bounds on the two discordance probabilities; Bonferroni union bound",
                "paired_normal_approximation": "paired_mcse and normal_approx_mc_lower/upper are descriptive plug-in approximations, which can degenerate when no discordances are observed",
                "prediction_interval_method": "approximate pointwise 95% log-delta Monte Carlo intervals for finite-variance fixed-coefficient SDs and independent-arm SD ratios",
                "same_schedules": "conditional benchmarks share schedules with observed errors; no independent validation claim",
                "coverage_diagnostics": {"A": "fixed-rho center and known population width",
                                         "B": "fitted PFML center and known population width",
                                         "C": "fitted PFML center and actual reported Hessian width"},
                "marginal_normal_comparator": "uses mean variance and RMS population width; replaces both the mixture and random widths",
                "source_module_sha256": file_digest(Path(__file__)), **evidence,
                "files": {name: file_digest(output/name) for name in tables}}
    if _sealed_inputs(package) != evidence["original_seal"]:
        raise ValueError("original seal changed while producing review")
    atomic_json(output/"manifest.json", manifest, immutable=True)
    return manifest


def verify(package, output):
    """Recompute every review table and recheck original sealed inputs."""
    output = Path(output)
    manifest = json.loads((output/"manifest.json").read_text())
    if file_digest(Path(__file__)) != manifest["source_module_sha256"]:
        raise ValueError("review source module changed")
    tables, evidence = _calculate(package)
    for key, value in evidence.items():
        if manifest[key] != value:
            raise ValueError("bound review evidence changed: " + key)
    if set(p.name for p in output.iterdir()) != set(tables) | {"manifest.json"}:
        raise ValueError("review output inventory changed")
    for name, table in tables.items():
        if (output/name).read_bytes() != csv_bytes(table) or file_digest(output/name) != manifest["files"].get(name):
            raise ValueError("review arithmetic or digest changed: " + name)
    return {"status": "passed", "rederived_tables": len(tables), **evidence["original_seal"],
            **evidence["raw_readback"]}
