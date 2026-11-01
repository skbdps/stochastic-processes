"""Saved-result reporting and quantitative covariance checks for Week 10."""
from pathlib import Path
import gzip
import json
import math

import numpy as np
import pandas as pd
from scipy.special import ndtr

from ..experiments.analysis import csv_bytes
from ..experiments.comparisons import _interval, _ratio
from ..experiments.provenance import file_digest
from ..experiments.storage import atomic_bytes, atomic_json
from ..experiments.summaries import ExperimentSummarizer
from ..week8.design import E4Cell, full_design as e4_design
from .execution import load_records, readback
from .source import HorizonSource
from .theory import (conditional_mean_variance, conditional_pfml_variance,
                     pfml_mean_weights, population_predictions, population_rho)


def summaries(records, design):
    """Use original denominator/uncertainty rules, explicitly labeling random N."""
    data = pd.DataFrame(records)
    result = ExperimentSummarizer().summarize(data, design)
    cells = {c.identity: c for c in design.cells}
    result["mechanism"] = result.cell_id.map(lambda cid: cells[cid].mechanism)
    result["ratio"] = result.cell_id.map(lambda cid: cells[cid].ratio)
    result["horizon"] = result.cell_id.map(lambda cid: cells[cid].horizon)
    return result.rename(columns={"n": "nominal_n"})


def contrasts(summary):
    """Independent MMPP minus IID contrasts, with each arm's own eligibility."""
    rows = []
    for a in summary[summary.mechanism == "mmpp"].itertuples(index=False):
        selected = summary[(summary.mechanism == "iid") & (summary.ratio == a.ratio)
                           & (summary.estimator == a.estimator) & (summary.parameter == a.parameter)]
        if len(selected) != 1:
            raise ValueError("missing unique matched IID summary")
        b = selected.iloc[0]
        row = {"ratio": a.ratio, "estimator": a.estimator, "parameter": a.parameter,
               "cell_id": a.cell_id, "reference_cell_id": b.cell_id,
               "comparison": "independent_cells_equal_T_random_N"}
        for field in ("n_planned", "n_point_valid", "n_interval_valid", "n_covered",
                      "bias", "rmse", "coverage_operational", "coverage_conditional"):
            row["mmpp_" + field], row["iid_" + field] = getattr(a, field), b[field]
        for name, value, se in (("bias_difference", "bias", "mcse_bias"),
                                ("coverage_difference", "coverage_operational", "operational_mcse")):
            _interval(row, name, getattr(a, value) - b[value], math.hypot(getattr(a, se), b[se]))
        logse = (math.hypot(a.mcse_rmse / a.rmse, b.mcse_rmse / b.rmse)
                 if a.rmse > 0 and b.rmse > 0 else float("nan"))
        _ratio(row, a.rmse, b.rmse, logse)
        rows.append(row)
    return pd.DataFrame(rows)


def source_summary(records, design):
    data = pd.DataFrame(records).drop_duplicates(["cell_id", "rep"])
    rows = []
    for cell in design.cells:
        group = data[data.cell_id == cell.identity]
        row = {"cell_id": cell.identity, "mechanism": cell.mechanism,
               "ratio": cell.ratio, "n_paths": len(group), "nominal_n": cell.n,
               "fixed_horizon": cell.horizon}
        for field in ("actual_n", "realized_span", "realized_mean_gap", "terminal_gap",
                      "clip_count", "min_gap"):
            if field not in group:
                continue
            values = pd.to_numeric(group[field], errors="raise").dropna().to_numpy()
            row.update({field + "_available": len(values),
                        field + "_mean": float(np.mean(values)) if len(values) else np.nan,
                        field + "_sd": float(np.std(values, ddof=1)) if len(values) > 1 else np.nan,
                        field + "_min": float(np.min(values)) if len(values) else np.nan,
                        field + "_max": float(np.max(values)) if len(values) else np.nan})
        rows.append(row)
    return pd.DataFrame(rows)


def population_table(e4_summary):
    rows = []
    for cell in e4_design().cells:
        prediction = population_predictions(cell)
        subset = e4_summary[(e4_summary.cell_id == cell.identity)
                            & (e4_summary.estimator == "pfml") & (e4_summary.parameter == "mu")]
        if len(subset) != 1:
            raise ValueError("original E4 PFML mu summary key is missing")
        observed = subset.iloc[0]
        rows.append({"cell_id": cell.identity, "mechanism": cell.mechanism,
                     "ratio": cell.ratio, **prediction,
                     "observed_bias": observed.bias, "observed_rmse": observed.rmse,
                     "observed_sd": observed.empirical_sd_point_valid,
                     "observed_rms_se": observed.rms_reported_se,
                     "observed_coverage": observed.coverage_operational})
    return pd.DataFrame(rows)


def schedule_predictions(records, design):
    """Reconstruct each new source for exact fixed-rho conditional diagnostics.

    No model is refitted. Original source array hashes must agree. The working
    Gaussian coverage uses population working variance, not fitted intervals.
    """
    data = pd.DataFrame(records)
    source = HorizonSource(design.seed)
    rows = []
    for cell in design.cells:
        law = E4Cell(n=cell.n, q=cell.q, theta=cell.theta, mu=cell.mu, sigma=cell.sigma,
                     mechanism=cell.mechanism, ratio=cell.ratio,
                     floor_fraction=cell.floor_fraction, floor_cv_threshold=cell.floor_cv_threshold)
        rho = population_rho(law)
        saved = data[(data.cell_id == cell.identity) & (data.estimator == "pfml")].set_index("rep")
        for rep in range(design.replications):
            sample, diag = source.sample(cell.to_dict(), rep)
            original = saved.loc[rep]
            for name in ("times_sha256", "observations_sha256"):
                if diag[name] != original[name]:
                    raise ValueError("source changed during mechanism reconstruction: " + name)
            n = len(sample.times)
            variance = conditional_pfml_variance(sample.times, rho, cell.theta, cell.sigma)
            working = cell.sigma**2 / (2*cell.theta) * (1+rho) / ((n-1)*(1-rho))
            linear_mu = float(pfml_mean_weights(n, rho) @ sample.values)
            rows.append({"cell_id": cell.identity, "rep": rep, "mechanism": cell.mechanism,
                         "ratio": cell.ratio, "actual_n": n, "horizon": cell.horizon,
                         "population_rho": rho, "fixed_rho_conditional_variance": variance,
                         "ordinary_mean_conditional_variance": conditional_mean_variance(
                             sample.times, cell.theta, cell.sigma),
                         "population_working_variance": working,
                         "normal_working_coverage_approximation": float(2*ndtr(1.959963984540054
                                                                             * math.sqrt(working/variance))-1),
                         "fixed_rho_mu": linear_mu,
                         "pfml_mu": original.mu, "pfml_mu_se": original.mu_se,
                         "point_valid": bool(original.point_valid),
                         "interval_valid": bool(original.interval_valid),
                         "times_sha256": diag["times_sha256"],
                         "observations_sha256": diag["observations_sha256"]})
    return pd.DataFrame(rows)


def mechanism_summary(paths, summary):
    rows = []
    for (cid, mechanism, ratio), group in paths.groupby(["cell_id", "mechanism", "ratio"], sort=True):
        fitted = summary[(summary.cell_id == cid) & (summary.estimator == "pfml")
                         & (summary.parameter == "mu")].iloc[0]
        variances = group.fixed_rho_conditional_variance.to_numpy()
        squares = group.fixed_rho_mu.to_numpy()**2
        difference = squares - variances
        rows.append({"cell_id": cid, "mechanism": mechanism, "ratio": ratio, "n_paths": len(group),
                     "predicted_fixed_rho_sd": math.sqrt(float(variances.mean())),
                     "predicted_variance_mcse": float(variances.std(ddof=1)/math.sqrt(len(group))),
                     "empirical_fixed_rho_rmse": math.sqrt(float(squares.mean())),
                     "fixed_rho_mse_minus_predicted_variance": float(difference.mean()),
                     "fixed_rho_mse_minus_predicted_variance_mcse": float(difference.std(ddof=1)/math.sqrt(len(group))),
                     "observed_pfml_rmse": fitted.rmse, "observed_pfml_sd": fitted.empirical_sd_point_valid,
                     "observed_pfml_rms_se": fitted.rms_reported_se,
                     "observed_pfml_coverage": fitted.coverage_operational,
                     "normal_working_coverage_approximation": float(group.normal_working_coverage_approximation.mean()),
                     "population_rms_working_se": math.sqrt(float(group.population_working_variance.mean())),
                     "fitted_minus_fixed_mu_rms": math.sqrt(float(np.mean(
                         (group.loc[group.point_valid, "pfml_mu"]-group.loc[group.point_valid, "fixed_rho_mu"])**2)))})
    return pd.DataFrame(rows)


def build(design, run_dir, e4_dir, output):
    output, run_dir, e4_dir = Path(output), Path(run_dir), Path(e4_dir)
    if output.exists():
        raise FileExistsError("report output must be new; preserve previous attempts")
    readback(design, run_dir)
    records = load_records(run_dir, design)
    summary = summaries(records, design)
    e4_summary = pd.read_csv(e4_dir / "parameter_summary.csv")
    paths = schedule_predictions(records, design)
    tables = {"parameter_summary.csv": summary, "matched_iid_comparisons.csv": contrasts(summary),
              "source_summary.csv": source_summary(records, design),
              "original_e4_mechanism.csv": population_table(e4_summary),
              "horizon_mechanism.csv": mechanism_summary(paths, summary)}
    raw = pd.DataFrame(records)
    failures = raw[(~raw.point_valid) | (~raw.interval_valid) | raw.pipeline_error.fillna("").ne("")]
    output.mkdir(parents=True)
    for name, table in tables.items():
        atomic_bytes(output/name, csv_bytes(table), immutable=True)
    for name, table in (("schedule_predictions.csv.gz", paths), ("failed_fits.csv.gz", failures)):
        atomic_bytes(output/name, gzip.compress(csv_bytes(table), mtime=0), immutable=True)
    manifest = {"schema_version": 1, "scope": "post_result_supplement",
                "design": design.to_dict(), "raw_record_count": len(raw),
                "source_paths": len(paths), "failed_record_count": len(failures),
                "source_inputs": {
                    "horizon": {str(p.relative_to(run_dir)): file_digest(p) for p in
                        sorted(run_dir.rglob("*.jsonl.gz")) + [run_dir/"resolved_design.json",
                        run_dir/"run_metadata.json"]},
                    "original_e4": {"parameter_summary.csv": file_digest(e4_dir/"parameter_summary.csv")}},
                "files": {p.name: file_digest(p) for p in sorted(output.iterdir())}}
    atomic_json(output/"table_manifest.json", manifest, immutable=True)
    return manifest


def verify_tables(design, run_dir, output, e4_dir, *, replay_mechanism=False):
    """Recompute every fitted-parameter summary and contrast from durable records."""
    output = Path(output)
    manifest = json.loads((output/"table_manifest.json").read_text())
    for name, digest in manifest["files"].items():
        if file_digest(output/name) != digest:
            raise ValueError("report file digest mismatch: " + name)
    for role, root in (("horizon", Path(run_dir)), ("original_e4", Path(e4_dir))):
        for name, digest in manifest["source_inputs"][role].items():
            if file_digest(root/name) != digest:
                raise ValueError("source input changed: " + name)
    records = load_records(run_dir, design)
    summary = summaries(records, design)
    for name, table in (("parameter_summary.csv", summary),
                        ("matched_iid_comparisons.csv", contrasts(summary)),
                        ("source_summary.csv", source_summary(records, design))):
        if (output/name).read_bytes() != csv_bytes(table):
            raise ValueError("raw-to-summary arithmetic changed: " + name)
    paths = pd.read_csv(output/"schedule_predictions.csv.gz", float_precision="round_trip")
    for name, table in (("original_e4_mechanism.csv", population_table(pd.read_csv(
            Path(e4_dir)/"parameter_summary.csv"))),
                        ("horizon_mechanism.csv", mechanism_summary(paths, summary))):
        if (output/name).read_bytes() != csv_bytes(table):
            raise ValueError("mechanism summary arithmetic changed: " + name)
    if replay_mechanism:
        regenerated = schedule_predictions(records, design)
        if gzip.decompress((output/"schedule_predictions.csv.gz").read_bytes()) != csv_bytes(regenerated):
            raise ValueError("schedule-based mechanism predictions changed")
    return {"status": "passed", "raw_records": len(records), "rederived_tables": 5,
            "mechanism_sources_reconstructed": len(paths) if replay_mechanism else 0,
            "hashed_files": len(manifest["files"])}
