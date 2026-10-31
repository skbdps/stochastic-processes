"""Publication views of the frozen E1--E4 evidence.

These functions read accepted saved outcomes. They neither fit a model nor
change a scientific eligibility rule. Classifications are checked against the
unrounded registered bounds before they enter T1 or F7. A study identifier is
part of every cross-study key, including when scientific cell IDs coincide.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd


SUMMARY_STUDIES = ("e1", "e2", "e3_followup", "e3_floor", "e4")
SCOPES = {"e1": "registered", "e2": "registered", "e3_followup": "exploratory_followup",
          "e3_floor": "exploratory_paired_floor", "e4": "exploratory_matched_marginal"}
METHODS = ("exact", "pfml", "euler")
PARAMETERS = ("theta", "mu", "sigma")
CLASSES = ("within_criteria", "supported_breach", "unresolved", "not_applicable_zero_truth")
KEYS = ["cell_id", "estimator", "parameter"]
CLASS_COLUMNS = ["operational_wilson_lower", "operational_wilson_upper",
                 "relative_bias_mc_lower", "relative_bias_mc_upper"]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def _plain(value):
    if isinstance(value, dict):
        return {str(k): _plain(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [_plain(v) for v in value]
    if isinstance(value, np.generic):
        return _plain(value.item())
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def _json(value):
    return json.dumps(_plain(value), sort_keys=True, separators=(",", ":"), allow_nan=False)


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def select_one(frame, **keys):
    """Select exactly one scientific row; ambiguous editorial selections fail."""
    selected = frame
    for key, value in keys.items():
        require(key in selected, "unknown publication selector: " + key)
        selected = selected[selected[key].eq(value)]
    require(len(selected) == 1, "publication selector must identify exactly one row: " + _json(keys))
    return selected.iloc[0]


def classify_row(row):
    """Independently apply the unchanged pointwise theta/sigma criteria.

    Return (classification, driver). Equality at the practical threshold is
    within the bound, while an interval crossing it is unresolved. Missing
    bias uncertainty cannot establish within-criteria status.
    """
    parameter, truth = row["parameter"], float(row["truth"])
    require(parameter in PARAMETERS, "unknown parameter in publication summary")
    if truth == 0:
        require(parameter == "mu", "unexpected zero-truth parameter")
        return "not_applicable_zero_truth", "not_applicable"
    require(parameter != "mu", "the archived mu truth must remain zero")
    cl, cu = float(row["operational_wilson_lower"]), float(row["operational_wilson_upper"])
    require(math.isfinite(cl) and math.isfinite(cu) and 0 <= cl <= cu <= 1,
            "invalid operational coverage uncertainty")
    bl, bu = float(row["relative_bias_mc_lower"]), float(row["relative_bias_mc_upper"])
    have_bias = math.isfinite(bl) and math.isfinite(bu)
    require(not have_bias or bl <= bu, "reversed bias Monte Carlo interval")
    coverage_breach = cu < .8
    bias_breach = have_bias and (bl > .1 or bu < -.1)
    if coverage_breach or bias_breach:
        driver = "bias_and_coverage" if coverage_breach and bias_breach else (
            "coverage" if coverage_breach else "bias")
        return "supported_breach", driver
    if cl >= .8 and have_bias and bl >= -.1 and bu <= .1:
        return "within_criteria", "none"
    return "unresolved", "uncertain"


def _bound_table(bundle, study, filename):
    frame = bundle.table(study, filename).copy()
    path = bundle.study_dir(study) / filename
    frame.insert(0, "study", study)
    frame.insert(1, "study_scope", SCOPES[study])
    frame["source_path"] = path.relative_to(bundle.root).as_posix()
    frame["source_sha256"] = _sha(path)
    return frame


def _validate_summary(frame):
    require(not frame.duplicated(["study", *KEYS]).any(), "duplicate cross-study summary keys")
    require(set(frame.estimator) == set(METHODS), "incomplete summary estimator set")
    require(set(frame.parameter) == set(PARAMETERS), "incomplete summary parameter set")
    for _, group in frame.groupby(["study", "cell_id", "estimator"], sort=True):
        require(len(group) == 3 and set(group.parameter) == set(PARAMETERS),
                "a summary group does not contain every parameter")
        for field in ("n_planned", "n_point_valid", "n_interval_valid"):
            require(group[field].nunique() == 1, "joint-fit counts disagree between parameters")
        for row in group.to_dict("records"):
            counts = [float(row[c]) for c in ("n_planned", "n_point_valid", "n_interval_valid", "n_covered")]
            require(all(math.isfinite(c) and c == int(c) for c in counts), "noninteger summary counts")
            planned, point, interval, covered = map(int, counts)
            require(planned > 0 and 0 <= covered <= interval <= point <= planned,
                    "invalid summary count ordering")
            require(float(row["coverage_operational"]) == covered / planned,
                    "operational coverage denominator differs from planned paths")
            if interval:
                require(float(row["coverage_conditional"]) == covered / interval,
                        "conditional coverage denominator differs from valid intervals")
            else:
                require(pd.isna(row["coverage_conditional"]), "conditional coverage must be missing without intervals")
            wanted, _ = classify_row(row)
            require(row["practical_classification"] == wanted,
                    "saved classification differs from unrounded registered criteria")
            if row["parameter"] == "mu":
                relative = [name for name in frame if name.startswith("relative_")]
                require(all(pd.isna(row[name]) for name in relative), "relative mu metric at zero truth")


def _availability(bundle, summary):
    rows = []
    for study in SUMMARY_STUDIES:
        source = _bound_table(bundle, study, "outcome_counts.csv")
        require(not source.duplicated(["cell_id", "estimator", "point_status", "interval_status"]).any(),
                "duplicate outcome composition")
        expected = set(summary[summary.study.eq(study)][["cell_id", "estimator"]].itertuples(index=False, name=None))
        actual = set(source[["cell_id", "estimator"]].itertuples(index=False, name=None))
        require(expected == actual, "outcomes omit or add a planned cell/method")
        for (cid, method), composition in source.groupby(["cell_id", "estimator"], sort=True):
            row = select_one(summary, study=study, cell_id=cid, estimator=method, parameter="theta")
            counts = composition["count"].to_numpy()
            require(np.isfinite(counts).all() and (counts > 0).all() and np.equal(counts, counts.astype(int)).all(),
                    "invalid outcome counts")
            planned, point, interval = (int(row[k]) for k in ("n_planned", "n_point_valid", "n_interval_valid"))
            require(int(counts.sum()) == planned and composition.n_planned.eq(planned).all(),
                    "outcome accounting does not contain all planned paths")
            require(int(composition.loc[composition.point_status.eq("valid_interior"), "count"].sum()) == point,
                    "point eligibility differs from saved outcomes")
            require(int(composition.loc[composition.interval_status.eq("valid"), "count"].sum()) == interval,
                    "interval eligibility differs from saved outcomes")
            fields = {k: row[k] for k in ("study", "study_scope", "cell_id", "estimator", "n", "theta_mean_gap", "cv", "floor_fraction")}
            for name in ("ratio", "mechanism"):
                if name in row and not pd.isna(row[name]):
                    fields[name] = row[name]
            rows.append({**fields, "n_planned": planned, "n_point_valid": point,
                         "n_interval_valid": interval, "n_point_invalid": planned-point,
                         "n_point_valid_without_interval": point-interval,
                         "n_interval_unavailable": planned-interval,
                         "interval_unavailable_includes_point_invalid": True,
                         "parameter": "theta", "count_scope": "one_joint_fit_per_cell_replication_estimator_not_one_per_parameter",
                         "outcome_composition_json": _json(composition[["point_status", "interval_status", "count"]].to_dict("records")),
                         "source_path": row.source_path, "source_sha256": row.source_sha256,
                         "outcome_source_path": composition.iloc[0].source_path,
                         "outcome_source_sha256": composition.iloc[0].source_sha256})
    return pd.DataFrame(rows)


def _f1_ecdf(bundle, summary):
    fields = ("cell_id", "rep", "estimator", "point_valid", "theta", "mu", "sigma", "_source_path", "_source_sha256")
    raw = pd.DataFrame([{name: row[name] for name in fields} for row in bundle.f1_records()])
    require(not raw.duplicated(["cell_id", "rep", "estimator"]).any(), "duplicate retained F1 fit key")
    e1 = summary[summary.study.eq("e1")]
    expected_cells = set(e1.cell_id)
    require(set(raw.cell_id) == expected_cells, "F1 source set is incomplete")
    records = []
    for (cid, method), group in raw.groupby(["cell_id", "estimator"], sort=True):
        report = select_one(e1, cell_id=cid, estimator=method, parameter="theta")
        planned = int(report.n_planned)
        require(len(group) == planned and set(group.rep) == set(range(planned)), "F1 planned replications incomplete")
        require(group.point_valid.map(lambda v: isinstance(v, (bool, np.bool_))).all(), "invalid F1 eligibility flag")
        valid = group[group.point_valid]
        require(len(valid) == int(report.n_point_valid), "F1 point validity differs from audited summary")
        for parameter in PARAMETERS:
            report_p = select_one(e1, cell_id=cid, estimator=method, parameter=parameter)
            ordered = valid.sort_values([parameter, "rep"], kind="mergesort")
            require(np.isfinite(ordered[parameter]).all(), "nonfinite accepted F1 estimate")
            n_valid = len(ordered)
            for rank, row in enumerate(ordered.to_dict("records"), 1):
                records.append({"study": "e1", "cell_id": cid, "estimator": method,
                    "parameter": parameter, "n": int(report.n), "theta_mean_gap": float(report.theta_mean_gap),
                    "cv": float(report.cv), "rep": int(row["rep"]), "estimate": row[parameter],
                    "truth": float(report_p.truth), "ecdf": rank/n_valid, "rank": rank,
                    "n_planned": planned, "n_point_valid": n_valid,
                    "source_path": row["_source_path"], "source_sha256": row["_source_sha256"],
                    "sample_selection": "all_point_valid_saved_fits_in_cell_no_refit"})
    return pd.DataFrame(records)


def build_tables(bundle):
    """Create deterministic publication tables from verified source bytes."""
    frames = [_bound_table(bundle, study, "parameter_summary.csv") for study in SUMMARY_STUDIES]
    summary = pd.concat(frames, ignore_index=True, sort=False)
    _validate_summary(summary)
    classifications = [classify_row(row) for row in summary.to_dict("records")]
    summary = summary.assign(
        breach_driver=[driver for _, driver in classifications],
        classification_scope="pointwise_registered_practical_criteria_not_nominal_calibration",
        source_row_key=[_json({k: row[k] for k in KEYS}) for row in summary.to_dict("records")])
    summary = summary.sort_values(["study", "n", "theta_mean_gap", "cv", "floor_fraction", "cell_id", "estimator", "parameter"], kind="mergesort").reset_index(drop=True)
    t1_e2 = summary[summary.study.eq("e2") & summary.parameter.ne("mu")].reset_index(drop=True)
    t1_e1 = summary[summary.study.eq("e1") & summary.parameter.ne("mu")].reset_index(drop=True)
    unresolved = summary[summary.study.isin(["e1", "e2"]) & summary.practical_classification.eq("unresolved")].reset_index(drop=True)
    counts = []
    for (n, method, parameter), group in t1_e2.groupby(["n", "estimator", "parameter"], sort=True):
        counts.append({"study": "e2", "n": int(n), "estimator": method, "parameter": parameter,
                       "n_cells": len(group), **{name: int(group.practical_classification.eq(name).sum()) for name in CLASSES[:3]},
                       "bias_only": int(group.breach_driver.eq("bias").sum()),
                       "coverage_only": int(group.breach_driver.eq("coverage").sum()),
                       "bias_and_coverage": int(group.breach_driver.eq("bias_and_coverage").sum())})
    class_counts = []
    for study, group in summary.groupby("study", sort=True):
        class_counts.append({"study": study, "n_parameter_rows": len(group),
                             **{name: int(group.practical_classification.eq(name).sum()) for name in CLASSES}})
    calibration_fields = ["study", "study_scope", "cell_id", "estimator", "parameter", "n", "theta_mean_gap", "cv", "floor_fraction",
                          "mechanism", "ratio", "n_planned", "n_point_valid", "n_interval_valid",
                          "empirical_sd_interval_valid", "rms_reported_se", "mean_reported_se", "source_path", "source_sha256"]
    calibration = summary[[c for c in calibration_fields if c in summary]].copy()
    calibration = calibration.assign(
        sd_over_rms_se=np.divide(calibration.empirical_sd_interval_valid, calibration.rms_reported_se,
                                out=np.full(len(calibration), np.nan), where=calibration.rms_reported_se.gt(0)),
        eligibility_subset="same_interval_valid_fits_for_empirical_sd_and_reported_se",
        interpretation="descriptive_empirical_calibration_not_long_run_covariance_theorem")
    tables = {"parameter_summary.csv": summary, "T1_e2.csv": t1_e2, "T1_e1.csv": t1_e1,
              "T1_counts.csv": pd.DataFrame(counts), "classification_counts.csv": pd.DataFrame(class_counts),
              "unresolved.csv": unresolved, "availability.csv": _availability(bundle, summary),
              "calibration.csv": calibration,
              "e4_matched.csv": _bound_table(bundle, "e4", "matched_iid_comparisons.csv"),
              "e3_floor_contrasts.csv": _bound_table(bundle, "e3_floor", "epsilon_comparisons.csv"),
              "e4_sources.csv": _bound_table(bundle, "e4", "source_summary.csv"),
              "f1_ecdf.csv": _f1_ecdf(bundle, summary)}
    require((len(t1_e2), len(t1_e1), len(unresolved), len(tables["T1_counts.csv"])) == (360, 42, 6, 18),
            "frozen E1/E2 publication accounting changed")
    return tables


def _binding(frame, columns):
    """Retain exact original selectors and values behind each prose claim."""
    result = []
    derivations = {
        "n_point_invalid": ("n_planned - n_point_valid", ["n_planned", "n_point_valid"]),
        "n_point_valid_without_interval": ("n_point_valid - n_interval_valid", ["n_point_valid", "n_interval_valid"]),
        "n_interval_unavailable": ("n_planned - n_interval_valid", ["n_planned", "n_interval_valid"]),
    }
    originals = []
    derived = {}
    for column in columns:
        if column in derivations:
            expression, inputs = derivations[column]
            derived[column] = {"expression": expression, "source_columns": inputs}
            originals.extend(inputs)
        else:
            originals.append(column)
    originals = list(dict.fromkeys(originals))
    for (path, sha), group in frame.groupby(["source_path", "source_sha256"], sort=True):
        selectors = [c for c in ("cell_id", "estimator", "parameter", "reference_cell_id", "point_status", "interval_status") if c in group]
        result.append({"path": path, "sha256": sha, "selector_columns": selectors,
                       "columns": originals, "rows": _plain(group[selectors + [c for c in originals if c not in selectors]].to_dict("records")),
                       "derived_columns": derived})
    return result


def _claim(cid, text, frame, columns, *, status="supported", values=None, figures=(), qualification="", uncertainty="saved pointwise Monte Carlo uncertainty", comparison="within_cell_summary", units="mixed; see source columns"):
    studies = sorted(frame.study.unique())
    counts = [c for c in ("n_planned", "n_point_valid", "n_interval_valid", "cell_n_planned", "cell_n_point_valid", "cell_n_interval_valid",
                          "reference_n_planned", "reference_n_point_valid", "reference_n_interval_valid", "n_common_point_valid", "n_common_interval_valid") if c in frame]
    design_columns = [c for c in ("study", "cell_id", "n", "theta_mean_gap", "cv", "floor_fraction", "mechanism", "ratio") if c in frame]
    return {"claim_id": cid, "text": text, "status": status, "studies": studies,
            "experiment_scope": sorted(frame.study_scope.unique()), "values": _plain(values or {}),
            "parameters": sorted(frame.parameter.unique()) if "parameter" in frame else [],
            "estimators": sorted(frame.estimator.unique()) if "estimator" in frame else [],
            "design_points": _plain(frame[design_columns].drop_duplicates().to_dict("records")),
            "source_bindings": _binding(frame, list(dict.fromkeys([*columns, *counts]))),
            "units": units, "uncertainty": uncertainty, "comparison": comparison,
            "independent_unit": "source_path; paired floor contrasts use independent coupling groups",
            "denominator_definitions": {"point_error": "point-valid fits; within-path contrasts use the stated common-valid subset",
                "operational_coverage": "all planned paths; unavailable intervals count as noncoverage",
                "conditional_coverage": "interval-valid fits", "calibration": "same interval-valid fits for both quantities"},
            "qualifications": qualification, "figure_ids": list(figures)}


def build_claims(bundle, tables):
    """Generate numerical prose and its exact source-row bindings together."""
    summary, available = tables["parameter_summary.csv"], tables["availability.csv"]
    claims = []
    for study in ("e1", "e2", "e4"):
        group = summary[summary.study.eq(study)]
        counts = {label: int(group.practical_classification.eq(label).sum()) for label in CLASSES}
        text = (f"{study.upper()} has {counts['within_criteria']} within-criteria, {counts['supported_breach']} supported-breach, "
                f"and {counts['unresolved']} unresolved theta/sigma classifications; {counts['not_applicable_zero_truth']} mu classifications are inapplicable at zero truth.")
        claims.append(_claim(study+"_classifications", text, group, ["practical_classification", *CLASS_COLUMNS], values=counts, figures=("F8",) if study == "e4" else ("F7", "T1"),
                             qualification="Pointwise grid statements under the broad 80%-coverage/10%-relative-bias criteria; within criteria does not establish nominal 95% calibration.", units="parameter/design classification counts"))
        a = available[available.study.eq(study)]
        values = {c: int(a[c].sum()) for c in ("n_planned", "n_point_valid", "n_interval_valid", "n_point_invalid", "n_point_valid_without_interval", "n_interval_unavailable")}
        paths = int(a[a.estimator.eq("exact")].n_planned.sum())
        values["planned_source_paths"] = paths
        text = (f"{study.upper()} retains {paths:,} source paths and {values['n_planned']:,} estimator records, "
                f"including {values['n_point_invalid']} invalid points and {values['n_point_valid_without_interval']} additional interval rejections "
                f"({values['n_interval_unavailable']} unavailable intervals in total).")
        claims.append(_claim(study+"_availability", text, a, ["n_point_invalid", "n_point_valid_without_interval", "n_interval_unavailable"], values=values,
                             figures=("F8",) if study == "e4" else ("F2", "F6"), qualification="Point-invalid records are already included in interval-unavailable records. Counts are estimator records, not independent source failures.", uncertainty="complete retained-record accounting", units="source paths and estimator records"))
    short = summary[summary.study.eq("e2") & summary.estimator.eq("exact") & summary.parameter.eq("theta") & summary.n.eq(250) & summary.theta_mean_gap.eq(.1)]
    require(len(short) == 5 and short.breach_driver.eq("bias").all(), "short-span exact-theta explanation changed")
    bias = [float(short.relative_bias.min()), float(short.relative_bias.max())]
    cov = [float(short.coverage_operational.min()), float(short.coverage_operational.max())]
    claims.append(_claim("e2_short_span_exact_theta", f"At n=250 and q=0.1, all five exact-theta breaches are bias-driven: relative bias ranges from {100*bias[0]:.1f}% to {100*bias[1]:.1f}%, while operational coverage ranges from {100*cov[0]:.1f}% to {100*cov[1]:.1f}%.",
        short, ["n", "theta_mean_gap", "cv", "nominal_dimensionless_span", "relative_bias", "relative_bias_mc_lower", "relative_bias_mc_upper", "coverage_operational", "operational_wilson_lower", "operational_wilson_upper"],
        values={"relative_bias_range": bias, "operational_coverage_range": cov, "n_cells": len(short)}, figures=("F5", "F7", "T1"), qualification="The regular reference is included. The anticipated coverage collapse did not occur; finite-span effects cannot all be attributed to irregularity.", units="bias and coverage proportions"))
    for method, parameter in (("pfml", "theta"), ("euler", "sigma")):
        group = summary[summary.study.eq("e2") & summary.n.eq(250) & summary.theta_mean_gap.eq(.1) & summary.estimator.eq(method) & summary.parameter.eq(parameter)].sort_values("cv")
        sequence = group.practical_classification.tolist()
        claims.append(_claim(f"e2_nonmonotone_{method}_{parameter}", f"At n=250 and q=0.1, {method.upper()} {parameter} classifications across ascending CV are " + ", ".join(sequence) + ".",
            group, ["cv", "practical_classification", *CLASS_COLUMNS], values={"cv": group.cv.tolist(), "classifications": sequence}, figures=("F7", "T1"),
            qualification="A single maximum acceptable CV or a monotone interpolated operating boundary would misrepresent these observed classifications.", uncertainty="saved pointwise classifications", units="classification labels"))
    unresolved = tables["unresolved.csv"]
    claims.append(_claim("headline_unresolved", f"All {len(unresolved)} unresolved E1/E2 parameter-specific classifications are retained separately from supported breaches and within-criteria results.",
        unresolved, ["n", "theta_mean_gap", "cv", "practical_classification", *CLASS_COLUMNS], values={"unresolved": len(unresolved)}, figures=("F7", "T1"),
        qualification="Overlapping E1/E2 nominal cells use independent seeds and are not pooled or selected between to obtain a preferred classification.", units="parameter/design classifications"))
    # An explicit denominator example makes operational versus conditional coverage tangible.
    r = select_one(summary, study="e2", n=250, theta_mean_gap=2., cv=0., estimator="exact", parameter="theta")
    f = r.to_frame().T
    claims.append(_claim("e2_operational_conditional_example", f"For regular exact theta at n=250,q=2, {int(r.n_covered)} covering intervals give {100*r.coverage_operational:.2f}% operational coverage over {int(r.n_planned)} planned paths versus {100*r.coverage_conditional:.2f}% conditional coverage over {int(r.n_interval_valid)} available intervals.",
        f, ["n_covered", "coverage_operational", "coverage_conditional", "operational_wilson_lower", "operational_wilson_upper", "conditional_wilson_lower", "conditional_wilson_upper"], figures=("F6",),
        values={"n_covered": int(r.n_covered), "operational": float(r.coverage_operational), "conditional": float(r.coverage_conditional)}, qualification="Unavailable intervals remain noncoverage operationally; bias/RMSE use point-valid fits.", units="coverage proportions"))
    # Registered data include favorable as well as adverse finite-sample comparisons.
    contrasts = _bound_table(bundle, "e2", "regular_comparisons.csv")
    favorable = contrasts[contrasts.n.eq(4000) & contrasts.theta_mean_gap.eq(2.) & contrasts.cv.eq(3.) & contrasts.estimator.eq("exact")].sort_values("parameter")
    require(len(favorable) == 3, "nominal-design example missing parameters")
    ratios = {row.parameter: float(row.rmse_ratio) for row in favorable.itertuples()}
    claims.append(_claim("e2_irregularity_mixed_rmse", f"At n=4000,q=2,CV=3, exact-estimator RMSE ratios to the regular reference are {ratios['theta']:.3f} for theta, {ratios['sigma']:.3f} for sigma and {ratios['mu']:.3f} for mu.",
        favorable, ["n", "theta_mean_gap", "cv", "rmse_ratio", "rmse_ratio_mc_lower", "rmse_ratio_mc_upper"], values=ratios, figures=("F4",), comparison="independent_cells_nominal_design_reference",
        qualification="The nominal design is matched; realized spans differ. This is evidence against a blanket assertion that irregularity worsens every parameter.", units="RMSE ratios"))
    ranking = summary[summary.study.eq("e2") & summary.n.eq(250) & summary.theta_mean_gap.eq(.1)
                      & summary.cv.eq(3.) & summary.parameter.eq("theta")].sort_values("estimator")
    require(len(ranking) == 3, "small-span RMSE ranking requires all three estimators")
    errors = {row.estimator: float(row.rmse) for row in ranking.itertuples()}
    claims.append(_claim("e2_exact_not_uniform_lowest_rmse", f"At n=250,q=0.1,CV=3, observed theta RMSE is {errors['exact']:.6f} for exact, {errors['pfml']:.6f} for PFML and {errors['euler']:.6f} for Euler.",
        ranking, ["n", "theta_mean_gap", "cv", "rmse", "mcse_rmse", "relative_bias", "relative_bias_mc_lower", "relative_bias_mc_upper"],
        values={"rmse_by_method": errors}, status="descriptive", figures=("F1", "F4"),
        comparison="three_estimators_on_shared_source_paths",
        qualification="This descriptive point-error ranking does not establish a general ordering. Within-source inferential method contrasts require the saved paired-comparison uncertainty; practical classification assesses bias and coverage separately.",
        uncertainty="saved per-method RMSE MCSE; this sentence asserts no paired significance test", units="theta RMSE in natural units"))
    # Exploratory follow-up availability is distinct from the registered grid.
    a = available[available.study.eq("e3_followup")]
    vals = {c: int(a[c].sum()) for c in ("n_planned", "n_point_invalid", "n_point_valid_without_interval", "n_interval_unavailable")}
    vals["planned_source_paths"] = int(a[a.estimator.eq("exact")].n_planned.sum())
    claims.append(_claim("e3_followup_accounting", f"The independent exploratory follow-up retains {vals['planned_source_paths']:,} paths and {vals['n_planned']:,} fits, with {vals['n_point_invalid']} invalid points and {vals['n_point_valid_without_interval']} additional interval rejections ({vals['n_interval_unavailable']} unavailable intervals).",
        a, ["n_point_invalid", "n_point_valid_without_interval", "n_interval_unavailable"], values=vals, figures=("S_E3",), qualification="These exploratory results do not overwrite E1/E2 classifications. One exact point remains boundary_suspected.", uncertainty="complete retained-record accounting", units="source paths and estimator records"))
    pending = _bound_table(bundle, "e3_followup", "outcome_counts.csv")
    pending = pending[pending.point_status.eq("boundary_suspected")]
    require(len(pending) == 1 and int(pending["count"].sum()) == 1 and pending.estimator.eq("exact").all(),
            "exploratory unresolved-boundary disposition changed")
    claims[-1]["source_bindings"].extend(_binding(pending, ["point_status", "interval_status", "count", "n_planned"]))
    selected_path = bundle.study_dir("e3_selected") / "audit.json"
    selected = json.loads(selected_path.read_text())
    require(selected["status"] == "passed" and selected["evaluated_sources"] == selected["planned_sources"]
            and selected["verified_fit_replays"] == selected["original_fit_records"], "selected diagnostic audit is incomplete")
    selected_columns = ["planned_sources", "evaluated_sources", "original_fit_records", "verified_fit_replays",
                        "independent_recalculation_of_all_diagnostic_math"]
    claims.append({"claim_id": "e3_selected_replay_scope",
        "text": f"The selected-path diagnostic archive evaluated {selected['evaluated_sources']} sources and verified {selected['verified_fit_replays']} original fit replays.",
        "status": "supported", "studies": ["e3_selected"], "experiment_scope": ["exploratory_selected_path_diagnostics"],
        "values": {k: selected[k] for k in selected_columns},
        "source_bindings": [{"path": selected_path.relative_to(bundle.root).as_posix(), "sha256": _sha(selected_path),
                             "columns": selected_columns, "values": {k: selected[k] for k in selected_columns}}],
        "units": "selected sources and original fit replays", "uncertainty": "complete fixed-selection accounting; no prevalence inference",
        "comparison": "original_fit_replay_at_fixed_selected_sources", "independent_unit": "selected source; selection is not a random frequency study",
        "qualifications": "This selected audit does not independently recalculate all diagnostic mathematics, prove global optimality, or replace headline outcomes.",
        "figure_ids": ["S_E3"]})
    audit_path = bundle.study_dir("e3_floor") / "audit.json"
    audit = json.loads(audit_path.read_text())
    floor = tables["e3_floor_contrasts.csv"]
    maximum = floor.loc[(floor.rmse_ratio-1).abs().idxmax()].to_frame().T
    m = maximum.iloc[0]
    max_cov = float(floor.operational_coverage_difference.abs().max())
    c = _claim("e3_floor_sensitivity", f"Across the specified paired floor contrasts, the largest estimated RMSE-ratio departure from one is {m.rmse_ratio:.8f}, with pointwise MC interval [{m.rmse_ratio_mc_lower:.8f}, {m.rmse_ratio_mc_upper:.8f}]; the largest absolute observed coverage difference is {100*max_cov:.2f} percentage points.",
        floor, ["cv", "floor_fraction", "reference_floor_fraction", "rmse_ratio", "rmse_ratio_mc_lower", "rmse_ratio_mc_upper", "operational_coverage_difference"],
        values={"largest_rmse_ratio_departure": float(m.rmse_ratio), "largest_absolute_coverage_difference": max_cov,
                "source_realizations": audit["planned_source_realizations"], "independent_coupling_groups": audit["planned_independent_coupling_groups"]},
        figures=("S_E3",), comparison="paired_floor_shared_draws_and_standardized_innovations", qualification="Small effects within this grid do not establish equivalence or universal floor invariance. The physical realizations are not all independent.", units="RMSE ratios and coverage proportions")
    c["source_bindings"].append({"path": audit_path.relative_to(bundle.root).as_posix(), "sha256": _sha(audit_path), "columns": ["planned_source_realizations", "planned_independent_coupling_groups"], "values": {k: audit[k] for k in ("planned_source_realizations", "planned_independent_coupling_groups")}})
    claims.append(c)
    matched = tables["e4_matched.csv"]
    for ratio in (5, 20):
        row = select_one(matched, ratio=ratio, estimator="pfml", parameter="mu")
        one = row.to_frame().T
        values = {name: float(row[name]) for name in ("rmse_ratio", "rmse_ratio_mc_lower", "rmse_ratio_mc_upper", "operational_coverage_difference", "operational_coverage_difference_mc_lower", "operational_coverage_difference_mc_upper", "excess_bias", "excess_bias_mc_lower", "excess_bias_mc_upper")}
        text = (f"At r={ratio}, PFML mu has MMPP/IID RMSE ratio {row.rmse_ratio:.6f} "
                f"[95% MC {row.rmse_ratio_mc_lower:.6f}, {row.rmse_ratio_mc_upper:.6f}] and an operational-coverage difference "
                f"of {100*row.operational_coverage_difference:+.2f} percentage points "
                f"[{100*row.operational_coverage_difference_mc_lower:+.3f}, {100*row.operational_coverage_difference_mc_upper:+.3f}].")
        claims.append(_claim(f"e4_pfml_mu_r{ratio}", text, one, ["ratio", *values], values=values, figures=("F8",),
            comparison="independent_MMPP_minus_full_marginal_IID_cells",
            uncertainty="pointwise normal MC interval for coverage difference; log-delta MC interval for RMSE ratio",
            qualification="Exploratory fixed design including realized-span consequences. This contrast does not establish a mean-bias shift or a universal clustering effect. PFML intervals use the registered working-model covariance.", units="RMSE ratios; coverage differences as proportions in values and percentage points in prose; mu bias in original state units"))
    bias_null = matched.excess_bias_mc_lower.le(0) & matched.excess_bias_mc_upper.ge(0)
    theta_sigma = matched[matched.parameter.ne("mu")]
    ratio_null = theta_sigma.rmse_ratio_mc_lower.le(1) & theta_sigma.rmse_ratio_mc_upper.ge(1)
    claims.append(_claim("e4_null_contrasts", f"All {int(bias_null.sum())} of {len(matched)} matched bias intervals contain zero; all {int(ratio_null.sum())} of {len(theta_sigma)} theta/sigma RMSE-ratio intervals contain one.",
        matched, ["ratio", "excess_bias_mc_lower", "excess_bias_mc_upper", "rmse_ratio_mc_lower", "rmse_ratio_mc_upper"], values={"bias_containing_zero": int(bias_null.sum()), "bias_contrasts": len(matched), "theta_sigma_ratio_containing_one": int(ratio_null.sum()), "theta_sigma_ratio_contrasts": len(theta_sigma)},
        figures=("F8",), comparison="independent_MMPP_minus_full_marginal_IID_cells", qualification="Intervals containing the null do not demonstrate equivalence. These comparisons are pointwise, not a familywise guarantee.", units="contrast counts"))
    isolated = select_one(matched, ratio=20, estimator="exact", parameter="sigma").to_frame().T
    i = isolated.iloc[0]
    claims.append(_claim("e4_isolated_exact_sigma_coverage", f"At r=20, the exact-sigma operational-coverage contrast is {100*i.operational_coverage_difference:+.2f} percentage points, with pointwise MC interval [{100*i.operational_coverage_difference_mc_lower:+.5f}, {100*i.operational_coverage_difference_mc_upper:+.5f}].",
        isolated, ["ratio", "operational_coverage_difference", "operational_coverage_difference_mc_lower", "operational_coverage_difference_mc_upper"],
        values={k: float(i[k]) for k in ("operational_coverage_difference", "operational_coverage_difference_mc_lower", "operational_coverage_difference_mc_upper")},
        figures=("F8",), comparison="independent_MMPP_minus_full_marginal_IID_cells",
        qualification="This isolated pointwise result among multiple comparisons does not establish a general beneficial clustering effect.", units="coverage proportions"))
    cal = tables["calibration.csv"]
    cal = cal[cal.study.eq("e4") & cal.parameter.eq("mu") & cal.estimator.eq("pfml") & cal.ratio.eq(20) & cal.mechanism.isin(["mmpp", "iid"])].sort_values("mechanism")
    require(len(cal) == 2, "E4 mu calibration example requires both independent arms")
    values = {row.mechanism: {"empirical_sd": float(row.empirical_sd_interval_valid), "rms_reported_se": float(row.rms_reported_se), "sd_over_rms_se": float(row.sd_over_rms_se)} for row in cal.itertuples()}
    claims.append(_claim("e4_pfml_mu_calibration", f"At r=20, PFML mu empirical SD / RMS reported SE is {values['mmpp']['sd_over_rms_se']:.3f} under MMPP and {values['iid']['sd_over_rms_se']:.3f} under IID sampling.",
        cal, ["ratio", "mechanism", "empirical_sd_interval_valid", "rms_reported_se"], values=values, figures=("S_calibration",), status="descriptive",
        qualification="Both quantities use the same interval-valid fits. This descriptive comparison is not a derived dependence-robust covariance benchmark.", uncertainty="descriptive saved variability estimates; no new MC interval claimed", units="dimensionless SD/SE ratios"))
    # Negative scope claims bind the prospective derivation rather than inventing a covariance result.
    theory_path = bundle.root / "docs/week8_population_targets.md"
    theory_sha = _sha(theory_path)
    for cid, text, status in (
        ("shared_expected_objective_targets", "Under the stated exogeneity and stationary initialization assumptions, matched full one-gap marginals share the derived expected-objective minimizers; this does not imply equal finite-sample variance or coverage.", "theoretical_under_stated_assumptions"),
        ("beyond_first_order_theory", "A demonstrated departure from a fully matched first-order coverage theory is not established by the current experiments or population-target calculation.", "unsupported_claim_excluded"),
        ("general_dependent_gap_covariance", "A general dependent-gap CLT and a dependence-robust PFML/Euler covariance formula are not established in this study.", "unsupported_claim_excluded")):
        claims.append({"claim_id": cid, "text": text, "status": status, "studies": ["e4"], "experiment_scope": ["prospective_mathematical_interpretation"],
            "values": {}, "source_bindings": [{"path": theory_path.relative_to(bundle.root).as_posix(), "sha256": theory_sha,
                "section": "Mean-gap PFML; What matched marginals do and do not imply"}], "units": "not applicable", "uncertainty": "no covariance or coverage theorem claimed",
            "comparison": "stationary_expected_objective_calculation", "independent_unit": "not applicable",
            "qualifications": "Population targets are interpretation benchmarks; all simulation error summaries continue to use the true OU parameters.", "figure_ids": []})
    require(len({claim["claim_id"] for claim in claims}) == len(claims), "duplicate claim identifier")
    return _plain(claims)
