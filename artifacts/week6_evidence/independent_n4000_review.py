"""Independent saved-data review; no simulation, fitting, or pipeline mutation.

This reviewer imports no project reporting/audit functions. The raw subset was
selected before E2 outcome tables existed: two n=250 designs plus the four
n=4000 (q, CV) corners. All classification labels and report totals are checked;
raw summary arithmetic is checked only in the stated six cells. The completed
project audit must already exist and pass. Visual inspection remains a separate
human/model review, not a claim made by this script.

Run from the repository: .venv/bin/python artifacts/week6_evidence/independent_n4000_review.py
Only the JSON report beside this script can be written. Existing unequal bytes
are never overwritten. The experiment directory is read-only throughout.
"""
import csv
from collections import Counter
import gzip
import hashlib
from itertools import product
import json
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
RUN = ROOT / "artifacts/week6_e2"
OUTPUT = Path(__file__).with_suffix(".json")
FREEZE = "960436615a710d21de3439d50257f50af72a45e4"
SELECTED = ((250, 2.0, 0.0), (250, .1, 3.0),
            (4000, .1, 0.0), (4000, .1, 3.0), (4000, 2.0, 0.0), (4000, 2.0, 3.0))
METHODS = ("exact", "pfml", "euler")
PARAMETERS = ("theta", "mu", "sigma")
CLASSES = ("within_criteria", "supported_breach", "unresolved", "not_applicable_zero_truth")
Z = 1.959963984540054
checked = 0
fingerprints = {}


def require(condition, message):
    global checked
    if not condition:
        raise ValueError(message)
    checked += 1


def number(value):
    return float(value) if value not in (None, "") else float("nan")


def equal(actual, expected, name):
    a, b = number(actual), number(expected)
    require((math.isnan(a) and math.isnan(b)) or math.isclose(a, b, rel_tol=2e-11, abs_tol=2e-12),
            f"{name}: {a!r} != {b!r}")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def remember(path):
    value = digest(path)
    fingerprints[str(path.relative_to(ROOT))] = value
    return value


def document(path):
    remember(path)
    return json.loads(path.read_text())


def rows(path):
    remember(path)
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def wilson(successes, total):
    if total == 0:
        return float("nan"), float("nan")
    p = successes/total
    denominator = 1 + Z*Z/total
    center = (p+Z*Z/(2*total))/denominator
    radius = Z*math.sqrt(p*(1-p)/total+Z*Z/(4*total*total))/denominator
    return max(0.0, center-radius), min(1.0, center+radius)


def classify(truth, bias, se_bias, coverage_lower, coverage_upper):
    if truth == 0:
        return "not_applicable_zero_truth"
    lo, hi = (bias-Z*se_bias)/abs(truth), (bias+Z*se_bias)/abs(truth)
    if coverage_upper < .8 or lo > .1 or hi < -.1:
        return "supported_breach"
    if coverage_lower >= .8 and lo >= -.1 and hi <= .1:
        return "within_criteria"
    return "unresolved"


def clean(value):
    if isinstance(value, dict):
        return {k: clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean(v) for v in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def main():
    audit = document(RUN/"audit.json")
    require(audit["status"] == "passed" and audit["full_completion_audit"] is True,
            "A passing full experiment audit is required before this independent review")
    require(audit["code_commit"] == FREEZE, "Unexpected E2 freeze identity")
    equal(audit["planned_paths"], 120000, "Audit planned source count")
    equal(audit["planned_fit_records"], 360000, "Audit planned estimator count")
    equal(audit["pipeline_errors"], 0, "Audit pipeline exceptions")
    design = document(RUN/"resolved_design.json")
    require(digest(RUN/"resolved_design.json") == audit["audited_artifact_sha256"]["resolved_design.json"],
            "Saved design differs from the passed audit")
    require((design["name"], design["seed"], design["replications"], design["scope"])
            == ("week6_e2", 20261002, 2000, "headline"), "Unexpected E2 scientific identity")
    cells = {c["cell_id"]: c for c in design["cells"]}
    require(len(cells) == 60 and {(c["n"], c["theta_mean_gap"], c["cv"]) for c in cells.values()}
            == set(product((250, 1000, 4000), (.1, .5, 1, 2), (0, .5, 1, 2, 3))), "Incomplete E2 grid")
    manifest = document(RUN/"table_manifest.json")
    require(digest(RUN/"table_manifest.json") == audit["audited_artifact_sha256"]["table_manifest.json"],
            "Core table manifest differs from the passed audit")
    for name, entry in manifest.items():
        require(remember(RUN/name) == entry["sha256"] == audit["audited_artifact_sha256"][name],
                f"Core table/audit digest: {name}")
    figures = document(RUN/"figures/figure_manifest.json")
    require(figures["design_sha256"] == audit["design_sha256"], "Figure/audit design identity")
    expected_images = {f"{stem}_n{n}.{ext}" for n in (250, 1000, 4000)
                       for stem in ("F5_operational", "F6_conditional", "F6_availability") for ext in ("png", "pdf")}
    reports = {"classification_table.csv", "classification_counts.csv", "failed_cells.csv", "report_manifest.json"}
    require(set(figures["figures"]) == expected_images | reports, "Incomplete derived artifact manifest")
    for name, sha in figures["figures"].items():
        require(remember(RUN/"figures"/name) == sha, f"Derived artifact digest: {name}")
    report_manifest = document(RUN/"figures/report_manifest.json")
    require(report_manifest["design_sha256"] == audit["design_sha256"], "Report/audit design identity")
    require(report_manifest["input_table_sha256"] == figures["input_table_sha256"], "Report/figure parent disagreement")
    for name, sha in figures["input_table_sha256"].items():
        require(sha == manifest[name]["sha256"], "Incorrect parent input hash")

    summary_rows = rows(RUN/"parameter_summary.csv")
    summary = {(r["cell_id"], r["estimator"], r["parameter"]): r for r in summary_rows}
    require(len(summary_rows) == 540 and set(summary) == set(product(cells, METHODS, PARAMETERS)), "Summary Cartesian keys")
    class_counts = Counter()
    for key, row in summary.items():
        planned, point, interval, covered = (int(row[k]) for k in ("n_planned", "n_point_valid", "n_interval_valid", "n_covered"))
        require(planned == 2000 and 0 <= covered <= interval <= point <= planned, f"Count nesting: {key}")
        for kind, valid in (("point", point), ("interval", interval)):
            equal(row[kind+"_availability"], valid/planned, f"{key} {kind} availability")
        for scope, denominator in (("operational", planned), ("conditional", interval)):
            estimate = covered/denominator if denominator else float("nan")
            lo, hi = wilson(covered, denominator)
            equal(row["coverage_"+scope], estimate, f"{key} {scope} coverage")
            equal(row[scope+"_wilson_lower"], lo, f"{key} {scope} Wilson lower")
            equal(row[scope+"_wilson_upper"], hi, f"{key} {scope} Wilson upper")
        lo, hi = wilson(covered, planned)
        truth, bias, se = (number(row[k]) for k in ("truth", "bias", "mcse_bias"))
        label = classify(truth, bias, se, lo, hi)
        require(label == row["practical_classification"], f"Independent classification: {key}")
        class_counts[(int(row["n"]), row["parameter"], row["estimator"], label)] += 1
        if truth == 0:
            require(all(math.isnan(number(row[k])) for k in ("relative_bias", "relative_rmse", "relative_mcse_bias",
                    "relative_mcse_rmse", "relative_bias_mc_lower", "relative_bias_mc_upper")), "Zero-truth relative metric")

    count_rows = rows(RUN/"figures/classification_counts.csv")
    expected_count_keys = set(product((250, 1000, 4000), PARAMETERS, METHODS, CLASSES))
    require(len(count_rows) == 108 and {(int(r["n"]), r["parameter"], r["estimator"], r["practical_classification"])
            for r in count_rows} == expected_count_keys, "Classification-count keys")
    for row in count_rows:
        key = (int(row["n"]), row["parameter"], row["estimator"], row["practical_classification"])
        equal(row["n_cells"], class_counts[key], "Independent class total")
        equal(row["n_planned_cells"], 20, "Planned per-n class denominator")
    classification_rows = rows(RUN/"figures/classification_table.csv")
    require(len(classification_rows) == 540 and {(r["cell_id"], r["estimator"], r["parameter"])
            for r in classification_rows} == set(summary), "Classification table complete keys")
    for row in classification_rows:
        key = (row["cell_id"], row["estimator"], row["parameter"])
        require(row["practical_classification"] == summary[key]["practical_classification"], "Exported class mismatch")

    failures = rows(RUN/"figures/failed_cells.csv")
    require(len(failures) == 180 and {(r["cell_id"], r["estimator"]) for r in failures}
            == set(product(cells, METHODS)), "Failure log complete keys")
    for row in failures:
        key = row["cell_id"], row["estimator"]
        source = summary[(*key, "theta")]
        p, s, v = (int(source[k]) for k in ("n_planned", "n_point_valid", "n_interval_valid"))
        for name, expected in (("n_point_invalid", p-s), ("n_interval_unavailable", p-v), ("n_point_valid_without_interval", s-v)):
            equal(row[name], expected, f"{key} {name}")
        require((row["has_failure"] == "True") == (v < p), "Failure flag")
        require(sum(int(x["count"]) for x in json.loads(row["status_counts_json"])) == p, "Failure status accounting")

    selected_stats = []
    selected_partitions = 0
    for selection in SELECTED:
        matching = [c for c in cells.values() if (c["n"], c["theta_mean_gap"], c["cv"]) == selection]
        require(len(matching) == 1, "Raw selection identity")
        cell = matching[0]
        raw = {m: {} for m in METHODS}
        for path in sorted((RUN/"cells"/cell["cell_id"]).glob("part-*.jsonl.gz")):
            require(remember(path) == audit["audited_artifact_sha256"][str(path.relative_to(RUN))], "Selected raw digest")
            selected_partitions += 1
            with gzip.open(path, "rt") as stream:
                for line in stream:
                    record = json.loads(line)
                    method, replication = record["estimator"], record["rep"]
                    require(method in raw and replication not in raw[method], "Duplicate selected raw key")
                    raw[method][replication] = record
        require(all(set(r) == set(range(2000)) for r in raw.values()), "Missing selected raw planned key")
        require(all(len({(raw[m][i]["times_sha256"], raw[m][i]["observations_sha256"]) for m in METHODS}) == 1
                    for i in range(2000)), "Selected estimators must share each source")
        for method in METHODS:
            records = [raw[method][i] for i in range(2000)]
            valid_points = [r for r in records if r["point_valid"]]
            valid_intervals = [r for r in records if r["interval_valid"]]
            require(all(math.isfinite(r["nll"]) for r in valid_points), "Selected valid point NLL")
            for parameter in PARAMETERS:
                truth = cell["true_"+parameter]
                errors = np.asarray([r[parameter]-truth for r in valid_points], dtype=float)
                s, v = len(errors), len(valid_intervals)
                bias = float(errors.mean()) if s else float("nan")
                rmse = float(np.sqrt(np.mean(errors*errors))) if s else float("nan")
                se_bias = float(errors.std(ddof=1)/math.sqrt(s)) if s >= 2 else float("nan")
                se_rmse = float((errors*errors).std(ddof=1)/(2*math.sqrt(s)*rmse)) if s >= 2 and rmse > 0 else float("nan")
                covered = sum(r[parameter+"_lower"] <= truth <= r[parameter+"_upper"] for r in valid_intervals)
                negative = sum(r[parameter+"_lower"] < 0 for r in valid_intervals)
                values = {"n_point_valid": s, "n_interval_valid": v, "n_covered": covered, "negative_lower_count": negative,
                          "bias": bias, "rmse": rmse, "mcse_bias": se_bias, "mcse_rmse": se_rmse,
                          "coverage_operational": covered/2000,
                          "coverage_conditional": covered/v if v else float("nan")}
                row = summary[(cell["cell_id"], method, parameter)]
                for name, value in values.items():
                    equal(row[name], value, f"Saved raw {selection} {method} {parameter} {name}")
                selected_stats.append({"n": selection[0], "q": selection[1], "cv": selection[2], "estimator": method,
                                       "parameter": parameter, **values})

    total_point = sum(int(r["n_point_valid"]) for r in summary_rows if r["parameter"] == "theta")
    total_interval = sum(int(r["n_interval_valid"]) for r in summary_rows if r["parameter"] == "theta")
    equal(audit["point_valid_records"], total_point, "All-cell point total against completed audit")
    equal(audit["interval_valid_records"], total_interval, "All-cell interval total against completed audit")
    for name, sha in fingerprints.items():
        require(digest(ROOT/name) == sha, "Input changed during read-only review: "+name)
    result = {
        "status": "passed", "review_kind": "independent saved-data arithmetic and reporting validation; no scientific rerun",
        "scope": "all540 summary classifications and totals; selected raw summaries only; visual review separate",
        "code_freeze": FREEZE, "run_identity_sha256": audit["run_identity_sha256"],
        "full_experiment_audit_status": audit["status"], "checks_passed": checked,
        "planned_paths": 120000, "planned_fit_records": 360000,
        "point_valid_records": total_point, "interval_valid_records": total_interval,
        "selected_raw_cells": SELECTED, "selected_raw_paths": 2000*len(SELECTED),
        "selected_raw_fit_records": 6000*len(SELECTED), "selected_raw_partitions": selected_partitions,
        "selected_raw_summaries": selected_stats, "classification_counts": count_rows,
        "failure_rows": [r for r in failures if r["has_failure"] == "True"],
        "n4000_summary": [r for r in summary_rows if int(r["n"]) == 4000],
        "reviewer_script_sha256": digest(Path(__file__)), "input_sha256": fingerprints,
        "limitations": ["Does not independently replay paths or estimators; the full project audit records its own replay scope.",
                        "Classification checks use saved bias and MCSE except where raw summary arithmetic is explicitly checked.",
                        "All uncertainty and classifications remain pointwise; no simultaneous or continuous-boundary claim.",
                        "A passing review does not identify the cause of an earlier transient validation metadata-status reset."]}
    data = (json.dumps(clean(result), indent=2, sort_keys=True, allow_nan=False)+"\n").encode()
    if OUTPUT.exists():
        require(OUTPUT.read_bytes() == data, "Refusing to overwrite a different independent review")
    else:
        with OUTPUT.open("xb") as stream:
            stream.write(data)
    print(json.dumps({k: result[k] for k in ("status", "checks_passed", "point_valid_records", "interval_valid_records",
                                            "selected_raw_paths", "selected_raw_fit_records")}, sort_keys=True))
    print(str(OUTPUT.relative_to(ROOT)))


if __name__ == "__main__":
    main()
