"""Read-only E2 provenance/failure review; no simulation or fitting imports.

Stream all saved raw records to reconcile keys, joint validity, coverage counts,
status composition and the exact failed-fit export. This is narrower than the
project's full mathematical/replay audit and does not replace it. The three
n=1000 PNGs were separately inspected by the reviewer before this script ran.
Only the adjacent review JSON is written; existing different evidence is kept.
"""
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
from itertools import product
import json
import math
from pathlib import Path
import subprocess


ROOT = Path(__file__).resolve().parents[2]
RUN = ROOT / "artifacts/week6_e2"
OUTPUT = Path(__file__).with_suffix(".json")
FREEZE = "960436615a710d21de3439d50257f50af72a45e4"
METHODS = ("exact", "pfml", "euler")
PARAMETERS = ("theta", "mu", "sigma")
R = 2000
LIFECYCLE = {"status", "audit_sha256", "analysis_finished_utc", "completed_utc",
             "visual_review_sha256", "artifact_sha256", "execution_attempts", "resume_events",
             "raw_finished_utc", "postprocessing_error_type", "postprocessing_error",
             "termination_type", "termination_reason"}
GATES = {"hessian_condition_number": 1e10, "hessian_relative_difference": .005,
         "hessian_whitened_difference": .01, "standardized_score": .01}
checks = 0
inputs = {}


def need(condition, label):
    global checks
    if not condition:
        raise ValueError(label)
    checks += 1


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def remember(path):
    path = Path(path)
    value = digest(path)
    inputs[str(path.relative_to(ROOT))] = value
    return value


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def semantic(metadata):
    return hashlib.sha256(canonical({k: v for k, v in metadata.items() if k not in LIFECYCLE}).encode()).hexdigest()


def read_json(path):
    remember(path)
    return json.loads(Path(path).read_text())


def read_csv(path):
    remember(path)
    with Path(path).open(newline="") as stream:
        return list(csv.DictReader(stream))


def equal(a, b, label):
    a = float(a) if a not in ("", None) else float("nan")
    b = float(b) if b not in ("", None) else float("nan")
    need((math.isnan(a) and math.isnan(b)) or math.isclose(a, b, rel_tol=1e-12, abs_tol=1e-13), label)


def main():
    metadata = json.loads((RUN / "run_metadata.json").read_text())
    initial_metadata_hash = digest(RUN / "run_metadata.json")
    provenance = metadata["provenance"]
    audit = read_json(RUN / "audit.json")
    need(metadata["status"] in {"audited_pending_visual_review", "completed"}, "Run not ready for review")
    need(audit["status"] == "passed" and audit["full_completion_audit"] is True, "Full project audit required")
    need(semantic(metadata) == audit["audited_run_metadata_sha256"], "Audited scientific metadata mismatch")
    need(provenance["code_commit"] == audit["code_commit"] == FREEZE, "Recorded freeze mismatch")
    need(not provenance["uncommitted_dependency_files"] and not metadata["scientific_overrides"], "Unfrozen or overridden run")
    need(metadata["pipeline_errors"] == 0, "Pipeline exceptions require separate review")
    identity = {k: provenance[k] for k in ("design", "design_sha256", "preregistration_commit", "source_files", "environment")}
    need(hashlib.sha256(canonical(identity).encode()).hexdigest() == provenance["identity_sha256"] == audit["run_identity_sha256"], "Run identity digest")
    expected_sources = {
        "requirements.txt", "requirements-lock.txt", "preregistration.md", "docs/plan_amendment_week3_validation.md",
        "ou_irregular/__init__.py", "ou_irregular/config.py", "ou_irregular/estimators.py", "ou_irregular/ou.py", "ou_irregular/spacing.py",
        *(f"ou_irregular/inference/{x}.py" for x in ("__init__", "adapters", "contracts", "coordinates", "hessian", "wald")),
        *(f"ou_irregular/week4/{x}.py" for x in ("__init__", "config", "simulation", "summaries")),
        *(f"ou_irregular/experiments/{x}.py" for x in ("__init__", "design", "summaries", "comparisons", "provenance", "storage", "execution", "analysis", "audit", "workflow")),
        *(f"ou_irregular/week6/{x}.py" for x in ("__init__", "config", "reporting", "plots", "runner")),
        "docs/week5_mathematics.md", "docs/week6_mathematics.md", "configs/week6_e2.yaml"}
    need(set(provenance["source_files"]) == expected_sources and len(expected_sources) == 37, "Explicit dependency set")
    for name, expected in provenance["source_files"].items():
        need(remember(ROOT / name) == expected, "Current dependency bytes: " + name)
        committed = subprocess.check_output(["git", "show", f"{FREEZE}:{name}"], cwd=ROOT)
        need(hashlib.sha256(committed).hexdigest() == expected, "Freeze dependency bytes: " + name)

    design = read_json(RUN / "resolved_design.json")
    need(design == provenance["design"], "Saved/provenance design mismatch")
    need((design["name"], design["seed"], design["replications"], design["scope"], tuple(design["estimators"]))
         == ("week6_e2", 20261002, R, "headline", METHODS), "Registered E2 identity")
    cells = {c["cell_id"]: c for c in design["cells"]}
    need(len(design["cells"]) == len(cells) == 60, "Unique 60-cell design")
    need({(c["n"], c["theta_mean_gap"], c["cv"]) for c in cells.values()}
         == set(product((250, 1000, 4000), (.1, .5, 1., 2.), (0., .5, 1., 2., 3.))), "Complete E2 grid")
    for cell in cells.values():
        need((cell["true_theta"], cell["true_mu"], cell["true_sigma"], cell["floor_fraction"], cell["floor_cv_threshold"])
             == (1., 0., .5, 1e-6, 1.5), "Truth/floor settings")
    need(set(metadata["completed_cells"]) == set(cells), "Root completion cell set")
    need({p.name for p in (RUN / "cells").iterdir()} == set(cells), "Saved cell directory set")
    audited = audit["audited_artifact_sha256"]
    manifest = read_json(RUN / "table_manifest.json")
    for name in ("resolved_design.json", "table_manifest.json", *manifest):
        need(remember(RUN / name) == audited[name], "Saved evidence differs from audited bytes: " + name)
        if name in manifest:
            need(inputs[str((RUN/name).relative_to(ROOT))] == manifest[name]["sha256"], "Table manifest digest")
    summary_rows = read_csv(RUN / "parameter_summary.csv")
    summary = {(r["cell_id"], r["estimator"], r["parameter"]): r for r in summary_rows}
    need(len(summary_rows) == 540 and set(summary) == set(product(cells, METHODS, PARAMETERS)), "540 summary keys")
    outcome_rows = read_csv(RUN / "outcome_counts.csv")
    outcomes = {(r["cell_id"], r["estimator"], r["point_status"], r["interval_status"]): int(r["count"]) for r in outcome_rows}
    need(len(outcomes) == len(outcome_rows), "Duplicate exported outcome keys")
    figures = read_json(RUN / "figures/figure_manifest.json")
    need(figures["design_sha256"] == audit["design_sha256"] == provenance["design_sha256"], "Figure design binding")
    for name, expected in figures["figures"].items():
        need(remember(RUN / "figures" / name) == expected, "Derived file digest: " + name)
    failures_report_rows = read_csv(RUN / "figures/failed_cells.csv")
    failures_report = {(r["cell_id"], r["estimator"]): r for r in failures_report_rows}
    need(len(failures_report_rows) == 180 and set(failures_report) == set(product(cells, METHODS)), "180 failure-log keys")

    statuses = Counter()
    raw_outcomes = Counter()
    raw_failures = {}
    affected_sources = set()
    totals = {method: Counter() for method in METHODS}
    numerical_maxima = {scope: {} for scope in ("all_finite", "interval_valid")}
    point_maxima = {}
    failure_groups = Counter()
    partitions = 0
    for cell_number, (cid, cell) in enumerate(cells.items(), start=1):
        directory = RUN / "cells" / cid
        cell_manifest = read_json(directory / "cell_manifest.json")
        need(inputs[str((directory/"cell_manifest.json").relative_to(ROOT))] == audited[f"cells/{cid}/cell_manifest.json"], "Cell manifest audit binding")
        need(cell_manifest["cell"] == cell and cell_manifest["run_identity"] == provenance["identity_sha256"], "Cell identity")
        need(cell_manifest["status"] == "complete" and cell_manifest["completed_replications"] == R and cell_manifest["replications"] == R, "Cell incomplete")
        seen, source_hashes = set(), {}
        counts = {method: Counter() for method in METHODS}
        cell_statuses = Counter()
        next_rep = 0
        for partition in cell_manifest["partitions"]:
            start, stop = partition["start"], partition["stop"]
            need(start == next_rep and start < stop <= R, "Partition continuity")
            path = directory / partition["filename"]
            need(remember(path) == partition["sha256"] == audited[str(path.relative_to(RUN))], "Raw partition digest")
            n_rows = 0
            with gzip.open(path, "rt") as stream:
                for line in stream:
                    row = json.loads(line)
                    method, rep = row["estimator"], row["rep"]
                    need(row["cell_id"] == cid and method in METHODS and type(rep) is int and start <= rep < stop, "Raw scientific key")
                    need((rep, method) not in seen, "Duplicate raw key")
                    seen.add((rep, method)); n_rows += 1
                    pv, iv = row["point_valid"], row["interval_valid"]
                    need(type(pv) is bool and type(iv) is bool and (not iv or pv), "Joint validity nesting")
                    need(pv == (row["point_status"] == "valid_interior") and iv == (row["interval_status"] == "valid"), "Validity/status agreement")
                    need(row["source_valid"] is True and not row["source_error"] and not row["pipeline_error"], "Source/pipeline exception")
                    hashes = (row["times_sha256"], row["observations_sha256"], row.get("gap_seed"), row.get("path_seed"))
                    need(rep not in source_hashes or source_hashes[rep] == hashes, "Source sharing across methods")
                    source_hashes[rep] = hashes
                    if pv:
                        need(row["optimizer_success"] is True and all(math.isfinite(row[p]) for p in (*PARAMETERS, "nll"))
                             and row["theta"] > 0 and row["sigma"] > 0, "Malformed valid point")
                        for parameter in PARAMETERS:
                            key = (method, parameter)
                            value = abs(row[parameter])
                            if key not in point_maxima or value > point_maxima[key]["absolute_value"]:
                                point_maxima[key] = {"absolute_value": value, "value": row[parameter], "cell_id": cid, "rep": rep,
                                                     "n": cell["n"], "q": cell["theta_mean_gap"], "cv": cell["cv"]}
                    for field, threshold in GATES.items():
                        value = row[field]
                        if value is not None:
                            need(math.isfinite(value) and value >= 0, "Malformed numerical diagnostic")
                            for scope in (("all_finite", "interval_valid") if iv else ("all_finite",)):
                                key = (method, field)
                                if key not in numerical_maxima[scope] or value > numerical_maxima[scope][key]["value"]:
                                    numerical_maxima[scope][key] = {"value": value, "threshold": threshold, "cell_id": cid, "rep": rep,
                                                                  "n": cell["n"], "q": cell["theta_mean_gap"], "cv": cell["cv"]}
                        if iv:
                            need(value is not None and value <= threshold, "Available interval violates recorded gate")
                    counts[method].update(planned=1, point_valid=int(pv), interval_valid=int(iv))
                    status = (method, row["point_status"], row["interval_status"])
                    statuses[status] += 1
                    raw_outcomes[(cid, *status)] += 1
                    cell_statuses[status] += 1
                    for parameter in PARAMETERS:
                        covered = False
                        if iv:
                            lo, hi = row[parameter+"_lower"], row[parameter+"_upper"]
                            need(math.isfinite(lo) and math.isfinite(hi) and lo <= hi and math.isfinite(row[parameter+"_se"])
                                 and row[parameter+"_se"] > 0, "Malformed available interval")
                            covered = lo <= cell["true_"+parameter] <= hi
                            counts[method][parameter+"_negative_lower"] += int(lo < 0)
                        counts[method][parameter+"_covered"] += int(covered)
                    if not pv or not iv:
                        raw_failures[(cid, rep, method)] = row
                        affected_sources.add((cid, rep))
                        failure_groups[(cell["n"], cell["theta_mean_gap"], cell["cv"], *status)] += 1
            need(n_rows == partition["rows"] == 3*(stop-start), "Partition planned row count")
            next_rep = stop; partitions += 1
        need(next_rep == R and seen == set(product(range(R), METHODS)), "Complete per-cell planned keys")
        need(metadata["completed_cells"][cid] == {"cell_id": cid, "paths": R, "rows": 3*R, "pipeline_errors": 0}, "Root cell accounting")
        for method, count in counts.items():
            totals[method].update(count)
            for parameter in PARAMETERS:
                s = summary[(cid, method, parameter)]
                for field, expected in (("n_planned", R), ("n_point_valid", count["point_valid"]), ("n_interval_valid", count["interval_valid"]),
                                        ("n_covered", count[parameter+"_covered"]), ("negative_lower_count", count[parameter+"_negative_lower"])):
                    equal(s[field], expected, "Raw/summary count: " + field)
                equal(s["point_availability"], count["point_valid"]/R, "Point availability")
                equal(s["interval_availability"], count["interval_valid"]/R, "Interval availability")
                equal(s["coverage_operational"], count[parameter+"_covered"]/R, "Operational coverage")
                equal(s["coverage_conditional"], count[parameter+"_covered"]/count["interval_valid"] if count["interval_valid"] else None, "Conditional coverage")
            f = failures_report[(cid, method)]
            for field, expected in (("n_planned", R), ("n_point_valid", count["point_valid"]), ("n_interval_valid", count["interval_valid"]),
                                    ("n_point_invalid", R-count["point_valid"]), ("n_interval_unavailable", R-count["interval_valid"]),
                                    ("n_point_valid_without_interval", count["point_valid"]-count["interval_valid"]),
                                    ("failure_fraction", (R-count["interval_valid"])/R)):
                equal(f[field], expected, "Failure report count: " + field)
            need(f["has_failure"] == str(count["interval_valid"] < R), "Failure flag")
            composition = {(method, x["point_status"], x["interval_status"]): x["count"] for x in json.loads(f["status_counts_json"])}
            need(composition == {k: v for k, v in cell_statuses.items() if k[0] == method}, "Complete failure status composition")
        if cell_number % 10 == 0:
            print(f"Reviewed {cell_number}/60 saved cells; {partitions} immutable partitions", flush=True)

    need(raw_outcomes == outcomes, "Raw/exported outcome status counts")
    with gzip.open(RUN / "failed_fits.csv.gz", "rt", newline="") as stream:
        exported_failures = list(csv.DictReader(stream))
    exported_keys = [(x["cell_id"], int(x["rep"]), x["estimator"]) for x in exported_failures]
    need(len(exported_keys) == len(set(exported_keys)) and set(exported_keys) == set(raw_failures), "Exact failed-fit export key set")
    for exported, key in zip(exported_failures, exported_keys):
        raw = raw_failures[key]
        need(set(exported) == set(raw), "Failed-fit export schema")
        for field, value in raw.items():
            if value is None:
                need(exported[field] == "", "Missing candidate field changed")
            elif isinstance(value, bool):
                need(exported[field] == str(value), "Failed-fit Boolean changed")
            elif isinstance(value, (float, int)):
                need(float(exported[field]) == value, "Failed-fit numeric candidate changed: " + field)
            else:
                need(exported[field] == value, "Failed-fit diagnostic/source field changed: " + field)
    need(sum(c["planned"] for c in totals.values()) == audit["planned_fit_records"] == 360000, "Full fit count")
    need(sum(c["point_valid"] for c in totals.values()) == audit["point_valid_records"], "Point-valid count")
    need(sum(c["interval_valid"] for c in totals.values()) == audit["interval_valid_records"], "Interval-valid count")
    need(partitions == audit["partitions_verified"] == 240, "Partition count")
    need(audit["planned_paths"] == 120000 and audit["replayed_source_paths"] == 120, "Project audit scope")
    for name, expected in inputs.items():
        need(digest(ROOT/name) == expected, "Review input changed: " + name)
    final_metadata = json.loads((RUN / "run_metadata.json").read_text())
    need(semantic(final_metadata) == semantic(metadata), "Scientific metadata changed during review")

    reviewed_pngs = {name: figures["figures"][name] for name in (
        "F5_operational_n1000.png", "F6_conditional_n1000.png", "F6_availability_n1000.png")}
    anomalous_valid_points = [{k: row[k] for k in (
        "cell_id", "n", "theta_mean_gap", "cv", "rep", "estimator", "theta", "mu", "sigma", "nll", "point_status",
        "point_reason", "interval_status", "interval_reason", "hessian_relative_difference", "hessian_whitened_difference",
        "hessian_condition_number", "standardized_score", "point_diagnostics_json", "interval_diagnostics_json")}
        for row in raw_failures.values() if row["point_valid"]]
    result = {"status": "passed", "reviewer": "independent provenance/failure/n1000 reviewer",
        "review_kind": "read-only saved-data verification; no source generation, refitting or interval substitution",
        "freeze_commit": FREEZE, "run_identity_sha256": provenance["identity_sha256"], "source_dependency_files_verified": 37,
        "cells_verified": 60, "planned_paths": 120000, "planned_fit_records": 360000, "partitions_verified": partitions,
        "checks_passed": checks, "method_accounting": {m: dict(c) for m, c in totals.items()},
        "status_counts": [{"estimator": k[0], "point_status": k[1], "interval_status": k[2], "count": v} for k, v in sorted(statuses.items())],
        "failed_fit_rows": len(raw_failures), "affected_source_paths": len(affected_sources),
        "failure_groups": [{"n": k[0], "q": k[1], "cv": k[2], "estimator": k[3], "point_status": k[4], "interval_status": k[5], "count": v} for k, v in sorted(failure_groups.items())],
        "all_failed_candidate_fields_preserved": True,
        "valid_point_interval_failures": anomalous_valid_points,
        "numerical_maxima": {scope: [{"estimator": key[0], "field": key[1], **value} for key, value in sorted(values.items())] for scope, values in numerical_maxima.items()},
        "point_absolute_maxima": [{"estimator": key[0], "parameter": key[1], **value} for key, value in sorted(point_maxima.items())],
        "visual_review": {"status": "passed", "reviewed_png_sha256": reviewed_pngs,
            "findings": "Original-resolution n1000 figures inspected: complete4x5 grids, clearREGISTERED E2/R2000 labels, legible annotations/captions, no clipping/overlap; mu omits relative-bias classes. AllS2000 here, so these PNGs do not exercise unavailable-interval NA."},
        "week7_follow_up": [
            "Preserve and investigate the n250/q2 regular AR-slope<=0 boundary paths; check paired exact/PFML disposition without replacing these outcomes.",
            "Investigate the additional PFML boundaries at n250/q2/CV0.5 and exact rep465's weak-curvature valid point/whitened-Hessian rejection under a separately dated exploratory protocol.",
            "Retain point-valid outliers in point metrics and unavailable intervals as operational noncoverage; neither is an automatic software-defect diagnosis.",
            "The earlier validation root-status reset remains unexplained; this passing E2 review does not establish its cause."],
        "limitations": ["Counts and coverage were independently recomputed from all raw rows; bias/RMSE and covariance matrices were not recomputed in this review.",
                        "No independent source replay was performed here; the project audit separately records120 replayed source paths.",
                        "Only the three n1000 PNGs were visually inspected by this reviewer; no claim about other PNGs or personal theory mastery."],
        "root_metadata_initial_sha256": initial_metadata_hash, "root_metadata_final_sha256": digest(RUN/"run_metadata.json"),
        "scientific_metadata_unchanged": True, "input_sha256": inputs, "reviewer_script_sha256": digest(Path(__file__))}
    payload = (json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+"\n").encode()
    if OUTPUT.exists():
        need(OUTPUT.read_bytes() == payload, "Refusing to replace different review evidence")
    else:
        with OUTPUT.open("xb") as stream:
            stream.write(payload)
    print(json.dumps({k: result[k] for k in ("status", "checks_passed", "planned_fit_records", "failed_fit_rows", "affected_source_paths", "partitions_verified")}, sort_keys=True))
    print(str(OUTPUT.relative_to(ROOT)))


if __name__ == "__main__":
    main()
