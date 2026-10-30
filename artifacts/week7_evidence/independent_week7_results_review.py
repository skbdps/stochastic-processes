#!/usr/bin/env python3
"""Read-only, stdlib-only review of the two fixed full Week 7 archives.

Prepared before inspecting the full-study results. Run only after the study
owner confirms that both archives are complete. No project imports, random
numbers, source reconstruction, optimization, inference, or archive writes.

This review targets denominator closure, failure retention, independent scalar
arithmetic, and covariance-aware epsilon contrasts. It complements the existing
full numerical audits; it does not independently certify their source/fit replay,
all Hessian identities, visual quality, or historical runtime-incident causes.
"""

import argparse
from collections import Counter
from datetime import datetime, timezone
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import sys


METHODS = ("exact", "pfml", "euler")
PARAMETERS = ("theta", "mu", "sigma")
TRUTHS = (1.0, 0.0, 0.5)
EPSILONS = (1e-8, 1e-6, 1e-4)
Z = 1.959963984540054
RTOL, ATOL = 2e-10, 2e-12
COUPLING = "week7_paired_floor_v1_shared_gamma_stationary_x0_and_standard_normal_innovations"
REPO = Path("/tmp/ou_stochastic_recovery_20261030")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def digest(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def unique_object(pairs):
    result = {}
    for key, value in pairs:
        require(key not in result, "duplicate JSON key: " + key)
        result[key] = value
    return result


def reject_constant(value):
    raise ValueError("nonstandard JSON numeric constant: " + value)


def loads(value):
    return json.loads(value, object_pairs_hook=unique_object, parse_constant=reject_constant)


def read_json(path):
    return loads(Path(path).read_bytes())


def snapshot(root):
    require(root.is_dir() and not root.is_symlink(), "archive is missing or a symlink")
    result = {}
    for path in sorted(root.rglob("*")):
        require(not path.is_symlink(), "archive contains a symlink")
        if path.is_file() and path.name != ".execution.lock":
            result[str(path.relative_to(root))] = digest(path)
    return result


def cell(n, q, cv, epsilon):
    value = dict(n=n, theta_mean_gap=float(q), mean_gap=float(q), cv=float(cv),
                 true_theta=1.0, true_mu=0.0, true_sigma=0.5,
                 floor_fraction=float(epsilon), floor_cv_threshold=1.5)
    return dict(cell_id=hashlib.sha256(canonical(value).encode()).hexdigest(), **value)


def design_for(kind):
    if kind == "followup":
        cells = [cell(250, 2, cv, 1e-6) for cv in (0, 0.5)]
        name, seed, replications = "week7_e3_followup", 20261007, 10000
    else:
        cells = [cell(1000, 0.5, cv, epsilon)
                 for cv in (0, 1.5, 2, 3) for epsilon in EPSILONS]
        name, seed, replications = "week7_floor_sensitivity", 20261009, 2000
    return dict(name=name, seed=seed, replications=replications, scope="exploratory",
                estimators=list(METHODS), cells=cells)


def csv_rows(path):
    opener = gzip.open if Path(path).suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as stream:
        reader = csv.reader(stream)
        header = next(reader, [])
        require(header and len(set(header)) == len(header), "missing or duplicate CSV header")
        rows = []
        for values in reader:
            require(len(values) == len(header), "CSV row width changed")
            rows.append(dict(zip(header, values)))
    return header, rows


def indexed(rows, keys):
    result = {}
    for row in rows:
        key = tuple(row[name] for name in keys)
        require(key not in result, "duplicate table key: " + repr(key))
        result[key] = row
    return result


def finite(value, label):
    require(type(value) in (int, float) and math.isfinite(value), "invalid numeric value: " + label)
    return float(value)


def same(saved, expected, label):
    if expected is None:
        require(saved == "", "missing value changed: " + label)
    elif type(expected) is bool:
        require(saved == str(expected), "Boolean changed: " + label)
    elif isinstance(expected, str):
        require(saved == expected, "text changed: " + label)
    elif type(expected) is int:
        require(saved == str(expected), "integer count changed: " + label)
    else:
        require(saved != "", "numeric value missing: " + label)
        actual = float(saved)
        require(math.isfinite(actual) and math.isfinite(expected)
                and math.isclose(actual, expected, rel_tol=RTOL, abs_tol=ATOL),
                f"arithmetic differs: {label}: saved={saved}; independent={expected!r}")


def compare_fields(saved, expected, label):
    for key, value in expected.items():
        require(key in saved, "missing table field: " + label + "." + key)
        same(saved[key], value, label + "." + key)


def mean(values):
    return math.fsum(values) / len(values) if values else None


def variance(values):
    if len(values) < 2:
        return None
    center = mean(values)
    return math.fsum((value - center) ** 2 for value in values) / (len(values) - 1)


def mcse(values):
    value = variance(values)
    return math.sqrt(value / len(values)) if value is not None else None


def scaled(value, denominator):
    return value / denominator if value is not None and denominator else None


def interval(value, standard_error):
    if value is None or standard_error is None:
        return None, None
    return value - Z * standard_error, value + Z * standard_error


def wilson(count, total):
    if not total:
        return None, None
    p = count / total
    weight = 1 + Z * Z / total
    center = (p + Z * Z / (2 * total)) / weight
    radius = Z * math.sqrt((p * (1 - p) + Z * Z / (4 * total)) / total) / weight
    return max(0.0, min(p, center - radius)), min(1.0, max(p, center + radius))


def binomial(count, total):
    p = count / total if total else None
    se = math.sqrt(p * (1 - p) / total) if total else None
    return p, se, *wilson(count, total)


def contrast_fields(name, values):
    value, se = mean(values), mcse(values)
    lower, upper = interval(value, se)
    return {name: value, name + "_mcse": se,
            name + "_mc_lower": lower, name + "_mc_upper": upper}


def point_fields(errors, truth):
    squares = [value * value for value in errors]
    bias, mse = mean(errors), mean(squares)
    rmse = math.sqrt(mse) if mse is not None else None
    se_bias, se_squares = mcse(errors), mcse(squares)
    se_rmse = scaled(se_squares, 2 * rmse) if rmse is not None else None
    lower, upper = interval(bias, se_bias)
    return dict(bias=bias, rmse=rmse, mcse_bias=se_bias, mcse_rmse=se_rmse,
                bias_mc_lower=lower, bias_mc_upper=upper,
                relative_bias=scaled(bias, abs(truth)),
                relative_mcse_bias=scaled(se_bias, abs(truth)),
                relative_rmse=scaled(rmse, abs(truth)),
                relative_mcse_rmse=scaled(se_rmse, abs(truth)),
                relative_bias_mc_lower=scaled(lower, abs(truth)),
                relative_bias_mc_upper=scaled(upper, abs(truth)))


def summary_fields(rows, parameter_index, total):
    truth = TRUTHS[parameter_index]
    points = [row for row in rows if row["point"]]
    intervals = [row for row in rows if row["interval"]]
    count = sum(row["hits"][parameter_index] for row in rows)
    expected = dict(truth=truth, n_planned=total, n_point_valid=len(points),
                    n_interval_valid=len(intervals), n_covered=count,
                    negative_lower_count=sum(row["negative"][parameter_index] for row in intervals))
    for name, count_ in (("point", len(points)), ("interval", len(intervals))):
        p, se, lo, hi = binomial(count_, total)
        expected.update({name + "_availability": p, name + "_availability_mcse": se,
                         name + "_availability_wilson_lower": lo,
                         name + "_availability_wilson_upper": hi})
    for name, denominator in (("conditional", len(intervals)), ("operational", total)):
        p, se, lo, hi = binomial(count, denominator)
        expected.update({"coverage_" + name: p, name + "_mcse": se,
                         name + "_wilson_lower": lo, name + "_wilson_upper": hi})
    expected.update(point_fields([row["errors"][parameter_index] for row in points], truth))
    expected.update(point_metric_scope="conditional_on_valid_point_fit",
                    primary_coverage_scope="all_planned_paths_unavailable_counts_as_not_covering")
    blo, bhi = expected["relative_bias_mc_lower"], expected["relative_bias_mc_upper"]
    clo, chi = expected["operational_wilson_lower"], expected["operational_wilson_upper"]
    if truth == 0:
        classification = "not_applicable_zero_truth"
    elif chi < 0.80 or (blo is not None and (blo > 0.10 or bhi < -0.10)):
        classification = "supported_breach"
    elif clo >= 0.80 and blo is not None and blo >= -0.10 and bhi <= 0.10:
        classification = "within_criteria"
    else:
        classification = "unresolved"
    expected["practical_classification"] = classification
    return expected


def paired_fields(a, b, parameter_index):
    total, truth = len(a), TRUTHS[parameter_index]
    common = [(x, y) for x, y in zip(a, b) if x["point"] and y["point"]]
    ea = [x["errors"][parameter_index] for x, y in common]
    eb = [y["errors"][parameter_index] for x, y in common]
    sa, sb = [v * v for v in ea], [v * v for v in eb]
    ma, mb, ba, bb = mean(sa), mean(sb), mean(ea), mean(eb)
    ra = math.sqrt(ma) if ma is not None else None
    rb = math.sqrt(mb) if mb is not None else None
    expected = dict(n_planned=total, n_point_valid_a=sum(x["point"] for x in a),
                    n_point_valid_b=sum(x["point"] for x in b),
                    n_interval_valid_a=sum(x["interval"] for x in a),
                    n_interval_valid_b=sum(x["interval"] for x in b),
                    n_common_point_valid=len(common),
                    n_common_interval_valid=sum(x["interval"] and y["interval"] for x, y in zip(a, b)),
                    bias_a_common=ba, bias_b_common=bb, mse_a_common=ma, mse_b_common=mb,
                    rmse_a_common=ra, rmse_b_common=rb,
                    relative_bias_a_common=scaled(ba, abs(truth)),
                    relative_bias_b_common=scaled(bb, abs(truth)),
                    relative_rmse_a_common=scaled(ra, abs(truth)),
                    relative_rmse_b_common=scaled(rb, abs(truth)))
    covariance = (math.fsum((x - ma) * (y - mb) for x, y in zip(sa, sb)) / (len(sa) - 1)
                  if len(sa) > 1 else None)
    expected["squared_error_covariance_common"] = covariance
    expected.update(contrast_fields("bias_difference", [x - y for x, y in zip(ea, eb)]))
    expected.update(contrast_fields("mse_difference", [x - y for x, y in zip(sa, sb)]))
    for suffix in ("", "_mcse", "_mc_lower", "_mc_upper"):
        expected["relative_bias_difference" + suffix] = scaled(expected["bias_difference" + suffix], abs(truth))
    difference = ra - rb if ra is not None and rb is not None else None
    se_difference = mcse([x / (2 * ra) - y / (2 * rb) for x, y in zip(sa, sb)]) if ra and rb else None
    lo, hi = interval(difference, se_difference)
    expected.update(rmse_difference=difference, rmse_difference_mcse=se_difference,
                    rmse_difference_mc_lower=lo, rmse_difference_mc_upper=hi)
    ratio = scaled(ra, rb)
    log_se = mcse([0.5 * (x / ma - y / mb) for x, y in zip(sa, sb)]) if ma and mb else None
    expected.update(rmse_ratio=ratio, rmse_ratio_log_mcse=log_se,
                    rmse_ratio_mcse=ratio * log_se if ratio is not None and log_se is not None else None,
                    rmse_ratio_mc_lower=ratio * math.exp(-Z * log_se) if ratio is not None and log_se is not None else None,
                    rmse_ratio_mc_upper=ratio * math.exp(Z * log_se) if ratio is not None and log_se is not None else None,
                    rmse_ratio_uncertainty="delta_method_log_ratio_95pct")
    hits = [(x["hits"][parameter_index], y["hits"][parameter_index]) for x, y in zip(a, b)]
    ha, hb = sum(x for x, y in hits), sum(y for x, y in hits)
    both = sum(x and y for x, y in hits)
    expected.update(n_covered_a=ha, n_covered_b=hb, n_both_cover=both,
                    n_a_only_cover=ha - both, n_b_only_cover=hb - both,
                    n_neither_cover=total - ha - hb + both,
                    operational_coverage_a=ha / total, operational_coverage_b=hb / total)
    expected.update(contrast_fields("operational_coverage_difference", [int(x) - int(y) for x, y in hits]))
    return expected


def review(kind, root):
    expected_design = design_for(kind)
    total = expected_design["replications"]
    cells = expected_design["cells"]
    before = snapshot(root)
    metadata, audit = read_json(root / "run_metadata.json"), read_json(root / "audit.json")
    require(metadata["status"] in {"audited_pending_visual_review", "completed"}, "archive is not fully audited")
    require(metadata["execution_kind"] == kind and metadata["require_frozen"] is True, "wrong execution identity")
    require(canonical(read_json(root / "resolved_design.json")) == canonical(expected_design), "full design differs")
    require(audit["status"] == "passed" and audit["full_completion_audit"] is True
            and audit["replay_performed"] is True, "saved full audit did not pass")
    require(metadata["audit_sha256"] == before["audit.json"], "root audit binding differs")
    identity = metadata["provenance"]["identity_sha256"]
    require(audit["run_identity_sha256"] == identity, "audit run identity differs")
    design_hash = hashlib.sha256(canonical(expected_design).encode()).hexdigest()
    require(audit["design_sha256"] == design_hash, "audit design identity differs")
    paths, independent = len(cells) * total, (4 if kind == "floor" else len(cells)) * total
    policy = metadata["source_policy"]
    require(policy["planned_source_realizations"] == paths
            and policy["planned_independent_coupling_groups"] == independent
            and policy["planned_fit_records"] == 3 * paths
            and policy["root_seed"] == expected_design["seed"]
            and policy["stationary_initialization"] is True
            and policy["renormalize_floored_gaps"] is False, "source policy or planned counts differ")
    require(policy["name"] == (COUPLING if kind == "floor" else "registered_ou_source_v1"), "source policy name differs")
    manifest = read_json(root / "table_manifest.json")
    for name, entry in manifest.items():
        require(before[name] == entry["sha256"], "saved table digest differs: " + name)
    _, summaries = csv_rows(root / "parameter_summary.csv")
    summary_index = indexed(summaries, ("cell_id", "estimator", "parameter"))
    expected_keys = {(c["cell_id"], method, p) for c in cells for method in METHODS for p in PARAMETERS}
    require(set(summary_index) == expected_keys, "parameter-summary Cartesian keys differ")
    require(manifest["parameter_summary.csv"]["rows"] == len(summaries) == 9 * len(cells), "parameter-summary counts differ")
    require({p.name for p in (root / "cells").iterdir()} == {c["cell_id"] for c in cells}, "cell set differs")
    failures, outcomes, failure_reasons = {}, Counter(), Counter()
    groups, reviewed_summaries, counts_by_method = {}, [], []
    raw_schema, partitions, point_count, interval_count = None, 0, 0, 0
    regular_control = {}
    for c in cells:
        cid, cv, epsilon = c["cell_id"], c["cv"], c["floor_fraction"]
        directory = root / "cells" / cid
        state = read_json(directory / "cell_manifest.json")
        require(state["status"] == "complete" and state["completed_replications"] == total
                and state["replications"] == total and state["run_identity"] == identity
                and canonical(state["cell"]) == canonical(c), "incomplete or changed cell manifest")
        anchor = cell(c["n"], c["theta_mean_gap"], cv, 1e-6)["cell_id"] if kind == "floor" else cid
        require(policy["actual_to_seed_anchor"][cid] == anchor, "source anchor differs")
        words = [int(anchor[i:i + 8], 16) for i in range(0, 64, 8)]
        by_method = {method: [] for method in METHODS}
        next_rep, filenames = 0, {"cell_manifest.json"}
        for part in state["partitions"]:
            start, stop = part["start"], part["stop"]
            filename = f"part-{start:06d}-{stop:06d}.jsonl.gz"
            require(type(start) is int and type(stop) is int and start == next_rep
                    and start < stop <= total and part["filename"] == filename
                    and part["rows"] == 3 * (stop - start), "partition range/count differs")
            require(before[f"cells/{cid}/{filename}"] == part["sha256"], "partition hash differs")
            filenames.add(filename)
            count, source_signature = 0, None
            with gzip.open(directory / filename, "rt", encoding="utf-8") as stream:
                for line in stream:
                    row = loads(line)
                    rep, method = start + count // 3, METHODS[count % 3]
                    require(row["rep"] == rep and type(row["rep"]) is int and row["estimator"] == method,
                            "raw records do not equal ordered complete replication/method grid")
                    require(all(type(row[k]) is type(v) and row[k] == v for k, v in c.items()), "raw cell identity differs")
                    if raw_schema is None:
                        raw_schema = sorted(row)
                    require(sorted(row) == raw_schema, "raw schema differs within study")
                    require(row["source_valid"] is True and row["source_error"] == ""
                            and row["pipeline_error"] == "", "source/pipeline error needs implementation review")
                    require(row["source_seed"] == expected_design["seed"], "raw seed differs")
                    address = [4, *words, rep]
                    for field, value in (("source_spawn_key", address), ("gap_spawn_key", address + [0]),
                                         ("path_spawn_key", address + [1])):
                        require(loads(row[field]) == value, "saved source address differs")
                    signature = tuple(row[field] for field in ("times_sha256", "observations_sha256", "source_spawn_key",
                                                                "gap_spawn_key", "path_spawn_key"))
                    for hashed in signature[:2]:
                        require(type(hashed) is str and len(hashed) == 64
                                and all(ch in "0123456789abcdef" for ch in hashed), "invalid source hash")
                    if method == "exact":
                        source_signature = signature
                    require(signature == source_signature, "within-source estimators do not share arrays/address")
                    if kind == "floor":
                        require(row["actual_cell_id"] == cid and row["seed_anchor_cell_id"] == anchor
                                and row["seed_anchor_floor_fraction"] == 1e-6
                                and row["source_coupling_policy"] == COUPLING, "paired floor evidence differs")
                    point, available = row["point_valid"], row["interval_valid"]
                    require(type(point) is bool and type(available) is bool and type(row["optimizer_success"]) is bool,
                            "non-Boolean validity flag")
                    require(point == (row["point_status"] == "valid_interior")
                            and available == (row["interval_status"] == "valid")
                            and (not point or row["optimizer_success"])
                            and (not available or point), "inconsistent inference status/eligibility")
                    errors, hits, negative = [], [], []
                    if point:
                        finite(row["nll"], "valid NLL")
                    for p, truth in zip(PARAMETERS, TRUTHS):
                        estimate = finite(row[p], p) if point else None
                        errors.append(estimate - truth if point else None)
                        if available:
                            lower, upper = finite(row[p + "_lower"], p + " lower"), finite(row[p + "_upper"], p + " upper")
                            require(lower <= upper and finite(row[p + "_se"], p + " SE") > 0, "invalid eligible interval")
                            hits.append(lower <= truth <= upper)
                            negative.append(lower < 0)
                        else:
                            hits.append(False)
                            negative.append(False)
                    compact = dict(point=point, interval=available, errors=tuple(errors), hits=tuple(hits), negative=tuple(negative))
                    by_method[method].append(compact)
                    if kind == "floor" and cv == 0:
                        control = (signature, compact, row["point_status"], row["interval_status"],
                                   tuple(row[p + suffix] for p in PARAMETERS for suffix in ("", "_se", "_lower", "_upper")))
                        key = (rep, method)
                        if epsilon == EPSILONS[0]:
                            regular_control[key] = control
                        else:
                            require(control == regular_control[key], "CV=0 shared negative control differs across floors")
                    outcomes[(cid, method, row["point_status"], row["interval_status"])] += 1
                    if not point or not available:
                        failures[(cid, str(rep), method)] = row
                        failure_reasons[(cid, method, row["point_status"], row["interval_status"],
                                         row["point_reason"], row["interval_reason"])] += 1
                    point_count += point
                    interval_count += available
                    count += 1
            require(count == part["rows"], "partition row closure differs")
            next_rep, partitions = stop, partitions + 1
        require(next_rep == total and {p.name for p in directory.iterdir()} == filenames, "cell range/layout closure differs")
        for method, rows in by_method.items():
            require(len(rows) == total, "method denominator differs")
            if kind == "floor":
                groups[(cv, epsilon, method)] = rows
            counts_by_method.append(dict(cell_id=cid, cv=cv, floor_fraction=epsilon, estimator=method,
                                         n_planned=total, n_point_valid=sum(r["point"] for r in rows),
                                         n_interval_valid=sum(r["interval"] for r in rows)))
            for i, p in enumerate(PARAMETERS):
                expected = summary_fields(rows, i, total)
                key = (cid, method, p)
                compare_fields(summary_index[key], expected, "summary" + repr(key))
                reviewed_summaries.append(dict(cell_id=cid, cv=cv, floor_fraction=epsilon,
                                               estimator=method, parameter=p, **expected))
    header, exported = csv_rows(root / "failed_fits.csv.gz")
    require(header == raw_schema, "failure export does not retain complete exact ordered raw header")
    failed_index = indexed(exported, ("cell_id", "rep", "estimator"))
    require(set(failed_index) == set(failures) and len(failures) == manifest["failed_fits.csv.gz"]["rows"],
            "failure export membership/count differs")
    for key, original in failures.items():
        compare_fields(failed_index[key], original, "failure" + repr(key))
    _, outcome_rows = csv_rows(root / "outcome_counts.csv")
    outcome_index = indexed(outcome_rows, ("cell_id", "estimator", "point_status", "interval_status"))
    require(set(outcome_index) == set(outcomes), "outcome-type table membership differs")
    for key, count in outcomes.items():
        lo, hi = wilson(count, total)
        compare_fields(outcome_index[key], dict(count=count, n_planned=total, fraction=count / total,
                                               wilson_lower=lo, wilson_upper=hi), "outcome" + repr(key))
    require(manifest["outcome_counts.csv"]["rows"] == len(outcome_rows), "outcome manifest count differs")
    paired_report = []
    if kind == "floor":
        _, paired = csv_rows(root / "epsilon_comparisons.csv")
        paired_index = indexed(paired, ("cell_id", "estimator", "parameter"))
        wanted = {(c["cell_id"], method, p) for c in cells if c["floor_fraction"] != 1e-6
                  for method in METHODS for p in PARAMETERS}
        require(set(paired_index) == wanted and len(paired) == 72, "epsilon contrast grid differs")
        epsilon_manifest = read_json(root / "epsilon_table_manifest.json")
        require(epsilon_manifest["table"]["sha256"] == before["epsilon_comparisons.csv"]
                and epsilon_manifest["table"]["rows"] == 72, "epsilon manifest binding differs")
        for c in cells:
            if c["floor_fraction"] == 1e-6:
                continue
            cv, epsilon = c["cv"], c["floor_fraction"]
            reference = cell(1000, 0.5, cv, 1e-6)["cell_id"]
            for method in METHODS:
                a, b = groups[(cv, epsilon, method)], groups[(cv, 1e-6, method)]
                for i, p in enumerate(PARAMETERS):
                    key = (c["cell_id"], method, p)
                    expected = paired_fields(a, b, i)
                    expected.update(reference_cell_id=reference, seed_anchor_cell_id=reference,
                                    floor_fraction=epsilon, reference_floor_fraction=1e-6,
                                    comparison_type="paired_floor_shared_standardized_innovations",
                                    difference_direction="epsilon_minus_registered_floor",
                                    ratio_direction="epsilon_over_registered_floor",
                                    point_subset="intersection_of_valid_point_fits",
                                    coverage_subset="all_planned_paths_unavailable_counts_as_not_covering")
                    compare_fields(paired_index[key], expected, "epsilon" + repr(key))
                    paired_report.append(dict(cell_id=c["cell_id"], cv=cv, estimator=method, parameter=p, **expected))
        require(audit["planned_source_realizations"] == paths
                and audit["planned_independent_coupling_groups"] == independent
                and audit["source_replay_performed"] is True
                and audit["source_paths_verified"] == paths
                and audit["source_array_hashes_verified"] == 2 * paths, "saved full-source audit counts differ")
        require(audit["tables_verified"]["epsilon_comparisons.csv"] == 72, "audit epsilon count differs")
    observed = dict(cells_verified=len(cells), partitions_verified=partitions, planned_paths=paths,
                    planned_fit_records=3 * paths, point_valid_records=point_count,
                    interval_valid_records=interval_count, retained_failure_records=len(failures), pipeline_errors=0)
    for name, value in observed.items():
        require(audit[name] == value, "saved audit accounting differs: " + name)
    require(audit["tables_verified"]["parameter_summary.csv"] == len(summaries)
            and audit["tables_verified"]["outcome_counts.csv"] == len(outcome_rows), "saved audit table counts differ")
    require(snapshot(root) == before, "archive changed during independent read-only review")
    return dict(status="passed", archive=str(root), archive_status=metadata["status"],
                script_scope="saved_scalar_arithmetic_and_accounting_no_source_or_fit_replay",
                code_commit=metadata["provenance"]["code_commit"], run_identity=identity,
                design_sha256=design_hash, artifact_sha256_before_and_after=before,
                planned_independent_coupling_groups=independent, **observed,
                parameter_summary_rows_checked=len(reviewed_summaries), epsilon_rows_checked=len(paired_report),
                availability_by_cell_and_method=counts_by_method,
                failure_types=[dict(cell_id=k[0], estimator=k[1], point_status=k[2], interval_status=k[3],
                                    point_reason=k[4], interval_reason=k[5], count=v)
                               for k, v in sorted(failure_reasons.items())],
                parameter_summary_independent=reviewed_summaries, epsilon_independent=paired_report)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--followup", type=Path, default=REPO / "artifacts/week7_followup")
    parser.add_argument("--floor", type=Path, default=REPO / "artifacts/week7_floor")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    roots = {"followup": args.followup.resolve(), "floor": args.floor.resolve()}
    output = args.output.resolve()
    require(not output.exists(), "review output already exists; choose a fresh evidence filename")
    require(all(output != root and not output.is_relative_to(root) for root in roots.values()),
            "review output must be outside both immutable archives")
    require(output.parent.is_dir(), "review output parent must already exist")
    script = Path(__file__).resolve()
    report = dict(schema_version=1, review="independent_week7_scalar_accounting_v1",
                  script=str(script), script_sha256=digest(script),
                  started_utc=datetime.now(timezone.utc).isoformat(), python=sys.version,
                  numerical_comparison=dict(relative_tolerance=RTOL, absolute_tolerance=ATOL),
                  new_sources=0, new_fits=0, archive_writes=0,
                  limitations=["Saved-array identity/address checks do not regenerate observations.",
                               "Saved full-audit flags/counts are checked but its source and fit replay is not repeated.",
                               "Hessian identities, all other table fields, figures, and historical runtime causes remain outside this review.",
                               "Conditional point estimates and common-valid paired contrasts do not erase retained failures.",
                               "Monte Carlo intervals are pointwise simulation uncertainty, not replacement parameter intervals.",
                               "Paired floor variants share innovations; they are not one continuous path on three grids.",
                               "A passed read-only review does not seal an archive awaiting visual review."], studies={})
    try:
        for kind, root in roots.items():
            report["studies"][kind] = review(kind, root)
        require(digest(script) == report["script_sha256"], "review script changed during execution")
        report["status"] = "passed"
    except Exception as exc:
        report.update(status="failed", error_type=type(exc).__name__, error=str(exc))
    report["finished_utc"] = datetime.now(timezone.utc).isoformat()
    with output.open("x", encoding="utf-8") as stream:
        json.dump(report, stream, sort_keys=True, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps(dict(status=report["status"], script_sha256=report["script_sha256"],
                          output=str(output), output_sha256=digest(output),
                          studies={kind: {key: value[key] for key in (
                              "planned_paths", "planned_fit_records", "point_valid_records", "interval_valid_records",
                              "retained_failure_records", "parameter_summary_rows_checked", "epsilon_rows_checked")}
                                   for kind, value in report["studies"].items()},
                          error=report.get("error")), sort_keys=True))
    return 0 if report["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
