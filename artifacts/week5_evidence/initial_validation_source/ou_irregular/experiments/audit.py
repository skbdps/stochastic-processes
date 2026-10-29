"""Independent audit of complete Monte Carlo artifacts and fixed-address replay.

The auditor does not call the reporting summarizers or comparison routines.
It independently derives every published table from one saved cell at a time.
The replay subset is fixed before seeing results: replication 0 and R-1 in
EVERY cell, with all three estimators. Replay is read-only and does not replace
an unavailable fit or produce new headline observations.

Every saved fine-step Hessian is checked against its covariance, coordinate
map, score and endpoints. Coarse Hessians and full paths are not stored in raw
partitions: their independent numerical reproduction is therefore checked on
the predetermined replay subset, not falsely claimed for every saved row.
"""
from collections import Counter
from dataclasses import replace
import csv
import gzip
import hashlib
from itertools import combinations
import json
from pathlib import Path
import re
import subprocess

import numpy as np
import pandas as pd
from scipy.special import gammainc, gammaincc

from ..inference.adapters import OuLikelihood, Week3Estimator
from ..inference.hessian import FiniteDifferenceHessian
from ..inference.wald import WaldInference
from ..week4.simulation import OuSampleSource
from .provenance import CORE_FILES, ROOT, canonical, file_digest, verify


PARAMETERS = ("theta", "mu", "sigma")
METHODS = ("exact", "pfml", "euler")
Z = 1.959963984540054
ARITHMETIC_RTOL = 2e-10
ARITHMETIC_ATOL = 2e-12
REPLAY_RTOL = 2e-10
REPLAY_ATOL = 2e-12
EXPECTED_PROCEDURE = {
    "level": .95, "relative_step": .001, "max_relative_difference": .005,
    "max_whitened_difference": .01, "max_condition_number": 1e10,
    "max_standardized_score": .01, "interval_scale": "natural_symmetric_unclipped",
    "covariance": "full_conditional_nll_inverse_hessian_model_based",
    "source_address": "SeedSequence(root,spawn_key=(4,*8_cell_words,rep)).spawn(2)",
}
AUDIT_FILES = (
    "ou_irregular/experiments/storage.py", "ou_irregular/experiments/execution.py",
    "ou_irregular/experiments/analysis.py", "ou_irregular/experiments/audit.py",
)
TABLE_KEYS = {
    "parameter_summary.csv": ("cell_id", "estimator", "parameter"),
    "regular_comparisons.csv": ("cell_id", "estimator", "parameter"),
    "paired_comparisons.csv": ("cell_id", "estimator_a", "estimator_b", "parameter"),
    "spacing_summary.csv": ("cell_id",),
    "outcome_counts.csv": ("cell_id", "n", "theta_mean_gap", "cv", "estimator", "point_status", "interval_status"),
    "numerical_diagnostics.csv": ("cell_id", "estimator"),
}
SOURCE_FIELDS = (
    "source_seed", "source_spawn_key", "gap_spawn_key", "path_spawn_key",
    "times_sha256", "observations_sha256", "realized_mean_gap", "realized_cv",
    "realized_span", "clip_fraction", "clip_count", "min_gap", "raw_mean_gap",
    "raw_cv", "family", "target_mean_gap", "requested_cv", "floor_applied",
    "floor_value", "expected_mean_gap_after_floor", "expected_clip_fraction",
    "source_valid", "source_error",
)
NUMERICAL_FIELDS = (
    "hessian_condition_number", "hessian_relative_difference",
    "hessian_whitened_difference", "standardized_score",
)
INTERVAL_STATUSES = {
    "valid", "invalid_point_fit", "invalid_coordinates", "nonfinite_derivatives",
    "unstable_hessian", "nonpositive_hessian", "invalid_hessian",
    "ill_conditioned_hessian", "nonstationary_point_fit", "invalid_covariance",
    "invalid_interval", "numerical_failure",
}


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _json(text, label):
    def invalid(value):
        raise ValueError(f"nonstandard nonfinite JSON value {value} in {label}")
    return json.loads(text, parse_constant=invalid)


def _read_json(path):
    return _json(Path(path).read_text(), str(path))


def _finite_number(value, label, *, missing=False):
    if value is None and missing:
        return np.nan
    _require(isinstance(value, (int, float, np.integer, np.floating))
             and not isinstance(value, (bool, np.bool_)), f"nonnumeric {label}")
    value = float(value)
    _require(np.isfinite(value), f"nonfinite {label}")
    return value


def _integer(value, label, *, minimum=0):
    _require(isinstance(value, (int, np.integer)) and not isinstance(value, (bool, np.bool_))
             and value >= minimum, f"invalid integer {label}")
    return int(value)


def _close(actual, expected, label, *, rtol=ARITHMETIC_RTOL, atol=ARITHMETIC_ATOL):
    a, b = np.asarray(actual, dtype=float), np.asarray(expected, dtype=float)
    _require(a.shape == b.shape and np.allclose(a, b, rtol=rtol, atol=atol, equal_nan=True),
             f"arithmetic mismatch: {label}")


def _scalar(actual, expected, label):
    """Compare a parsed CSV scalar with an independently derived value."""
    if isinstance(expected, str):
        _require(actual == expected or (expected == "" and pd.isna(actual)),
                 f"text mismatch: {label}")
    elif isinstance(expected, (bool, np.bool_)):
        _require(actual in (bool(expected), str(bool(expected))), f"Boolean mismatch: {label}")
    elif expected is None or pd.isna(expected):
        _require(actual is None or actual == "" or pd.isna(actual), f"missing-value mismatch: {label}")
    elif isinstance(expected, (int, np.integer)):
        _require(float(actual) == int(expected), f"count mismatch: {label}")
    else:
        _close(float(actual), float(expected), label)


def _mean_se(values):
    values = np.asarray(values, dtype=float)
    if len(values) < 2:
        return np.nan
    centered = values - values.mean()
    return float(np.sqrt(np.dot(centered, centered) / (len(values) * (len(values) - 1))))


def _sd(values):
    return _mean_se(values) * np.sqrt(len(values)) if len(values) > 1 else np.nan


def _wilson(successes, total):
    if not total:
        return np.nan, np.nan
    p = successes / total
    # Algebraic form obtained by multiplying the conventional center/width
    # formula by N. This is deliberately independent of the reporter helper.
    center = (successes + Z * Z / 2) / (total + Z * Z)
    half = Z * np.sqrt(total * p * (1-p) + Z * Z / 4) / (total + Z * Z)
    return max(0., min(p, center-half)), min(1., max(p, center+half))


def _classification(truth, coverage_lo, coverage_hi, bias_lo, bias_hi):
    if truth == 0:
        return "not_applicable_zero_truth"
    bias_known = np.isfinite(bias_lo) and np.isfinite(bias_hi)
    if (np.isfinite(coverage_hi) and coverage_hi < .80) or (
            bias_known and (bias_lo > .10 or bias_hi < -.10)):
        return "supported_breach"
    if (np.isfinite(coverage_lo) and coverage_lo >= .80 and bias_known
            and bias_lo >= -.10 and bias_hi <= .10):
        return "within_criteria"
    return "unresolved"


def _matrix(row, kind):
    matrix = np.empty((3, 3))
    for i, p in enumerate(PARAMETERS):
        for j in range(i, 3):
            key = f"cov_{p}_{PARAMETERS[j]}" if kind == "cov" else f"{kind}_{i}_{j}"
            matrix[i, j] = matrix[j, i] = _finite_number(row[key], key, missing=True)
    return matrix


def _source(row, cell, design):
    for key in SOURCE_FIELDS:
        _require(key in row, f"missing source diagnostic: {key}")
    _require(type(row["source_valid"]) is bool and row["source_valid"], "unavailable source requires implementation review")
    _require(row["source_error"] == "" and row["pipeline_error"] == "", "pipeline exception requires implementation review")
    _require(_integer(row["source_seed"], "source seed") == design.seed, "source seed mismatch")
    words = [int(cell.identity[i:i+8], 16) for i in range(0, 64, 8)]
    prefix = [4, *words, row["rep"]]
    for field, expected in (("source_spawn_key", prefix), ("gap_spawn_key", [*prefix, 0]),
                            ("path_spawn_key", [*prefix, 1])):
        value = _json(row[field], field)
        _require(isinstance(value, list) and all(type(x) is int for x in value)
                 and value == expected, f"registered source addressing mismatch: {field}")
    for field in ("times_sha256", "observations_sha256"):
        _require(isinstance(row[field], str) and re.fullmatch(r"[0-9a-f]{64}", row[field]), f"invalid source array hash: {field}")
    k, nominal = cell.n-1, cell.q/cell.theta
    floored = cell.cv >= cell.floor_cv_threshold
    _require(type(row["floor_applied"]) is bool and row["floor_applied"] == floored, "floor rule mismatch")
    _require(row["family"] == ("equidistant" if cell.cv == 0 else "exponential" if cell.cv == 1 else "gamma"), "spacing family mismatch")
    _close(row["target_mean_gap"], nominal, "nominal mean gap")
    _close(row["requested_cv"], cell.cv, "nominal CV")
    floor = cell.floor_fraction * nominal if floored else 0.
    _close(row["floor_value"], floor, "literal floor")
    count = _integer(row["clip_count"], "clip count")
    _require(count <= k, "clip count exceeds transitions")
    _close(row["clip_fraction"], count/k, "clip fraction")
    for field in ("raw_mean_gap", "realized_mean_gap", "realized_span", "min_gap"):
        _require(_finite_number(row[field], field) > 0, f"nonpositive {field}")
    for field in ("raw_cv", "realized_cv", "expected_clip_fraction"):
        _require(_finite_number(row[field], field) >= 0, f"negative {field}")
    _close(row["realized_span"], k*row["realized_mean_gap"], "span and realized mean")
    _require(row["min_gap"] <= row["realized_mean_gap"] * (1+1e-12), "minimum gap exceeds mean")
    # With K positive summands, sequential floating-point accumulation has a
    # relative error bounded by gamma_K = K*eps/(1-K*eps), conservatively using
    # machine epsilon rather than half-epsilon unit roundoff. Both recorded
    # means also use scaling, division, summation and a final multiplication.
    # Eight gamma_(K+8) terms safely cover these operations and their difference;
    # a few ulps cover the final scalar operations. A fixed number of span ulps
    # alone would NOT bound accumulation for K up to 3999. This is an artifact
    # arithmetic check, not a change to any registered inference gate.
    eps = np.finfo(float).eps
    operation_error = (k+8)*eps
    _require(operation_error < 1, "path length exceeds the accumulation-error bound domain")
    gamma = operation_error/(1-operation_error)
    span_scale = max(row["realized_span"], k*nominal)
    _finite_number(span_scale, "span arithmetic scale")
    rounding = 8*gamma*span_scale + 8*np.spacing(span_scale)
    delta_span = k*(row["realized_mean_gap"] - row["raw_mean_gap"])
    _require(-rounding <= delta_span <= count*floor+rounding, "gap changes violate literal flooring without renormalization")
    if floored:
        shape, scale = cell.cv**-2, nominal*cell.cv**2
        clip = float(gammainc(shape, floor/scale))
        expected_mean = floor*clip + nominal*float(gammaincc(shape+1, floor/scale))
        # An individual gap is the subtraction of consecutive accumulated
        # timestamps. Its local error is a few ulps at the timestamp scale;
        # the much larger whole-sum bound above must not loosen this check.
        gap_rounding = 4*np.spacing(row["realized_span"])
        _require(row["min_gap"] >= floor-gap_rounding, "realized gap violates floor beyond timestamp roundoff")
    else:
        clip, expected_mean = 0., nominal
        _require(count == 0, "clipping recorded where floor is inactive")
    _close(row["expected_clip_fraction"], clip, "population clip probability")
    _close(row["expected_mean_gap_after_floor"], expected_mean, "population floored mean")
    if cell.cv == 0:
        for field in ("raw_mean_gap", "realized_mean_gap", "min_gap"):
            _close(row[field], nominal, f"regular {field}")
        _require(row["raw_cv"] == 0 and row["realized_cv"] == 0, "regular CV diagnostic mismatch")


def _inference(row):
    """Check all retained parameters and every available inference identity."""
    for flag in ("point_valid", "interval_valid", "optimizer_success"):
        _require(type(row[flag]) is bool, f"non-Boolean {flag}")
    pv, iv = row["point_valid"], row["interval_valid"]
    _require(pv == (row["point_status"] == "valid_interior"), "point validity/status mismatch")
    _require(iv == (row["interval_status"] == "valid"), "interval validity/status mismatch")
    _require(row["interval_status"] in INTERVAL_STATUSES, "unknown registered interval status")
    _require(not iv or pv, "valid interval without valid point fit")
    _require(not pv or row["optimizer_success"], "valid point without optimizer success")
    _require(isinstance(row["point_reason"], str) and row["point_reason"], "missing point reason")
    _require(isinstance(row["interval_reason"], str) and row["interval_reason"], "missing interval reason")
    _require(row["covariance_type"] == "inverse_observed_hessian_model_based", "covariance interpretation mismatch")
    point = _json(row["point_diagnostics_json"], "point diagnostics")
    diag = _json(row["interval_diagnostics_json"], "interval diagnostics")
    _require(isinstance(point, dict) and isinstance(diag, dict), "diagnostic objects required")
    estimates = np.array([_finite_number(row[p], p, missing=not pv) for p in PARAMETERS])
    if pv:
        _require(estimates[0] > 0 and estimates[2] > 0, "valid point has nonpositive theta/sigma")
        _finite_number(row["nll"], "valid point NLL")
    for field, expected in (("estimator", row["estimator"]), ("fit_status", row["point_status"]),
                            ("optimizer_success", row["optimizer_success"]), ("valid_for_point_summary", pv)):
        _require(point.get(field) == expected, f"point diagnostic mismatch: {field}")
    for p in PARAMETERS:
        _scalar(point[p], row[p], f"point diagnostics {p}")
    _scalar(point["nll"], row["nll"], "point diagnostic NLL")
    for key, expected in (("confidence_level", .95), ("interval_scale", "natural_symmetric"),
                          ("conditions_on_initial_value", True), ("covariance_type", "inverse_observed_hessian"),
                          ("max_condition_number", 1e10), ("max_standardized_score", .01)):
        _require(diag.get(key) == expected, f"inference procedure mismatch: {key}")
    covariance = _matrix(row, "cov")
    working_covariance = _matrix(row, "working_cov")
    hessian = _matrix(row, "working_hessian")
    gradient = np.array([_finite_number(row[f"gradient_{i}"], f"gradient {i}", missing=True) for i in range(3)])
    lower = np.array([_finite_number(row[f"{p}_lower"], f"{p} lower", missing=True) for p in PARAMETERS])
    upper = np.array([_finite_number(row[f"{p}_upper"], f"{p} upper", missing=True) for p in PARAMETERS])
    se = np.array([_finite_number(row[f"{p}_se"], f"{p} SE", missing=True) for p in PARAMETERS])
    if not iv:
        _require(np.isnan(np.concatenate((covariance.ravel(), working_covariance.ravel(), lower, upper, se))).all(),
                 "unavailable interval retains published covariance/endpoints")
    if not pv:
        _require(row["interval_status"] == "invalid_point_fit", "invalid point received an inappropriate interval outcome")
        _require(np.isnan(hessian).all() and np.isnan(gradient).all(), "invalid point has a stale Hessian trace")
        return
    if row["interval_status"] == "invalid_coordinates":
        _require(not iv, "invalid coordinates marked valid")
        return
    offset = _finite_number(diag["state_offset"], "observed offset")
    scale = _finite_number(diag["state_scale"], "observed state scale")
    time_scale = _finite_number(diag["time_scale"], "observed time scale")
    _require(scale > 0 and time_scale > 0, "nonpositive numerical scale")
    _close(time_scale, row["realized_mean_gap"], "observed mean gap normalization")
    for key, value in (("state_center", offset), ("state_scale", scale), ("time_scale", time_scale)):
        _close(point[key], value, f"point/interval coordinate agreement {key}")
    eta = np.asarray(diag["working_parameters"], dtype=float)
    expected_eta = np.array([np.log(estimates[0])+np.log(time_scale),
                             float((np.longdouble(estimates[1])-offset)/scale),
                             np.log(estimates[2])+.5*np.log(time_scale)-np.log(scale)])
    _close(eta, expected_eta, "normalized working parameters")
    _close(point["nll_state_jacobian"], (row["n"]-1)*np.log(scale), "conditional state density Jacobian")
    _close(point["normalized_nll"]+point["nll_state_jacobian"], row["nll"], "normalized conditional NLL")
    hd = diag["hessian_diagnostics"]
    for key, expected in (("relative_step", .001), ("max_relative_difference", .005),
                          ("max_whitened_difference", .01), ("selected_step_fraction", .5)):
        _require(hd.get(key) == expected, f"Hessian procedure mismatch: {key}")
    _close(hd["steps"], .001*np.maximum(1., np.abs(eta)), "registered finite-difference steps")
    for field, expected in (("hessian_relative_difference", diag.get("hessian_relative_difference")),
                            ("hessian_whitened_difference", hd.get("whitened_relative_difference")),
                            ("hessian_condition_number", diag.get("hessian_condition_number")),
                            ("standardized_score", diag.get("standardized_score"))):
        _scalar(row[field], expected, f"nested/raw {field}")
    hessian_status = diag["hessian_status"]
    for field in NUMERICAL_FIELDS:
        if row[field] is not None:
            _require(_finite_number(row[field], field) >= 0, f"negative numerical diagnostic: {field}")
    if hessian_status != "stable":
        _require(not iv and row["interval_status"] == hessian_status,
                 "derivative rejection and interval status disagree")
        if hessian_status == "unstable_hessian":
            rel, white = row["hessian_relative_difference"], row["hessian_whitened_difference"]
            failed = rel is None or rel > .005 or (
                "whitened_relative_difference" in hd and (white is None or white > .01))
            _require(failed, "unstable Hessian without a recorded failed stability gate")
        elif hessian_status == "nonpositive_hessian":
            _require(np.isfinite(hessian).all(), "nonpositive Hessian rejection lacks a finite matrix")
            positive = True
            try:
                np.linalg.cholesky(hessian)
            except np.linalg.LinAlgError:
                positive = False
            _require(not positive, "nonpositive Hessian rejection has a positive-definite matrix")
        elif hessian_status == "nonfinite_derivatives":
            _require(not (np.isfinite(hessian).all() and np.isfinite(gradient).all()),
                     "nonfinite-derivative rejection has finite saved derivatives")
        else:
            _require(hessian_status == "invalid_coordinates", "unknown derivative rejection")
        return
    if hessian_status == "stable":
        for field, maximum in (("hessian_relative_difference", .005), ("hessian_whitened_difference", .01)):
            value = _finite_number(row[field], field)
            _require(0 <= value <= maximum, f"registered numerical gate violated: {field}")
        _require(np.isfinite(hessian).all() and np.isfinite(gradient).all(), "accepted derivative provider has nonfinite derivatives")
    if not np.isfinite(hessian).all() or not np.isfinite(gradient).all():
        _require(not iv, "nonfinite accepted derivatives")
        return
    eigenvalues = np.linalg.eigvalsh(hessian)
    if "hessian_eigenvalues" in diag:
        _close(diag["hessian_eigenvalues"], eigenvalues, "Hessian eigenvalues")
    if iv:
        _require(eigenvalues[0] > 0, "valid interval Hessian is not positive definite")
    if eigenvalues[0] <= 0:
        _require(row["interval_status"] == "nonpositive_hessian", "nonpositive curvature has an inconsistent rejection status")
        return
    condition = float(eigenvalues[-1]/eigenvalues[0])
    if row["hessian_condition_number"] is not None:
        _close(row["hessian_condition_number"], condition, "Hessian condition number")
    if iv:
        _require(np.isfinite(condition) and condition <= 1e10, "registered condition-number gate violated")
    elif condition > 1e10:
        _require(row["interval_status"] == "ill_conditioned_hessian", "failed conditioning gate has inconsistent rejection status")
        return
    inverse = np.linalg.inv(hessian)
    score = float(np.sqrt(max(0., gradient @ inverse @ gradient)))
    if row["standardized_score"] is not None:
        _close(row["standardized_score"], score, "standardized score")
    if iv:
        _require(score <= .01, "registered standardized-score gate violated")
        _close(working_covariance, inverse, "full inverse-Hessian covariance", rtol=2e-9,
               atol=2e-12*max(1., np.linalg.norm(inverse, ord="fro")))
        # The middle Jacobian entry is the OBSERVED state scale, not one and
        # not the simulation's sigma or stationary standard deviation.
        jacobian = np.diag([estimates[0], scale, estimates[2]])
        natural = jacobian @ inverse @ jacobian.T
        _close(covariance, natural, "full natural covariance transformation", rtol=2e-9,
               atol=2e-12*max(1., np.linalg.norm(natural, ord="fro")))
        _require(np.isfinite(covariance).all() and (np.diag(covariance) > 0).all(), "invalid accepted natural covariance")
        expected_se = np.sqrt(np.diag(covariance))
        _close(se, expected_se, "reported standard errors")
        _close(diag["standard_errors"], expected_se, "diagnostic standard errors")
        _close(diag["critical_value"], Z, "95 percent critical value")
        _close(lower, estimates-Z*expected_se, "unclipped natural lower endpoints")
        _close(upper, estimates+Z*expected_se, "natural upper endpoints")
    elif row["interval_status"] == "ill_conditioned_hessian":
        _require(condition > 1e10, "invalid condition-number rejection")
    elif row["interval_status"] == "nonstationary_point_fit":
        _require(score > .01, "invalid score rejection")
    else:
        # Remaining registered numerical rejections occur after accepted
        # derivatives. They must not be used to conceal an ordinary available
        # interval. A genuine nonfinite covariance/endpoint is retained.
        _require(score <= .01, "failed score gate has inconsistent rejection status")
        with np.errstate(over="ignore", invalid="ignore"):
            jacobian = np.diag([estimates[0], scale, estimates[2]])
            natural = jacobian @ inverse @ jacobian.T
            good_covariance = np.isfinite(natural).all() and (np.diag(natural) > 0).all()
            if good_covariance:
                errors = Z*np.sqrt(np.diag(natural))
                good_endpoints = np.isfinite(estimates-errors).all() and np.isfinite(estimates+errors).all()
            else:
                good_endpoints = False
        _require((not good_covariance and row["interval_status"] == "invalid_covariance")
                 or (good_covariance and not good_endpoints and row["interval_status"] == "invalid_interval"),
                 "interval rejected although all recorded numerical gates passed; implementation review required")


def _load_cell(output, cell, design, identity):
    directory = output/"cells"/cell.identity
    manifest = _read_json(directory/"cell_manifest.json")
    _require(manifest["schema_version"] == 1 and manifest["run_identity"] == identity, "cell provenance identity mismatch")
    _require(canonical(manifest["cell"]) == canonical(cell.to_dict()), "cell scientific identity mismatch")
    _require(_integer(manifest["replications"], "cell replications", minimum=1) == design.replications,
             "cell replication plan mismatch")
    _require(manifest["status"] == "complete"
             and _integer(manifest["completed_replications"], "completed replications") == design.replications,
             "incomplete cell cannot pass artifact audit")
    rows, names, start = [], set(), 0
    for part in manifest["partitions"]:
        begin, stop = _integer(part["start"], "partition start"), _integer(part["stop"], "partition stop")
        _require(begin == start and begin < stop <= design.replications, "incomplete or overlapping partition ranges")
        name = f"part-{begin:06d}-{stop:06d}.jsonl.gz"
        _require(part["filename"] == name and name not in names, "duplicate or unsafe partition filename")
        names.add(name)
        _require(file_digest(directory/name) == part["sha256"], "partition digest mismatch")
        with gzip.open(directory/name, "rt", encoding="utf-8") as stream:
            current = [_json(line, name) for line in stream if line.strip()]
        _require(len(current) == _integer(part["rows"], "partition row count") == 3*(stop-begin), "partition row count mismatch")
        keys = [(r.get("rep"), r.get("estimator")) for r in current]
        expected = {(r, method) for r in range(begin, stop) for method in METHODS}
        _require(len(set(keys)) == len(keys) and set(keys) == expected, "partition is missing complete planned estimator keys")
        for row in current:
            _integer(row["rep"], "replication index")
            for key, value in cell.to_dict().items():
                _require(key in row and type(row[key]) is type(value) and row[key] == value,
                         f"raw scientific metadata mismatch: {key}")
        rows.extend(current)
        start = stop
    _require(start == design.replications, "missing final planned partition")
    actual = {p.name for p in directory.iterdir()}
    _require(actual == names | {"cell_manifest.json"}, "unexpected unreferenced cell evidence")
    rows.sort(key=lambda r: (r["rep"], METHODS.index(r["estimator"])))
    for rep in range(design.replications):
        triplet = rows[3*rep:3*rep+3]
        _source(triplet[0], cell, design)
        for row in triplet:
            for field in SOURCE_FIELDS:
                _require(field in row and type(row[field]) is type(triplet[0][field])
                         and row[field] == triplet[0][field], f"shared source diagnostic mismatch: {field}")
            _require(row["pipeline_error"] == "", "pipeline exception requires implementation review")
            _inference(row)
    return rows, len(names)


def _hits(rows, parameter):
    return np.array([bool(row["interval_valid"] and row[f"{parameter}_lower"] <= row[f"true_{parameter}"]
                          <= row[f"{parameter}_upper"]) for row in rows])


def _parameter_summaries(rows, cell, design):
    result = []
    for method in METHODS:
        group = [r for r in rows if r["estimator"] == method]
        point = [r for r in group if r["point_valid"]]
        interval = [r for r in group if r["interval_valid"]]
        n, s, v = design.replications, len(point), len(interval)
        for p in PARAMETERS:
            truth = cell.to_dict()[f"true_{p}"]
            errors = np.array([r[p]-truth for r in point])
            with np.errstate(over="ignore", invalid="ignore"):
                _require(np.isfinite(errors**2).all(), "point-error moments overflow")
                bias = float(errors.mean()) if s else np.nan
                rmse = float(np.sqrt(np.dot(errors, errors)/s)) if s else np.nan
                se_bias = _mean_se(errors)
                se_rmse = _mean_se(errors**2)/(2*rmse) if rmse > 0 else np.nan
            if s:
                _require(np.isfinite(bias) and np.isfinite(rmse), "valid point metrics overflow")
            if s > 1:
                _require(np.isfinite(se_bias), "bias Monte Carlo uncertainty overflow")
                if rmse > 0:
                    _require(np.isfinite(se_rmse), "RMSE Monte Carlo uncertainty overflow")
            relative = bias/abs(truth) if truth else np.nan
            relative_se = se_bias/abs(truth) if truth else np.nan
            blo, bhi = relative-Z*relative_se, relative+Z*relative_se
            hits = _hits(group, p)
            count = int(hits.sum())
            conditional, operational = count/v if v else np.nan, count/n
            clo, chi = _wilson(count, v)
            olo, ohi = _wilson(count, n)
            estimates = np.array([r[p] for r in interval])
            ses = np.array([r[f"{p}_se"] for r in interval])
            row = {"cell_id": cell.identity, "estimator": method, "parameter": p, "truth": truth,
                   "n_planned": n, "n_point_valid": s, "n_interval_valid": v, "n_covered": count,
                   "point_availability": s/n, "interval_availability": v/n,
                   "bias": bias, "rmse": rmse, "mcse_bias": se_bias,
                   "relative_bias": relative, "relative_mcse_bias": relative_se,
                   "relative_bias_mc_lower": blo, "relative_bias_mc_upper": bhi,
                   "coverage_conditional": conditional, "coverage_operational": operational,
                   "conditional_mcse": np.sqrt(conditional*(1-conditional)/v) if v else np.nan,
                   "operational_mcse": np.sqrt(operational*(1-operational)/n),
                   "conditional_wilson_lower": clo, "conditional_wilson_upper": chi,
                   "operational_wilson_lower": olo, "operational_wilson_upper": ohi,
                   "negative_lower_count": sum(r[f"{p}_lower"] < 0 for r in interval),
                   "practical_classification": _classification(truth, olo, ohi, blo, bhi),
                   "point_metric_scope": "conditional_on_valid_point_fit",
                   "primary_coverage_scope": "all_planned_paths_unavailable_counts_as_not_covering",
                   "mcse_rmse": se_rmse, "relative_rmse": rmse/abs(truth) if truth else np.nan,
                   "relative_mcse_rmse": se_rmse/abs(truth) if truth else np.nan,
                   "empirical_sd_point_valid": _sd(errors),
                   "bias_mc_lower": bias-Z*se_bias, "bias_mc_upper": bias+Z*se_bias,
                   "empirical_sd_interval_valid": _sd(estimates),
                   "mean_reported_se": float(ses.mean()) if v else np.nan,
                   "rms_reported_se": float(np.sqrt(np.dot(ses, ses)/v)) if v else np.nan,
                   "n": cell.n, "theta_mean_gap": cell.q, "cv": cell.cv,
                   "mean_gap": cell.q/cell.theta, "floor_fraction": cell.floor_fraction,
                   "nominal_dimensionless_span": (cell.n-1)*cell.q}
            for name, number in (("point", s), ("interval", v)):
                lo, hi = _wilson(number, n)
                row.update({f"{name}_availability_mcse": np.sqrt(number/n*(1-number/n)/n),
                            f"{name}_availability_wilson_lower": lo,
                            f"{name}_availability_wilson_upper": hi})
            result.append(row)
    return result


def _difference(row, key, value, standard_error):
    row.update({key: value, f"{key}_mcse": standard_error,
                f"{key}_mc_lower": value-Z*standard_error,
                f"{key}_mc_upper": value+Z*standard_error})


def _relative_difference(row, name, truth):
    for suffix in ("", "_mcse", "_mc_lower", "_mc_upper"):
        row[f"relative_{name}{suffix}"] = row[f"{name}{suffix}"]/abs(truth) if truth else np.nan


def _ratio(row, numerator, denominator, log_se):
    ratio = numerator/denominator if np.isfinite(denominator) and denominator > 0 else np.nan
    row.update(rmse_ratio=ratio, rmse_ratio_mcse=ratio*log_se,
               rmse_ratio_log_mcse=log_se,
               rmse_ratio_mc_lower=ratio*np.exp(-Z*log_se),
               rmse_ratio_mc_upper=ratio*np.exp(Z*log_se),
               rmse_ratio_uncertainty="delta_method_log_ratio_95pct")


def _regular_comparisons(summaries, design):
    index = {(r["cell_id"], r["estimator"], r["parameter"]): r for r in summaries}
    cells = {c.identity: c for c in design.cells}
    result = []
    for a in summaries:
        cell = cells[a["cell_id"]]
        reference = replace(cell, cv=0).identity
        key = (reference, a["estimator"], a["parameter"])
        _require(key in index, "missing complete regular reference")
        b = index[key]
        same = reference == cell.identity
        row = {"cell_id": cell.identity, "reference_cell_id": reference,
               "estimator": a["estimator"], "parameter": a["parameter"], "truth": a["truth"],
               "n": cell.n, "theta_mean_gap": cell.q, "cv": cell.cv,
               "comparison_type": "self_reference" if same else "independent_cells",
               "point_subset": "each_cell_own_valid_point_fits",
               "matching": "same_n_truth_nominal_mean_gap_and_floor_settings_not_realized_span"}
        for name in ("n_planned", "n_point_valid", "n_interval_valid", "n_covered", "bias", "rmse",
                     "mcse_bias", "mcse_rmse", "coverage_operational", "coverage_conditional"):
            row[f"cell_{name}"], row[f"reference_{name}"] = a[name], b[name]
        for name, metric, uncertainty in (("excess_bias", "bias", "mcse_bias"),
                ("operational_coverage_difference", "coverage_operational", "operational_mcse"),
                ("conditional_coverage_difference", "coverage_conditional", "conditional_mcse")):
            value = a[metric]-b[metric]
            error = (0. if np.isfinite(value) else np.nan) if same else np.sqrt(a[uncertainty]**2+b[uncertainty]**2)
            _difference(row, name, value, error)
        _relative_difference(row, "excess_bias", a["truth"])
        log_se = np.nan
        if a["rmse"] > 0 and b["rmse"] > 0:
            log_se = 0. if same else np.sqrt((a["mcse_rmse"]/a["rmse"])**2+(b["mcse_rmse"]/b["rmse"])**2)
        _ratio(row, a["rmse"], b["rmse"], log_se)
        result.append(row)
    return result


def _paired_comparisons(rows, cell, design):
    groups = {name: sorted((r for r in rows if r["estimator"] == name), key=lambda r: r["rep"])
              for name in METHODS}
    result = []
    for name_a, name_b in combinations(METHODS, 2):
        a, b = groups[name_a], groups[name_b]
        common = [(x, y) for x, y in zip(a, b) if x["point_valid"] and y["point_valid"]]
        s = len(common)
        for p in PARAMETERS:
            truth = cell.to_dict()[f"true_{p}"]
            ea = np.array([x[p]-truth for x, _ in common])
            eb = np.array([y[p]-truth for _, y in common])
            u, v = ea**2, eb**2
            ma, mb = (float(u.mean()), float(v.mean())) if s else (np.nan, np.nan)
            ra, rb = np.sqrt(ma), np.sqrt(mb)
            covariance = float(np.dot(u-ma, v-mb)/(s-1)) if s > 1 else np.nan
            row = {"cell_id": cell.identity, "n": cell.n, "theta_mean_gap": cell.q, "cv": cell.cv,
                   "estimator_a": name_a, "estimator_b": name_b, "parameter": p, "truth": truth,
                   "comparison_type": "paired_shared_source_paths", "direction": "a_minus_b_and_a_over_b",
                   "n_planned": design.replications, "n_point_valid_a": sum(x["point_valid"] for x in a),
                   "n_point_valid_b": sum(x["point_valid"] for x in b), "n_common_point_valid": s,
                   "n_common_interval_valid": sum(x["interval_valid"] and y["interval_valid"] for x, y in zip(a, b)),
                   "point_subset": "intersection_of_valid_point_fits",
                   "coverage_subset": "all_planned_paths_unavailable_counts_as_not_covering",
                   "bias_a_common": float(ea.mean()) if s else np.nan,
                   "bias_b_common": float(eb.mean()) if s else np.nan,
                   "rmse_a_common": ra, "rmse_b_common": rb, "squared_error_covariance_common": covariance}
            for name, differences in (("bias_difference", ea-eb), ("mse_difference", u-v)):
                _difference(row, name, float(differences.mean()) if s else np.nan, _mean_se(differences))
            _relative_difference(row, "bias_difference", truth)
            # Centered influences independently reproduce the reporting
            # module's uncentered-influence sample-variance calculation.
            rmse_se = _mean_se((u-ma)/(2*ra)-(v-mb)/(2*rb)) if ra > 0 and rb > 0 else np.nan
            _difference(row, "rmse_difference", ra-rb, rmse_se)
            log_se = _mean_se(.5*((u-ma)/ma-(v-mb)/mb)) if ma > 0 and mb > 0 else np.nan
            _ratio(row, ra, rb, log_se)
            ha, hb = _hits(a, p), _hits(b, p)
            differences = ha.astype(int)-hb.astype(int)
            row.update(n_covered_a=int(ha.sum()), n_covered_b=int(hb.sum()),
                       n_both_cover=int((ha & hb).sum()), n_a_only_cover=int((ha & ~hb).sum()),
                       n_b_only_cover=int((~ha & hb).sum()), n_neither_cover=int((~ha & ~hb).sum()),
                       operational_coverage_a=float(ha.mean()), operational_coverage_b=float(hb.mean()))
            _difference(row, "operational_coverage_difference", float(differences.mean()), _mean_se(differences))
            result.append(row)
    return result


def _other_summaries(rows, cell, design):
    sources = [r for r in rows if r["estimator"] == "exact"]
    spacing = {"cell_id": cell.identity, "n_planned": design.replications,
               "n": cell.n, "theta_mean_gap": cell.q, "cv": cell.cv}
    for field in ("realized_mean_gap", "realized_cv", "realized_span", "clip_fraction", "min_gap"):
        values = np.array([r[field] for r in sources])
        spacing.update({f"{field}_n_available": len(values), f"{field}_mean": float(values.mean()),
                        f"{field}_mcse": _mean_se(values), f"{field}_min": float(values.min()),
                        f"{field}_max": float(values.max())})
    counts = Counter((r["estimator"], r["point_status"], r["interval_status"]) for r in rows)
    outcomes = []
    for (method, point_status, interval_status), number in counts.items():
        lo, hi = _wilson(number, design.replications)
        outcomes.append({"cell_id": cell.identity, "n": cell.n, "theta_mean_gap": cell.q, "cv": cell.cv,
                         "estimator": method, "point_status": point_status, "interval_status": interval_status,
                         "count": number, "n_planned": design.replications, "fraction": number/design.replications,
                         "wilson_lower": lo, "wilson_upper": hi})
    diagnostics = []
    for method in METHODS:
        group = [r for r in rows if r["estimator"] == method]
        diagnostic = {"cell_id": cell.identity, "estimator": method, "n_planned": design.replications,
                      "n_interval_valid": sum(r["interval_valid"] for r in group),
                      "n": cell.n, "theta_mean_gap": cell.q, "cv": cell.cv}
        for field in NUMERICAL_FIELDS:
            finite = np.array([r[field] for r in group if r[field] is not None])
            _require(np.isfinite(finite).all(), "nonfinite JSON numerical diagnostics")
            diagnostic.update({f"{field}_n_finite": len(finite), f"{field}_n_nonfinite": len(group)-len(finite),
                               f"{field}_max": float(finite.max()) if len(finite) else np.nan,
                               f"{field}_median": float(np.percentile(finite, 50)) if len(finite) else np.nan,
                               f"{field}_p95": float(np.percentile(finite, 95)) if len(finite) else np.nan})
        diagnostics.append(diagnostic)
    return [spacing], outcomes, diagnostics


def _check_table(table, expected, name):
    keys = TABLE_KEYS[name]
    _require(expected and set(table.columns) == set(expected[0]), f"table schema mismatch: {name}")
    _require(len(table) == len(expected) and not table.duplicated(list(keys)).any(), f"table row/key count mismatch: {name}")
    actual = {tuple(row[k] for k in keys): row for row in table.to_dict("records")}
    desired = {tuple(row[k] for k in keys): row for row in expected}
    _require(actual.keys() == desired.keys(), f"table complete key set mismatch: {name}")
    for key, wanted in desired.items():
        for column, value in wanted.items():
            _scalar(actual[key][column], value, f"{name} {key} {column}")


class _ReplayHessian:
    def __init__(self):
        self.provider = FiniteDifferenceHessian(relative_step=.001, max_relative_difference=.005,
                                               max_whitened_difference=.01)
        self.last = None

    def evaluate(self, objective, point):
        self.last = self.provider.evaluate(objective, point)
        return self.last


def _replay(rows, cell, design):
    source = OuSampleSource(design.seed)
    records = {(r["rep"], r["estimator"]): r for r in rows}
    report = []
    for rep in sorted({0, design.replications-1}):
        sample, source_diagnostics = source.sample(cell.to_dict(), rep)
        hashes = {"times_sha256": hashlib.sha256(np.asarray(sample.times, dtype="<f8").tobytes()).hexdigest(),
                  "observations_sha256": hashlib.sha256(np.asarray(sample.values, dtype="<f8").tobytes()).hexdigest()}
        reference = records[(rep, "exact")]
        for key, value in source_diagnostics.items():
            _scalar(reference[key], value, f"source replay {cell.identity} rep={rep} {key}")
        for key, value in hashes.items():
            _require(reference[key] == value, f"actual source array hash replay mismatch: {key}")
        for method in METHODS:
            row = records[(rep, method)]
            trace = _ReplayHessian()
            service = WaldInference(Week3Estimator(method), OuLikelihood(method), hessian=trace,
                                    level=.95, max_condition_number=1e10, max_standardized_score=.01)
            fitted = service.infer(sample)
            point = fitted.point_fit
            _require((point.valid, point.optimizer_success, point.status, fitted.valid, fitted.status)
                     == (row["point_valid"], row["optimizer_success"], row["point_status"], row["interval_valid"], row["interval_status"]),
                     "replay point/interval classification mismatch")
            actual = np.array([row[p] if row[p] is not None else np.nan for p in PARAMETERS])
            _close(actual, point.parameters, "replayed point parameters", rtol=REPLAY_RTOL, atol=REPLAY_ATOL)
            _close(np.nan if row["nll"] is None else row["nll"], point.nll, "replayed NLL", rtol=REPLAY_RTOL, atol=REPLAY_ATOL)
            for suffix, expected in (("lower", fitted.lower), ("upper", fitted.upper)):
                saved = [row[f"{p}_{suffix}"] if row[f"{p}_{suffix}"] is not None else np.nan for p in PARAMETERS]
                _close(saved, expected, f"replayed {suffix}", rtol=REPLAY_RTOL, atol=REPLAY_ATOL)
            _close(_matrix(row, "cov"), fitted.covariance, "replayed covariance", rtol=REPLAY_RTOL, atol=REPLAY_ATOL)
            _close(_matrix(row, "working_cov"), fitted.working_covariance, "replayed working covariance", rtol=REPLAY_RTOL, atol=REPLAY_ATOL)
            if trace.last is not None:
                _close(_matrix(row, "working_hessian"), trace.last.hessian, "replayed Hessian", rtol=REPLAY_RTOL, atol=REPLAY_ATOL)
                _close([np.nan if row[f"gradient_{i}"] is None else row[f"gradient_{i}"] for i in range(3)],
                       trace.last.gradient, "replayed gradient", rtol=REPLAY_RTOL, atol=REPLAY_ATOL)
                for field, value in (("hessian_relative_difference", trace.last.relative_difference),
                                     ("hessian_whitened_difference", trace.last.diagnostics.get("whitened_relative_difference", np.nan))):
                    _close(np.nan if row[field] is None else row[field], value, f"replayed {field}", rtol=REPLAY_RTOL, atol=REPLAY_ATOL)
            diag = fitted.diagnostics
            if "state_scale" in diag:
                # Recompute the numerical constants from the observed sample;
                # no parameter truth enters the inference normalization.
                values = sample.values.astype(np.longdouble)
                center = float(values.mean())
                centered = values-np.longdouble(center)
                amplitude = np.max(np.abs(centered))
                scale = float(amplitude*np.sqrt(np.mean((centered/amplitude)**2)))
                _close(diag["state_offset"], center, "replayed observed offset")
                _close(diag["state_scale"], scale, "replayed observed state scale")
                _close(diag["time_scale"], float(np.diff(sample.times).astype(np.longdouble).mean()), "replayed observed time scale")
                saved_diag = _json(row["interval_diagnostics_json"], "saved inference diagnostics")
                for field in ("state_offset", "state_scale", "time_scale"):
                    _close(saved_diag[field], diag[field], f"replayed {field}")
        report.append({"cell_id": cell.identity, "replication": rep, **hashes, "estimators_checked": list(METHODS)})
    return report


def audit_run(output_dir, design, *, replay=True):
    """Raise ValueError on any defect; return a read-only JSON-safe audit report.

    ``replay=False`` is useful for fault-injection tests or arithmetic checks,
    but does not constitute the full scientific completion audit. The report
    explicitly records that distinction. Raw memory is bounded by one cell;
    the failure export is checked as a stream in canonical saved-row order.
    """
    try:
        return _audit_run(Path(output_dir), design, replay=replay)
    except ValueError:
        raise
    except (KeyError, TypeError, OSError, OverflowError, IndexError, np.linalg.LinAlgError) as exc:
        raise ValueError(f"artifact audit failed: {type(exc).__name__}: {exc}") from exc


def artifact_fingerprint(output):
    """Hash audited numerical evidence, excluding later review/figure products."""
    output = Path(output)
    excluded = {"run_metadata.json", ".execution.lock", "audit.json", "visual_review.json"}
    return {str(p.relative_to(output)): file_digest(p) for p in sorted(output.rglob("*"))
            if p.is_file() and p.name not in excluded and "figures" not in p.relative_to(output).parts}


def _audit_run(output, design, *, replay):
    _require(type(replay) is bool, "replay option must be Boolean")
    before = artifact_fingerprint(output)
    metadata_before = file_digest(output/"run_metadata.json")
    metadata = _read_json(output/"run_metadata.json")
    _require(metadata["schema_version"] == 1 and metadata["project"] == "ou-irregular-experiments", "unknown run schema")
    _require(metadata["status"] in {"raw_complete_pending_audit", "audited_pending_visual_review", "completed"}, "run is incomplete or requires implementation review")
    _require(metadata["scientific_overrides"] == [], "overridden scientific services cannot pass registered-procedure audit")
    _require(canonical(metadata["procedure"]) == canonical(EXPECTED_PROCEDURE), "registered procedure changed")
    _require(_integer(metadata["pipeline_errors"], "pipeline error count") == 0, "pipeline exceptions require implementation review")
    provenance = metadata["provenance"]
    _require(set(CORE_FILES) | set(AUDIT_FILES) <= set(provenance["source_files"]), "explicit dependency manifest omits an audited source file")
    verify(provenance, expected_design=design)
    _require(canonical(_read_json(output/"resolved_design.json")) == canonical(design.to_dict()), "resolved design mismatch")
    if design.scope == "headline":
        _require(provenance["uncommitted_dependency_files"] == [], "headline code was not frozen")
        for name, expected in provenance["source_files"].items():
            committed = subprocess.run(["git", "show", f"{provenance['code_commit']}:{name}"], cwd=ROOT, capture_output=True)
            _require(committed.returncode == 0 and hashlib.sha256(committed.stdout).hexdigest() == expected,
                     f"recorded headline freeze does not contain dependency: {name}")
    expected_cells = {c.identity for c in design.cells}
    _require({p.name for p in (output/"cells").iterdir()} == expected_cells, "cell directories do not equal the full planned cell set")
    _require(set(metadata["completed_cells"]) == expected_cells, "completed-cell manifest omits planned cells")
    for cid, completed in metadata["completed_cells"].items():
        _require(completed == {"cell_id": cid, "paths": design.replications, "rows": 3*design.replications, "pipeline_errors": 0},
                 "root completed-cell accounting mismatch")
    manifest = _read_json(output/"table_manifest.json")
    _require(set(manifest) == set(TABLE_KEYS) | {"failed_fits.csv.gz"}, "table manifest omits required outputs")
    tables = {}
    for name, entry in manifest.items():
        _require(file_digest(output/name) == entry["sha256"], f"table digest mismatch: {name}")
        _integer(entry["rows"], f"{name} row count")
        if name in TABLE_KEYS:
            table = pd.read_csv(output/name, float_precision="round_trip")
            _require(len(table) == entry["rows"], f"table manifest count mismatch: {name}")
            tables[name] = table
    expected_tables = {name: [] for name in TABLE_KEYS}
    replay_records, partitions, failure_count, point_count, interval_count = [], 0, 0, 0, 0
    with gzip.open(output/"failed_fits.csv.gz", "rt", newline="", encoding="utf-8") as stream:
        failures = csv.DictReader(stream)
        for cell in design.cells:
            rows, count = _load_cell(output, cell, design, provenance["identity_sha256"])
            partitions += count
            point_count += sum(r["point_valid"] for r in rows)
            interval_count += sum(r["interval_valid"] for r in rows)
            _require(set(failures.fieldnames or []) == set(rows[0]), "failure log schema mismatch")
            for row in rows:
                if not row["point_valid"] or not row["interval_valid"]:
                    exported = next(failures, None)
                    _require(exported is not None, "failure log omits a planned failure")
                    for column, value in row.items():
                        _scalar(exported[column], value, f"failure log {row['cell_id']} {row['rep']} {row['estimator']} {column}")
                    failure_count += 1
            expected_tables["parameter_summary.csv"].extend(_parameter_summaries(rows, cell, design))
            expected_tables["paired_comparisons.csv"].extend(_paired_comparisons(rows, cell, design))
            spacing, outcomes, numerical = _other_summaries(rows, cell, design)
            expected_tables["spacing_summary.csv"].extend(spacing)
            expected_tables["outcome_counts.csv"].extend(outcomes)
            expected_tables["numerical_diagnostics.csv"].extend(numerical)
            if replay:
                replay_records.extend(_replay(rows, cell, design))
        _require(next(failures, None) is None, "failure log contains unexpected or duplicated rows")
    _require(failure_count == manifest["failed_fits.csv.gz"]["rows"], "failure log row count mismatch")
    expected_tables["regular_comparisons.csv"] = _regular_comparisons(expected_tables["parameter_summary.csv"], design)
    for name, expected in expected_tables.items():
        _check_table(tables[name], expected, name)
    verify(provenance, expected_design=design)
    _require(before == artifact_fingerprint(output)
             and metadata_before == file_digest(output/"run_metadata.json"), "evidence changed during read-only audit")
    report = {
        "schema_version": 1, "status": "passed", "full_completion_audit": replay,
        "design_sha256": design.digest, "run_identity_sha256": provenance["identity_sha256"],
        "code_commit": provenance["code_commit"], "source_dependency_files_verified": len(provenance["source_files"]),
        "audited_artifact_sha256": before,
        "cells_verified": len(design.cells), "partitions_verified": partitions,
        "planned_paths": len(design.cells)*design.replications,
        "planned_fit_records": len(design.cells)*design.replications*3,
        "point_valid_records": point_count, "interval_valid_records": interval_count,
        "retained_failure_records": failure_count, "pipeline_errors": 0,
        "tables_verified": {name: len(rows) for name, rows in expected_tables.items()},
        "replay_performed": replay, "replay_selection": "replications 0 and R-1 in every cell, all estimators",
        "replayed_source_paths": len(replay_records), "replayed_fits": 3*len(replay_records),
        "replay_records": replay_records,
        "tolerances": {"arithmetic_rtol": ARITHMETIC_RTOL, "arithmetic_atol": ARITHMETIC_ATOL,
                       "replay_rtol": REPLAY_RTOL, "replay_atol": REPLAY_ATOL,
                       "covariance_rtol": 2e-9, "covariance_atol": "2e-12 * max(1, expected matrix Frobenius norm)",
                       "span_roundoff": "8*gamma_(K+8)*max(span,K*nominal_gap)+8 ulps; gamma_m=m*eps/(1-m*eps)",
                       "minimum_gap_roundoff": "4 ulps at the realized-span scale"},
        "scope": {
            "all_records": "complete keys, source address/shared diagnostics, point/interval semantics, saved fine-Hessian and covariance identities, recorded numerical gates, endpoints",
            "all_tables": "independent arithmetic, denominators, Monte Carlo uncertainty, classification and paired/independent contrasts",
            "replay_subset": "actual source array hashes, source diagnostics, point fits, full finite-difference inference and observed-data scaling",
            "limitation": "coarse Hessians and full source arrays are not stored; their full numerical reproduction is checked on the fixed replay subset",
        },
    }
    # Enforce the public contract even if a future field accidentally contains
    # a NumPy scalar or a nonfinite numerical value.
    json.dumps(report, allow_nan=False)
    return report
