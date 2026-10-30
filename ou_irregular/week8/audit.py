"""Read-only E4 audit: independent source addressing/races and table arithmetic.

No producer raw loader or reporting calculator is used for expected results.
Neutral strict-parsing/inference checks and independently implemented standard
summary formulas are composed from the reviewed core auditor. The calibrated
cell model, analytical gap-moment routines, exact OU primitive and original
endpoint fit/inference services remain shared dependencies; the report states
these limits. Source hashes are exact; diagnostic arithmetic is tolerance-bound.
"""
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import subprocess

import numpy as np

from ..experiments import audit as core
from ..experiments.execution import registered_wald, source_records
from ..experiments.provenance import CORE_FILES, ROOT, canonical, file_digest, verify
from ..inference.contracts import Sample
from ..ou import simulate_ou
from ..week7 import floor_audit as strict
from . import reporting, source
from .design import full_design, parse_cell, parse_design, source_policy, validate_design
from .execution import PROCEDURE

FILES = tuple(sorted(set(reporting.FILES) | set(strict.FILES) | {
    "ou_irregular/week8/audit.py", "docs/week8_audit.md",
}))
METHODS = ("exact", "pfml", "euler")
PARAMETERS = ("theta", "mu", "sigma")
RTOL, ATOL = 2e-10, 2e-12
AUDIT_CONTRACT = "week8_independent_e4_audit_v1"
TABLE_KEYS = {
    "parameter_summary.csv": ("cell_id", "estimator", "parameter"),
    "paired_comparisons.csv": ("cell_id", "estimator_a", "estimator_b", "parameter"),
    "regular_comparisons.csv": ("cell_id", "estimator", "parameter"),
    "matched_iid_comparisons.csv": ("cell_id", "estimator", "parameter"),
    "spacing_summary.csv": ("cell_id",), "source_summary.csv": ("cell_id",),
    "outcome_counts.csv": ("cell_id", "estimator", "point_status", "interval_status"),
    "numerical_diagnostics.csv": ("cell_id", "estimator"),
}
SOURCE_FIELDS = frozenset(source.SOURCE_FIELDS)
INFERENCE_FIELDS = frozenset(strict.INFERENCE_FIELDS)
_CELL_FIELDS = frozenset(full_design().cells[0].to_dict())
RAW_FIELDS = _CELL_FIELDS | SOURCE_FIELDS | INFERENCE_FIELDS | {
    "rep", "source_valid", "source_error", "pipeline_error"}
_TEXT_INFERENCE = strict._TEXT_INFERENCE
_FLAGS = strict._FLAGS
_LIFECYCLE = strict._LIFECYCLE
_require, _int = strict._require, strict._int
strict_json, _json_file = strict.strict_json, strict._json_file
_compare, _nested_json = strict._compare, strict._nested_json
_inference = core._inference
check_observed_coordinates = strict.check_observed_coordinates
_csv_scalar = strict._csv_scalar


def _digest(value):
    return hashlib.sha256(canonical(value).encode()).hexdigest()


def _array_hash(array, dtype="<f8"):
    return hashlib.sha256(np.asarray(array, dtype=dtype).tobytes()).hexdigest()


def reconstruct_schedule(cell, rep, seed):
    """Independent address and arrival/switch race; no OU draw or model fit.

    Retain the published numeric calibration and RNG call order so byte-level
    replay is meaningful. There is no phase advancement through floor padding,
    burn-in, cap, resampling or renormalization.
    """
    cell = parse_cell(cell.to_dict())
    _require(cell.identity in {c.identity for c in full_design().cells}, "unplanned E4 source cell")
    _int(rep, "source replication")
    _int(seed, "source root seed")
    data = cell.to_dict()
    scientific = {k: v for k, v in data.items() if k != "cell_id"}
    identity = hashlib.sha256(json.dumps(scientific, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()
    _require(identity == cell.identity, "source cell hash does not bind its complete mapping")
    address = (8, *(int(identity[i:i+8], 16) for i in range(0, 64, 8)), rep)
    children = np.random.SeedSequence(seed, spawn_key=address).spawn(2)
    rng = np.random.default_rng(children[0])
    count = cell.n-1
    raw = np.empty(count, dtype=np.float64)
    phases = np.full((3, count), 0, dtype=np.int64)
    if cell.mechanism == "regular":
        raw[:] = cell.q/cell.theta
        phases[:2] = -1
    else:
        arrivals = (data["calm_arrival_rate"], data["burst_arrival_rate"])
        switches = (data["calm_switch_rate"], data["burst_switch_rate"])
        state = None
        for i in range(count):
            if i == 0 or cell.mechanism == "iid":
                state = int(float(rng.random()) >= data["calm_event_probability"])
            phases[0, i] = state
            duration = 0.0
            arriving = False
            while not arriving:
                total_rate = arrivals[state]+switches[state]
                wait = float(rng.exponential(1.0/total_rate))
                _require(math.isfinite(wait) and wait > 0, "invalid source wait; no replacement permitted")
                duration += wait
                _require(math.isfinite(duration), "source duration overflow; no replacement permitted")
                arriving = float(rng.random()) < arrivals[state]/total_rate
                if not arriving:
                    phases[2, i] += 1
                    state = 1-state
            raw[i], phases[1, i] = duration, state
    floor = cell.floor_fraction*(cell.q/cell.theta) if cell.cv >= cell.floor_cv_threshold else 0.0
    transformed = np.maximum(raw, floor)
    times = np.concatenate((np.zeros(1), np.cumsum(transformed)))
    _require(np.isfinite(times).all() and np.all(np.diff(times) > 0), "reconstructed source times invalid")
    diagnostic = {
        "family": {"mmpp": "event_stationary_mmpp", "iid": "iid_full_palm_marginal", "regular": "equidistant"}[cell.mechanism],
        "target_mean_gap": cell.q/cell.theta, "requested_cv": cell.cv,
        "floor_applied": cell.cv >= cell.floor_cv_threshold, "floor_value": float(floor),
        "source_seed": seed, "source_spawn_key": canonical(list(address)),
        "gap_spawn_key": canonical([*address, 0]), "path_spawn_key": canonical([*address, 1]),
        "raw_gaps_sha256": _array_hash(raw), "transformed_gaps_sha256": _array_hash(transformed),
        "phase_trace_sha256": hashlib.sha256(np.asarray(phases, dtype="<i8").tobytes()).hexdigest(),
        "times_sha256": _array_hash(times), "source_initial_phase": int(phases[0, 0]),
        "source_switch_count": int(sum(int(v) for v in phases[2])),
        "source_max_switches_per_gap": int(max(phases[2])),
        "source_phase_reset_count": count if cell.mechanism == "iid" else 1 if cell.mechanism == "mmpp" else 0,
    }
    for label, value in (("calm", 0), ("burst", 1)):
        for side, index in (("start", 0), ("end", 1)):
            diagnostic[f"source_{label}_{side}_count"] = int(np.count_nonzero(phases[index] == value))
    clipped = int(np.count_nonzero(raw < floor))
    diagnostic.update(clip_count=clipped, clip_fraction=clipped/count, raw_below_floor_count=clipped,
                      transformed_below_floor_count=int(np.count_nonzero(transformed < floor)),
                      transformed_at_floor_count=int(np.count_nonzero(transformed == floor)))
    for kind, vector in (("raw", raw), ("transformed", transformed)):
        extended = vector.astype(np.longdouble)
        mean = extended.sum(dtype=np.longdouble)/count
        diagnostic.update({f"{kind}_mean_gap": float(mean),
                           f"{kind}_cv": float(np.sqrt(np.sum((extended-mean)**2)/count)/mean),
                           f"{kind}_second_moment": float(np.sum(extended**2)/count),
                           f"{kind}_adjacent_product": float(np.sum(extended[:-1]*extended[1:])/(count-1)),
                           f"{kind}_min_gap": float(min(extended)), f"{kind}_span": float(extended.sum())})
        for fraction, label in ((.1, "0p1"), (.5, "0p5"), (1., "1"), (2., "2"), (5., "5")):
            diagnostic[f"{kind}_cdf_t_over_m_{label}"] = float(np.count_nonzero(vector <= fraction*cell.q/cell.theta)/count)
        # The published analytical law is a shared dependency, not a claimed
        # second derivation. Its independent rational/spectral tests precede E4.
        theory = source.moments(cell, transformed=kind == "transformed")
        diagnostic.update({f"expected_{kind}_mean_gap": theory["mean"],
                           f"expected_{kind}_second_moment": theory["second_moment"],
                           f"expected_{kind}_cv": theory["cv"],
                           f"expected_{kind}_adjacent_product": source.adjacent_product(cell, transformed=kind == "transformed")})
    actual = np.diff(times).astype(np.longdouble)
    average = np.sum(actual)/count
    diagnostic.update(realized_mean_gap=float(average),
                      realized_cv=float(np.sqrt(np.sum((actual-average)**2)/count)/average),
                      realized_span=float(times[-1]-times[0]), min_gap=float(np.min(actual)),
                      expected_mean_gap_after_floor=diagnostic["expected_transformed_mean_gap"],
                      expected_clip_fraction=source.cdf(cell, floor))
    _require(set(diagnostic) == SOURCE_FIELDS-{"observations_sha256"}, "reconstructed source schema changed")
    return times, diagnostic, children[1]


def independent_source(cell, rep, seed):
    times, diagnostic, child = reconstruct_schedule(cell, rep, seed)
    values = simulate_ou(cell.theta, cell.mu, cell.sigma, times, np.random.default_rng(child))
    sample = Sample(values, times)
    diagnostic["observations_sha256"] = _array_hash(sample.values)
    return sample, diagnostic


def _validate_row(row, cell, rep, method):
    _require(type(row) is dict and set(row) == RAW_FIELDS, "raw record fields differ from the complete schema")
    _require(type(row["rep"]) is int and row["rep"] == rep and type(row["estimator"]) is str
             and row["estimator"] == method, "raw Cartesian key or ordering changed")
    scientific = {key: row[key] for key in _CELL_FIELDS}
    _compare(scientific, cell.to_dict(), "canonical actual cell", exact=True)
    _require(canonical(scientific) == canonical(cell.to_dict()), "canonical Cell representation changed")
    for key in _FLAGS | {"source_valid", "floor_applied"}:
        _require(type(row[key]) is bool, "raw flag must be Boolean: " + key)
    _require(row["source_valid"] and row["source_error"] == "" and row["pipeline_error"] == ""
             and type(row["source_error"]) is str and type(row["pipeline_error"]) is str,
             "saved source or pipeline failure requires implementation review")
    for key in _TEXT_INFERENCE:
        _require(type(row[key]) is str and bool(row[key]), "missing raw text diagnostic: " + key)
    _require(row["point_status"] in {"valid_interior", "boundary", "boundary_suspected", "degenerate", "numerical_failure"},
             "unknown registered point-fit status")
    for key in INFERENCE_FIELDS - _TEXT_INFERENCE - _FLAGS:
        _require(row[key] is None or type(row[key]) is float and math.isfinite(row[key]),
                 "inference numbers must be finite floating values or JSON null: " + key)
    _nested_json(row)
    for field in ("point_diagnostics_json", "interval_diagnostics_json"):
        _require(type(strict_json(row[field], field)) is dict, "diagnostic object required")
    point = strict_json(row["point_diagnostics_json"])
    for key in ("optimizer_success", "valid_for_point_summary"):
        _require(type(point.get(key)) is bool, "point diagnostic validity flag type changed")
    interval = strict_json(row["interval_diagnostics_json"])
    _require(type(interval.get("conditions_on_initial_value")) is bool, "conditional-likelihood flag type changed")
    _inference(row)


def load_cell(output, cell, design, identity):
    """Own strict loader; no producer CellStore, reporting or audit loader."""
    directory = Path(output)/"cells"/cell.identity
    state, payload = _json_file(directory/"cell_manifest.json")
    expected = {"schema_version", "run_identity", "cell", "replications", "status",
                "completed_replications", "partitions", "recovered_orphans"}
    _require(set(state) == expected and _int(state["schema_version"], "cell schema", 1) == 1,
             "cell manifest schema changed")
    _require(state["run_identity"] == identity and canonical(state["cell"]) == canonical(cell.to_dict()),
             "cell manifest identity changed")
    _require(_int(state["replications"], "cell plan", 1) == design.replications
             and _int(state["completed_replications"], "completed sources", 1) == design.replications
             and state["status"] == "complete", "incomplete raw cell")
    _require(type(state["partitions"]) is list and type(state["recovered_orphans"]) is list,
             "cell partition/recovery history must be a list")
    rows, parts, start = [], [], 0
    for part in state["partitions"]:
        _require(type(part) is dict and set(part) == {"filename", "start", "stop", "rows", "sha256"},
                 "partition descriptor schema changed")
        begin, stop = _int(part["start"], "start"), _int(part["stop"], "stop", 1)
        _require(begin == start and begin < stop <= design.replications, "partition range or contiguous plan changed")
        name = f"part-{begin:06d}-{stop:06d}.jsonl.gz"
        _require(part["filename"] == name and file_digest(directory/name) == part["sha256"], "partition name/hash changed")
        with gzip.open(directory/name, "rt", encoding="utf-8", newline="") as stream:
            current = []
            for line in stream:
                row = strict_json(line, name)
                _require(line == canonical(row) + "\n", "noncanonical or empty raw JSONL row")
                current.append(row)
        expected_keys = [(rep, method) for rep in range(begin, stop) for method in METHODS]
        _require(len(current) == _int(part["rows"], "partition rows", 1) == len(expected_keys), "partition rows changed")
        for row, (rep, method) in zip(current, expected_keys):
            _validate_row(row, cell, rep, method)
        rows.extend(current)
        parts.append(dict(part))
        start = stop
    _require(start == design.replications, "missing planned sources")
    _require({p.name for p in directory.iterdir()} == {"cell_manifest.json", *(p["filename"] for p in parts)},
             "unexpected cell evidence")
    names = {p["filename"] for p in parts}
    recovered = set()
    for entry in state["recovered_orphans"]:
        _require(type(entry) is dict and set(entry) == {"filename", "verified_utc", "method"}
                 and entry["filename"] in names and entry["filename"] not in recovered
                 and type(entry["verified_utc"]) is str and bool(entry["verified_utc"])
                 and entry["method"] == "exact_byte_equal_deterministic_replay", "invalid orphan-adoption history")
        recovered.add(entry["filename"])
    return rows, {"manifest_sha256": hashlib.sha256(payload).hexdigest(), "partitions": parts}


class FailureExport:
    """Streaming exact-schema failure verification in canonical raw row order."""

    def __init__(self, path):
        self.path = Path(path)
        self.stream = None
        self.reader = None
        self.count = 0

    def __enter__(self):
        self.stream = gzip.open(self.path, "rt", encoding="utf-8", newline="")
        self.reader = csv.reader(self.stream)
        try:
            header = next(self.reader, None)
            _require(header == sorted(RAW_FIELDS) and len(set(header)) == len(header),
                     "failure export must retain unique complete canonical headers")
            self.header = header
        except BaseException:
            self.stream.close()
            raise
        return self

    def check_rows(self, rows):
        for row in rows:
            if row["point_valid"] and row["interval_valid"] and not row["pipeline_error"]:
                continue
            exported = next(self.reader, None)
            _require(exported is not None and len(exported) == len(self.header), "failure export omitted or malformed a row")
            for column, value in zip(self.header, exported):
                _csv_scalar(value, row[column], f"{row['cell_id']}/{row['rep']}/{row['estimator']}/{column}")
            self.count += 1

    def finish(self, expected_count):
        _require(self.count == _int(expected_count, "failure manifest rows") and next(self.reader, None) is None,
                 "failure export count, order or complete membership changed")
        return self.count

    def __exit__(self, *_):
        self.stream.close()


def verify_provenance(metadata, design, extra_files=()):
    """Read-only exact dependency/runtime/freeze checks, valid for sealed runs."""
    validate_design(design)
    expected_files = set(CORE_FILES) | set(FILES) | set(extra_files)
    _require(all(type(name) is str and not Path(name).is_absolute() and ".." not in Path(name).parts
                 for name in expected_files), "unsafe dependency path")
    provenance = metadata["provenance"]
    _require(type(provenance) is dict and type(provenance.get("source_files")) is dict
             and set(provenance["source_files"]) == expected_files, "explicit audit dependency membership changed")
    _require(canonical(provenance.get("design")) == canonical(design.to_dict()), "provenance design representation changed")
    verify(provenance, expected_design=design)
    frozen = design.scope == "exploratory"
    _require(metadata.get("require_frozen") is frozen, "recorded freeze requirement changed")
    dirty = provenance.get("uncommitted_dependency_files")
    _require(type(dirty) is list and all(type(name) is str and name in expected_files for name in dirty)
             and len(set(dirty)) == len(dirty), "invalid recorded dependency-dirty list")
    commit = provenance.get("code_commit")
    _require(type(commit) is str and re.fullmatch(r"[0-9a-f]{40}", commit) is not None, "invalid recorded Git commit")
    if frozen:
        _require(dirty == [], "full exploratory code was not frozen before execution")
        for name, expected in provenance["source_files"].items():
            completed = subprocess.run(["git", "show", f"{commit}:{name}"], cwd=ROOT, capture_output=True)
            _require(completed.returncode == 0 and hashlib.sha256(completed.stdout).hexdigest() == expected,
                     "recorded Git freeze does not contain dependency: " + name)
    return provenance


def _enrich(rows, cell):
    metadata = cell.to_dict()
    for row in rows:
        for key in set(row) & set(metadata):
            _require(canonical(row[key]) == canonical(metadata[key]), "independent summary cell metadata differs")
        yield {**row, **metadata}


def independent_contrasts(summary, design, kind):
    """Independent references and own-eligible-subset contrasts, not paired rows."""
    validate_design(design)
    _require(kind in {"regular", "matched_iid"}, "unknown independent comparison kind")
    index = {(r["cell_id"], r["estimator"], r["parameter"]): r for r in summary}
    _require(len(index) == len(summary) == 9*len(design.cells), "incomplete or duplicate independent summary keys")
    expected = {(c.identity, m, p) for c in design.cells for m in METHODS for p in PARAMETERS}
    _require(set(index) == expected, "unplanned independent summary keys")
    regular = next(c for c in design.cells if c.mechanism == "regular")
    iid = {c.ratio: c for c in design.cells if c.mechanism == "iid"}
    result = []
    for cell in design.cells:
        if cell.mechanism == "regular" or kind == "matched_iid" and cell.mechanism != "mmpp":
            continue
        reference = regular if kind == "regular" else iid[cell.ratio]
        for method in METHODS:
            for parameter in PARAMETERS:
                a, b = (index[(c.identity, method, parameter)] for c in (cell, reference))
                row = {**cell.to_dict(), **{"reference_"+k: v for k, v in reference.to_dict().items()},
                       "estimator": method, "parameter": parameter, "truth": a["truth"],
                       "comparison_kind": kind, "comparison_type": "independent_cells",
                       "point_subset": "each_cell_own_valid_point_fits",
                       "coverage_subset": "all_planned_paths_unavailable_counts_as_not_covering",
                       "difference_direction": "cell_minus_reference", "ratio_direction": "cell_over_reference",
                       "matching": "same_full_event_gap_marginal_and_floor_transform_independent_paths" if kind == "matched_iid"
                       else "same_n_truth_nominal_mean_gap_and_floor_settings_not_realized_span"}
                for field in reporting.COMPONENT_FIELDS:  # declarative output contract, not arithmetic
                    row["cell_"+field], row["reference_"+field] = a[field], b[field]
                for name, value, se in (("excess_bias", "bias", "mcse_bias"),
                                        ("operational_coverage_difference", "coverage_operational", "operational_mcse"),
                                        ("conditional_coverage_difference", "coverage_conditional", "conditional_mcse")):
                    core._difference(row, name, a[value]-b[value], math.hypot(a[se], b[se]))
                core._relative_difference(row, "excess_bias", a["truth"])
                log_se = (math.hypot(a["mcse_rmse"]/a["rmse"], b["mcse_rmse"]/b["rmse"])
                          if a["rmse"] > 0 and b["rmse"] > 0 else math.nan)
                core._ratio(row, a["rmse"], b["rmse"], log_se)
                result.append(row)
    return result


def independent_source_summary(rows, cell, design):
    """One source per path; arithmetic uses independent paths, not 3R fit rows."""
    sources = [r for r in rows if r["estimator"] == "exact"]
    _require(len(sources) == design.replications, "source summary is not complete")
    constants = ("family", "target_mean_gap", "requested_cv", "floor_applied", "floor_value",
                 *sorted(k for k in SOURCE_FIELDS if k.startswith("expected_")))
    for row in sources:
        for key in constants:
            _compare(row[key], sources[0][key], "within-cell source constant", exact=True)
    result = {**cell.to_dict(), **{k: sources[0][k] for k in constants},
              "n_planned": design.replications, "n_unique_sources": len(sources),
              "n_source_valid": len(sources), "n_fit_records": len(rows),
              "source_metric_mcse": "across_independent_paths_sample_sd_ddof1_over_sqrtR",
              "raw_metric_scope": "generated_raw_gap_vector_before_cumsum",
              "transformed_metric_scope": "posthoc_floored_gap_vector_before_cumsum",
              "realized_metric_scope": "differences_of_stored_cumulative_times"}
    for key in reporting.SOURCE_METRICS:  # shared names, independently recomputed reductions
        values = np.array([r[key] for r in sources], dtype=float)
        _require(np.isfinite(values).all(), "nonfinite source metric")
        result.update({key+"_n_available": len(values), key+"_mean": math.fsum(values)/len(values),
                       key+"_mcse": core._mean_se(values), key+"_min": float(min(values)), key+"_max": float(max(values))})
        if key.endswith("_count") or key == "source_max_switches_per_gap":
            result[key+"_total"] = sum(_int(r[key], key) for r in sources)
    for label, phase in (("calm", 0), ("burst", 1), ("regular", -1)):
        indicators = np.array([int(r["source_initial_phase"] == phase) for r in sources], dtype=float)
        result.update({f"initial_phase_{label}_count": int(indicators.sum()),
                       f"initial_phase_{label}_fraction": float(indicators.mean()),
                       f"initial_phase_{label}_mcse": core._mean_se(indicators)})
        result["expected_initial_phase_"+label] = (float(phase == -1) if cell.mechanism == "regular" else
                                                   0. if phase == -1 else cell.to_dict()[label+"_event_probability"])
    return [result]


def _layout(output, design):
    _require(output.is_dir() and not any(p.is_symlink() for p in (output, *output.parents)), "audit path is missing or linked")
    required = {"run_metadata.json", "resolved_design.json", "cells", "table_manifest.json", "table_binding.json",
                "failed_fits.csv.gz", *TABLE_KEYS}
    allowed = required | {".execution.lock", "audit.json", "visual_review.json", "figures"}
    _require(required <= {p.name for p in output.iterdir()} <= allowed, "missing or unexpected E4 evidence")
    for path in output.rglob("*"):
        _require(not path.is_symlink() and (path.is_file() or path.is_dir()), "nonregular or linked evidence")
    for path in output.iterdir():
        _require(path.is_dir() if path.name in {"cells", "figures"} else path.is_file(), "root evidence type changed")
    _require({p.name for p in (output/"cells").iterdir()} == {c.identity for c in design.cells}, "incomplete cell set")
    for path in (output/"cells").iterdir():
        _require(path.is_dir(), "cell path must be a directory")
        for child in path.iterdir():
            _require(child.is_file() and (child.name == "cell_manifest.json" or
                     re.fullmatch(r"part-[0-9]{6}-[0-9]{6}\.jsonl\.gz", child.name)), "unknown cell artifact")
    figures = output/"figures"
    if figures.exists() and any(figures.iterdir()):
        figure, _ = _json_file(figures/"figure_manifest.json")
        mapping = figure.get("figures")
        _require(type(mapping) is dict and set(mapping) | {"figure_manifest.json"} == {p.name for p in figures.iterdir()}, "figure file set changed")
        for name, digest in mapping.items():
            _require(type(name) is str and Path(name).name == name and Path(name).suffix in {".png", ".pdf"}
                     and (figures/name).is_file() and file_digest(figures/name) == digest, "figure artifact changed")


def _snapshot(output, design):
    _layout(output, design)
    return {"directories": sorted(str(p.relative_to(output)) for p in output.rglob("*") if p.is_dir()),
            "files": {str(p.relative_to(output)): file_digest(p) for p in sorted(output.rglob("*")) if p.is_file()}}


def _artifact_map(snapshot):
    return {name: digest for name, digest in snapshot["files"].items()
            if name not in {"run_metadata.json", ".execution.lock", "audit.json", "visual_review.json"}
            and not name.startswith("figures/")}


def semantic_metadata_digest(metadata):
    return _digest({k: v for k, v in metadata.items() if k not in _LIFECYCLE})


def _metadata(output, design, extra_files):
    metadata, _ = _json_file(output/"run_metadata.json")
    expected = {"schema_version": 1, "project": "ou-week8-e4", "execution_kind": "e4",
                "execution_contract": "week8_guarded_cell_execution_v1", "scientific_overrides": [],
                "procedure": PROCEDURE, "source_policy": source_policy(design)}
    for key, value in expected.items():
        _compare(metadata.get(key), value, "root/"+key, exact=True)
    _require(metadata.get("status") in {"raw_complete_pending_audit", "audited_pending_visual_review", "completed"}, "incomplete or review-required run")
    _require(_int(metadata.get("pipeline_errors"), "pipeline errors") == 0, "pipeline errors block audit")
    _require(all(type(metadata.get(k)) is list for k in ("execution_attempts", "resume_events")), "missing lifecycle history")
    completed = {c.identity: {"cell_id": c.identity, "paths": design.replications,
                             "rows": 3*design.replications, "pipeline_errors": 0} for c in design.cells}
    _compare(metadata.get("completed_cells"), completed, "complete planned cell counts", exact=True)
    provenance = verify_provenance(metadata, design, extra_files)
    resolved, payload = _json_file(output/"resolved_design.json")
    _require(canonical(parse_design(resolved).to_dict()) == canonical(design.to_dict()), "resolved design changed")
    binding = {"run_identity": provenance["identity_sha256"], "design_sha256": design.digest,
               "resolved_design_sha256": hashlib.sha256(payload).hexdigest(), "source_policy": source_policy(design), "cells": {}}
    return metadata, binding


def _table_scalar(actual, expected, label):
    if expected is None or isinstance(expected, (float, np.floating)) and math.isnan(expected):
        _require(actual == "", "table missing value changed: "+label)
    elif isinstance(expected, (bool, np.bool_)):
        _require(actual == str(bool(expected)), "table Boolean changed: "+label)
    elif isinstance(expected, (int, np.integer)):
        _require(actual == str(int(expected)), "table integer changed: "+label)
    elif isinstance(expected, str):
        _require(actual == expected, "table text changed: "+label)
    else:
        _require(actual != "" and math.isfinite(float(actual)) and math.isfinite(float(expected))
                 and math.isclose(float(actual), float(expected), rel_tol=RTOL, abs_tol=ATOL), "table arithmetic changed: "+label)


def _verify_table(path, name, expected):
    _require(bool(expected) and all(set(row) == set(expected[0]) for row in expected), "independent table schema inconsistent")
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.reader(stream)
        header = next(reader, None)
        _require(header == sorted(expected[0]) and tuple(header) == reporting.TABLE_COLUMNS[name], "table exact schema changed: "+name)
        actual = []
        for values in reader:
            _require(len(values) == len(header), "malformed table row")
            actual.append(dict(zip(header, values)))
    keys = TABLE_KEYS[name]
    received = {tuple(r[k] for k in keys): r for r in actual}
    desired = {tuple(r[k] for k in keys): r for r in expected}
    _require(len(received) == len(actual) == len(expected) == len(desired) and received.keys() == desired.keys(), "table complete key set changed: "+name)
    for key, row in desired.items():
        for column, value in row.items():
            _table_scalar(received[key][column], value, f"{name}/{key}/{column}")
    return len(actual)


def _table_manifest(output):
    manifest, _ = _json_file(output/"table_manifest.json")
    _require(set(manifest) == set(TABLE_KEYS) | {"failed_fits.csv.gz"}, "table manifest must declare exactly nine artifacts")
    for name, entry in manifest.items():
        _require(type(entry) is dict and set(entry) == {"rows", "sha256"}, "table manifest descriptor changed")
        _int(entry["rows"], "table rows")
        _require(file_digest(output/name) == entry["sha256"], "table hash changed")
    return manifest


class _CachedSource:
    def __init__(self, cell, rep, sample, diagnostics):
        self.cell, self.rep, self.value, self.diagnostics = cell, rep, sample, diagnostics

    def sample(self, mapping, rep):
        _require(canonical(mapping) == canonical(self.cell.to_dict()) and type(rep) is int and rep == self.rep,
                 "endpoint replay requested a different source")
        return self.value, dict(self.diagnostics)


def _verify_sources(rows, cell, design, replay):
    coordinates = {"point_coordinate_rows": 0, "interval_coordinate_rows": 0, "coordinate_scalars_checked": 0}
    endpoints = []
    services = {method: registered_wald(method) for method in METHODS} if replay else None
    for rep in range(design.replications):
        sample, diagnostic = independent_source(cell, rep, design.seed)
        triplet = rows[3*rep:3*rep+3]
        first = {key: triplet[0][key] for key in SOURCE_FIELDS}
        for row in triplet:
            saved = {key: row[key] for key in SOURCE_FIELDS}
            _compare(saved, first, "shared source triplet", exact=True)
            _compare(saved, diagnostic, "independent E4 source reconstruction")
            for key, count in check_observed_coordinates(row, sample).items():
                coordinates[key] += count
        if rep in {0, design.replications-1}:
            descriptor = {"cell_id": cell.identity, "replication": rep,
                          "times_sha256": diagnostic["times_sha256"], "observations_sha256": diagnostic["observations_sha256"],
                          "estimators_checked": list(METHODS)}
            endpoints.append(descriptor)
            if replay:
                fitted = strict_json(canonical(source_records(cell, rep, _CachedSource(cell, rep, sample, diagnostic), services)))
                _require(len(fitted) == 3, "endpoint replay did not retain all estimators")
                for actual, wanted in zip(triplet, fitted):
                    _require(wanted["pipeline_error"] == "", "endpoint fit replay raised a pipeline exception")
                    _compare(actual, wanted, f"endpoint/{cell.identity}/{rep}/{actual['estimator']}")
    return coordinates, endpoints


def _closed_bindings(output, metadata, snapshot):
    closed = metadata["status"] in {"audited_pending_visual_review", "completed"}
    saved = None
    if (output/"audit.json").exists():
        saved, _ = _json_file(output/"audit.json")
        _require(saved.get("audit_contract") == AUDIT_CONTRACT and saved.get("status") == "passed", "saved E4 audit contract changed")
        if closed or "audit_sha256" in metadata:
            _require(metadata.get("audit_sha256") == file_digest(output/"audit.json"), "root audit seal changed")
    _require(not closed or saved is not None and saved.get("full_completion_audit") is True, "closed state requires full endpoint audit")
    if metadata["status"] == "completed":
        complete = {name: digest for name, digest in snapshot["files"].items() if name not in {"run_metadata.json", ".execution.lock"}}
        _require(metadata.get("artifact_sha256") == complete, "completed artifact seal changed")
        review, _ = _json_file(output/"visual_review.json")
        figure, _ = _json_file(output/"figures"/"figure_manifest.json")
        pngs = {name for name in figure["figures"] if name.endswith(".png")}
        _require(metadata.get("visual_review_sha256") == file_digest(output/"visual_review.json")
                 and bool(pngs) and review.get("all_panels_legible") is True
                 and type(review.get("reviewed_pngs")) is list
                 and all(type(name) is str for name in review["reviewed_pngs"])
                 and len(review["reviewed_pngs"]) == len(set(review["reviewed_pngs"]))
                 and set(review["reviewed_pngs"]) == pngs
                 and review.get("figure_manifest_sha256") == file_digest(output/"figures"/"figure_manifest.json")
                 and type(review.get("reviewer")) is str and bool(review["reviewer"].strip()), "completed visual review binding changed")
    return saved


def _audit(output, design, extra_files, replay, snapshot):
    metadata, binding = _metadata(output, design, extra_files)
    saved_audit = _closed_bindings(output, metadata, snapshot)
    manifest = _table_manifest(output)
    desired = {name: [] for name in TABLE_KEYS}
    coordinates = {"point_coordinate_rows": 0, "interval_coordinate_rows": 0, "coordinate_scalars_checked": 0}
    endpoints = []
    partition_count = point_count = interval_count = 0
    with FailureExport(output/"failed_fits.csv.gz") as failures:
        for cell in sorted(design.cells, key=lambda c: c.identity):
            rows, cell_binding = load_cell(output, cell, design, binding["run_identity"])
            binding["cells"][cell.identity] = cell_binding
            partition_count += len(cell_binding["partitions"])
            counts, cell_endpoints = _verify_sources(rows, cell, design, replay)
            endpoints.extend(cell_endpoints)
            for key, value in counts.items():
                coordinates[key] += value
            failures.check_rows(sorted(rows, key=lambda r: (r["cell_id"], r["rep"], r["estimator"])))
            point_count += sum(r["point_valid"] for r in rows)
            interval_count += sum(r["interval_valid"] for r in rows)
            desired["parameter_summary.csv"].extend(_enrich(core._parameter_summaries(rows, cell, design), cell))
            desired["paired_comparisons.csv"].extend(_enrich(core._paired_comparisons(rows, cell, design), cell))
            spacing, outcomes, numerical = core._other_summaries(rows, cell, design)
            for name, values in (("spacing_summary.csv", spacing), ("outcome_counts.csv", outcomes), ("numerical_diagnostics.csv", numerical)):
                desired[name].extend(_enrich(values, cell))
            desired["source_summary.csv"].extend(independent_source_summary(rows, cell, design))
        failure_count = failures.finish(manifest["failed_fits.csv.gz"]["rows"])
    for kind in ("regular", "matched_iid"):
        desired[kind+"_comparisons.csv"] = independent_contrasts(desired["parameter_summary.csv"], design, kind)
    for name, expected in desired.items():
        _require(_verify_table(output/name, name, expected) == manifest[name]["rows"], "table row count disagrees with manifest")
    saved_binding, _ = _json_file(output/"table_binding.json")
    expected_binding = {"schema_version": 1, "project": "ou-week8-e4-reporting", "reporting_contract": "week8_saved_raw_tables_v1",
                        "run_identity": binding["run_identity"], "design_sha256": design.digest,
                        "raw_inputs": binding, "raw_input_sha256": _digest(binding)}
    _compare(saved_binding, expected_binding, "complete raw-to-table binding", exact=True)
    verify_provenance(metadata, design, extra_files)
    paths = len(design.cells)*design.replications
    replayed = endpoints if replay else []
    report = {"schema_version": 1, "audit_contract": AUDIT_CONTRACT, "status": "passed", "full_completion_audit": replay,
              "design_sha256": design.digest, "run_identity_sha256": binding["run_identity"],
              "code_commit": metadata["provenance"]["code_commit"],
              "source_dependency_files_verified": len(metadata["provenance"]["source_files"]),
              "audited_artifact_sha256": _artifact_map(snapshot), "audited_run_metadata_sha256": semantic_metadata_digest(metadata),
              "raw_input_sha256": _digest(binding), "scientific_overrides": [], "source_policy": source_policy(design),
              "cells_verified": len(design.cells), "partitions_verified": partition_count,
              "planned_paths": paths, "planned_source_realizations": paths, "planned_independent_source_paths": paths,
              "planned_fit_records": 3*paths, "source_replay_performed": True, "source_paths_verified": paths,
              "source_array_hashes_verified": 2*paths, "schedule_trace_hashes_verified": 3*paths,
              "point_valid_records": point_count, "interval_valid_records": interval_count,
              "retained_failure_records": failure_count, "pipeline_errors": 0, **coordinates,
              "tables_verified": {**{name: len(values) for name, values in desired.items()}, "failed_fits.csv.gz": failure_count},
              "replay_performed": replay, "replay_selection": "replications 0 and R-1 in every E4 cell, all estimators",
              "replayed_source_paths": len(replayed), "replayed_fits": 3*len(replayed), "replay_records": replayed,
              "tolerances": {"arithmetic_rtol": RTOL, "arithmetic_atol": ATOL, "replay_rtol": RTOL, "replay_atol": ATOL,
                             "source_hashes": "exact SHA-256 of both arrays, both gap vectors and full phase/switch trace",
                             "covariance_rtol": 2e-9, "covariance_atol": "2e-12 * max(1, expected matrix Frobenius norm)"},
              "scope": {"all_sources": "independent namespace-8 addressing and competing arrival/switch race, every path and all diagnostics",
                        "all_records": "exact keys/types/strict nested JSON, recorded coordinates, fine-Hessian/covariance identities, validity gates and unclipped endpoints",
                        "all_tables": "independent core summary formulas, E4 own-subset contrasts and source reductions; full ordered failure export",
                        "endpoint_fits": "complete original-service replay at fixed endpoints, including nested diagnostics and missingness",
                        "limitations": "calibrated cell definitions, analytical gap-moment routines and exact OU simulator are shared; endpoint fits share original inference and serializer; full fit reproduction/coarse-Hessian reconstruction is endpoint-only; actual visual judgment is outside this numerical audit"}}
    strict._verify_saved_report(saved_audit, report, endpoints, replay=replay)
    json.dumps(report, allow_nan=False)
    return report


def audit_run(output, design, *, extra_files=(), replay=True):
    """Audit complete E4 evidence without writing files, locks or lifecycle state.

    replay=False still reconstructs every source; it skips endpoint model fits
    and cannot produce a full-completion report. Changed evidence is never
    regenerated or replaced. Existing sealed reports are checked, not trusted.
    """
    _require(type(replay) is bool, "replay must be Boolean")
    validate_design(design)
    _require(not isinstance(extra_files, (str, bytes)), "extra_files must be repository-relative names")
    extra_files = tuple(extra_files)
    output = Path(output).absolute()
    before = _snapshot(output, design)
    try:
        return _audit(output, design, extra_files, replay, before)
    except ValueError:
        raise
    except (KeyError, TypeError, AttributeError, OSError, OverflowError, IndexError, np.linalg.LinAlgError) as exc:
        raise ValueError(f"E4 artifact audit failed: {type(exc).__name__}: {exc}") from exc
    finally:
        _require(_snapshot(output, design) == before, "artifact/root bytes changed during read-only E4 audit")
