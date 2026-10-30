"""Guarded Week 7 orchestration around the unchanged durable cell engine.

This module chooses only the two frozen exploratory designs or their separate
validation namespaces.  It does not alter any source, estimator, inference gate
or partition format.  Raw completion remains pending independent reporting,
audit and review.  Read-only preflight never creates output or lock files.
"""
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
import gzip
import json
import multiprocessing
from pathlib import Path
import re
import time

from ..experiments.audit import TABLE_KEYS
from ..experiments.design import ESTIMATORS, RunDesign, integer
from ..experiments.execution import PROCEDURE, _execute_cell
from ..experiments.provenance import canonical, capture, file_digest, utc_now
from ..experiments.storage import CellStore, atomic_json, run_lock
from .floor_source import ANCHOR_FLOOR_FRACTION, COUPLING_POLICY, PairedFloorSource, anchor_cell
from .path_diagnostics import (FILES as GUARD_FILES, _ManifestGuard, _check_frozen,
                               _check_incident_block, _reject_archive_ancestors)
from .selection import floor_design, followup_design


FILES = tuple(sorted(set(GUARD_FILES) | {
    "ou_irregular/week7/floor_source.py", "ou_irregular/week7/execution.py",
    "docs/week7_execution.md"}))
SYNC_SCRIPT = Path("/usr/local/scripts/sync_share.sh")
SYNC_ROOT = Path("/workspace/scratch")
EXECUTION_CONTRACT = "week7_guarded_cell_execution_v1"
FLOOR_SOURCE_ADDRESS = "SeedSequence(root,spawn_key=(4,*8_registered_floor_anchor_cell_words,rep)).spawn(2)"
_OPEN = {"running_incomplete", "incomplete"}
_CACHED = {"raw_complete_pending_audit", "implementation_review_required"}
_CLOSED = {"audited_pending_visual_review", "completed"}
_RAW_FILES = {".execution.lock", "run_metadata.json", "resolved_design.json"}
_DERIVED_FILES = set(TABLE_KEYS) | {"failed_fits.csv.gz", "table_manifest.json", "audit.json",
                                   "visual_review.json", "epsilon_comparisons.csv", "epsilon_table_manifest.json"}
_PARTITION = re.compile(r"part-(\d{6})-(\d{6})\.jsonl\.gz")


def require_design(kind, design):
    """Exact full design or exact separate-namespace validation(R), including order."""
    if kind not in {"followup", "floor"} or not isinstance(design, RunDesign):
        raise ValueError("kind must be followup/floor with a canonical RunDesign")
    registered = followup_design() if kind == "followup" else floor_design()
    expected = registered.validation(design.replications) if design.scope == "validation" else registered
    if canonical(design.to_dict()) != canonical(expected.to_dict()):
        raise ValueError("Week 7 design differs from its fixed full or separate validation design")
    return design


def source_policy(kind, design):
    """Declared planned addressing/counts, separate from completed-cell accounting."""
    require_design(kind, design)
    anchors = {cell.identity: (anchor_cell(cell).identity if kind == "floor" else cell.identity)
               for cell in design.cells}
    return {"name": COUPLING_POLICY if kind == "floor" else "registered_ou_source_v1",
            "root_seed": design.seed, "actual_to_seed_anchor": anchors,
            "seed_anchor_floor_fraction": ANCHOR_FLOOR_FRACTION if kind == "floor" else None,
            "planned_source_realizations": len(design.cells) * design.replications,
            "planned_independent_coupling_groups": len(set(anchors.values())) * design.replications,
            "planned_fit_records": len(design.cells) * design.replications * len(ESTIMATORS),
            "stationary_initialization": True, "renormalize_floored_gaps": False}


def _procedure(kind):
    return {**PROCEDURE, **({"source_address": FLOOR_SOURCE_ADDRESS} if kind == "floor" else {})}


def _unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError("duplicate JSON object key")
        result[key] = value
    return result


def _loads(raw):
    def nonstandard(value):
        raise ValueError("nonstandard nonfinite JSON constant: " + value)
    return json.loads(raw, object_pairs_hook=_unique_object, parse_constant=nonstandard)


def _json_file(path, *, canonical_required=True):
    raw = Path(path).read_bytes()
    value = _loads(raw)
    if not isinstance(value, dict):
        raise ValueError("manifest must be a JSON object")
    if canonical_required and raw != (canonical(value) + "\n").encode():
        raise ValueError("manifest is not the canonical saved JSON representation")
    return value, raw


def _output_path(output_dir):
    supplied = Path(output_dir).absolute()
    if any(path.is_symlink() for path in (supplied, *supplied.parents)):
        raise ValueError("symbolic links are not permitted in mutable output paths")
    output = supplied.resolve()
    if SYNC_SCRIPT.exists() and output.is_relative_to(SYNC_ROOT.resolve()):
        raise ValueError("mutable runs in synchronized /workspace/scratch are blocked; use isolated /tmp output")
    _reject_archive_ancestors(output)
    _check_incident_block(output)
    return output


def _root_layout(output, design):
    if not output.exists():
        return
    if not output.is_dir():
        raise ValueError("execution output must be a directory")
    allowed = _RAW_FILES | _DERIVED_FILES | {"cells", "figures"}
    for path in output.iterdir():
        if path.name not in allowed or path.is_symlink():
            raise ValueError("unknown or symbolic-link root artifact")
        if path.name in {"cells", "figures"}:
            if not path.is_dir():
                raise ValueError("artifact directory has an invalid object type")
        elif not path.is_file():
            raise ValueError("root artifact is not a regular file")
    cells = output / "cells"
    if cells.exists():
        expected = {cell.identity for cell in design.cells}
        for path in cells.iterdir():
            if path.name not in expected or path.is_symlink() or not path.is_dir():
                raise ValueError("unknown cell directory or symbolic link")
            for child in path.iterdir():
                if child.is_symlink() or not child.is_file() or (child.name != "cell_manifest.json" and not _PARTITION.fullmatch(child.name)):
                    raise ValueError("unknown cell artifact or symbolic link")
    figures = output / "figures"
    if figures.exists():
        for path in figures.iterdir():
            if path.is_symlink() or not path.is_file():
                raise ValueError("unknown nonregular figure artifact")
        if (figures / "figure_manifest.json").exists():
            manifest, _ = _json_file(figures / "figure_manifest.json")
            files = manifest.get("figures")
            if not isinstance(files, dict) or set(files) | {"figure_manifest.json"} != {p.name for p in figures.iterdir()}:
                raise ValueError("figure layout does not equal its recorded manifest")
            for name, digest in files.items():
                if Path(name).name != name or file_digest(figures / name) != digest:
                    raise ValueError("figure identity or digest changed")
        elif any(path.suffix not in {".png", ".pdf"} for path in figures.iterdir()):
            # Atomic figure publication can leave a partial set of completed
            # PNG/PDF files before its manifest.  Raw preflight preserves those
            # bytes without claiming figure completion.  The later plot/audit
            # stage must validate exact planned names and byte-equal adoption.
            raise ValueError("unknown partial figure artifact without a manifest")
    table_manifest = output / "table_manifest.json"
    if table_manifest.exists():
        manifest, _ = _json_file(table_manifest)
        for name, entry in manifest.items():
            if name not in _DERIVED_FILES or not isinstance(entry, dict):
                raise ValueError("unknown derived table manifest entry")
            integer(entry.get("rows"), "table row count", 0)
            if file_digest(output / name) != entry.get("sha256"):
                raise ValueError("derived table digest changed")


class SavedPipelineError(ValueError):
    """Saved implementation failures require review, never an automatic retry."""


def _strict_cell_manifest(path, cell, design, identity):
    state, _ = _json_file(path)
    expected = {"schema_version", "run_identity", "cell", "replications", "status",
                "completed_replications", "partitions", "recovered_orphans"}
    if set(state) != expected:
        raise ValueError("cell manifest fields differ from the storage contract")
    integer(state["schema_version"], "cell schema version", 1)
    integer(state["replications"], "cell replication count", 1)
    integer(state["completed_replications"], "cell committed prefix", 0)
    if (state["schema_version"] != 1 or state["run_identity"] != identity
            or canonical(state["cell"]) != canonical(cell.to_dict())
            or state["replications"] != design.replications
            or state["completed_replications"] > design.replications
            or state["status"] not in {"running_incomplete", "complete"}
            or not isinstance(state["partitions"], list) or not isinstance(state["recovered_orphans"], list)):
        raise ValueError("cell manifest identity or prefix is invalid")
    for part in state["partitions"]:
        if not isinstance(part, dict) or set(part) != {"filename", "start", "stop", "rows", "sha256"}:
            raise ValueError("partition manifest has unexpected fields")
        for name in ("start", "stop", "rows"):
            integer(part[name], "partition " + name, 0)
        if not 0 <= part["start"] < part["stop"] <= design.replications:
            raise ValueError("partition range is outside the fixed design")
    return state


def _scan_cells(output, design, identity, *, block_pipeline_errors):
    """Pure prefix/orphan inspection of every durable cell, including root lag."""
    reports, failures = {}, []
    for cell in design.cells:
        directory = output / "cells" / cell.identity
        manifest = directory / "cell_manifest.json"
        if not manifest.exists():
            if directory.exists() and any(directory.iterdir()):
                raise ValueError("cell evidence exists without its manifest")
            reports[cell.identity] = {"cell_id": cell.identity, "paths": 0, "rows": 0,
                                      "pipeline_errors": 0, "complete": False, "orphan_stop": None}
            continue
        _strict_cell_manifest(manifest, cell, design, identity)
        # resume=True is essential: this constructor only reads/verifies and
        # must not create a missing manifest during pure preflight.
        store = CellStore(directory, cell, design, identity, resume=True)
        committed = {part["filename"] for part in store.state["partitions"]}
        actual = {path.name for path in directory.iterdir()} - {"cell_manifest.json"}
        orphans = actual - committed
        if len(orphans) > 1 or (bool(orphans) != (store.orphan_stop is not None)):
            raise ValueError("cell has an unidentified orphan partition")
        errors = 0
        for filename in sorted(actual):
            with gzip.open(directory / filename, "rt", encoding="utf-8") as stream:
                for line in stream:
                    if not line.strip():
                        raise ValueError("blank line in durable partition")
                    row = _loads(line)
                    if not isinstance(row, dict) or type(row.get("pipeline_error")) is not str:
                        raise ValueError("every durable fit must retain its pipeline-error field")
                    if row["pipeline_error"]:
                        errors += 1
                        failures.append((cell.identity, row.get("rep"), row.get("estimator")))
        reports[cell.identity] = {"cell_id": cell.identity, "paths": store.next_replication,
                                  "rows": store.next_replication * len(ESTIMATORS),
                                  "pipeline_errors": errors, "complete": store.next_replication == design.replications,
                                  "orphan_stop": store.orphan_stop}
    if failures and block_pipeline_errors:
        raise SavedPipelineError(f"{len(failures)} durable pipeline-error records require review; no source is retried")
    return reports


def _cell_result(report):
    return {key: report[key] for key in ("cell_id", "paths", "rows", "pipeline_errors")}


@dataclass(frozen=True)
class Preflight:
    output: Path
    metadata: dict | None
    expected_metadata_bytes: bytes | None
    provenance: dict
    cells: dict
    cached_completion: bool
    require_frozen: bool


def _prepare(kind, design, output_dir, *, extra_files, resume, allow_lock_only=False):
    require_design(kind, design)
    if type(resume) is not bool:
        raise ValueError("resume must be Boolean")
    output = _output_path(output_dir)
    _root_layout(output, design)
    require_frozen = design.scope == "exploratory"
    fresh = capture(design, (*FILES, *tuple(extra_files)), require_frozen=require_frozen)
    _check_frozen(fresh, design, fresh, require_frozen)
    root = output / "run_metadata.json"
    if not root.exists():
        names = {p.name for p in output.iterdir()} if output.exists() else set()
        allowed = {".execution.lock"} if allow_lock_only else set()
        if resume or names - allowed:
            raise ValueError("cannot initialize or resume unidentified evidence without a root manifest")
        return Preflight(output, None, None, fresh, {}, False, require_frozen)
    if not resume:
        raise FileExistsError("existing experiment evidence requires explicit resume")
    metadata, expected_bytes = _json_file(root)
    integer(metadata.get("schema_version"), "root schema version", 1)
    expected_project = "ou-irregular-experiments" if kind == "followup" else "ou-week7-paired-floor"
    if (metadata["schema_version"] != 1 or metadata.get("project") != expected_project
            or metadata.get("execution_kind") != kind or metadata.get("execution_contract") != EXECUTION_CONTRACT
            or metadata.get("require_frozen") is not require_frozen
            or canonical(metadata.get("source_policy")) != canonical(source_policy(kind, design))
            or canonical(metadata.get("procedure")) != canonical(_procedure(kind))
            or metadata.get("scientific_overrides") != ([] if kind == "followup" else ["paired_floor_source"])):
        raise ValueError("Week 7 root procedure, source policy or freeze identity changed")
    status = metadata.get("status")
    if status in _CLOSED:
        raise ValueError("audited or completed evidence cannot enter raw execution resume")
    if status not in _OPEN | _CACHED:
        raise ValueError("unrecognized raw execution state")
    for field in ("execution_attempts", "resume_events"):
        if not isinstance(metadata.get(field), list):
            raise ValueError("execution lifecycle history is missing")
    _check_frozen(metadata["provenance"], design, fresh, require_frozen)
    resolved, _ = _json_file(output / "resolved_design.json")
    if canonical(resolved) != canonical(design.to_dict()):
        raise ValueError("resolved design differs from the complete canonical plan")
    cells = _scan_cells(output, design, metadata["provenance"]["identity_sha256"], block_pipeline_errors=True)
    completed = metadata.get("completed_cells")
    if not isinstance(completed, dict) or not set(completed) <= set(cells):
        raise ValueError("root completed-cell map is invalid")
    for cid, value in completed.items():
        if not cells[cid]["complete"] or canonical(value) != canonical(_cell_result(cells[cid])):
            raise ValueError("root completed-cell accounting disagrees with durable evidence")
        for name in ("paths", "rows", "pipeline_errors"):
            integer(value[name], "completed-cell " + name, 0)
    if "pipeline_errors" in metadata:
        integer(metadata["pipeline_errors"], "root pipeline-error count", 0)
        if metadata["pipeline_errors"] != 0:
            raise SavedPipelineError("recorded root pipeline errors require review")
    cached = status in _CACHED
    if cached and (not all(report["complete"] for report in cells.values())
                   or set(completed) != set(cells) or metadata.get("pipeline_errors") != 0):
        raise ValueError("cached/postprocessing resume requires all cells complete with zero pipeline errors")
    # Catch root changes during the potentially long pure cell scan, without
    # creating incident or lock directories during preflight itself.
    if root.read_bytes() != expected_bytes:
        raise ValueError("root metadata changed during read-only preflight")
    return Preflight(output, metadata, expected_bytes, fresh, cells, cached, require_frozen)


def preflight(kind, design, output_dir, *, extra_files=(), resume=False):
    """Read-only preflight of raw execution; no mkdir, lock or lifecycle reset.

    This is not a generic read-only archive auditor: mutation policy deliberately
    rejects closed runs and synchronized output locations even during preflight.
    """
    return _prepare(kind, design, output_dir, extra_files=extra_files, resume=resume)


def run_raw(kind, design, output_dir, *, extra_files=(), workers=1, partition_size=500,
            resume=False, progress=True):
    """Run/resume raw fixed cells, preserving statistical failures and raw evidence."""
    integer(workers, "workers", 1)
    integer(partition_size, "partition_size", 1)
    if type(progress) is not bool:
        raise ValueError("progress must be Boolean")
    extra_files = tuple(extra_files)
    prepared = preflight(kind, design, output_dir, extra_files=extra_files, resume=resume)
    if prepared.cached_completion:
        return prepared.metadata
    output = prepared.output
    output.mkdir(parents=True, exist_ok=True)
    with run_lock(output):
        # Recheck after acquiring the exclusive lock; another process cannot
        # turn the earlier read-only snapshot into permission to overwrite it.
        prepared = _prepare(kind, design, output, extra_files=extra_files, resume=resume, allow_lock_only=True)
        if prepared.cached_completion:
            return prepared.metadata
        path = output / "run_metadata.json"
        guard = _ManifestGuard(path, prepared.expected_metadata_bytes)
        if prepared.metadata is None:
            metadata = {"schema_version": 1,
                        "project": "ou-irregular-experiments" if kind == "followup" else "ou-week7-paired-floor",
                        "execution_kind": kind, "execution_contract": EXECUTION_CONTRACT,
                        "source_policy": source_policy(kind, design), "require_frozen": prepared.require_frozen,
                        "status": "running_incomplete", "started_utc": utc_now(), "provenance": prepared.provenance,
                        "procedure": _procedure(kind),
                        "scientific_overrides": [] if kind == "followup" else ["paired_floor_source"],
                        "resume_events": [], "completed_cells": {}, "execution_attempts": []}
            guard.check("before_resolved_design")
            atomic_json(output / "resolved_design.json", design.to_dict(), immutable=True)
        else:
            metadata = prepared.metadata
            metadata["resume_events"].append({"utc": utc_now(), "workers": workers, "partition_size": partition_size})
        metadata["status"] = "running_incomplete"
        attempt = {"started_utc": utc_now(), "workers": workers, "partition_size": partition_size}
        metadata["execution_attempts"].append(attempt)
        guard.write(metadata, "attempt_start")
        started = time.perf_counter()
        try:
            source = PairedFloorSource(design.seed) if kind == "floor" else None
            args = [(str(output), cell, design, metadata["provenance"]["identity_sha256"],
                     partition_size, source, None, progress) for cell in design.cells]

            def commit(result):
                cid = result.get("cell_id")
                if cid not in {cell.identity for cell in design.cells}:
                    raise ValueError("worker returned an unplanned cell")
                for field in ("paths", "rows", "pipeline_errors"):
                    integer(result.get(field), "worker " + field, 0)
                if set(result) != {"cell_id", "paths", "rows", "pipeline_errors"} or result["paths"] != design.replications or result["rows"] != 3 * design.replications:
                    raise ValueError("worker returned incomplete or unexpected cell accounting")
                metadata["completed_cells"][cid] = result
                guard.write(metadata, "cell_commit:" + cid)

            if workers == 1:
                for arguments in args:
                    guard.check("before_cell_execution:" + arguments[1].identity)
                    commit(_execute_cell(*arguments))
            else:
                # Spawn avoids inheriting the held root flock.  Only the parent
                # writes its root manifest; each worker owns one disjoint cell.
                with ProcessPoolExecutor(max_workers=workers, mp_context=multiprocessing.get_context("spawn")) as pool:
                    guard.check("before_worker_submission")
                    futures = [pool.submit(_execute_cell, *arguments) for arguments in args]
                    try:
                        for future in as_completed(futures):
                            commit(future.result())
                    except BaseException:
                        for future in futures:
                            future.cancel()
                        raise
            guard.check("before_raw_completion_verification")
            _root_layout(output, design)
            _check_frozen(metadata["provenance"], design, prepared.provenance, prepared.require_frozen)
            cells = _scan_cells(output, design, metadata["provenance"]["identity_sha256"], block_pipeline_errors=False)
            if (not all(report["complete"] for report in cells.values())
                    or canonical(metadata["completed_cells"]) != canonical({cid: _cell_result(report) for cid, report in cells.items()})):
                raise ValueError("raw completion disagrees with readback cell accounting")
            errors = sum(report["pipeline_errors"] for report in cells.values())
            metadata.update(status="implementation_review_required" if errors else "raw_complete_pending_audit",
                            pipeline_errors=errors, raw_finished_utc=utc_now())
        except BaseException as exc:
            if not guard.drifted:
                metadata.update(status="incomplete", termination_type=type(exc).__name__, termination_reason=str(exc))
                attempt.update(finished_utc=utc_now(), elapsed_seconds=time.perf_counter() - started)
                guard.write(metadata, "attempt_failure")
            raise
        else:
            attempt.update(finished_utc=utc_now(), elapsed_seconds=time.perf_counter() - started)
            guard.write(metadata, "attempt_finish")
        return metadata
