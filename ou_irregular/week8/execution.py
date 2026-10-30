"""Guarded execution of the two fixed E4 fitted-run plans.

``preflight(design, output, extra_files=(), resume=False)`` performs a pure
raw-execution preflight. ``run_raw`` executes the unchanged durable cell engine
with E4Source and the registered inference services, stopping at raw completion.
``readback(design, output, extra_files=())`` checks saved identity, provenance,
layout and cell accounting without mutation restrictions or scientific replay;
it permits closed states but is not a source audit or completion certification.

The imported Week 7 JSON, cell-scan, manifest-guard and freeze helpers are source
neutral. Their complete declared dependencies are pinned below. E4 layouts,
source policy and root authority are defined here, never borrowed from a Gamma
or paired-floor scientific design. No source/inference overrides are exposed.
"""
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
import multiprocessing
from pathlib import Path
import time

from ..experiments.design import ESTIMATORS, integer
from ..experiments.execution import PROCEDURE as REGISTERED_PROCEDURE, _execute_cell
from ..experiments.provenance import CORE_FILES, canonical, capture, file_digest, utc_now
from ..experiments.storage import atomic_json, run_lock
from ..week7.execution import (FILES as SHARED_FILES, SavedPipelineError, _PARTITION,
                               _cell_result, _json_file, _loads, _scan_cells)
from ..week7.path_diagnostics import (ManifestDriftError, _ManifestGuard, _check_frozen,
                                      _check_incident_block, _reject_archive_ancestors)
from .design import SOURCE_ADDRESS, parse_design, source_policy, validate_design
from .source import E4Source, FILES as SOURCE_FILES


FILES = tuple(sorted(set(SHARED_FILES) | set(SOURCE_FILES) | {
    "ou_irregular/week5/__init__.py", "ou_irregular/week6/__init__.py",
    "ou_irregular/week8/execution.py", "docs/week8_execution.md",
}))
EXECUTION_CONTRACT = "week8_guarded_cell_execution_v1"
PROCEDURE = {**REGISTERED_PROCEDURE, "source_address": SOURCE_ADDRESS}
SYNC_SCRIPT = Path("/usr/local/scripts/sync_share.sh")
SYNC_ROOT = Path("/workspace/scratch")
TABLE_FILES = frozenset({
    "parameter_summary.csv", "paired_comparisons.csv", "regular_comparisons.csv",
    "matched_iid_comparisons.csv", "spacing_summary.csv", "source_summary.csv",
    "outcome_counts.csv", "numerical_diagnostics.csv", "failed_fits.csv.gz",
})
_RAW_FILES = {".execution.lock", "run_metadata.json", "resolved_design.json"}
_DERIVED_FILES = TABLE_FILES | {"table_manifest.json", "table_binding.json", "audit.json", "visual_review.json"}
_OPEN = {"running_incomplete", "incomplete"}
_CACHED = {"raw_complete_pending_audit", "implementation_review_required"}
_CLOSED = {"audited_pending_visual_review", "completed"}


def _ordinary_path(output):
    """Resolve an ordinary path for readback, without authorizing mutation."""
    supplied = Path(output).absolute()
    if any(path.is_symlink() for path in (supplied, *supplied.parents)):
        raise ValueError("symbolic links are not permitted in E4 evidence paths")
    return supplied.resolve()


def _output_path(output):
    """Mutation policy; intentionally separate from read-only archive identity."""
    output = _ordinary_path(output)
    if SYNC_SCRIPT.exists() and output.is_relative_to(SYNC_ROOT.resolve()):
        raise ValueError("mutable runs in synchronized /workspace/scratch are blocked; use isolated /tmp output")
    _reject_archive_ancestors(output)
    for ancestor in output.parents:
        if any((ancestor / name).is_file() for name in ("validation_plan.json", "validation.json", "table_binding.json")):
            raise ValueError("execution output is nested in saved scientific evidence")
    _check_incident_block(output)
    return output


def _root_layout(output, design):
    """Validate known E4 artifacts, retaining partial derived publication."""
    output = Path(output)
    if not output.exists():
        return
    if not output.is_dir():
        raise ValueError("execution output must be a directory")
    for path in output.iterdir():
        if path.name not in _RAW_FILES | _DERIVED_FILES | {"cells", "figures"} or path.is_symlink():
            raise ValueError("unknown or symbolic-link root artifact")
        if path.name in {"cells", "figures"}:
            if not path.is_dir():
                raise ValueError("artifact directory has an invalid object type")
        elif not path.is_file():
            raise ValueError("root artifact is not a regular file")
    cells = output / "cells"
    if cells.exists():
        planned = {cell.identity for cell in design.cells}
        for path in cells.iterdir():
            if path.name not in planned or path.is_symlink() or not path.is_dir():
                raise ValueError("unknown cell directory or symbolic link")
            for child in path.iterdir():
                if (child.is_symlink() or not child.is_file()
                        or (child.name != "cell_manifest.json" and not _PARTITION.fullmatch(child.name))):
                    raise ValueError("unknown cell artifact or symbolic link")
    for name in ("table_binding.json", "audit.json", "visual_review.json"):
        if (output / name).exists():
            _json_file(output / name)
    table_manifest = output / "table_manifest.json"
    if table_manifest.exists():
        manifest, _ = _json_file(table_manifest)
        for name, entry in manifest.items():
            if name not in TABLE_FILES or not isinstance(entry, dict) or set(entry) != {"rows", "sha256"}:
                raise ValueError("unknown or malformed derived table manifest entry")
            integer(entry["rows"], "table row count", 0)
            if file_digest(output / name) != entry["sha256"]:
                raise ValueError("derived table digest changed")
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
                if Path(name).name != name or name == "figure_manifest.json" or file_digest(figures / name) != digest:
                    raise ValueError("figure identity or digest changed")
        elif any(path.suffix not in {".png", ".pdf"} for path in figures.iterdir()):
            raise ValueError("unknown partial figure artifact without a manifest")


def _dependencies(extra_files):
    if isinstance(extra_files, (str, bytes)):
        raise ValueError("extra_files must be an iterable of repository-relative file names")
    extra_files = tuple(extra_files)
    for name in extra_files:
        if (type(name) is not str or not name or Path(name).is_absolute()
                or ".." in Path(name).parts or str(Path(name)) != name):
            raise ValueError("dependency paths must be canonical repository-relative file names")
    return tuple(sorted(set(FILES) | set(extra_files)))


def _fresh_provenance(design, extra_files):
    dependencies = _dependencies(extra_files)
    frozen = design.scope == "exploratory"
    fresh = capture(design, dependencies, require_frozen=frozen)
    if set(fresh["source_files"]) != set(CORE_FILES) | set(dependencies):
        raise ValueError("E4 dependency manifest membership differs from the exact requested union")
    _check_frozen(fresh, design, fresh, frozen)
    return fresh


@dataclass(frozen=True)
class Preflight:
    output: Path
    metadata: dict | None
    expected_metadata_bytes: bytes | None
    provenance: dict
    cells: dict
    cached_completion: bool
    require_frozen: bool


def _prepare(design, output, *, extra_files, resume, mutable, allow_lock_only=False):
    validate_design(design)
    if type(resume) is not bool:
        raise ValueError("resume must be Boolean")
    output = _output_path(output) if mutable else _ordinary_path(output)
    _root_layout(output, design)
    fresh = _fresh_provenance(design, extra_files)
    frozen = design.scope == "exploratory"
    root = output / "run_metadata.json"
    if not root.exists():
        if not mutable:
            raise FileNotFoundError("E4 readback requires an existing root manifest")
        names = {path.name for path in output.iterdir()} if output.exists() else set()
        allowed = {".execution.lock"} if allow_lock_only else set()
        if resume or names - allowed:
            raise ValueError("cannot initialize or resume unidentified evidence without a root manifest")
        return Preflight(output, None, None, fresh, {}, False, frozen)
    if not resume:
        raise FileExistsError("existing experiment evidence requires explicit resume")
    metadata, payload = _json_file(root)
    expected = {"schema_version": 1, "project": "ou-week8-e4", "execution_kind": "e4",
                "execution_contract": EXECUTION_CONTRACT, "require_frozen": frozen,
                "source_policy": source_policy(design), "procedure": PROCEDURE, "scientific_overrides": []}
    if any(key not in metadata or type(metadata[key]) is not type(value)
           or canonical(metadata[key]) != canonical(value) for key, value in expected.items()):
        raise ValueError("E4 root procedure, source policy or freeze identity changed")
    status = metadata.get("status")
    if type(status) is not str or status not in _OPEN | _CACHED | _CLOSED:
        raise ValueError("unrecognized E4 execution state")
    if mutable and status in _CLOSED:
        raise ValueError("audited or completed evidence cannot enter raw execution resume")
    for field in ("execution_attempts", "resume_events"):
        if not isinstance(metadata.get(field), list):
            raise ValueError("execution lifecycle history is missing")
    if not isinstance(metadata.get("provenance"), dict):
        raise ValueError("saved E4 provenance is missing")
    _check_frozen(metadata["provenance"], design, fresh, frozen)
    resolved, _ = _json_file(output / "resolved_design.json")
    if canonical(parse_design(resolved).to_dict()) != canonical(design.to_dict()):
        raise ValueError("resolved design differs from the complete canonical E4 plan")
    cells = _scan_cells(output, design, metadata["provenance"]["identity_sha256"], block_pipeline_errors=mutable)
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
        if mutable and metadata["pipeline_errors"]:
            raise SavedPipelineError("recorded root pipeline errors require review")
    cached = status in _CACHED | _CLOSED
    if cached:
        errors = sum(report["pipeline_errors"] for report in cells.values())
        if (not all(report["complete"] for report in cells.values()) or set(completed) != set(cells)
                or metadata.get("pipeline_errors") != errors):
            raise ValueError("saved raw-complete/closed state disagrees with durable cell accounting")
        if status != "implementation_review_required" and errors:
            raise SavedPipelineError("raw-complete or closed state cannot certify saved pipeline errors")
    if root.read_bytes() != payload:
        raise ValueError("root metadata changed during read-only preflight")
    return Preflight(output, metadata, payload, fresh, cells, cached, frozen)


def readback(design, output, *, extra_files=()):
    """Pure saved identity/cell inspection, including closed or review-required runs.

    It does not generate sources, replay fits, acquire locks or certify derived
    scientific results. Mutation restrictions remain in preflight/run_raw.
    """
    return _prepare(design, output, extra_files=extra_files, resume=True, mutable=False)


def preflight(design, output, *, extra_files=(), resume=False):
    """Read-only raw-execution permission check; no files or locks are created."""
    return _prepare(design, output, extra_files=extra_files, resume=resume, mutable=True)


def run_raw(design, output, *, extra_files=(), workers=1, partition_size=250,
            resume=False, progress=True):
    """Run fixed E4 sources through the unchanged registered durable cell engine."""
    integer(workers, "workers", 1)
    integer(partition_size, "partition_size", 1)
    if type(progress) is not bool:
        raise ValueError("progress must be Boolean")
    extra_files = _dependencies(extra_files)
    prepared = preflight(design, output, extra_files=extra_files, resume=resume)
    if prepared.cached_completion:
        return prepared.metadata
    output = prepared.output
    output.mkdir(parents=True, exist_ok=True)
    with run_lock(output):
        prepared = _prepare(design, output, extra_files=extra_files, resume=resume,
                            mutable=True, allow_lock_only=True)
        if prepared.cached_completion:
            return prepared.metadata
        guard = _ManifestGuard(output / "run_metadata.json", prepared.expected_metadata_bytes)
        if prepared.metadata is None:
            metadata = {"schema_version": 1, "project": "ou-week8-e4", "execution_kind": "e4",
                        "execution_contract": EXECUTION_CONTRACT, "source_policy": source_policy(design),
                        "require_frozen": prepared.require_frozen, "status": "running_incomplete",
                        "started_utc": utc_now(), "provenance": prepared.provenance, "procedure": PROCEDURE,
                        "scientific_overrides": [], "resume_events": [], "completed_cells": {},
                        "execution_attempts": []}
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
            source = E4Source(design.seed)
            args = [(str(output), cell, design, metadata["provenance"]["identity_sha256"],
                     partition_size, source, None, progress) for cell in design.cells]

            def commit(result, cell):
                if not isinstance(result, dict) or set(result) != {"cell_id", "paths", "rows", "pipeline_errors"}:
                    raise ValueError("worker returned unexpected cell accounting fields")
                if result["cell_id"] != cell.identity:
                    raise ValueError("worker returned the wrong planned cell identity")
                for field in ("paths", "rows", "pipeline_errors"):
                    integer(result[field], "worker " + field, 0)
                if result["paths"] != design.replications or result["rows"] != len(ESTIMATORS) * design.replications:
                    raise ValueError("worker returned incomplete cell accounting")
                metadata["completed_cells"][cell.identity] = result
                guard.write(metadata, "cell_commit:" + cell.identity)

            if workers == 1:
                for arguments in args:
                    guard.check("before_cell_execution:" + arguments[1].identity)
                    commit(_execute_cell(*arguments), arguments[1])
            else:
                with ProcessPoolExecutor(max_workers=workers, mp_context=multiprocessing.get_context("spawn")) as pool:
                    guard.check("before_worker_submission")
                    futures = {pool.submit(_execute_cell, *arguments): arguments[1] for arguments in args}
                    try:
                        for future in as_completed(futures):
                            commit(future.result(), futures[future])
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
