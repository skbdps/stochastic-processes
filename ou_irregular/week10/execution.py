"""Durable fixed-horizon execution and independent saved-source readback."""
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
import multiprocessing
from pathlib import Path
import time

import numpy as np
import pandas as pd

from ..experiments.design import ESTIMATORS, integer
from ..experiments.execution import (ENGINE_FILES, PROCEDURE, _execute_cell,
                                     registered_wald, source_records)
from ..experiments.provenance import capture, canonical, file_digest, utc_now, verify
from ..experiments.storage import CellStore, atomic_json, run_lock
from ..experiments.summaries import validate_records
from .design import parse_design, source_policy, validate_design
from .source import HorizonSource

FILES = tuple(sorted(set(ENGINE_FILES) | {
    "ou_irregular/week8/__init__.py", "ou_irregular/week8/design.py", "ou_irregular/week8/source.py",
    "ou_irregular/week10/__init__.py", "ou_irregular/week10/design.py",
    "ou_irregular/week10/source.py", "ou_irregular/week10/execution.py",
    "configs/week10_supplement.json", "docs/week10_supplement_protocol.md",
    "tests/week10/test_source.py", "tests/week10/test_execution.py",
}))
CONTRACT = "week10_frozen_horizon_execution_v1"
HORIZON_PROCEDURE = {**PROCEDURE,
    "source_address": "SeedSequence(seed,spawn_key=(10,*eight_cell_id_words,rep)).spawn(2)"}


def freeze(design, *, extra_files=()):
    """Capture explicit pre-execution source/design/runtime hashes.

    This supplement records uncommitted files honestly; a Git commit is not
    substituted for the actual source byte inventory. The returned mapping
    must be saved and passed explicitly before a full study can execute.
    """
    validate_design(design)
    return capture(design, (*FILES, *extra_files), require_frozen=False)


def _path(output):
    path = Path(output).absolute()
    if any(p.is_symlink() for p in (path, *path.parents)):
        raise ValueError("symbolic links are not permitted in experiment output paths")
    for ancestor in path.parents:
        if any((ancestor / name).exists() for name in
               ("run_metadata.json", "package_state.json", "table_binding.json", "validation.json")):
            raise ValueError("new experiment output cannot be nested within saved evidence")
    return path.resolve()


def _stores(output, design, identity):
    directory = output / "cells"
    if directory.exists():
        expected = {c.identity for c in design.cells}
        if any(p.is_symlink() or not p.is_dir() or p.name not in expected for p in directory.iterdir()):
            raise ValueError("unknown or nonregular cell directory")
    stores = {}
    for cell in design.cells:
        path = directory / cell.identity
        if not path.exists():
            continue
        if any(p.is_symlink() or not p.is_file() or
               (p.name != "cell_manifest.json" and not p.name.startswith("part-")) for p in path.iterdir()):
            raise ValueError("unknown or nonregular cell artifact")
        stores[cell.identity] = CellStore(path, cell, design, identity, resume=True)
    return stores


def readback(design, output, *, check_environment=True):
    """Verify saved bytes and complete planned records; do not regenerate data."""
    validate_design(design)
    output = _path(output)
    metadata = json.loads((output / "run_metadata.json").read_text())
    if (metadata.get("contract") != CONTRACT or metadata.get("procedure") != HORIZON_PROCEDURE
            or metadata.get("source_policy") != source_policy(design)):
        raise ValueError("horizon procedure/source identity mismatch")
    resolved = parse_design((output / "resolved_design.json").read_text())
    if resolved.digest != design.digest:
        raise ValueError("resolved design mismatch")
    verify(metadata["provenance"], expected_design=design, check_environment=check_environment)
    saved_freeze = json.loads((output / "execution_freeze.json").read_text())
    if saved_freeze != metadata["provenance"]:
        raise ValueError("saved pre-execution freeze changed")
    stores = _stores(output, design, metadata["provenance"]["identity_sha256"])
    if set(stores) != {c.identity for c in design.cells}:
        raise ValueError("not all planned cells are present")
    rows = []
    for cell in design.cells:
        store = stores[cell.identity]
        if store.next_replication != design.replications or store.orphan_stop is not None:
            raise ValueError("cell is incomplete or has uncommitted evidence")
        rows.extend(store.records())
    data = validate_records(rows, design)
    source_fields = [key for key in data if key.startswith("source_") or key.endswith("_sha256")]
    source_fields += ["actual_n", "actual_transitions", "realized_span", "terminal_gap", "terminal_parent_overshoot"]
    available = data.source_valid.to_numpy(bool)
    if not data.loc[available, "realized_span"].eq(data.loc[available, "horizon"]).all():
        raise ValueError("fixed-horizon source span changed")
    if data.groupby(["cell_id", "rep"])[source_fields].nunique(dropna=False).gt(1).any().any():
        raise ValueError("estimators have mismatched source evidence")
    if not data.loc[available, "actual_n"].eq(data.loc[available, "actual_transitions"] + 1).all():
        raise ValueError("source observation-count mismatch")
    errors = int(data.pipeline_error.fillna("").astype(bool).sum())
    expected = {c.identity: {"cell_id": c.identity, "paths": design.replications,
                            "rows": 3 * design.replications,
                            "pipeline_errors": int(data.loc[data.cell_id == c.identity,
                                                           "pipeline_error"].fillna("").astype(bool).sum())}
                for c in design.cells}
    if metadata.get("completed_cells") != expected or metadata.get("pipeline_errors") != errors:
        raise ValueError("root accounting disagrees with immutable partitions")
    if metadata.get("status") not in {"raw_complete_pending_audit", "implementation_review_required", "completed"}:
        raise ValueError("root does not record raw completion")
    return metadata, data


def load_records(output, design, *, check_environment=True):
    return readback(design, output, check_environment=check_environment)[1].to_dict("records")


def run_raw(design, output, *, workers=1, partition_size=250, resume=False,
            frozen_provenance=None, extra_files=(), progress=True):
    validate_design(design)
    integer(workers, "workers", 1)
    integer(partition_size, "partition_size", 1)
    if type(resume) is not bool or type(progress) is not bool:
        raise ValueError("resume and progress must be Boolean")
    if frozen_provenance is None:
        if design.scope != "validation":
            raise ValueError("full execution requires explicit pre-execution frozen_provenance")
        frozen_provenance = freeze(design, extra_files=extra_files)
    verify(frozen_provenance, expected_design=design)
    if not set(FILES).issubset(frozen_provenance["source_files"]):
        raise ValueError("execution freeze omits required dependencies")
    output = _path(output)
    if output.exists() and any(output.iterdir()) and not resume:
        raise FileExistsError("existing experiment evidence requires explicit resume")
    output.mkdir(parents=True, exist_ok=True)
    with run_lock(output):
        manifest = output / "run_metadata.json"
        if manifest.exists():
            metadata = json.loads(manifest.read_text())
            if (metadata.get("contract") != CONTRACT or metadata.get("procedure") != HORIZON_PROCEDURE
                    or metadata.get("source_policy") != source_policy(design)
                    or metadata.get("provenance") != frozen_provenance):
                raise ValueError("resume source, runtime or frozen identity changed")
            if metadata.get("status") == "completed" or (output / "audit.json").exists():
                raise ValueError("audited or completed evidence cannot enter execution resume")
            if parse_design((output / "resolved_design.json").read_text()).digest != design.digest:
                raise ValueError("resume design changed")
            if json.loads((output / "execution_freeze.json").read_text()) != frozen_provenance:
                raise ValueError("pre-execution freeze changed")
            stores = _stores(output, design, frozen_provenance["identity_sha256"])
            for store in stores.values():
                if any(row.get("pipeline_error") for row in store.records()):
                    raise ValueError("saved pipeline error requires review before resumption")
            if metadata.get("status") == "raw_complete_pending_audit":
                readback(design, output)
                return metadata
            if metadata.get("status") not in {"incomplete", "running_incomplete"}:
                raise ValueError("unrecognized or blocked execution status")
            metadata["resume_events"].append({"utc": utc_now(), "workers": workers,
                                               "partition_size": partition_size})
        else:
            if resume or any(p.name != ".execution.lock" for p in output.iterdir()):
                raise ValueError("cannot start or resume unidentified output evidence")
            atomic_json(output / "execution_freeze.json", frozen_provenance, immutable=True)
            atomic_json(output / "resolved_design.json", design.to_dict(), immutable=True)
            metadata = {"schema_version": 1, "contract": CONTRACT, "source_policy": source_policy(design),
                        "procedure": HORIZON_PROCEDURE, "provenance": frozen_provenance,
                        "started_utc": utc_now(), "completed_cells": {}, "resume_events": [],
                        "execution_attempts": []}
        metadata["status"] = "running_incomplete"
        attempt = {"started_utc": utc_now(), "workers": workers, "partition_size": partition_size}
        metadata["execution_attempts"].append(attempt)
        atomic_json(manifest, metadata)
        started = time.perf_counter()
        try:
            source = HorizonSource(design.seed)
            args = [(str(output), cell, design, frozen_provenance["identity_sha256"],
                     partition_size, source, None, progress) for cell in design.cells]

            def commit(result):
                metadata["completed_cells"][result["cell_id"]] = result
                atomic_json(manifest, metadata)

            if workers == 1:
                for argument in args:
                    commit(_execute_cell(*argument))
            else:
                with ProcessPoolExecutor(max_workers=workers,
                                         mp_context=multiprocessing.get_context("spawn")) as pool:
                    futures = [pool.submit(_execute_cell, *argument) for argument in args]
                    for future in as_completed(futures):
                        commit(future.result())
            verify(frozen_provenance, expected_design=design)
            errors = sum(row["pipeline_errors"] for row in metadata["completed_cells"].values())
            metadata.update(status="implementation_review_required" if errors else "raw_complete_pending_audit",
                            pipeline_errors=errors, raw_finished_utc=utc_now())
        except BaseException as exc:
            metadata.update(status="incomplete", termination_type=type(exc).__name__, termination_reason=str(exc))
            raise
        finally:
            attempt.update(finished_utc=utc_now(), elapsed_seconds=time.perf_counter() - started)
            atomic_json(manifest, metadata)
        readback(design, output)
    return metadata


def audit(design, output, *, progress=True):
    """Reconstruct every schedule/path and replay fixed endpoint fits.

    PFML mean is independently checked by unscaled OLS on every valid source.
    This audit does not declare figure review or the complete supplement done.
    """
    metadata, data = readback(design, output)
    if metadata["pipeline_errors"]:
        raise ValueError("pipeline errors block a successful scientific audit")
    source = HorizonSource(design.seed)
    services = {name: registered_wald(name) for name in ESTIMATORS}
    indexed = data.set_index(["cell_id", "rep", "estimator"])
    sources = identity_checks = endpoint_fits = 0
    max_mu_difference = 0.
    for cell in design.cells:
        for rep in range(design.replications):
            sample, diagnostics = source.sample(cell.to_dict(), rep)
            for name in ESTIMATORS:
                saved = indexed.loc[cell.identity, rep, name]
                if any(canonical(saved[key]) != canonical(value) for key, value in diagnostics.items()):
                    raise ValueError(f"source replay mismatch at {cell.identity}/{rep}/{name}")
            saved = indexed.loc[cell.identity, rep, "pfml"]
            if saved.point_valid:
                x, y = sample.values[:-1], sample.values[1:]
                xc, yc = x - x.mean(), y - y.mean()
                rho = float(np.dot(xc, yc) / np.dot(xc, xc))
                if not 0 < rho < 1:
                    raise ValueError("valid PFML record has no interior independent OLS solution")
                mu = float((y.mean() - rho * x.mean()) / (1 - rho))
                difference = abs(mu - float(saved.mu))
                if difference > 1e-10 * max(1., abs(mu), float(np.std(sample.values))):
                    raise ValueError("PFML mean differs from independent OLS identity")
                max_mu_difference = max(max_mu_difference, difference)
                identity_checks += 1
            if rep in {0, design.replications - 1}:
                replays = source_records(cell, rep, source, services)
                for replay in replays:
                    original = indexed.loc[cell.identity, rep, replay["estimator"]].to_dict()
                    original.update(cell_id=cell.identity, rep=rep, estimator=replay["estimator"])
                    if canonical(original) != canonical(replay):
                        differing = [k for k in replay if canonical(original.get(k)) != canonical(replay[k])]
                        raise ValueError(f"fit replay mismatch: {cell.identity}/{rep}/{replay['estimator']}: {differing}")
                    endpoint_fits += 1
            sources += 1
        if progress:
            print(f"horizon audit: {sources}/{len(design.cells) * design.replications} sources", flush=True)
    partitions = {str(path.relative_to(output)): file_digest(path)
                  for path in sorted(Path(output).glob("cells/*/*.jsonl.gz"))}
    result = {"schema_version": 1, "status": "passed", "audited_utc": utc_now(),
              "design_sha256": design.digest, "execution_identity_sha256": metadata["provenance"]["identity_sha256"],
              "sources_reconstructed": sources, "fit_records": len(data),
              "endpoint_fit_replays": endpoint_fits, "independent_pfml_mean_checks": identity_checks,
              "maximum_absolute_pfml_mean_difference": max_mu_difference,
              "pipeline_errors": metadata["pipeline_errors"],
              "point_valid": int(data.point_valid.sum()), "interval_valid": int(data.interval_valid.sum()),
              "all_realized_spans_equal_horizon": True, "partitions": partitions,
              "scope": "source hashes, complete planned keys, independent PFML mean and endpoint fits; no visual review"}
    path = Path(output) / "audit.json"
    if path.exists():
        previous = json.loads(path.read_text())
        if {k: v for k, v in previous.items() if k != "audited_utc"} != {k: v for k, v in result.items() if k != "audited_utc"}:
            raise ValueError("existing audit disagrees with fresh scientific readback")
        return previous
    atomic_json(path, result, immutable=True)
    return result
