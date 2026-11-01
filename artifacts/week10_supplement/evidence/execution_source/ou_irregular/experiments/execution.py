"""Execute fixed source addresses through the unchanged Week 3/4 services.

The unit of scheduling is a cell. Workers write disjoint, durable partitions;
worker count and partition size do not enter random streams. This module only
finishes raw execution. Scientific completion additionally requires derived
tables, figures, independent arithmetic/source replay and documented review.
"""
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
import multiprocessing
from pathlib import Path
import time

import numpy as np

from ..inference.adapters import OuLikelihood, Week3Estimator
from ..inference.contracts import PointFit
from ..inference.hessian import FiniteDifferenceHessian
from ..inference.wald import WaldInference, WaldResult
from ..week4.simulation import OuSampleSource
from ..week4.summaries import PARAMETERS
from .design import ESTIMATORS, integer
from .provenance import capture, canonical, utc_now, verify
from .storage import CellStore, atomic_json, run_lock


ENGINE_FILES = ("ou_irregular/experiments/storage.py", "ou_irregular/experiments/execution.py")
PROCEDURE = {"level": .95, "relative_step": .001, "max_relative_difference": .005,
             "max_whitened_difference": .01, "max_condition_number": 1e10,
             "max_standardized_score": .01, "interval_scale": "natural_symmetric_unclipped",
             "covariance": "full_conditional_nll_inverse_hessian_model_based",
             "source_address": "SeedSequence(root,spawn_key=(4,*8_cell_words,rep)).spawn(2)"}


class HessianTrace:
    """Record the exact returned derivative object without modifying it."""

    def __init__(self):
        self.provider = FiniteDifferenceHessian(relative_step=.001, max_relative_difference=.005,
                                               max_whitened_difference=.01)
        self.last = None

    def evaluate(self, objective, point):
        self.last = self.provider.evaluate(objective, point)
        return self.last


def registered_wald(name):
    return WaldInference(Week3Estimator(name), OuLikelihood(name), hessian=HessianTrace(),
                         level=.95, max_condition_number=1e10, max_standardized_score=.01)


def failed_result(name, status, reason, point=None):
    point = point if point is not None else PointFit(
        name, np.nan, np.nan, np.nan, np.nan, False, False, status, reason)
    return WaldResult(point, False, status, reason, np.full(3, np.nan), np.full(3, np.nan),
                      np.full((3, 3), np.nan), np.full((3, 3), np.nan))


def inference_record(result, trace=None):
    point, diagnostic = result.point_fit, result.diagnostics
    hd = diagnostic.get("hessian_diagnostics", {})
    row = {"estimator": point.estimator, "point_valid": bool(point.valid),
           "optimizer_success": bool(point.optimizer_success), "point_status": point.status,
           "point_reason": point.reason, "nll": point.nll,
           "point_diagnostics_json": canonical(point.diagnostics),
           "interval_valid": bool(result.valid), "interval_status": result.status,
           "interval_reason": result.reason, "interval_diagnostics_json": canonical(diagnostic),
           "covariance_type": "inverse_observed_hessian_model_based",
           "hessian_condition_number": diagnostic.get("hessian_condition_number", np.nan),
           "hessian_relative_difference": diagnostic.get("hessian_relative_difference", np.nan),
           "hessian_whitened_difference": hd.get("whitened_relative_difference", np.nan),
           "standardized_score": diagnostic.get("standardized_score", np.nan)}
    for i, p in enumerate(PARAMETERS):
        row.update({p: point.parameters[i], f"{p}_lower": result.lower[i],
                    f"{p}_upper": result.upper[i], f"{p}_se": np.sqrt(result.covariance[i, i]),
                    f"gradient_{i}": trace.gradient[i] if trace is not None else np.nan})
        for j in range(i, 3):
            row[f"cov_{p}_{PARAMETERS[j]}"] = result.covariance[i, j]
            row[f"working_cov_{i}_{j}"] = result.working_covariance[i, j]
            row[f"working_hessian_{i}_{j}"] = trace.hessian[i, j] if trace is not None else np.nan
    return row


def source_records(cell, rep, source, services):
    """Every planned estimator survives source, fit or inference exceptions."""
    try:
        sample, diagnostics = source.sample(cell.to_dict(), rep)
        source_error = ""
    except Exception as exc:
        sample, diagnostics = None, {}
        source_error = f"{type(exc).__name__}: {exc}"
    rows = []
    for name in ESTIMATORS:
        service, point, error = services[name], None, source_error
        tracer = getattr(service, "hessian", None)
        if isinstance(tracer, HessianTrace):
            tracer.last = None
        if sample is None:
            result = failed_result(name, "simulation_failure", source_error)
        else:
            try:
                point = service.estimator.fit(sample)
                result = service.infer(sample, point_fit=point)
            except Exception as exc:
                error = f"{type(exc).__name__}: {exc}"
                result = failed_result(name, "pipeline_exception", error, point)
        trace = tracer.last if isinstance(tracer, HessianTrace) else None
        rows.append({**cell.to_dict(), "rep": rep, **diagnostics,
                     "source_valid": sample is not None, "source_error": source_error,
                     **inference_record(result, trace), "pipeline_error": error})
    return rows


def _execute_cell(output, cell, design, identity, partition_size, source, wald_factory, progress):
    directory = Path(output) / "cells" / cell.identity
    store = CellStore(directory, cell, design, identity,
                      resume=(directory / "cell_manifest.json").exists())
    source = source if source is not None else OuSampleSource(design.seed)
    services = {name: (wald_factory or registered_wald)(name) for name in ESTIMATORS}
    while store.next_replication < design.replications:
        start = store.next_replication
        stop = store.orphan_stop or min(start + partition_size, design.replications)
        rows = [row for rep in range(start, stop) for row in source_records(cell, rep, source, services)]
        store.append(rows, stop)
        if progress:
            print(f"{design.name} n={cell.n} q={cell.q:g} CV={cell.cv:g}: {stop}/{design.replications} paths saved", flush=True)
    errors = sum(bool(row["pipeline_error"]) for row in store.records())
    return {"cell_id": cell.identity, "paths": store.next_replication,
            "rows": store.next_replication * len(ESTIMATORS), "pipeline_errors": errors}


def run_raw(design, output_dir, *, extra_files=(), workers=1, partition_size=250,
            resume=False, source=None, wald_factory=None, headline_design_digest=None, progress=True):
    integer(workers, "workers", 1)
    integer(partition_size, "partition_size", 1)
    overrides = [name for name, value in (("source", source), ("wald_factory", wald_factory)) if value is not None]
    if design.scope == "headline" and (overrides or headline_design_digest != design.digest):
        raise ValueError("headline execution requires a frozen entry-point design and registered services")
    if overrides and workers != 1:
        raise ValueError("injected diagnostic services require a single worker")
    output = Path(output_dir)
    if output.exists() and any(output.iterdir()) and not resume:
        raise FileExistsError("existing experiment evidence requires explicit --resume")
    output.mkdir(parents=True, exist_ok=True)
    with run_lock(output):
        path = output / "run_metadata.json"
        if path.exists():
            metadata = json.loads(path.read_text())
            verify(metadata["provenance"], expected_design=design)
            if metadata["scientific_overrides"] != overrides or metadata["procedure"] != PROCEDURE:
                raise ValueError("resume numerical-procedure identity mismatch")
            if overrides:
                raise ValueError("injected diagnostic services cannot claim verified automatic resume")
            metadata["resume_events"].append({"utc": utc_now(), "workers": workers,
                                                "partition_size": partition_size})
        else:
            if resume and any(p.name != ".execution.lock" for p in output.iterdir()):
                raise ValueError("cannot resume evidence without a root manifest")
            provenance = capture(design, (*ENGINE_FILES, *extra_files), require_frozen=design.scope == "headline")
            metadata = {"schema_version": 1, "project": "ou-irregular-experiments",
                        "status": "running_incomplete", "started_utc": utc_now(),
                        "provenance": provenance, "procedure": PROCEDURE,
                        "scientific_overrides": overrides, "resume_events": [], "completed_cells": {},
                        "execution_attempts": []}
            atomic_json(output / "resolved_design.json", design.to_dict(), immutable=True)
        metadata["status"] = "running_incomplete"
        attempt = {"started_utc": utc_now(), "workers": workers, "partition_size": partition_size}
        metadata["execution_attempts"].append(attempt)
        atomic_json(path, metadata)
        started = time.perf_counter()
        try:
            args = [(str(output), cell, design, metadata["provenance"]["identity_sha256"],
                     partition_size, source, wald_factory, progress) for cell in design.cells]
            if workers == 1:
                for arguments in args:
                    result = _execute_cell(*arguments)
                    metadata["completed_cells"][result["cell_id"]] = result
                    atomic_json(path, metadata)
            else:
                # Spawn avoids inheriting a held run lock or BLAS thread state.
                with ProcessPoolExecutor(max_workers=workers, mp_context=multiprocessing.get_context("spawn")) as pool:
                    futures = [pool.submit(_execute_cell, *arguments) for arguments in args]
                    for future in as_completed(futures):
                        result = future.result()
                        metadata["completed_cells"][result["cell_id"]] = result
                        atomic_json(path, metadata)
            verify(metadata["provenance"], expected_design=design)
            total_errors = sum(c["pipeline_errors"] for c in metadata["completed_cells"].values())
            metadata.update(status="implementation_review_required" if total_errors else "raw_complete_pending_audit",
                            pipeline_errors=total_errors, raw_finished_utc=utc_now())
        except BaseException as exc:
            metadata.update(status="incomplete", termination_type=type(exc).__name__, termination_reason=str(exc))
            raise
        finally:
            attempt.update(finished_utc=utc_now(), elapsed_seconds=time.perf_counter()-started)
            atomic_json(path, metadata)
    return metadata


def load_records(output_dir, design):
    output = Path(output_dir)
    metadata = json.loads((output / "run_metadata.json").read_text())
    verify(metadata["provenance"], expected_design=design)
    for cell in design.cells:
        store = CellStore(output / "cells" / cell.identity, cell, design,
                          metadata["provenance"]["identity_sha256"], resume=True)
        if store.next_replication != design.replications:
            raise ValueError("cannot summarize an incomplete experiment")
        yield from store.records()
