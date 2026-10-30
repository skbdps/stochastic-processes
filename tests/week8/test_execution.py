"""E4 execution contracts using literal records and real durable stores.

No schedules, OU observations, or fitted results are generated. The simulation/
inference boundary is replaced by deterministic records. Provenance capture and
Git objects are explicitly mocked, while real dependency, identity, cell-store,
root-guard, and filesystem checks remain active. Process-pool tests use Futures
without starting processes.
"""
from concurrent.futures import Future
from copy import deepcopy
from dataclasses import replace
import gzip
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from ou_irregular.experiments import execution as engine
from ou_irregular.experiments import provenance as provenance_module
from ou_irregular.experiments.design import ESTIMATORS
from ou_irregular.experiments.storage import read_partition, run_lock
from ou_irregular.week7 import path_diagnostics as guard_module
from ou_irregular.week8 import execution
from ou_irregular.week8 import source as source_module
from ou_irregular.week8.design import SOURCE_ADDRESS, full_design, validation_design


ROOT = Path(__file__).resolve().parents[2]
MOCK_COMMIT = "8" * 40


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def write_json(path, value):
    Path(path).write_text(canonical(value) + "\n")


def state(output):
    return json.loads((Path(output) / "run_metadata.json").read_bytes())


def snapshot(root):
    """Include directories and links so a nominally pure call cannot create them."""
    result = {}
    for path in Path(root).rglob("*"):
        name = str(path.relative_to(root))
        if path.is_symlink():
            result[name] = ("symlink", str(path.readlink()))
        elif path.is_dir():
            result[name] = ("directory", None)
        elif path.is_file():
            result[name] = ("file", path.read_bytes())
    return result


def resign(provenance):
    identity = {name: provenance[name] for name in (
        "design", "design_sha256", "preregistration_commit", "source_files", "environment")}
    provenance["identity_sha256"] = hashlib.sha256(canonical(identity).encode()).hexdigest()


def literal_records(cell, rep, source, services):
    """Literal eligible candidates; these hashes label fixtures, not observations."""
    assert isinstance(source, source_module.E4Source), "old Gamma source must not be used"
    assert tuple(services) == ("exact", "pfml", "euler")
    words = [int(cell.identity[i:i + 8], 16) for i in range(0, 64, 8)]
    address = [8, *words, rep]
    common = {**cell.to_dict(), "rep": rep, "source_seed": source.seed,
              "source_spawn_key": canonical(address), "gap_spawn_key": canonical(address + [0]),
              "path_spawn_key": canonical(address + [1]), "source_valid": True,
              "source_error": "", "pipeline_error": "",
              "times_sha256": hashlib.sha256(f"literal-times:{cell.identity}:{rep}".encode()).hexdigest(),
              "observations_sha256": hashlib.sha256(f"literal-values:{cell.identity}:{rep}".encode()).hexdigest()}
    rows = []
    for method in ESTIMATORS:
        row = {**common, "estimator": method, "point_valid": True, "point_status": "valid_interior",
               "point_reason": "literal fixture", "optimizer_success": True,
               "interval_valid": True, "interval_status": "valid", "interval_reason": "literal fixture",
               "nll": 1., "point_diagnostics_json": '{"fixture":"literal"}',
               "interval_diagnostics_json": '{"fixture":"literal"}'}
        for parameter in ("theta", "mu", "sigma"):
            truth = row["true_" + parameter]
            row.update({parameter: truth + .01 * rep, parameter + "_se": .1,
                        parameter + "_lower": truth - .25, parameter + "_upper": truth + .25})
        rows.append(row)
    return rows


def test_declared_dependency_closure_includes_transitively_executed_package_initializers():
    # The reused Week7 helpers import Week5/6 configuration modules, which
    # execute both package initializers even though they are presently small.
    assert {"ou_irregular/week5/__init__.py", "ou_irregular/week6/__init__.py"} <= set(execution.FILES)


@pytest.fixture
def harness(tmp_path, monkeypatch):
    calls, allocations, captures = [], [], []
    environment = provenance_module.environment()

    def forbidden(*args, **kwargs):
        raise AssertionError("execution unit tests must not draw sources or fit models")

    monkeypatch.setattr(source_module, "generate_gaps", forbidden)
    monkeypatch.setattr(source_module, "simulate_ou", forbidden)
    monkeypatch.setattr(source_module.E4Source, "sample", forbidden)
    monkeypatch.setattr(source_module.E4Source, "gaps", forbidden)
    monkeypatch.setattr(engine.OuSampleSource, "sample", forbidden)
    monkeypatch.setattr(np.random, "default_rng", forbidden)
    monkeypatch.setattr(engine, "registered_wald", lambda name: SimpleNamespace(
        estimator=SimpleNamespace(fit=forbidden), infer=forbidden))

    def source_records(cell, rep, source, services):
        calls.append((cell.identity, rep, source.seed))
        return literal_records(cell, rep, source, services)

    def allocate(seed):
        allocations.append(seed)
        return source_module.E4Source(seed)

    def capture(design, extra_files=(), *, require_frozen=False):
        # Mock the recorded commit, not the runtime identity checks. Real
        # current dependency hashes/environment are verified by shared helpers.
        captures.append((design.scope, require_frozen, tuple(extra_files)))
        paths = set(provenance_module.CORE_FILES) | set(extra_files)
        record = {"design": design.to_dict(), "design_sha256": design.digest,
                  "preregistration_commit": provenance_module.PREREG_COMMIT,
                  "source_files": {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
                                   for name in sorted(paths)},
                  "environment": deepcopy(environment), "code_commit": MOCK_COMMIT,
                  "uncommitted_dependency_files": [], "captured_utc": "MOCK_ONLY"}
        resign(record)
        return record

    original_git = guard_module.subprocess.run

    def mock_git(arguments, *args, **kwargs):
        if (isinstance(arguments, (list, tuple)) and arguments[:2] == ["git", "show"]
                and arguments[2].startswith(MOCK_COMMIT + ":")):
            relative = arguments[2].split(":", 1)[1]
            path = ROOT / relative
            return SimpleNamespace(returncode=0 if path.is_file() else 1,
                                   stdout=path.read_bytes() if path.is_file() else b"")
        return original_git(arguments, *args, **kwargs)

    monkeypatch.setattr(guard_module.subprocess, "run", mock_git)
    monkeypatch.setattr(execution, "capture", capture)
    monkeypatch.setattr(execution, "E4Source", allocate)
    monkeypatch.setattr(engine, "source_records", source_records)

    def run(*, design=None, output=None, **kwargs):
        return execution.run_raw(design or validation_design(), output or tmp_path / "run",
                                 progress=False, **kwargs)

    return SimpleNamespace(root=tmp_path, output=tmp_path / "run", calls=calls,
                           allocations=allocations, captures=captures, capture=capture,
                           source_records=source_records, run=run, monkeypatch=monkeypatch)


@pytest.mark.parametrize("factory,frozen,paths,fits", [
    (full_design, True, 14000, 42000), (validation_design, False, 28, 84),
])
def test_exact_fixed_root_and_source_policy_with_pure_preflight(harness, factory, frozen, paths, fits):
    design = factory()
    before = snapshot(harness.root)
    checked = execution.preflight(design, harness.root / "missing" / "new")
    assert checked.metadata is None and not checked.cached_completion
    assert checked.require_frozen is frozen and checked.cells == {}
    assert snapshot(harness.root) == before and not harness.allocations and not harness.calls
    policy = execution.source_policy(design)
    assert policy["root_seed"] == (20261011 if frozen else 21261011)
    assert policy["source_namespace"] == 8 and policy["source_address"] == SOURCE_ADDRESS
    assert policy["planned_source_realizations"] == policy["planned_independent_source_paths"] == paths
    assert policy["planned_fit_records"] == fits
    assert policy["actual_to_seed_anchor"] == {cell.identity: cell.identity for cell in design.cells}
    assert all(entry[1] is frozen for entry in harness.captures)


@pytest.mark.parametrize("change", ["r3", "seed", "name", "scope", "subset", "order", "cell", "validation_full_seed"])
def test_rejects_any_nonfixed_execution_design_before_capture(harness, change):
    full, validation = full_design(), validation_design()
    changed = {
        "r3": full.validation(3), "seed": replace(full, seed=full.seed + 1),
        "name": replace(full, name="unregistered"), "scope": replace(full, scope="headline"),
        "subset": replace(full, cells=full.cells[:-1]), "order": replace(full, cells=full.cells[::-1]),
        "cell": replace(full, cells=(replace(full.cells[0], q=1.), *full.cells[1:])),
        "validation_full_seed": replace(validation, seed=full.seed),
    }[change]
    with pytest.raises(ValueError):
        harness.run(design=changed)
    assert not harness.output.exists() and not harness.captures and not harness.calls


def test_full_freeze_refusal_precedes_output_lock_or_source(harness, monkeypatch):
    def not_frozen(design, extra_files=(), **kwargs):
        assert design.replications == 2000 and kwargs["require_frozen"] is True
        raise ValueError("literal uncommitted dependency refusal")
    monkeypatch.setattr(execution, "capture", not_frozen)
    with pytest.raises(ValueError, match="uncommitted"):
        harness.run(design=full_design())
    assert not harness.output.exists() and not harness.allocations and not harness.calls


@pytest.mark.parametrize("defect", ["missing_dependency", "uncommitted_full", "wrong_commit_bytes"])
def test_fresh_freeze_requires_exact_membership_and_recorded_git_bytes(harness, monkeypatch, defect):
    if defect == "wrong_commit_bytes":
        original = guard_module.subprocess.run
        def wrong_git(arguments, *args, **kwargs):
            if arguments[:2] == ["git", "show"] and arguments[2] == MOCK_COMMIT + ":ou_irregular/week8/execution.py":
                return SimpleNamespace(returncode=0, stdout=b"different recorded module bytes")
            return original(arguments, *args, **kwargs)
        monkeypatch.setattr(guard_module.subprocess, "run", wrong_git)
    else:
        def changed_capture(*args, **kwargs):
            record = harness.capture(*args, **kwargs)
            if defect == "missing_dependency":
                record["source_files"].pop("ou_irregular/week8/source.py")
                resign(record)
            else:
                record["uncommitted_dependency_files"] = ["ou_irregular/week8/execution.py"]
            return record
        monkeypatch.setattr(execution, "capture", changed_capture)
    with pytest.raises(ValueError):
        execution.preflight(full_design(), harness.output)
    assert not harness.output.exists() and not harness.calls and not harness.allocations


def test_raw_completion_binds_e4_source_and_unchanged_inference_then_noop_resume(harness, monkeypatch):
    design = validation_design()
    result = harness.run(partition_size=2)
    assert result["project"] == "ou-week8-e4" and result["execution_kind"] == "e4"
    assert result["execution_contract"] == "week8_guarded_cell_execution_v1"
    assert result["require_frozen"] is False and result["scientific_overrides"] == []
    assert result["procedure"] == {**engine.PROCEDURE, "source_address": SOURCE_ADDRESS}
    assert result["status"] == "raw_complete_pending_audit" and result["pipeline_errors"] == 0
    assert set(result["completed_cells"]) == {cell.identity for cell in design.cells}
    assert sum(cell["paths"] for cell in result["completed_cells"].values()) == 28
    assert sum(cell["rows"] for cell in result["completed_cells"].values()) == 84
    assert harness.calls == [(cell.identity, rep, 21261011) for cell in design.cells for rep in range(4)]
    assert harness.allocations == [21261011]
    before, calls, allocations = snapshot(harness.output), list(harness.calls), list(harness.allocations)
    monkeypatch.setattr(execution, "run_lock", lambda *args: pytest.fail("cached result acquired a lock"))
    monkeypatch.setattr(execution, "_execute_cell", lambda *args: pytest.fail("cached result entered worker"))
    checked = execution.preflight(design, harness.output, resume=True)
    assert checked.cached_completion
    assert harness.run(resume=True, workers=2, partition_size=1) == result
    assert snapshot(harness.output) == before and harness.calls == calls and harness.allocations == allocations


@pytest.mark.parametrize("status", ["completed", "audited_pending_visual_review", "unknown"])
def test_closed_or_unknown_states_refuse_raw_entry_but_closed_readback_is_pure(harness, status):
    harness.run()
    metadata = state(harness.output)
    metadata["status"] = status
    write_json(harness.output / "run_metadata.json", metadata)
    before, calls, allocations = snapshot(harness.output), list(harness.calls), list(harness.allocations)
    with pytest.raises(ValueError):
        harness.run(resume=True)
    if status != "unknown":
        # This checks identity and cell accounting only; the literal archive
        # has no scientific audit or seal, and this call does not claim one.
        checked = execution.readback(validation_design(), harness.output)
        assert checked.metadata == metadata and checked.cached_completion
    else:
        with pytest.raises(ValueError):
            execution.readback(validation_design(), harness.output)
    assert snapshot(harness.output) == before and harness.calls == calls and harness.allocations == allocations


def test_review_required_clean_raw_is_noop_and_preserves_original_postprocessing_error(harness):
    harness.run()
    metadata = state(harness.output)
    metadata.update(status="implementation_review_required", postprocessing_error_type="LiteralFailure",
                    postprocessing_error="retained later-stage failure")
    write_json(harness.output / "run_metadata.json", metadata)
    before, calls = snapshot(harness.output), list(harness.calls)
    assert harness.run(resume=True) == metadata
    assert snapshot(harness.output) == before and harness.calls == calls


def test_fresh_pipeline_error_retains_all_records_and_readback_but_blocks_resume(harness, monkeypatch):
    first = validation_design().cells[0].identity
    def failure(cell, rep, source, services):
        rows = harness.source_records(cell, rep, source, services)
        if cell.identity == first and rep == 0:
            rows[0].update(pipeline_error="RuntimeError: literal retained error",
                           interval_valid=False, interval_status="pipeline_exception")
        return rows
    monkeypatch.setattr(engine, "source_records", failure)
    result = harness.run(partition_size=2)
    assert result["status"] == "implementation_review_required" and result["pipeline_errors"] == 1
    assert len(harness.calls) == 28 and sum(c["rows"] for c in result["completed_cells"].values()) == 84
    before, calls = snapshot(harness.output), list(harness.calls)
    inspected = execution.readback(validation_design(), harness.output)
    assert inspected.metadata == result and sum(c["pipeline_errors"] for c in inspected.cells.values()) == 1
    with pytest.raises(execution.SavedPipelineError):
        harness.run(resume=True)
    assert snapshot(harness.output) == before and harness.calls == calls


def first_partition(output):
    cell = validation_design().cells[0]
    directory = output / "cells" / cell.identity
    manifest = json.loads((directory / "cell_manifest.json").read_text())
    part = manifest["partitions"][0]
    return directory, manifest, part, directory / part["filename"]


def rewrite_partition(path, mutator):
    rows = read_partition(path)
    mutator(rows)
    payload = gzip.compress(("\n".join(canonical(row) for row in rows) + "\n").encode(), compresslevel=6, mtime=0)
    path.write_bytes(payload)
    return hashlib.sha256(payload).hexdigest()


@pytest.mark.parametrize("defect", ["nonempty", "missing", "nonstring"])
def test_saved_pipeline_evidence_is_checked_even_when_root_lags(harness, defect):
    harness.run(partition_size=2)
    directory, manifest, part, path = first_partition(harness.output)
    def change(rows):
        if defect == "nonempty":
            rows[0].update(pipeline_error="ArithmeticError: saved failure", interval_valid=False,
                           interval_status="pipeline_exception")
        elif defect == "missing":
            rows[0].pop("pipeline_error")
        else:
            rows[0]["pipeline_error"] = None
    part["sha256"] = rewrite_partition(path, change)
    write_json(directory / "cell_manifest.json", manifest)
    metadata = state(harness.output)
    metadata.update(status="incomplete", completed_cells={})
    write_json(harness.output / "run_metadata.json", metadata)
    before, calls = snapshot(harness.output), list(harness.calls)
    with pytest.raises(ValueError, match="(?i)(pipeline|error)"):
        harness.run(resume=True)
    assert snapshot(harness.output) == before and harness.calls == calls


def test_finite_statistical_failure_is_retained_without_pipeline_reclassification(harness):
    harness.run(partition_size=2)
    directory, manifest, part, path = first_partition(harness.output)
    def change(rows):
        rows[0].update(point_valid=False, point_status="boundary", optimizer_success=True,
                       interval_valid=False, interval_status="invalid_point_fit", theta=1e100)
    part["sha256"] = rewrite_partition(path, change)
    write_json(directory / "cell_manifest.json", manifest)
    before, calls = snapshot(harness.output), list(harness.calls)
    assert harness.run(resume=True)["status"] == "raw_complete_pending_audit"
    assert read_partition(path)[0]["theta"] == 1e100
    assert snapshot(harness.output) == before and harness.calls == calls


@pytest.mark.parametrize("defect", [
    "resolved_bytes", "identity_digest", "dependency_omission", "environment", "source_policy",
    "procedure", "schema_bool", "completed_count", "partition_hash", "unknown_root",
    "floor_specific_table", "unknown_cell", "unknown_cell_file",
])
def test_corruption_and_foreign_layout_block_before_mutation(harness, defect):
    harness.run(partition_size=2)
    metadata = state(harness.output)
    changed_metadata = True
    if defect == "resolved_bytes":
        path = harness.output / "resolved_design.json"
        path.write_bytes(path.read_bytes() + b" ")
        changed_metadata = False
    elif defect == "identity_digest":
        metadata["provenance"]["identity_sha256"] = "0" * 64
    elif defect == "dependency_omission":
        metadata["provenance"]["source_files"].pop("ou_irregular/week8/source.py")
        resign(metadata["provenance"])
    elif defect == "environment":
        metadata["provenance"]["environment"]["python"] = "0.0.0"
        resign(metadata["provenance"])
    elif defect == "source_policy":
        metadata["source_policy"]["source_namespace"] = 4
    elif defect == "procedure":
        metadata["procedure"]["max_standardized_score"] = .1
    elif defect == "schema_bool":
        metadata["schema_version"] = True
    elif defect == "completed_count":
        metadata["completed_cells"][validation_design().cells[0].identity]["paths"] = 3
    else:
        changed_metadata = False
        if defect == "partition_hash":
            _, _, _, path = first_partition(harness.output)
            path.write_bytes(path.read_bytes() + b"tamper")
        elif defect in {"unknown_root", "floor_specific_table"}:
            name = "unidentified.txt" if defect == "unknown_root" else "epsilon_comparisons.csv"
            (harness.output / name).write_text("unidentified evidence\n")
        elif defect == "unknown_cell":
            (harness.output / "cells" / ("0" * 64)).mkdir()
        else:
            directory, _, _, _ = first_partition(harness.output)
            (directory / "unidentified.txt").write_text("preserve this evidence\n")
    if changed_metadata:
        write_json(harness.output / "run_metadata.json", metadata)
    before, calls = snapshot(harness.output), list(harness.calls)
    with pytest.raises((ValueError, FileNotFoundError)):
        harness.run(resume=True)
    assert snapshot(harness.output) == before and harness.calls == calls


def test_dependency_extension_cannot_be_silently_removed_on_resume(harness):
    extra = ("README.md",)
    harness.run(extra_files=extra)
    before, calls = snapshot(harness.output), list(harness.calls)
    with pytest.raises(ValueError, match="(?i)(dependenc|manifest)"):
        harness.run(resume=True)
    assert snapshot(harness.output) == before and harness.calls == calls
    assert harness.run(resume=True, extra_files=extra)["status"] == "raw_complete_pending_audit"
    assert snapshot(harness.output) == before and harness.calls == calls


def make_orphan(harness):
    harness.run(partition_size=2)
    directory, manifest, _, _ = first_partition(harness.output)
    orphan = manifest["partitions"].pop()
    assert (orphan["start"], orphan["stop"]) == (2, 4)
    manifest.update(completed_replications=2, status="running_incomplete")
    write_json(directory / "cell_manifest.json", manifest)
    metadata = state(harness.output)
    metadata["completed_cells"].pop(validation_design().cells[0].identity)
    metadata["status"] = "incomplete"
    write_json(harness.output / "run_metadata.json", metadata)
    return directory / orphan["filename"]


def test_orphan_adoption_keeps_original_boundary_despite_new_partition_size(harness):
    orphan = make_orphan(harness)
    before, calls = orphan.read_bytes(), len(harness.calls)
    checked = execution.preflight(validation_design(), harness.output, resume=True)
    first = validation_design().cells[0].identity
    assert checked.cells[first]["paths"] == 2 and checked.cells[first]["orphan_stop"] == 4
    result = harness.run(resume=True, partition_size=1)
    assert result["status"] == "raw_complete_pending_audit"
    assert harness.calls[calls:] == [(first, 2, 21261011), (first, 3, 21261011)]
    assert orphan.read_bytes() == before
    manifest = json.loads((orphan.parent / "cell_manifest.json").read_text())
    assert [(p["start"], p["stop"]) for p in manifest["partitions"]] == [(0, 2), (2, 4)]
    assert manifest["recovered_orphans"][0]["method"] == "exact_byte_equal_deterministic_replay"


def test_byte_different_orphan_is_never_overwritten(harness):
    orphan = make_orphan(harness)
    with gzip.open(orphan, "rb") as stream:
        raw = stream.read()
    orphan.write_bytes(gzip.compress(raw, compresslevel=6, mtime=1))
    before = orphan.read_bytes()
    with pytest.raises(FileExistsError, match="immutable"):
        harness.run(resume=True, partition_size=1)
    assert orphan.read_bytes() == before and state(harness.output)["status"] == "incomplete"


def test_pipeline_error_in_orphan_blocks_before_deterministic_replay(harness):
    orphan = make_orphan(harness)
    def change(rows):
        rows[0].update(pipeline_error="RuntimeError: orphan failure", interval_valid=False,
                       interval_status="pipeline_exception")
    rewrite_partition(orphan, change)
    before, calls = snapshot(harness.output), list(harness.calls)
    with pytest.raises(execution.SavedPipelineError):
        harness.run(resume=True)
    assert snapshot(harness.output) == before and harness.calls == calls


def test_completed_durable_cells_ahead_of_root_are_adopted_without_new_records(harness):
    harness.run(partition_size=2)
    metadata = state(harness.output)
    metadata.update(status="incomplete", completed_cells={})
    write_json(harness.output / "run_metadata.json", metadata)
    parts = {name: value for name, value in snapshot(harness.output).items() if name.endswith(".jsonl.gz")}
    calls = list(harness.calls)
    result = harness.run(resume=True, partition_size=1)
    assert result["status"] == "raw_complete_pending_audit" and len(result["completed_cells"]) == 7
    assert harness.calls == calls
    assert all(snapshot(harness.output)[name] == value for name, value in parts.items())


def test_root_drift_preserves_observed_bytes_and_incident_blocks_automatic_resume(harness, monkeypatch):
    injected = b'{"status":"literal_external_change","retain":"observed bytes"}\n'
    def drift(cell, rep, source, services):
        rows = harness.source_records(cell, rep, source, services)
        (harness.output / "run_metadata.json").write_bytes(injected)
        return rows
    monkeypatch.setattr(engine, "source_records", drift)
    with pytest.raises(execution.ManifestDriftError, match="manifest drift"):
        harness.run(partition_size=1)
    assert (harness.output / "run_metadata.json").read_bytes() == injected
    incidents = list(harness.root.glob("run-manifest-incident-*"))
    assert len(incidents) == 1
    incident = json.loads((incidents[0] / "incident.json").read_text())
    assert incident["observed_manifest_overwritten_by_handler"] is False
    assert incident["observed_sha256"] == hashlib.sha256(injected).hexdigest()
    assert (incidents[0] / "observed_run_metadata.json").read_bytes() == injected
    restored = (incidents[0] / "expected_run_metadata.json").read_bytes()
    (harness.output / "run_metadata.json").write_bytes(restored)
    monkeypatch.setattr(engine, "source_records", harness.source_records)
    before, calls = snapshot(harness.root), list(harness.calls)
    with pytest.raises(ValueError, match="incident"):
        harness.run(resume=True)
    assert snapshot(harness.root) == before and harness.calls == calls


@pytest.mark.parametrize("marker", ["run_metadata.json", "validation_plan.json", "validation.json", "table_binding.json"])
def test_output_cannot_nest_in_saved_scientific_evidence(harness, marker):
    archive = harness.root / "archive"
    archive.mkdir()
    (archive / marker).write_text("{}\n")
    before = snapshot(archive)
    with pytest.raises(ValueError, match="(?i)(nested|saved)"):
        harness.run(output=archive / "new")
    assert snapshot(archive) == before and not harness.calls


def test_sync_containment_and_symlinks_block_mutation_before_creation(harness, monkeypatch):
    marker = harness.root / "sync_presence_marker"
    marker.write_text("never execute this marker\n")
    watched = harness.root / "watched"
    monkeypatch.setattr(execution, "SYNC_SCRIPT", marker)
    monkeypatch.setattr(execution, "SYNC_ROOT", watched)
    with pytest.raises(ValueError, match="synchronized"):
        harness.run(output=watched / "new")
    assert not watched.exists()
    target, alias = harness.root / "target", harness.root / "alias"
    target.mkdir()
    alias.symlink_to(target, target_is_directory=True)
    with pytest.raises(ValueError, match="symbolic"):
        harness.run(output=alias / "new")
    assert not list(target.iterdir()) and not harness.calls and not harness.allocations


def test_spawn_pool_exception_preserves_completed_work_and_resumes_only_missing_cells(harness, monkeypatch):
    futures, contexts = [], []
    class MockPool:
        def __init__(self, *, max_workers, mp_context):
            assert max_workers == 2
            contexts.append(mp_context.get_start_method())
        def __enter__(self):
            return self
        def __exit__(self, *args):
            return False
        def submit(self, function, *args):
            future = Future()
            futures.append(future)
            if len(futures) == 1:
                future.set_result(function(*args))
            else:
                future.set_exception(RuntimeError("literal worker unavailable before execution"))
            return future
    monkeypatch.setattr(execution, "ProcessPoolExecutor", MockPool)
    monkeypatch.setattr(execution, "as_completed", lambda items: list(items))
    with pytest.raises(RuntimeError, match="worker unavailable"):
        harness.run(workers=2, partition_size=2)
    metadata = state(harness.output)
    assert contexts == ["spawn"] and metadata["status"] == "incomplete"
    assert metadata["termination_type"] == "RuntimeError"
    assert len(metadata["completed_cells"]) == 1 and len(harness.calls) == 4
    result = harness.run(resume=True, workers=1, partition_size=1)
    assert result["status"] == "raw_complete_pending_audit"
    assert len(harness.calls) == len(set(harness.calls)) == 28


@pytest.mark.parametrize("defect", ["wrong_cell", "wrong_count", "boolean_count", "invented_completion"])
def test_worker_reports_cannot_replace_exact_identity_and_durable_readback(harness, monkeypatch, defect):
    def wrong(output, cell, design, identity, partition_size, source, factory, progress):
        assert factory is None and isinstance(source, source_module.E4Source)
        report = {"cell_id": cell.identity, "paths": 4, "rows": 12, "pipeline_errors": 0}
        if defect == "wrong_cell":
            report["cell_id"] = "0" * 64
        elif defect == "wrong_count":
            report["rows"] = 9
        elif defect == "boolean_count":
            report["pipeline_errors"] = False
        # Even perfect counters must not certify cells that were never saved.
        return report
    monkeypatch.setattr(execution, "_execute_cell", wrong)
    with pytest.raises(ValueError):
        harness.run()
    assert state(harness.output)["status"] == "incomplete" and not harness.calls


@pytest.mark.parametrize("arguments", [
    {"workers": 0}, {"workers": True}, {"workers": 1.0},
    {"partition_size": 0}, {"partition_size": False}, {"partition_size": 2.5},
    {"progress": 0}, {"resume": 1},
])
def test_invalid_execution_controls_rejected_before_output_or_capture(harness, arguments):
    with pytest.raises(ValueError):
        execution.run_raw(validation_design(), harness.output, **arguments)
    assert not harness.output.exists() and not harness.captures and not harness.calls


@pytest.mark.parametrize("keyword", ["source", "wald_factory", "source_factory", "require_frozen", "headline_design_digest"])
def test_no_undocumented_scientific_injection_or_freeze_override(harness, keyword):
    with pytest.raises(TypeError):
        execution.run_raw(validation_design(), harness.output, **{keyword: object()})
    assert not harness.output.exists() and not harness.allocations and not harness.calls


@pytest.mark.parametrize("extra", ["README.md", ("/tmp/absolute.py",), ("../escape.py",), ("docs//week8_gate.md",), (Path("README.md"),)])
def test_dependency_names_must_be_canonical_repository_relative_strings(harness, extra):
    with pytest.raises(ValueError):
        execution.preflight(validation_design(), harness.output, extra_files=extra)
    assert not harness.output.exists() and not harness.captures and not harness.calls


def test_lock_refusal_preserves_incomplete_evidence(harness):
    harness.run()
    metadata = state(harness.output)
    metadata.update(status="incomplete", completed_cells={})
    write_json(harness.output / "run_metadata.json", metadata)
    before, calls = snapshot(harness.output), list(harness.calls)
    with run_lock(harness.output):
        with pytest.raises(RuntimeError, match="another process"):
            harness.run(resume=True)
    assert snapshot(harness.output) == before and harness.calls == calls
