"""Independent lifecycle checks using literal synthetic observations only.

These tests deliberately replace archive discovery and scientific replay. They
test durable-state behavior and exact binding, not any E1/E2 source outcome.
"""
from copy import deepcopy
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from ou_irregular.experiments.design import Cell, ESTIMATORS, RunDesign
from ou_irregular.experiments.provenance import canonical, digest_bytes, file_digest
from ou_irregular.experiments.storage import run_lock
from ou_irregular.inference.contracts import Sample
from ou_irregular.week4.simulation import array_digest
from ou_irregular.week7 import path_diagnostics as paths


def write_json(path, value):
    Path(path).write_text(canonical(value) + "\n")


def snapshot(directory):
    return {str(p.relative_to(directory)): p.read_bytes()
            for p in Path(directory).rglob("*") if p.is_file()}


@pytest.fixture
def run_fixture(tmp_path, monkeypatch):
    """Two addressed records, but no randomized source or headline archive."""
    cell = Cell(n=8, q=.5, cv=0)
    design = RunDesign("literal_lifecycle_fixture", 0, 2, (cell,), "validation")
    selected = {"scope": "synthetic_lifecycle_fixture", "sources": [
        {"experiment": "E1", "cell_id": cell.identity, "rep": rep,
         "reasons": ["synthetic:retained_failure"]} for rep in (0, 1)]}
    sample = Sample(values=np.array([.2, -.1, .35, .4, -.2, .15, .05, -.3]),
                    times=np.arange(8.) / 2)
    originals = {}
    for item in selected["sources"]:
        originals[paths._source_key(item)] = {
            method: {"cell_id": cell.identity, "rep": item["rep"], "estimator": method,
                     "point_valid": False, "interval_valid": False,
                     "point_status": "synthetic_boundary", "interval_status": "point_invalid",
                     "candidate_theta": 100., "theta": None, "mu": None, "sigma": None,
                     "point_diagnostics_json": '{"candidate":100.0,"accepted":false}',
                     "times_sha256": array_digest(sample.times),
                     "observations_sha256": array_digest(sample.values)}
            for method in ESTIMATORS}
    diagnostics = paths.analyze_sample(sample, originals[paths._source_key(selected["sources"][0])])
    calls = []

    def one(item, original, original_design, identity):
        calls.append(item["rep"])
        return {"schema_version": 1, "run_identity": identity, "selection": deepcopy(item),
                "status": "evaluated", "headline_results_replaced": False,
                "original_records": deepcopy(original),
                "original_replay": {"status": "passed", "records": deepcopy(original),
                    "source_diagnostics": {"times_sha256": array_digest(sample.times),
                                           "observations_sha256": array_digest(sample.values)},
                    "array_hashes": {"times_sha256": array_digest(sample.times),
                                     "observations_sha256": array_digest(sample.values)},
                    "rtol": paths.SETTINGS["replay_rtol"], "atol": paths.SETTINGS["replay_atol"]},
                "diagnostic": deepcopy(diagnostics)}

    def capture(descriptor, extra_files, require_frozen):
        return {"design": descriptor.to_dict(), "design_sha256": descriptor.digest,
                "source_files": {"fixture_dependency.py": "1" * 64},
                "identity_sha256": digest_bytes(descriptor.serialized.encode()),
                "code_commit": "1" * 40, "uncommitted_dependency_files": [],
                "environment": {"fixture_only": True}}

    monkeypatch.setattr(paths, "load_saved_design", lambda directory: design)
    monkeypatch.setattr(paths, "build_selection", lambda *args: deepcopy(selected))
    monkeypatch.setattr(paths, "_selected_originals", lambda *args: deepcopy(originals))
    monkeypatch.setattr(paths, "capture", capture)
    monkeypatch.setattr(paths, "verify", lambda *args, **kwargs: True)
    monkeypatch.setattr(paths, "_one_selected", one)

    def run(**kwargs):
        return paths.run_selected(selected, tmp_path / "E1", tmp_path / "E2", tmp_path / "out",
                                  require_frozen=False, progress=False, **kwargs)

    return SimpleNamespace(selected=selected, originals=originals, design=design, sample=sample,
                           calls=calls, one=one, run=run, output=tmp_path / "out", root=tmp_path)


def interrupt_after_first(fixture, monkeypatch):
    def interrupt(item, *args):
        if item["rep"] == 1:
            raise KeyboardInterrupt("fixed synthetic interruption")
        return fixture.one(item, *args)
    monkeypatch.setattr(paths, "_one_selected", interrupt)
    with pytest.raises(KeyboardInterrupt):
        fixture.run()
    assert fixture.calls == [0]
    monkeypatch.setattr(paths, "_one_selected", fixture.one)


def test_completed_output_is_immutable_and_audit_scope_is_honest(run_fixture):
    result = run_fixture.run()
    assert result["status"] == "completed_diagnostic_evidence"
    assert run_fixture.calls == [0, 1]
    report = json.loads((run_fixture.output / "audit.json").read_text())
    assert report["verified_fit_replays"] == 6
    assert report["retained_sources"] == report["planned_sources"] == 2
    assert report["independent_recalculation_of_all_diagnostic_math"] is False
    before = snapshot(run_fixture.output)
    with pytest.raises(ValueError, match="closed"):
        run_fixture.run(resume=True)
    assert snapshot(run_fixture.output) == before
    assert run_fixture.calls == [0, 1]


def test_verified_resume_reuses_committed_source_without_replay(run_fixture, monkeypatch):
    interrupt_after_first(run_fixture, monkeypatch)
    part = run_fixture.output / "sources" / paths._filename(run_fixture.selected["sources"][0])
    before = part.read_bytes()
    result = run_fixture.run(resume=True)
    assert result["status"] == "completed_diagnostic_evidence"
    assert run_fixture.calls == [0, 1]
    assert part.read_bytes() == before
    assert len(result["resume_events"]) == 1


def test_committed_diagnostic_failure_is_retained_never_retried(run_fixture, monkeypatch):
    def failed_first(item, *args):
        if item["rep"] == 1:
            raise KeyboardInterrupt("after committed failure")
        result = run_fixture.one(item, *args)
        result.update(status="diagnostic_failure", diagnostic={"status": "failed", "error": "SyntheticError: fixed"})
        return result
    monkeypatch.setattr(paths, "_one_selected", failed_first)
    with pytest.raises(KeyboardInterrupt):
        run_fixture.run()
    monkeypatch.setattr(paths, "_one_selected", run_fixture.one)
    result = run_fixture.run(resume=True)
    assert result["status"] == "diagnostic_review_required"
    assert run_fixture.calls == [0, 1]
    report = json.loads((run_fixture.output / "audit.json").read_text())
    assert len(report["failures"]) == 1
    before = snapshot(run_fixture.output)
    with pytest.raises(ValueError, match="closed"):
        run_fixture.run(resume=True)
    assert snapshot(run_fixture.output) == before


def make_orphan(fixture, monkeypatch):
    """Simulate loss between an immutable source flush and its root commit."""
    original_atomic = paths.atomic_bytes
    fired = False

    def crash(path, payload, **kwargs):
        nonlocal fired
        value = json.loads(payload) if Path(path).name == "run_metadata.json" else {}
        if Path(path).name == "run_metadata.json" and value.get("source_files") and not fired:
            fired = True
            raise KeyboardInterrupt("after source flush")
        return original_atomic(path, payload, **kwargs)

    monkeypatch.setattr(paths, "atomic_bytes", crash)
    with pytest.raises(KeyboardInterrupt):
        fixture.run()
    # A graceful finally can commit the in-memory prefix. Model a hard process
    # loss by restoring precisely the preceding initial manifest prefix.
    metadata_path = fixture.output / "run_metadata.json"
    metadata = json.loads(metadata_path.read_text())
    metadata["source_files"] = {}
    metadata["status"] = "incomplete"
    write_json(metadata_path, metadata)
    monkeypatch.setattr(paths, "atomic_bytes", original_atomic)


def test_only_exact_byte_equal_orphan_may_be_adopted(run_fixture, monkeypatch):
    make_orphan(run_fixture, monkeypatch)
    path = run_fixture.output / "sources" / paths._filename(run_fixture.selected["sources"][0])
    before = path.read_bytes()
    result = run_fixture.run(resume=True)
    assert result["status"] == "completed_diagnostic_evidence"
    assert run_fixture.calls == [0, 0, 1]
    assert path.read_bytes() == before
    assert [event["method"] for event in result["recovered_orphans"]] == [
        "exact_byte_equal_deterministic_replay"]


def test_semantically_equal_but_byte_changed_orphan_is_not_overwritten(run_fixture, monkeypatch):
    make_orphan(run_fixture, monkeypatch)
    path = run_fixture.output / "sources" / paths._filename(run_fixture.selected["sources"][0])
    path.write_text(json.dumps(json.loads(path.read_text()), indent=2) + "\n")
    before = path.read_bytes()
    with pytest.raises(FileExistsError, match="immutable"):
        run_fixture.run(resume=True)
    assert path.read_bytes() == before
    assert run_fixture.calls == [0, 0]


@pytest.mark.parametrize("change", ["selection_bytes", "descriptor_bytes", "dependencies", "payload_digest",
                                   "source_prefix", "unknown_artifact", "missing_source", "closed_state"])
def test_invalid_resume_is_rejected_before_any_mutation_or_replay(run_fixture, monkeypatch, change):
    interrupt_after_first(run_fixture, monkeypatch)
    out = run_fixture.output
    metadata_path = out / "run_metadata.json"
    metadata = json.loads(metadata_path.read_text())
    if change == "selection_bytes":
        (out / "selection.json").write_bytes((out / "selection.json").read_bytes() + b" ")
    elif change == "descriptor_bytes":
        (out / "diagnostic_design.json").write_bytes((out / "diagnostic_design.json").read_bytes() + b" ")
    elif change == "dependencies":
        metadata["provenance"]["source_files"] = {}
        write_json(metadata_path, metadata)
    elif change == "payload_digest":
        part = next((out / "sources").iterdir())
        part.write_bytes(part.read_bytes() + b" ")
    elif change == "source_prefix":
        old_name = paths._filename(run_fixture.selected["sources"][0])
        other_name = paths._filename(run_fixture.selected["sources"][1])
        metadata["source_files"][other_name] = metadata["source_files"].pop(old_name)
        write_json(metadata_path, metadata)
    elif change == "unknown_artifact":
        (out / "extra.txt").write_text("unidentified evidence")
    elif change == "missing_source":
        next((out / "sources").iterdir()).unlink()
    elif change == "closed_state":
        metadata["status"] = "diagnostic_review_required"
        write_json(metadata_path, metadata)
    before = snapshot(out)
    with pytest.raises((ValueError, FileNotFoundError)):
        run_fixture.run(resume=True)
    assert snapshot(out) == before
    assert run_fixture.calls == [0]


def test_root_lock_excludes_another_orchestrator(run_fixture, monkeypatch):
    interrupt_after_first(run_fixture, monkeypatch)
    before = snapshot(run_fixture.output)
    with run_lock(run_fixture.output):
        with pytest.raises(RuntimeError, match="another process"):
            run_fixture.run(resume=True)
    assert run_fixture.calls == [0]
    assert snapshot(run_fixture.output) == before


@pytest.mark.parametrize("label", ["missing_diagnostic", "empty_evaluated", "missing_method", "false_schema_type",
                                   "false_rep_type", "missing_precision", "missing_profile", "replacement_claim",
                                   "missing_boundary_evidence"])
def test_payload_readback_rejects_missing_or_misrepresented_diagnostics(run_fixture, label):
    item = run_fixture.selected["sources"][1]
    originals = run_fixture.originals[paths._source_key(item)]
    value = run_fixture.one(item, originals, run_fixture.design, "fixture_identity")
    if label == "missing_diagnostic":
        value.pop("diagnostic")
    elif label == "empty_evaluated":
        value["diagnostic"] = {"status": "evaluated"}
    elif label == "missing_method":
        value["diagnostic"]["fits"].pop("euler")
    elif label == "false_schema_type":
        value["schema_version"] = True
    elif label == "false_rep_type":
        value["selection"]["rep"] = True
    elif label == "missing_precision":
        value["diagnostic"].pop("precision")
    elif label == "missing_profile":
        value["diagnostic"].pop("exact_profile")
    elif label == "replacement_claim":
        value["diagnostic"]["headline_results_replaced"] = True
    else:
        value["diagnostic"]["ar1"] = {"status": "evaluated"}
    with pytest.raises(ValueError):
        paths._payload_identity(value, item, originals, "fixture_identity")


@pytest.mark.parametrize("status", ["original_replay_failed", "diagnostic_failure"])
def test_failed_payload_requires_explicit_retained_failure(run_fixture, status):
    item = run_fixture.selected["sources"][0]
    originals = run_fixture.originals[paths._source_key(item)]
    value = run_fixture.one(item, originals, run_fixture.design, "fixture_identity")
    value["status"] = status
    if status == "original_replay_failed":
        value["original_replay"] = {"status": "failed"}
        value.pop("diagnostic")
    else:
        value["diagnostic"] = {"status": "failed"}
    with pytest.raises(ValueError):
        paths._payload_identity(value, item, originals, "fixture_identity")


def test_frozen_run_requires_actual_recorded_git_bytes(monkeypatch):
    monkeypatch.setattr(paths, "verify", lambda *args, **kwargs: True)
    actual = b"frozen dependency bytes\n"
    provenance = {"source_files": {"x.py": digest_bytes(actual)}, "code_commit": "a" * 40,
                  "uncommitted_dependency_files": []}
    calls = []

    def git(args, **kwargs):
        calls.append(args)
        return SimpleNamespace(returncode=0, stdout=b"different bytes\n")

    monkeypatch.setattr(paths.subprocess, "run", git)
    with pytest.raises(ValueError, match="recorded freeze"):
        paths._check_frozen(provenance, object(), deepcopy(provenance), True)
    assert calls == [["git", "show", "a" * 40 + ":x.py"]]


def test_full_replay_rejects_external_selection_before_any_output(run_fixture):
    external = run_fixture.root / "external_selection.json"
    write_json(external, run_fixture.selected)
    with pytest.raises(ValueError, match="fixed repository selection"):
        paths.run_selected(external, run_fixture.root / "E1", run_fixture.root / "E2",
                           run_fixture.output, require_frozen=True, progress=False)
    assert not run_fixture.output.exists()
    assert not run_fixture.calls


def test_new_output_cannot_be_nested_in_an_unrelated_preserved_archive(run_fixture):
    archive = run_fixture.root / "unrelated_archived_experiment"
    archive.mkdir()
    (archive / "run_metadata.json").write_text('{"status":"completed","historical":true}\n')
    (archive / "evidence.txt").write_text("preserve these original bytes\n")
    before = snapshot(archive)
    with pytest.raises(ValueError, match="(?i)(archive|ancestor|nested|preserv)"):
        paths.run_selected(run_fixture.selected, run_fixture.root / "E1", run_fixture.root / "E2",
                           archive / "new_diagnostics", require_frozen=False, progress=False)
    assert snapshot(archive) == before
    assert not run_fixture.calls


def test_source_replay_checks_nested_candidates_and_discrete_types():
    original = {"point_valid": False, "rep": 1, "candidate": 100.,
                "point_diagnostics_json": '{"candidate":[100.0,null],"status":"boundary"}'}
    within = deepcopy(original)
    within["candidate"] += 1e-9
    paths._compare_replay(original, within)
    for changed in ({**original, "point_valid": 0}, {**original, "rep": True},
                    {**original, "candidate": None},
                    {**original, "point_diagnostics_json": '{"candidate":[100.0,0],"status":"boundary"}'},
                    {**original, "point_diagnostics_json": '{"candidate":[100.0,null],"status":"valid"}'}):
        with pytest.raises(ValueError):
            paths._compare_replay(original, changed)


def test_external_manifest_drift_is_preserved_without_stale_overwrite(run_fixture, monkeypatch):
    injected = b'{"status":"external_unattributed_state","evidence":"retain exact bytes"}\n'

    def interference(item, *args):
        value = run_fixture.one(item, *args)
        (run_fixture.output / "run_metadata.json").write_bytes(injected)
        return value

    monkeypatch.setattr(paths, "_one_selected", interference)
    with pytest.raises((ValueError, RuntimeError), match="(?i)(manifest|changed|drift|state)"):
        run_fixture.run()
    assert (run_fixture.output / "run_metadata.json").read_bytes() == injected
    assert run_fixture.calls == [0]
    siblings = [p for p in run_fixture.root.iterdir() if p.is_dir() and p != run_fixture.output]
    assert any(p.read_bytes() == injected for d in siblings for p in d.rglob("*") if p.is_file())


def test_drift_during_root_write_is_detected_and_not_overwritten(run_fixture, monkeypatch):
    injected = b'{"status":"synthetic_postwrite_drift"}\n'
    original_atomic = paths.atomic_bytes

    def interference(path, payload, **kwargs):
        result = original_atomic(path, payload, **kwargs)
        if Path(path) == run_fixture.output / "run_metadata.json":
            Path(path).write_bytes(injected)
        return result

    monkeypatch.setattr(paths, "atomic_bytes", interference)
    with pytest.raises(RuntimeError, match="manifest drift"):
        run_fixture.run()
    assert (run_fixture.output / "run_metadata.json").read_bytes() == injected
    assert not run_fixture.calls


def test_prior_incident_blocks_mutation_even_after_external_root_restoration(run_fixture, monkeypatch):
    injected = b'{"status":"external_state_for_incident_test"}\n'

    def interference(item, *args):
        value = run_fixture.one(item, *args)
        (run_fixture.output / "run_metadata.json").write_bytes(injected)
        return value

    monkeypatch.setattr(paths, "_one_selected", interference)
    with pytest.raises(RuntimeError, match="manifest drift"):
        run_fixture.run()
    incident = next(p for p in run_fixture.root.iterdir()
                    if p.is_dir() and p.name.startswith("out-manifest-incident-"))
    # Fault injection models a later external restoration, not an authorized
    # recovery procedure. The preserved incident must still block the runner.
    (run_fixture.output / "run_metadata.json").write_bytes((incident / "expected_run_metadata.json").read_bytes())
    monkeypatch.setattr(paths, "_one_selected", run_fixture.one)
    before = snapshot(run_fixture.root)
    with pytest.raises((ValueError, RuntimeError), match="incident"):
        run_fixture.run(resume=True)
    assert run_fixture.calls == [0]
    assert snapshot(run_fixture.root) == before


@pytest.fixture(scope="module")
def valid_point_payload():
    """An arbitrary finite point on literal data exercises stored schema only.

    No claim is made that this point maximizes a likelihood. Its interval is
    deliberately marked unavailable; the diagnostic must still retain ladders.
    """
    sample = Sample(values=[.2, .1, .35, .4, .2, .15, .05, -.1], times=np.arange(8.) / 2)
    item = {"experiment": "E1", "cell_id": "a" * 64, "rep": 0,
            "reasons": ["synthetic:schema_only"]}
    hashes = {"times_sha256": array_digest(sample.times),
              "observations_sha256": array_digest(sample.values)}
    originals = {method: {"estimator": method, "point_valid": True,
                          "point_status": "synthetic_finite_point", "interval_valid": False,
                          "interval_status": "synthetic_no_interval", "theta": .8,
                          "mu": .1, "sigma": .4, "nll": 10., **hashes}
                 for method in ESTIMATORS}
    diagnostic = paths.analyze_sample(sample, originals)
    value = {"schema_version": 1, "run_identity": "schema_identity", "selection": item,
             "status": "evaluated", "headline_results_replaced": False,
             "original_records": originals,
             "original_replay": {"status": "passed", "records": deepcopy(originals),
                                  "source_diagnostics": deepcopy(hashes), "array_hashes": deepcopy(hashes),
                                  "rtol": paths.SETTINGS["replay_rtol"],
                                  "atol": paths.SETTINGS["replay_atol"]},
             "diagnostic": diagnostic}
    return value, item, originals


@pytest.mark.parametrize("change", ["missing_ladder", "lost_step", "lost_adjacent", "changed_step",
                                   "missing_extended_ladder", "matrix_wrong_shape", "missing_profile_start",
                                   "truncated_profile_grid", "replacement_interval", "missing_best_candidate",
                                   "missing_limiting_nll", "missing_stored_full_nll"])
def test_valid_point_diagnostics_retain_all_planned_steps_and_searches(valid_point_payload, change):
    value, item, originals = deepcopy(valid_point_payload)
    paths._payload_identity(value, item, originals, "schema_identity")
    exact = value["diagnostic"]["fits"]["exact"]
    ladder = exact["float64_ladder"]
    profile = value["diagnostic"]["exact_profile"]
    if change == "missing_ladder":
        exact.pop("float64_ladder")
    elif change == "lost_step":
        ladder["entries"].pop()
    elif change == "lost_adjacent":
        ladder["adjacent"].pop()
    elif change == "changed_step":
        ladder["entries"][0]["relative_step"] = .002
    elif change == "missing_extended_ladder":
        exact.pop("independent_longdouble_ladder")
    elif change == "matrix_wrong_shape":
        ladder["entries"][0]["hessian"] = [[1., 0.], [0., 1.]]
    elif change == "missing_profile_start":
        profile["local_searches"].pop()
    elif change == "truncated_profile_grid":
        profile["grid_nll"].pop()
    elif change == "replacement_interval":
        profile["replacement_fit_or_interval"] = True
    elif change == "missing_best_candidate":
        profile.pop("best_finite")
    elif change == "missing_limiting_nll":
        profile["limits"]["brownian"] = {}
    elif change == "missing_stored_full_nll":
        profile["stored_full"] = {}
    with pytest.raises(ValueError):
        paths._payload_identity(value, item, originals, "schema_identity")
