"""Independent reporting lifecycle checks on literal saved raw records only."""
from copy import deepcopy
import csv
import gzip
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from ou_irregular.experiments import analysis as historical_analysis
from ou_irregular.experiments import execution as historical_execution
from ou_irregular.experiments.provenance import canonical, digest_bytes, file_digest
from ou_irregular.experiments.storage import read_partition, run_lock
from ou_irregular.week4.simulation import OuSampleSource
from ou_irregular.week7 import execution, reporting
from ou_irregular.week7.floor_comparisons import COLUMNS
from ou_irregular.week7.floor_source import PairedFloorSource
from ou_irregular.week7.selection import floor_design, followup_design
from .test_execution_lifecycle import literal_rows


def write_json(path, value):
    Path(path).write_text(canonical(value) + "\n")


def snapshot(path):
    return {str(p.relative_to(path)): p.read_bytes() for p in Path(path).rglob("*") if p.is_file()}


@pytest.fixture
def make_raw(tmp_path, monkeypatch):
    calls = []

    def source_records(cell, rep, source, services):
        calls.append((cell.identity, rep))
        rows = literal_rows(cell, rep, source, services)
        # These are typed literal diagnostics for lifecycle fixtures, not
        # simulated schedules or a claim about OU fit/information identities.
        for row in rows:
            row.update(realized_mean_gap=cell.q / cell.theta, realized_cv=cell.cv,
                       realized_span=(cell.n - 1) * cell.q / cell.theta, clip_fraction=0.,
                       clip_count=0, min_gap=cell.q / cell.theta, raw_mean_gap=cell.q / cell.theta,
                       raw_cv=cell.cv, family="equidistant" if cell.cv == 0 else "gamma",
                       target_mean_gap=cell.q / cell.theta, requested_cv=cell.cv,
                       floor_applied=cell.cv >= 1.5,
                       floor_value=cell.floor_fraction * cell.q / cell.theta if cell.cv >= 1.5 else 0.,
                       expected_mean_gap_after_floor=cell.q / cell.theta, expected_clip_fraction=0.,
                       hessian_condition_number=2., hessian_relative_difference=1e-8,
                       hessian_whitened_difference=1e-7, standardized_score=0.)
        return rows

    def forbidden(*args, **kwargs):
        raise AssertionError("reporting must not generate or fit a random source")

    monkeypatch.setattr(historical_execution, "source_records", source_records)
    monkeypatch.setattr(OuSampleSource, "sample", forbidden)
    monkeypatch.setattr(PairedFloorSource, "sample", forbidden)

    def make(kind="followup", *, with_reporting_dependencies=True):
        design = (followup_design() if kind == "followup" else floor_design()).validation(2)
        output = tmp_path / kind
        extra = reporting.FILES if with_reporting_dependencies else ()
        metadata = execution.run_raw(kind, design, output, extra_files=extra,
                                     partition_size=1, workers=1, progress=False)
        return SimpleNamespace(kind=kind, design=design, output=output, metadata=metadata,
                               calls=calls, root=tmp_path,
                               write=lambda: reporting.write_tables(kind, design, output))

    return make


@pytest.mark.parametrize("kind", ["followup", "floor"])
def test_tables_are_read_from_saved_records_and_root_bytes_never_change(make_raw, kind):
    case = make_raw(kind)
    root_bytes = (case.output / "run_metadata.json").read_bytes()
    original_calls = list(case.calls)
    before = snapshot(case.output)
    prepared = reporting.preflight_reporting(kind, case.design, case.output)
    assert prepared.standard is None and prepared.epsilon is None
    assert snapshot(case.output) == before
    result = case.write()
    assert set(result["standard"]) == set(historical_analysis.TABLE_NAMES) | {"failed_fits.csv.gz"}
    assert (result["epsilon"] is None) == (kind == "followup")
    assert (case.output / "run_metadata.json").read_bytes() == root_bytes
    assert case.calls == original_calls
    written = snapshot(case.output)
    again = case.write()
    assert again == result and snapshot(case.output) == written
    assert case.calls == original_calls


def test_reporting_dependencies_must_have_been_frozen_in_raw_identity_before_outputs(make_raw):
    case = make_raw(with_reporting_dependencies=False)
    before = snapshot(case.output)
    with pytest.raises(ValueError, match="(?i)(dependency|manifest)"):
        case.write()
    assert snapshot(case.output) == before
    assert not (case.output / "table_manifest.json").exists()


@pytest.mark.parametrize("status", ["running_incomplete", "incomplete", "audited_pending_visual_review", "completed"])
def test_nonwritable_lifecycle_states_do_not_create_or_recreate_artifacts(make_raw, status):
    case = make_raw()
    metadata = deepcopy(case.metadata)
    metadata["status"] = status
    write_json(case.output / "run_metadata.json", metadata)
    before = snapshot(case.output)
    with pytest.raises(ValueError):
        case.write()
    assert snapshot(case.output) == before
    assert not (case.output / "parameter_summary.csv").exists()


def test_reporting_preserves_original_postprocessing_review_status_and_error(make_raw):
    case = make_raw()
    metadata = deepcopy(case.metadata)
    metadata.update(status="implementation_review_required", postprocessing_error_type="EarlierSyntheticError",
                    postprocessing_error="retain this original failure")
    write_json(case.output / "run_metadata.json", metadata)
    before = (case.output / "run_metadata.json").read_bytes()
    result = case.write()
    assert result["standard"]
    assert (case.output / "run_metadata.json").read_bytes() == before


def test_reporting_owns_nonblocking_root_lock(make_raw):
    case = make_raw()
    before = snapshot(case.output)
    with run_lock(case.output):
        with pytest.raises(RuntimeError, match="another process"):
            case.write()
    assert snapshot(case.output) == before


def test_identical_partial_standard_output_is_adopted_without_replacing_it(make_raw, monkeypatch):
    case = make_raw()
    root_bytes = (case.output / "run_metadata.json").read_bytes()
    original = historical_analysis.atomic_bytes
    calls = []

    def interrupt(path, payload, **kwargs):
        calls.append(Path(path).name)
        if len(calls) == 2:
            raise KeyboardInterrupt("synthetic interruption between table writes")
        return original(path, payload, **kwargs)

    monkeypatch.setattr(historical_analysis, "atomic_bytes", interrupt)
    with pytest.raises(KeyboardInterrupt):
        case.write()
    partial = case.output / historical_analysis.TABLE_NAMES[0]
    partial_bytes = partial.read_bytes()
    assert not (case.output / "table_manifest.json").exists()
    assert (case.output / "run_metadata.json").read_bytes() == root_bytes
    monkeypatch.setattr(historical_analysis, "atomic_bytes", original)
    case.write()
    assert partial.read_bytes() == partial_bytes
    assert (case.output / "run_metadata.json").read_bytes() == root_bytes


def test_different_undeclared_partial_table_is_not_overwritten(make_raw):
    case = make_raw()
    partial = case.output / historical_analysis.TABLE_NAMES[0]
    partial.write_bytes(b"different partial evidence; preserve bytes\n")
    before = snapshot(case.output)
    with pytest.raises(FileExistsError, match="immutable"):
        case.write()
    assert snapshot(case.output) == before


@pytest.mark.parametrize("corruption", ["missing_file", "changed_file", "missing_entry", "extra_entry", "wrong_rows"])
def test_declared_standard_tables_cannot_be_recreated_or_silently_rebound(make_raw, corruption):
    case = make_raw()
    case.write()
    path = case.output / "table_manifest.json"
    manifest = json.loads(path.read_text())
    name = "parameter_summary.csv"
    if corruption == "missing_file":
        (case.output / name).unlink()
    elif corruption == "changed_file":
        (case.output / name).write_bytes((case.output / name).read_bytes() + b"changed\n")
    elif corruption == "missing_entry":
        manifest.pop(name)
        write_json(path, manifest)
    elif corruption == "extra_entry":
        manifest["epsilon_comparisons.csv"] = {"rows": 0, "sha256": "f" * 64}
        write_json(path, manifest)
    else:
        manifest[name]["rows"] += 1
        write_json(path, manifest)
    before = snapshot(case.output)
    with pytest.raises((ValueError, FileNotFoundError)):
        case.write()
    assert snapshot(case.output) == before


@pytest.mark.parametrize("artifact", ["epsilon_comparisons.csv", "epsilon_table_manifest.json"])
def test_followup_rejects_floor_artifacts_before_writing(make_raw, artifact):
    case = make_raw()
    (case.output / artifact).write_text("{}\n")
    before = snapshot(case.output)
    with pytest.raises(ValueError, match="(?i)(follow.up|paired.floor)"):
        case.write()
    assert snapshot(case.output) == before


def test_epsilon_manifest_binds_every_actual_raw_partition_and_full_table_schema(make_raw):
    case = make_raw("floor")
    result = case.write()
    manifest = result["epsilon"]
    assert set(manifest) == {"schema_version", "project", "reporting_contract", "run_identity", "design_sha256",
                             "raw_inputs", "raw_input_sha256", "table", "difference_direction", "ratio_direction"}
    binding = manifest["raw_inputs"]
    assert binding["run_identity"] == case.metadata["provenance"]["identity_sha256"]
    assert binding["design_sha256"] == case.design.digest
    assert binding["resolved_design_sha256"] == file_digest(case.output / "resolved_design.json")
    assert binding["source_policy"] == case.metadata["source_policy"]
    assert manifest["raw_input_sha256"] == hashlib.sha256(canonical(binding).encode()).hexdigest()
    assert set(binding["cells"]) == {cell.identity for cell in case.design.cells}
    for cell in case.design.cells:
        directory = case.output / "cells" / cell.identity
        item = binding["cells"][cell.identity]
        actual = json.loads((directory / "cell_manifest.json").read_text())
        assert item["manifest_sha256"] == file_digest(directory / "cell_manifest.json")
        assert item["partitions"] == actual["partitions"]
        assert all(part["sha256"] == file_digest(directory / part["filename"]) for part in item["partitions"])
    table = manifest["table"]
    assert table["filename"] == "epsilon_comparisons.csv" and table["rows"] == 72
    assert table["columns"] == list(COLUMNS) and len(table["columns"]) == 71
    assert table["sha256"] == file_digest(case.output / table["filename"])
    with (case.output / table["filename"]).open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    expected = {(cell.identity, method, parameter) for cell in case.design.cells if cell.floor_fraction != 1e-6
                for method in ("exact", "pfml", "euler") for parameter in ("theta", "mu", "sigma")}
    assert len(rows) == 72
    assert {(r["cell_id"], r["estimator"], r["parameter"]) for r in rows} == expected


@pytest.mark.parametrize("corruption", ["missing_csv", "changed_csv", "wrong_input_hash", "wrong_policy",
                                       "wrong_columns", "wrong_direction"])
def test_declared_epsilon_evidence_is_verified_before_outputs(make_raw, corruption):
    case = make_raw("floor")
    case.write()
    path = case.output / "epsilon_table_manifest.json"
    manifest = json.loads(path.read_text())
    if corruption == "missing_csv":
        (case.output / "epsilon_comparisons.csv").unlink()
    elif corruption == "changed_csv":
        target = case.output / "epsilon_comparisons.csv"
        target.write_bytes(target.read_bytes() + b"changed\n")
    elif corruption == "wrong_input_hash":
        manifest["raw_input_sha256"] = "0" * 64
    elif corruption == "wrong_policy":
        manifest["raw_inputs"]["source_policy"]["name"] = "unverified_pairing"
    elif corruption == "wrong_columns":
        manifest["table"]["columns"] = manifest["table"]["columns"][:-1]
    else:
        manifest["difference_direction"] = "reversed"
    if corruption not in {"missing_csv", "changed_csv"}:
        write_json(path, manifest)
    before = snapshot(case.output)
    with pytest.raises((ValueError, FileNotFoundError)):
        case.write()
    assert snapshot(case.output) == before


def test_partial_epsilon_csv_is_adopted_only_by_exact_bytes(make_raw, monkeypatch):
    case = make_raw("floor")
    original = reporting.atomic_json

    def interrupt(path, value, **kwargs):
        if Path(path).name == "epsilon_table_manifest.json":
            raise KeyboardInterrupt("synthetic interruption before epsilon manifest")
        return original(path, value, **kwargs)

    monkeypatch.setattr(reporting, "atomic_json", interrupt)
    with pytest.raises(KeyboardInterrupt):
        case.write()
    path = case.output / "epsilon_comparisons.csv"
    before = path.read_bytes()
    assert not (case.output / "epsilon_table_manifest.json").exists()
    monkeypatch.setattr(reporting, "atomic_json", original)
    case.write()
    assert path.read_bytes() == before


def mutate_raw_partition(case):
    cell = next(cell for cell in case.design.cells if cell.cv != 0)
    directory = case.output / "cells" / cell.identity
    manifest_path = directory / "cell_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    part = manifest["partitions"][0]
    path = directory / part["filename"]
    rows = read_partition(path)
    rows[0]["theta"] += .125
    payload = gzip.compress(("\n".join(canonical(row) for row in rows) + "\n").encode(), compresslevel=6, mtime=0)
    path.write_bytes(payload)
    part["sha256"] = digest_bytes(payload)
    write_json(manifest_path, manifest)
    return path, payload, manifest_path, manifest_path.read_bytes()


def test_raw_input_change_with_unchanged_root_blocks_epsilon_commit(make_raw, monkeypatch):
    case = make_raw("floor")
    root_bytes = (case.output / "run_metadata.json").read_bytes()
    original = reporting.standard_write_tables
    changed = []

    def interfere(output, design):
        result = original(output, design)
        changed.append(mutate_raw_partition(case))
        return result

    monkeypatch.setattr(reporting, "standard_write_tables", interfere)
    with pytest.raises(ValueError, match="raw input binding"):
        case.write()
    path, payload, manifest_path, manifest_bytes = changed[0]
    assert path.read_bytes() == payload and manifest_path.read_bytes() == manifest_bytes
    assert (case.output / "run_metadata.json").read_bytes() == root_bytes
    assert (case.output / "table_manifest.json").is_file()
    assert not (case.output / "epsilon_comparisons.csv").exists()
    assert not (case.output / "epsilon_table_manifest.json").exists()


def test_root_drift_during_reporting_is_preserved_and_blocks_automatic_retry(make_raw, monkeypatch):
    case = make_raw()
    original = reporting.standard_write_tables
    injected = b'{"status":"unattributed_external_reporting_state"}\n'

    def interfere(output, design):
        result = original(output, design)
        (case.output / "run_metadata.json").write_bytes(injected)
        return result

    monkeypatch.setattr(reporting, "standard_write_tables", interfere)
    with pytest.raises(RuntimeError, match="manifest drift"):
        case.write()
    assert (case.output / "run_metadata.json").read_bytes() == injected
    incident = next(p for p in case.root.iterdir() if p.is_dir() and "manifest-incident-" in p.name)
    assert (incident / "observed_run_metadata.json").read_bytes() == injected
    (case.output / "run_metadata.json").write_bytes((incident / "expected_run_metadata.json").read_bytes())
    monkeypatch.setattr(reporting, "standard_write_tables", original)
    before = snapshot(case.root)
    with pytest.raises(ValueError, match="incident"):
        case.write()
    assert snapshot(case.root) == before
