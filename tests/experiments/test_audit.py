"""Independent artifact verification and adversarial corruption tests.

All generated datasets use a separate small validation design. Mutations affect
only temporary copies; neither historical evidence nor headline streams run.
"""
import gzip
import json
from pathlib import Path
import shutil

import numpy as np
import pandas as pd
import pytest

from ou_irregular.experiments.analysis import write_tables
from ou_irregular.experiments import audit
from ou_irregular.experiments.audit import AUDIT_FILES, audit_run
from ou_irregular.experiments.design import Cell, RunDesign
from ou_irregular.experiments.execution import run_raw
from ou_irregular.experiments.provenance import canonical, file_digest


@pytest.fixture(scope="module")
def completed_run(tmp_path_factory):
    output = tmp_path_factory.mktemp("auditor_validation")
    design = RunDesign("auditor_validation", 333333333, 4,
                       (Cell(40, .5, 0), Cell(40, .5, 2)), "validation")
    run_raw(design, output, extra_files=AUDIT_FILES, partition_size=2, progress=False)
    write_tables(output, design)
    return output, design


@pytest.fixture
def copied_run(completed_run, tmp_path):
    original, design = completed_run
    output = tmp_path/"run"
    shutil.copytree(original, output)
    return output, design


def _rewrite_json(path, value):
    Path(path).write_text(canonical(value)+"\n")


def _mutate_raw(output, cell, mutation):
    directory = output/"cells"/cell.identity
    manifest = json.loads((directory/"cell_manifest.json").read_text())
    part = manifest["partitions"][0]
    path = directory/part["filename"]
    rows = [json.loads(line) for line in gzip.decompress(path.read_bytes()).decode().splitlines()]
    mutation(rows)
    path.write_bytes(gzip.compress(("\n".join(canonical(r) for r in rows)+"\n").encode(), mtime=0))
    part["sha256"], part["rows"] = file_digest(path), len(rows)
    _rewrite_json(directory/"cell_manifest.json", manifest)


def _mutate_table(output, name, column):
    path = output/name
    data = pd.read_csv(path, float_precision="round_trip")
    data.loc[0, column] = data.loc[0, column]+1
    data.to_csv(path, index=False, float_format="%.17g")
    manifest = json.loads((output/"table_manifest.json").read_text())
    manifest[name]["sha256"] = file_digest(path)
    _rewrite_json(output/"table_manifest.json", manifest)


def test_true_small_run_passes_independent_arithmetic_and_every_cell_replay(completed_run, monkeypatch):
    output, design = completed_run
    before = {str(p.relative_to(output)): file_digest(p) for p in output.rglob("*") if p.is_file()}
    calls = []
    actual_source = audit.OuSampleSource

    class CountingSource(actual_source):
        def sample(self, cell, rep):
            calls.append((cell["cell_id"], rep))
            return super().sample(cell, rep)

    monkeypatch.setattr(audit, "OuSampleSource", CountingSource)
    report = audit_run(output, design)
    assert report["status"] == "passed" and report["full_completion_audit"]
    assert report["cells_verified"] == 2 and report["partitions_verified"] == 4
    assert report["planned_paths"] == 8 and report["planned_fit_records"] == 24
    assert report["replayed_source_paths"] == 4 and report["replayed_fits"] == 12
    assert set(calls) == {(c.identity, r) for c in design.cells for r in (0, 3)}
    assert len(calls) == 4
    assert report["tables_verified"]["parameter_summary.csv"] == 18
    json.dumps(report, allow_nan=False)
    assert before == {str(p.relative_to(output)): file_digest(p) for p in output.rglob("*") if p.is_file()}


def test_arithmetic_only_audit_does_not_claim_replay(completed_run):
    output, design = completed_run
    report = audit_run(output, design, replay=False)
    assert report["status"] == "passed" and not report["full_completion_audit"]
    assert not report["replay_performed"] and report["replayed_fits"] == 0


def test_missing_whole_cell_and_estimator_group_are_rejected(copied_run):
    output, design = copied_run
    directory = output/"cells"/design.cells[1].identity
    directory.rename(output/"removed_cell")
    with pytest.raises(ValueError, match="planned cell set"):
        audit_run(output, design, replay=False)
    (output/"removed_cell").rename(directory)
    _mutate_raw(output, design.cells[0], lambda rows: rows.__setitem__(slice(None), [r for r in rows if r["estimator"] != "euler"]))
    with pytest.raises(ValueError, match="row count|planned estimator"):
        audit_run(output, design, replay=False)


def test_changed_partition_digest_is_rejected(copied_run):
    output, design = copied_run
    path = next((output/"cells"/design.cells[0].identity).glob("part-*.jsonl.gz"))
    path.write_bytes(path.read_bytes()+b"changed")
    with pytest.raises(ValueError, match="partition digest"):
        audit_run(output, design, replay=False)


def test_changed_table_digest_is_rejected(copied_run):
    output, design = copied_run
    path = output/"parameter_summary.csv"
    path.write_bytes(path.read_bytes()+b"\n")
    with pytest.raises(ValueError, match="table digest"):
        audit_run(output, design, replay=False)


@pytest.mark.parametrize(("name", "column"), [
    ("parameter_summary.csv", "bias"),
    ("parameter_summary.csv", "coverage_operational"),
    ("regular_comparisons.csv", "excess_bias"),
    ("paired_comparisons.csv", "rmse_ratio"),
    ("spacing_summary.csv", "realized_span_mean"),
    ("outcome_counts.csv", "count"),
    ("numerical_diagnostics.csv", "hessian_condition_number_max"),
])
def test_rehashed_arithmetic_corruption_is_rejected(copied_run, name, column):
    output, design = copied_run
    _mutate_table(output, name, column)
    with pytest.raises(ValueError, match="arithmetic mismatch|count mismatch"):
        audit_run(output, design, replay=False)


def test_rehashed_source_hash_corruption_is_caught_by_fixed_replay(copied_run):
    output, design = copied_run

    def corrupt(rows):
        for row in rows:
            if row["rep"] == 0:
                row["observations_sha256"] = "0"*64

    _mutate_raw(output, design.cells[0], corrupt)
    with pytest.raises(ValueError, match="source replay|source array hash"):
        audit_run(output, design, replay=True)


def test_source_pairing_and_registered_seed_addresses_are_checked(copied_run):
    output, design = copied_run

    def corrupt(rows):
        rows[0]["path_spawn_key"] = "[1,2,3]"

    _mutate_raw(output, design.cells[0], corrupt)
    with pytest.raises(ValueError, match="source addressing|shared source"):
        audit_run(output, design, replay=False)


def test_rehashed_valid_interval_with_failed_gate_is_rejected(copied_run):
    output, design = copied_run

    def corrupt(rows):
        row = next(r for r in rows if r["interval_valid"])
        diag = json.loads(row["interval_diagnostics_json"])
        diag["hessian_diagnostics"]["whitened_relative_difference"] = .02
        row["hessian_whitened_difference"] = .02
        row["interval_diagnostics_json"] = canonical(diag)

    _mutate_raw(output, design.cells[0], corrupt)
    with pytest.raises(ValueError, match="numerical gate"):
        audit_run(output, design, replay=False)


@pytest.mark.parametrize("field", ["theta", "nll"])
def test_invalid_valid_point_candidates_are_rejected(copied_run, field):
    output, design = copied_run

    def corrupt(rows):
        row = next(r for r in rows if r["point_valid"])
        row[field] = True if field == "theta" else None

    _mutate_raw(output, design.cells[0], corrupt)
    with pytest.raises(ValueError, match="nonnumeric"):
        audit_run(output, design, replay=False)


def test_unclipped_endpoint_identity_is_audited(copied_run):
    output, design = copied_run

    def corrupt(rows):
        row = next(r for r in rows if r["interval_valid"])
        row["theta_lower"] += .01

    _mutate_raw(output, design.cells[0], corrupt)
    with pytest.raises(ValueError, match="natural lower endpoints"):
        audit_run(output, design, replay=False)


def test_explicit_dependency_manifest_cannot_omit_auditor(copied_run):
    output, design = copied_run
    path = output/"run_metadata.json"
    metadata = json.loads(path.read_text())
    metadata["provenance"]["source_files"].pop("ou_irregular/experiments/audit.py")
    _rewrite_json(path, metadata)
    with pytest.raises(ValueError, match="manifest omits"):
        audit_run(output, design, replay=False)


def test_incomplete_or_pipeline_review_required_run_cannot_pass(copied_run):
    output, design = copied_run
    path = output/"run_metadata.json"
    metadata = json.loads(path.read_text())
    metadata["status"] = "implementation_review_required"
    _rewrite_json(path, metadata)
    with pytest.raises(ValueError, match="incomplete or requires"):
        audit_run(output, design, replay=False)


def test_real_boundary_failures_are_retained_audited_and_replayed(tmp_path):
    design = RunDesign("boundary_auditor_validation", 333333333, 4,
                       (Cell(6, 2, 0),), "validation")
    run_raw(design, tmp_path, extra_files=AUDIT_FILES, progress=False)
    write_tables(tmp_path, design)
    report = audit_run(tmp_path, design)
    assert report["retained_failure_records"] == 6
    assert report["point_valid_records"] == report["interval_valid_records"] == 6
    assert report["replayed_fits"] == 6
    path = tmp_path/"failed_fits.csv.gz"
    failed = pd.read_csv(path)
    assert set(failed.point_status) == {"boundary"}
    assert set(failed.interval_status) == {"invalid_point_fit"}
    payload = failed.iloc[1:].to_csv(index=False, float_format="%.17g").encode()
    path.write_bytes(gzip.compress(payload, mtime=0))
    manifest = json.loads((tmp_path/"table_manifest.json").read_text())
    manifest["failed_fits.csv.gz"] = {"rows": 5, "sha256": file_digest(path)}
    _rewrite_json(tmp_path/"table_manifest.json", manifest)
    with pytest.raises(ValueError, match="failure log"):
        audit_run(tmp_path, design, replay=False)


def test_independent_summary_arithmetic_rejects_overflow():
    cell = Cell(6, .5, 0)
    design = RunDesign("overflow_validation", 333333333, 2, (cell,), "validation")
    rows = [{**cell.to_dict(), "rep": r, "estimator": name,
             "point_valid": True, "interval_valid": False,
             "theta": 1e200, "mu": 0., "sigma": .5}
            for r in range(2) for name in ("exact", "pfml", "euler")]
    with pytest.raises(ValueError, match="point-error moments overflow"):
        audit._parameter_summaries(rows, cell, design)


def test_long_path_accumulation_allowance_is_distinct_from_floor_ulp_check():
    # These are validation-source diagnostics, not headline paths or modified
    # artifacts. The edited copies model bounded arithmetic discrepancies and
    # confirm that the whole-sum allowance is not reused for individual gaps.
    for cv in (.25, 3.):
        cell = Cell(4000, .5, cv)
        design = RunDesign("long_path_arithmetic_validation", 333333333, 1, (cell,), "validation")
        _, diagnostics = audit.OuSampleSource(design.seed).sample(cell.to_dict(), 0)
        row = {**cell.to_dict(), **diagnostics, "rep": 0, "source_valid": True,
               "source_error": "", "pipeline_error": ""}
        audit._source(row, cell, design)
        if cv < cell.floor_cv_threshold:
            rounded = dict(row)
            # 128 span ulps exceed the previous inadequate 32-ulp rule, while
            # lying well inside the K-dependent bound for 3999 additions.
            rounded["raw_mean_gap"] = row["realized_mean_gap"]-128*np.spacing(row["realized_span"])/(cell.n-1)
            audit._source(rounded, cell, design)
            rounded["raw_mean_gap"] = row["realized_mean_gap"]*(1-1e-8)
            with pytest.raises(ValueError, match="literal flooring"):
                audit._source(rounded, cell, design)
        else:
            rounded = dict(row)
            rounded["min_gap"] = row["floor_value"]-128*np.spacing(row["realized_span"])
            assert rounded["min_gap"] > 0
            with pytest.raises(ValueError, match="floor beyond timestamp roundoff"):
                audit._source(rounded, cell, design)
