from pathlib import Path

import pytest

from ou_irregular.week10.design import full_design, validation_design
from ou_irregular.week10.execution import audit, freeze, readback, run_raw


def test_full_requires_explicit_pre_execution_freeze(tmp_path):
    with pytest.raises(ValueError, match="explicit pre-execution"):
        run_raw(full_design(), tmp_path / "forbidden", progress=False)
    assert not (tmp_path / "forbidden").exists()


def test_complete_development_execution_replay_and_no_overwrite(tmp_path):
    design = validation_design()
    frozen = freeze(design)
    output = tmp_path / "development"
    metadata = run_raw(design, output, frozen_provenance=frozen, partition_size=2, progress=False)
    assert metadata["status"] == "raw_complete_pending_audit"
    assert metadata["pipeline_errors"] == 0
    metadata, records = readback(design, output)
    assert len(records) == 84
    assert records.realized_span.eq(499.5).all()
    assert records.actual_n.nunique() > 1
    assert records.actual_n.eq(records.actual_transitions + 1).all()
    original = {str(p.relative_to(output)): p.read_bytes() for p in output.rglob("*") if p.is_file()}
    cached = run_raw(design, output, frozen_provenance=frozen, resume=True, progress=False)
    assert cached == metadata
    assert original == {str(p.relative_to(output)): p.read_bytes() for p in output.rglob("*") if p.is_file()}
    with pytest.raises(FileExistsError):
        run_raw(design, output, frozen_provenance=frozen, progress=False)
    report = audit(design, output, progress=False)
    assert report["sources_reconstructed"] == 28
    assert report["endpoint_fit_replays"] == 42
    assert report["independent_pfml_mean_checks"] == 28
    assert audit(design, output, progress=False) == report
    with pytest.raises(ValueError, match="audited or completed"):
        run_raw(design, output, frozen_provenance=frozen, resume=True, progress=False)
    partition = next(output.glob("cells/*/part-*.jsonl.gz"))
    partition.write_bytes(partition.read_bytes() + b"corrupt")
    with pytest.raises(ValueError, match="partition digest"):
        readback(design, output)


def test_saved_freeze_tamper_blocks_execution(tmp_path):
    design = validation_design()
    frozen = freeze(design)
    frozen["identity_sha256"] = "0" * 64
    with pytest.raises(ValueError, match="provenance identity"):
        run_raw(design, tmp_path / "bad", frozen_provenance=frozen, progress=False)
    assert not (tmp_path / "bad").exists()


def test_nested_saved_evidence_is_rejected(tmp_path):
    archive = tmp_path / "old"
    archive.mkdir()
    (archive / "run_metadata.json").write_text("{}")
    with pytest.raises(ValueError, match="nested within saved evidence"):
        run_raw(validation_design(), archive / "new", progress=False)


def test_symlink_output_is_rejected(tmp_path):
    target = tmp_path / "target"
    target.mkdir()
    alias = tmp_path / "alias"
    alias.symlink_to(target, target_is_directory=True)
    with pytest.raises(ValueError, match="symbolic links"):
        run_raw(validation_design(), alias / "new", progress=False)
