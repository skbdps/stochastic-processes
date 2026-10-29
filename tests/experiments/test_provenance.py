from copy import deepcopy
from pathlib import Path
import subprocess

import numpy as np
import pytest

from ou_irregular.experiments.design import Cell, RunDesign
from ou_irregular.experiments import provenance as p


def design():
    return RunDesign("test", 3, 2, (Cell(10, .5, 0),), "validation")


def test_explicit_provenance_roundtrip_and_design_guard():
    d = design()
    record = p.capture(d)
    assert p.verify(record, expected_design=d)
    assert record["preregistration_commit"] == p.PREREG_COMMIT
    with pytest.raises(ValueError, match="design"):
        p.verify(record, expected_design=d.validation())


def test_tampering_and_uncommitted_freeze_rejected():
    record = p.capture(design())
    modified = deepcopy(record)
    modified["design"]["seed"] = 99
    with pytest.raises(ValueError, match="identity digest"):
        p.verify(modified)
    # Before the implementation freeze at least this newly added module is
    # untracked. Testing an explicitly missing dependency is stable afterwards.
    with pytest.raises(FileNotFoundError):
        p.capture(design(), ["no_such_dependency.py"], require_frozen=True)


def test_changed_file_and_unrelated_addition(tmp_path):
    record = p.capture(design())
    for path in record["source_files"]:
        target = tmp_path / path
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes((p.ROOT / path).read_bytes())
    (tmp_path / "unrelated.py").write_text("unrelated later package")
    assert p.verify(record, root=tmp_path)
    (tmp_path / "ou_irregular/estimators.py").write_text("changed")
    with pytest.raises(ValueError, match="dependency changed"):
        p.verify(record, root=tmp_path)


def test_json_safety_preserves_missing_and_types():
    text = p.canonical({"missing": np.nan, "array": np.array([1, 2]), "valid": np.bool_(False)})
    assert text == '{"array":[1,2],"missing":null,"valid":false}'


def test_freeze_rejects_modified_and_untracked_dependencies(tmp_path, monkeypatch):
    (tmp_path / "preregistration.md").write_bytes((p.ROOT / "preregistration.md").read_bytes())
    (tmp_path / "tracked.py").write_text("original\n")
    monkeypatch.setattr(p, "CORE_FILES", ("preregistration.md", "tracked.py"))
    commands = [["git", "init"], ["git", "add", "."],
                ["git", "-c", "user.name=Test", "-c", "user.email=test@example.invalid",
                 "commit", "-m", "fixture"]]
    for command in commands:
        subprocess.run(command, cwd=tmp_path, check=True, capture_output=True)
    assert p.capture(design(), root=tmp_path, require_frozen=True)["uncommitted_dependency_files"] == []
    (tmp_path / "tracked.py").write_text("changed\n")
    with pytest.raises(ValueError, match="must be committed"):
        p.capture(design(), root=tmp_path, require_frozen=True)
    (tmp_path / "tracked.py").write_text("original\n")
    (tmp_path / "untracked.py").write_text("new\n")
    with pytest.raises(ValueError, match="must be committed"):
        p.capture(design(), ["untracked.py"], root=tmp_path, require_frozen=True)


def test_resume_detects_changed_numerical_thread_environment(monkeypatch):
    record = p.capture(design())
    before = p.environment()["thread_environment"]["OPENBLAS_NUM_THREADS"]
    monkeypatch.setenv("OPENBLAS_NUM_THREADS", "2" if before != "2" else "1")
    with pytest.raises(ValueError, match="environment"):
        p.verify(record)
