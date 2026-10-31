"""Read-only catalog guarantees and real saved-E1 accounting, without simulation."""
import hashlib
import json
from pathlib import Path

import pytest

from ou_irregular.week9 import inputs


def digest(data):
    return hashlib.sha256(data).hexdigest()


@pytest.fixture(scope="module")
def published():
    return inputs.load_sources()


@pytest.fixture
def tiny(tmp_path, monkeypatch):
    root = tmp_path / "repo"
    source = root / "artifacts" / "study"
    source.mkdir(parents=True)
    payload = b"saved numerical evidence\n"
    (source / "data.csv").write_bytes(payload)
    relative = "artifacts/study/data.csv"
    entry = {"sha256": digest(payload), "git_blob": hashlib.sha1(
        b"blob " + str(len(payload)).encode() + b"\0" + payload).hexdigest(),
        "size_bytes": len(payload), "git_mode": "100644"}
    config = {"studies": {"e1": {"path": "artifacts/study", "files": {relative: entry}}},
              "additional_sources": {}, "reference_files": {}}
    config_path = root / "catalog.json"
    config_path.write_text(json.dumps(config))
    monkeypatch.setattr(inputs, "CONFIG_SHA256", digest(config_path.read_bytes()))
    bundle = inputs.SourceBundle(root, config, config_path, {relative: entry}, inputs._canonical(config))
    return bundle, source / "data.csv"


def test_published_catalog_has_all_six_studies_and_baseline_git_objects(published):
    assert published.config["baseline_commit"] == "b6edb85d1e7a21ed19fb54b5c71213ed04703309"
    assert set(published.config["studies"]) == {
        "e1", "e2", "e3_selected", "e3_followup", "e3_floor", "e4"}
    snapshot = published.snapshot()
    assert len(snapshot) == 831
    assert all(len(sha) == 64 for sha in snapshot.values())
    assert published.config["studies"]["e3_floor"]["expected"]["independent_coupling_groups"] == 8000
    assert published.config["studies"]["e3_selected"]["registration_status"] == "exploratory_selected"
    assert "docs/week8_population_targets.md" in snapshot


def test_published_summary_complete_keys_and_counts(published):
    expected = {"e1": (63, 42000, 42000), "e2": (540, 359915, 359914),
                "e3_followup": (18, 59579, 59556), "e3_floor": (108, 72000, 72000),
                "e4": (63, 42000, 42000)}
    for study, (rows, points, intervals) in expected.items():
        summary = published.table(study, "parameter_summary.csv")
        assert len(summary) == rows
        assert not summary.duplicated(["cell_id", "estimator", "parameter"]).any()
        theta = summary[summary.parameter == "theta"]
        assert int(theta.n_point_valid.sum()) == points
        assert int(theta.n_interval_valid.sum()) == intervals
        assert (summary[summary.parameter == "mu"].practical_classification
                == "not_applicable_zero_truth").all()


def test_f1_contains_all_saved_fit_keys_and_partition_bindings(published):
    records = published.f1_records()
    keys = {(r["cell_id"], r["rep"], r["estimator"]) for r in records}
    assert len(keys) == len(records) == 42000
    assert len({r["cell_id"] for r in records}) == 7
    assert {r["cv"] for r in records} == {0, .25, .5, 1, 1.5, 2, 3}
    assert {r["source_seed"] for r in records} == {20261001}
    assert all(r["point_valid"] and r["interval_valid"] for r in records)
    assert len({r["_source_path"] for r in records}) == 56
    assert all(published._entries[r["_source_path"]]["sha256"] == r["_source_sha256"] for r in records)


@pytest.mark.parametrize("name", ["../file", "/file", "a/../file", "a\\file", "a//file", "./file", ".git/config", ""])
def test_unsafe_source_names_are_rejected(name):
    with pytest.raises(ValueError):
        inputs._relative(name)


def test_source_bytes_are_checked_again_after_loading(tiny):
    bundle, path = tiny
    assert bundle.snapshot() == {"artifacts/study/data.csv": digest(path.read_bytes())}
    path.write_bytes(b"changed numerical evidence\n")
    with pytest.raises(ValueError, match="source digest mismatch"):
        bundle.read_bytes("artifacts/study/data.csv")
    with pytest.raises(ValueError, match="source digest mismatch"):
        bundle.snapshot()


def test_missing_and_unexpected_archive_files_are_rejected(tiny):
    bundle, path = tiny
    (path.parent / "extra.csv").write_text("orphan")
    with pytest.raises(ValueError, match="membership changed"):
        bundle.snapshot()
    (path.parent / "extra.csv").unlink()
    path.unlink()
    with pytest.raises(ValueError, match="membership changed"):
        bundle.snapshot()


def test_source_symlink_is_rejected_even_with_matching_content(tiny):
    bundle, path = tiny
    payload = path.read_bytes()
    path.unlink()
    target = bundle.root / "external.csv"
    target.write_bytes(payload)
    path.symlink_to(target)
    with pytest.raises(ValueError, match="symlink"):
        bundle.snapshot()
    with pytest.raises(ValueError, match="symlink"):
        bundle.read_bytes("artifacts/study/data.csv")


def test_in_memory_and_on_disk_catalog_mutations_are_rejected(tiny):
    bundle, _ = tiny
    bundle.config["studies"]["e1"]["path"] = "somewhere/else"
    with pytest.raises(ValueError, match="in-memory"):
        bundle.snapshot()


def test_on_disk_catalog_mutation_is_rejected(tiny):
    bundle, _ = tiny
    bundle.config_path.write_text("{}")
    with pytest.raises(ValueError, match="catalog digest"):
        bundle.snapshot()


def test_full_loader_rejects_changed_or_missing_catalog_keys(tmp_path):
    source = Path(inputs.ROOT) / inputs.DEFAULT_CONFIG
    config = json.loads(source.read_text())
    config["studies"].pop("e2")
    altered = tmp_path / "altered.json"
    altered.write_text(json.dumps(config))
    with pytest.raises(ValueError, match="catalog digest"):
        inputs.load_sources(altered)


def test_table_access_cannot_escape_catalog(published):
    for filename in ["../parameter_summary.csv", "/tmp/parameter_summary.csv", "cells/x.csv", "audit.json"]:
        with pytest.raises(ValueError):
            published.table("e1", filename)


def test_guard_output_requires_fresh_destination_and_no_source_overlap(tiny, tmp_path):
    bundle, path = tiny
    safe = tmp_path / "new-publication"
    assert inputs.guard_output(bundle, safe) == safe.resolve()
    assert not safe.exists()
    with pytest.raises(ValueError, match="fresh"):
        inputs.guard_output(bundle, bundle.root)
    with pytest.raises(ValueError, match="overlaps"):
        inputs.guard_output(bundle, path.parent / "new-publication")
    with pytest.raises(ValueError, match="repository metadata"):
        inputs.guard_output(bundle, bundle.root / ".git" / "new-publication")
    with pytest.raises(ValueError, match="parent traversal"):
        inputs.guard_output(bundle, tmp_path / "other" / ".." / "new-publication")


def test_guard_output_resolves_symlink_ancestor_before_overlap_check(tiny, tmp_path):
    bundle, path = tiny
    linked = tmp_path / "linked-source"
    linked.symlink_to(path.parent, target_is_directory=True)
    with pytest.raises(ValueError, match="overlaps"):
        inputs.guard_output(bundle, linked / "new-publication")
