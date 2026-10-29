from copy import deepcopy
import json

import pytest

from ou_irregular.experiments.storage import CellStore, atomic_json, run_lock
from test_summaries import fixture_records


def test_durable_prefix_and_resume_retain_exact_outcomes(tmp_path):
    design, data = fixture_records()
    cell = design.cells[0]
    store = CellStore(tmp_path, cell, design, "identity")
    store.append(data[data.rep < 2].to_dict("records"), 2)
    resumed = CellStore(tmp_path, cell, design, "identity", resume=True)
    assert resumed.next_replication == 2
    resumed.append(data[data.rep == 2].to_dict("records"), 3)
    assert resumed.verify() and resumed.state["status"] == "complete"
    rows = list(resumed.records())
    assert len(rows) == 9
    assert sum(r["point_valid"] for r in rows) == 6
    assert sum(r["interval_valid"] for r in rows) == 3
    with pytest.raises(FileExistsError):
        CellStore(tmp_path, cell, design, "identity")
    with pytest.raises(ValueError, match="identity"):
        CellStore(tmp_path, cell, design, "other", resume=True)


def test_corrupt_and_incomplete_partitions_cannot_resume(tmp_path):
    design, data = fixture_records()
    store = CellStore(tmp_path, design.cells[0], design, "identity")
    with pytest.raises(ValueError):
        store.append(data.iloc[:2].to_dict("records"), 1)
    store.append(data[data.rep == 0].to_dict("records"), 1)
    part = tmp_path / store.state["partitions"][0]["filename"]
    part.write_bytes(part.read_bytes() + b"corruption")
    with pytest.raises(ValueError, match="digest"):
        CellStore(tmp_path, design.cells[0], design, "identity", resume=True)


def test_crash_between_file_flush_and_manifest_commit_reconciles_only_equal_bytes(tmp_path):
    design, data = fixture_records()
    store = CellStore(tmp_path, design.cells[0], design, "identity")
    before = deepcopy(store.state)
    rows = data[data.rep == 0].to_dict("records")
    store.append(rows, 1)
    atomic_json(store.path, before)
    replay = CellStore(tmp_path, design.cells[0], design, "identity", resume=True)
    bad = deepcopy(rows)
    bad[0]["theta"] += 1
    with pytest.raises(FileExistsError):
        replay.append(bad, 1)
    replay.append(rows, 1)
    assert len(replay.state["recovered_orphans"]) == 1
    assert replay.next_replication == 1


def test_run_lock_rejects_concurrent_writer(tmp_path):
    with run_lock(tmp_path):
        with pytest.raises(RuntimeError, match="execution lock"):
            with run_lock(tmp_path):
                pass


def test_changed_partition_size_cannot_bypass_orphan_evidence(tmp_path):
    design, data = fixture_records()
    store = CellStore(tmp_path, design.cells[0], design, "identity")
    before = deepcopy(store.state)
    store.append(data[data.rep < 2].to_dict("records"), 2)
    atomic_json(store.path, before)
    resumed = CellStore(tmp_path, design.cells[0], design, "identity", resume=True)
    assert resumed.orphan_stop == 2
    with pytest.raises(ValueError, match="original partition range"):
        resumed.append(data[data.rep == 0].to_dict("records"), 1)
    bad = data[data.rep < 2].to_dict("records")
    bad[0]["theta"] = 8
    with pytest.raises(FileExistsError):
        resumed.append(bad, 2)
    resumed.append(data[data.rep < 2].to_dict("records"), 2)
    resumed.append(data[data.rep == 2].to_dict("records"), 3)
    extra = tmp_path / "part-000000-000001.jsonl.gz"
    extra.write_bytes((tmp_path / resumed.state["partitions"][0]["filename"]).read_bytes())
    with pytest.raises(ValueError, match="unreferenced"):
        resumed.verify()


def test_missing_manifest_read_resume_does_not_create_evidence(tmp_path):
    design, _ = fixture_records()
    with pytest.raises(FileNotFoundError):
        CellStore(tmp_path, design.cells[0], design, "identity", resume=True)
    assert not list(tmp_path.iterdir())
