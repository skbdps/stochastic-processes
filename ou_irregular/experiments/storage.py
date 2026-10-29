"""Durable immutable per-cell partitions and fail-closed resume identities.

Only a fully flushed partition is committed to a cell manifest. A crash after
the file is flushed but before the manifest advances may leave an orphan. A
deterministic replay may adopt that file only when its exact bytes agree; it
never overwrites different evidence or replaces a failed statistical outcome.
"""
from contextlib import contextmanager
from dataclasses import replace
import fcntl
import gzip
import json
import os
from pathlib import Path
import re
import tempfile

import numpy as np
import pandas as pd

from .design import ESTIMATORS
from .provenance import canonical, digest_bytes, file_digest, utc_now
from .summaries import validate_records


def atomic_bytes(path, payload, *, immutable=False):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if immutable and path.exists():
        if path.read_bytes() != payload:
            raise FileExistsError(f"refusing to replace different immutable evidence: {path}")
        return False
    descriptor, temporary = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(descriptor, "wb") as stream:
            stream.write(payload)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
        directory = os.open(path.parent, os.O_DIRECTORY)
        try:
            os.fsync(directory)
        finally:
            os.close(directory)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)
    return True


def atomic_json(path, value, *, immutable=False):
    return atomic_bytes(path, (canonical(value) + "\n").encode(), immutable=immutable)


@contextmanager
def run_lock(output):
    """A second orchestrator cannot concurrently mutate the same run manifest."""
    path = Path(output) / ".execution.lock"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as stream:
        try:
            fcntl.flock(stream, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError("another process holds the experiment execution lock") from exc
        try:
            yield
        finally:
            fcntl.flock(stream, fcntl.LOCK_UN)


def read_partition(path):
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        return [json.loads(line) for line in stream if line.strip()]


class CellStore:
    """A contiguous prefix of scheduled replications, three rows per source."""

    def __init__(self, directory, cell, design, identity, *, resume=False):
        self.directory, self.cell, self.design = Path(directory), cell, design
        self.identity = identity
        self.orphan_stop = None
        self.path = self.directory / "cell_manifest.json"
        if self.path.exists():
            if not resume:
                raise FileExistsError("cell evidence already exists; explicit resume required")
            self.state = json.loads(self.path.read_text())
            if (self.state["run_identity"] != identity or self.state["cell"] != cell.to_dict()
                    or self.state["replications"] != design.replications):
                raise ValueError("cell resume identity mismatch")
            self.verify()
        else:
            if resume:
                raise FileNotFoundError("cannot resume a cell without its manifest")
            if self.directory.exists() and any(self.directory.iterdir()):
                raise FileExistsError("unidentified nonempty cell directory")
            self.directory.mkdir(parents=True, exist_ok=True)
            self.state = {"schema_version": 1, "run_identity": identity, "cell": cell.to_dict(),
                          "replications": design.replications, "status": "running_incomplete",
                          "completed_replications": 0, "partitions": [], "recovered_orphans": []}
            atomic_json(self.path, self.state)

    @property
    def next_replication(self):
        return self.state["completed_replications"]

    def _validate(self, records, start, stop):
        if not 0 <= start < stop <= self.design.replications:
            raise ValueError("partition range outside planned replications")
        data = pd.DataFrame(records)
        if "rep" not in data or not data.rep.map(lambda v: isinstance(v, (int, np.integer))
                                                and not isinstance(v, (bool, np.bool_))).all():
            raise ValueError("partition replication IDs must be integers")
        data["rep"] = data.rep - start
        partition_design = replace(self.design, cells=(self.cell,), replications=stop-start)
        validate_records(data, partition_design)

    def append(self, records, stop):
        start = self.next_replication
        if self.orphan_stop is not None and stop != self.orphan_stop:
            raise ValueError("orphan recovery must preserve the original partition range")
        self._validate(records, start, stop)
        filename = f"part-{start:06d}-{stop:06d}.jsonl.gz"
        ordered = sorted(records, key=lambda row: (row["rep"], ESTIMATORS.index(row["estimator"])))
        raw = ("\n".join(canonical(row) for row in ordered) + "\n").encode()
        payload = gzip.compress(raw, compresslevel=6, mtime=0)
        created = atomic_bytes(self.directory / filename, payload, immutable=True)
        if not created:
            self.state["recovered_orphans"].append({"filename": filename, "verified_utc": utc_now(),
                                                     "method": "exact_byte_equal_deterministic_replay"})
        self.state["partitions"].append({"filename": filename, "start": start, "stop": stop,
                                         "rows": len(ordered), "sha256": digest_bytes(payload)})
        self.state["completed_replications"] = stop
        self.state["status"] = "complete" if stop == self.design.replications else "running_incomplete"
        atomic_json(self.path, self.state)
        self.orphan_stop = None

    def verify(self):
        expected_start = 0
        for part in self.state["partitions"]:
            filename = f"part-{part['start']:06d}-{part['stop']:06d}.jsonl.gz"
            if part["filename"] != filename or part["start"] != expected_start:
                raise ValueError("noncontiguous or unsafe partition manifest")
            path = self.directory / filename
            if file_digest(path) != part["sha256"]:
                raise ValueError("partition digest mismatch")
            records = read_partition(path)
            if len(records) != part["rows"]:
                raise ValueError("partition row count mismatch")
            self._validate(records, part["start"], part["stop"])
            expected_start = part["stop"]
        if expected_start != self.next_replication:
            raise ValueError("cell manifest completion count mismatch")
        complete = expected_start == self.design.replications
        if (self.state["status"] == "complete") != complete:
            raise ValueError("cell completion status mismatch")
        known = {p["filename"] for p in self.state["partitions"]}
        extras = sorted(p for p in self.directory.glob("*.jsonl.gz") if p.name not in known)
        self.orphan_stop = None
        if extras:
            if complete or len(extras) != 1:
                raise ValueError("unexpected unreferenced partition evidence")
            match = re.fullmatch(r"part-(\d{6})-(\d{6})\.jsonl\.gz", extras[0].name)
            if match is None or int(match[1]) != expected_start:
                raise ValueError("orphan partition overlaps or is outside the committed prefix")
            stop = int(match[2])
            self._validate(read_partition(extras[0]), expected_start, stop)
            # Its bytes still must agree with deterministic replay before it
            # can become committed evidence. A new partition size cannot skip it.
            self.orphan_stop = stop
        return True

    def records(self):
        for part in self.state["partitions"]:
            yield from read_partition(self.directory / part["filename"])
