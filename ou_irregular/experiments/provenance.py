"""Explicit, versioned dependencies for deterministic execution and resumption.

Unlike the historical Week 4 recursive source-map check, a run records the
specific files it actually depends on. Later unrelated packages may be added
without changing this record. Existing dependencies must still match exactly.
"""
from collections.abc import Mapping
from datetime import datetime, timezone
import hashlib
from importlib.metadata import PackageNotFoundError, version
import json
import os
from pathlib import Path
import platform
import subprocess

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
PREREG_COMMIT = "1b1255b0650ad7efcf700316bd1cccf704af4d2d"
PREREG_SHA256 = "d7b1a7b6a6afa6351d3115ade91100402f846eb931952b394ecec4b68dd2813a"
DEPENDENCIES = ("numpy", "scipy", "pandas", "matplotlib", "PyYAML")
CORE_FILES = (
    "requirements.txt", "requirements-lock.txt", "preregistration.md",
    "docs/plan_amendment_week3_validation.md", "ou_irregular/__init__.py",
    "ou_irregular/config.py", "ou_irregular/estimators.py", "ou_irregular/ou.py",
    "ou_irregular/spacing.py", "ou_irregular/inference/__init__.py",
    "ou_irregular/inference/adapters.py", "ou_irregular/inference/contracts.py",
    "ou_irregular/inference/coordinates.py", "ou_irregular/inference/hessian.py",
    "ou_irregular/inference/wald.py", "ou_irregular/week4/__init__.py",
    "ou_irregular/week4/config.py", "ou_irregular/week4/simulation.py",
    "ou_irregular/week4/summaries.py", "ou_irregular/experiments/__init__.py",
    "ou_irregular/experiments/design.py", "ou_irregular/experiments/summaries.py",
    "ou_irregular/experiments/comparisons.py", "ou_irregular/experiments/provenance.py",
)


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def json_safe(value):
    if isinstance(value, Mapping):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (tuple, list, np.ndarray)):
        return [json_safe(v) for v in value]
    if isinstance(value, np.generic):
        return json_safe(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def canonical(value):
    return json.dumps(json_safe(value), sort_keys=True, separators=(",", ":"), allow_nan=False)


def digest_bytes(value):
    return hashlib.sha256(value).hexdigest()


def file_digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def environment():
    # Versions fingerprint the numerical runtime, not every installed binary.
    # Platform details are recorded separately because cross-platform bitwise
    # identity is not promised. Resume checks require this same runtime record.
    optional = {}
    for name in ("nbformat", "ipython", "ipykernel", "pytest"):
        try:
            optional[name] = version(name)
        except PackageNotFoundError:
            optional[name] = None
    return {"python": platform.python_version(), "implementation": platform.python_implementation(),
            "system": platform.system(), "machine": platform.machine(),
            "thread_environment": {name: os.environ.get(name) for name in
                                   ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")},
            "dependencies": {name: version(name) for name in DEPENDENCIES},
            "other_requirements": optional}


def capture(design, extra_files=(), *, root=ROOT, require_frozen=False):
    root = Path(root)
    paths = sorted(set(CORE_FILES) | set(extra_files))
    if any(Path(p).is_absolute() or ".." in Path(p).parts for p in paths):
        raise ValueError("dependency paths must stay within the repository")
    files = {p: file_digest(root / p) for p in paths}
    if files["preregistration.md"] != PREREG_SHA256:
        raise ValueError("the published preregistration digest does not match")
    head = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip()
    dirty = []
    for path in paths:
        original = subprocess.run(["git", "show", f"{head}:{path}"], cwd=root, capture_output=True)
        if original.returncode or digest_bytes(original.stdout) != files[path]:
            dirty.append(path)
    if require_frozen and dirty:
        raise ValueError(f"headline dependencies must be committed before execution: {dirty}")
    identity = {"design": design.to_dict(), "design_sha256": design.digest,
                "preregistration_commit": PREREG_COMMIT, "source_files": files,
                "environment": environment()}
    return {**identity, "identity_sha256": digest_bytes(canonical(identity).encode()),
            "code_commit": head, "uncommitted_dependency_files": dirty,
            "captured_utc": utc_now()}


def verify(provenance, *, root=ROOT, expected_design=None, check_environment=True):
    """Reject tampered provenance, changed dependencies, runtime or run identity."""
    identity = {k: provenance[k] for k in ("design", "design_sha256", "preregistration_commit",
                                          "source_files", "environment")}
    if digest_bytes(canonical(identity).encode()) != provenance["identity_sha256"]:
        raise ValueError("provenance identity digest mismatch")
    if provenance["preregistration_commit"] != PREREG_COMMIT:
        raise ValueError("unexpected registration freeze")
    if expected_design is not None and (provenance["design_sha256"] != expected_design.digest
                                        or provenance["design"] != expected_design.to_dict()):
        raise ValueError("resume design identity mismatch")
    for path, expected in provenance["source_files"].items():
        if Path(path).is_absolute() or ".." in Path(path).parts:
            raise ValueError("unsafe dependency path")
        if file_digest(Path(root) / path) != expected:
            raise ValueError(f"dependency changed: {path}")
    if check_environment and environment() != provenance["environment"]:
        raise ValueError("numerical environment changed")
    return True
