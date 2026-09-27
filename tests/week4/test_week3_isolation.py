"""Week 4 must preserve the archived Week 3 program and its historical replay.

These checks use the committed provenance manifest, not Git working-tree state.
An installed distribution without the research artifacts skips archival checks;
the ordinary API unit tests remain available in that environment. Exact replay
requires the archived dependency versions, as the Week 3 public API specifies.
"""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import ou_irregular
from ou_irregular import runner
from ou_irregular.config import load_config
# Import the new extension before exercising historical APIs. An extension
# must not patch the old runner, estimator, or likelihood at import time.
from ou_irregular.inference.wald import WaldInference  # noqa: F401
from ou_irregular.inference.bootstrap import ParametricBootstrap  # noqa: F401


ROOT = Path(__file__).resolve().parents[2]
ARCHIVE = ROOT / "artifacts" / "week3_validation_rerun"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


@pytest.fixture(scope="module")
def archive():
    if not (ARCHIVE / "run_metadata.json").is_file():
        pytest.skip("archived Week 3 research artifacts are not installed")
    return ARCHIVE, json.loads((ARCHIVE / "run_metadata.json").read_text())


def test_week3_module_bytes_match_the_original_validation_manifest(archive):
    _, metadata = archive
    package = Path(ou_irregular.__file__).resolve().parent
    assert set(metadata["source_sha256"]) == {
        "__init__.py", "config.py", "estimators.py", "metrics.py", "ou.py",
        "plots.py", "runner.py", "spacing.py",
    }
    for name, expected in metadata["source_sha256"].items():
        assert digest(package / name) == expected, f"Week 3 source changed: {name}"


def test_archived_week3_outputs_retain_their_recorded_bytes(archive):
    path, metadata = archive
    for name, expected in metadata["output_sha256"].items():
        assert digest(path / name) == expected, f"Week 3 artifact changed: {name}"
    # Week 3 canonicalizes the resolved configuration before hashing, so this
    # also checks scientific parameters without relying on YAML formatting.
    config = load_config(path / "config.yaml")
    assert runner._provenance(config)["config_sha256"] == metadata["config_sha256"]


def test_original_week3_smoke_outputs_remain_unchanged():
    # The original smoke used an earlier Week 3 solver, so it deliberately has
    # different source hashes. Its results remain historical evidence and
    # must not be regenerated in place by the Week 4 extension.
    path = ROOT / "artifacts" / "week3_smoke"
    if not (path / "run_metadata.json").is_file():
        pytest.skip("original Week 3 smoke artifacts are not installed")
    metadata = json.loads((path / "run_metadata.json").read_text())
    for name, expected in metadata["output_sha256"].items():
        assert digest(path / name) == expected, f"original Week 3 artifact changed: {name}"


def test_historical_cv2_replication_replays_exactly_without_touching_archive(archive, tmp_path):
    path, metadata = archive
    config = load_config(path / "config.yaml")
    current = runner._provenance(config)
    differences = runner.provenance_differences(metadata, current)
    # Changed source must always fail; a different Python/dependency version
    # may legitimately run the other tests, but cannot claim exact replay.
    assert current["source_sha256"] == metadata["source_sha256"]
    if differences:
        pytest.skip("strict historical replay requires the archived Python/dependency versions")
    cell = next(cell for cell in config.cells()
                if cell["cv"] == 2.0 and cell["theta_mean_gap"] == 0.5)
    before = {str(file.relative_to(path)): digest(file)
              for file in path.rglob("*") if file.is_file()}
    saved = pd.read_csv(path / "replications.csv", float_precision="round_trip")
    saved = saved[(saved.cell_id == cell["cell_id"]) & (saved.rep == 0)].sort_values("estimator")
    replayed = runner.replay(path, cell["cell_id"], 0,
                             output_dir=tmp_path / "historical_replay").sort_values("estimator")
    assert len(saved) == len(replayed) == 3
    assert list(replayed.estimator) == list(saved.estimator) == ["euler", "exact", "pfml"]
    # Bit-identical input digests AND fitted values protect the actual result,
    # beyond merely checking unchanged function names or random seed labels.
    for column in ("times_sha256", "observations_sha256", "fit_status", "valid_for_point_summary"):
        assert replayed[column].tolist() == saved[column].tolist()
    np.testing.assert_array_equal(replayed[["theta", "mu", "sigma", "nll"]].to_numpy(),
                                  saved[["theta", "mu", "sigma", "nll"]].to_numpy())
    evidence = json.loads((tmp_path / "historical_replay" / "replay_metadata.json").read_text())
    assert evidence["kind"] == "strict_replay"
    assert evidence["input_identity"] is True
    assert not evidence["provenance_differences"]
    after = {str(file.relative_to(path)): digest(file)
             for file in path.rglob("*") if file.is_file()}
    assert after == before
