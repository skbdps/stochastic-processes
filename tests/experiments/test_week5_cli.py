import subprocess
import sys

from ou_irregular.experiments.provenance import ROOT


def test_cli_rejects_headline_replication_change_before_execution(tmp_path):
    result = subprocess.run([sys.executable, "-m", "ou_irregular.week5.runner", "run",
                             "--output", str(tmp_path / "forbidden"), "--replications", "1"],
                            cwd=ROOT, capture_output=True, text=True)
    assert result.returncode != 0
    assert "headline R is frozen" in result.stderr
    assert not (tmp_path / "forbidden").exists()


def test_cli_validation_rejects_zero_replications(tmp_path):
    result = subprocess.run([sys.executable, "-m", "ou_irregular.week5.runner", "run",
                             "--output", str(tmp_path / "forbidden"), "--validation", "--replications", "0"],
                            cwd=ROOT, capture_output=True, text=True)
    assert result.returncode != 0 and "integer >= 1" in result.stderr
