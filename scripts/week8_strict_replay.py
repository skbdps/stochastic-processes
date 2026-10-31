"""Run unchanged four-case historical replay on an AVX-512-capable host.

Hardware is screened before the test, never by weakening its exact comparisons.
"""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import numpy as np


def main():
    from week8_release_regression import IDENTITY_PLUGIN, RUNTIME_CAPTURE
    out = Path(sys.argv[1]); out.mkdir(exist_ok=False)
    cpu = np._core._multiarray_umath.__cpu_features__
    (out / "cpu.txt").write_text(subprocess.check_output(["lscpu"], text=True))
    runtime = subprocess.check_output([sys.executable, "-c", RUNTIME_CAPTURE], text=True)
    (out / "runtime.json").write_text(runtime)
    (out / "dependencies.txt").write_text(subprocess.check_output([sys.executable, "-m", "pip", "freeze"], text=True))
    if not cpu.get("AVX512F"):
        (out / "hardware_unavailable.json").write_text(json.dumps({"status": "compatible hardware unavailable", "AVX512F": False}) + "\n")
        return 2
    assert subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip() == "0ff57919f5d8345f52d8689732ff94448ca3ca1d"
    assert not subprocess.check_output(["git", "status", "--porcelain"], text=True).strip()
    (out / "week8_release_identity_plugin.py").write_text(IDENTITY_PLUGIN)
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join([str(Path.cwd()), str(out)])
    env["OU_RELEASE_COLLECTION_OUTPUT"] = str(out / "week3-isolation-nodeids.json")
    env["OU_RELEASE_REPORT_OUTPUT"] = str(out / "week3-isolation-reports.json")
    command = [sys.executable, "-m", "pytest", "-q", "-ra", "tests/week4/test_week3_isolation.py",
               "-p", "week8_release_identity_plugin", "--junitxml="+str(out / "week3-isolation-junit.xml"),
               "--basetemp="+str(out.parent / (out.name+"-pytest"))]
    (out / "command.json").write_text(json.dumps(command, indent=2)+"\n")
    with (out / "pytest-week3-isolation.txt").open("w") as f:
        result = subprocess.run(command, env=env, stdout=f, stderr=subprocess.STDOUT)
    (out / "exit.json").write_text(json.dumps({"returncode": result.returncode})+"\n")
    assert not subprocess.check_output(["git", "status", "--porcelain"], text=True).strip()
    manifest = {str(p.relative_to(out)): {"bytes": p.stat().st_size, "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}
                for p in out.rglob("*") if p.is_file()}
    (out / "evidence_sha256.json").write_text(json.dumps(manifest, indent=2)+"\n")
    print((out / "pytest-week3-isolation.txt").read_text(), flush=True)
    return result.returncode


if __name__ == "__main__":
    raise SystemExit(main())
