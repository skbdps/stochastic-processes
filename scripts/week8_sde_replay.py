"""Run the unchanged four-case Week 3 replay in one Intel SDE process.

Invoke with ``-p -follow_execv -p 0 -icl -- python this_script``. Each ``-p``
forwards the next token to Pin rather than treating it as an SDE tool option.
The effective NumPy runtime is captured in this same process before pytest.main;
no scientific function, test, tolerance, or NumPy CPU flag is changed.
"""
import argparse
from collections import Counter
import contextlib
from datetime import datetime, timezone
import hashlib
import io
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import traceback
import types
import xml.etree.ElementTree as ET

FREEZE = "0ff57919f5d8345f52d8689732ff94448ca3ca1d"
ARCHIVE_SHA256 = "94e97d623fec54385686e1e7ba65ebc9941748c05ee451423948334892bf2b50"
TEST_PATH = "tests/week4/test_week3_isolation.py"
EXPECTED = {TEST_PATH + "::" + name for name in (
    "test_week3_module_bytes_match_the_original_validation_manifest",
    "test_archived_week3_outputs_retain_their_recorded_bytes",
    "test_original_week3_smoke_outputs_remain_unchanged",
    "test_historical_cv2_replication_replays_exactly_without_touching_archive",
)}


def write(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def git(*args):
    return subprocess.check_output(["git", *args], text=True).strip()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("output", "emulator", "emulator-archive", "version-file", "help-file", "command-file"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    root = Path.cwd().resolve()
    out = args.output.resolve()
    if out.exists() or out.is_relative_to(root) or root.is_relative_to(out):
        parser.error("output must be a fresh directory outside the scientific checkout")
    out.mkdir(parents=True)
    sys.path.insert(0, str(root))
    record = {"started_utc": datetime.now(timezone.utc).isoformat(), "status": "failed",
              "scientific_freeze": FREEZE, "control_commit": os.environ.get("CONTROL_COMMIT"),
              "workflow_run_id": os.environ.get("GITHUB_RUN_ID"),
              "execution_mode": "Intel SDE instruction emulation (-p -follow_execv -p 0 -icl), same-process pytest; Pin exec-follow configured off",
              "emulator_download_url": "https://downloadmirror.intel.com/924984/sde-external-10.13.1-2026-07-28-lin.tar.xz",
              "emulator_documentation_url": "https://www.intel.com/content/www/us/en/developer/articles/tool/software-development-emulator.html",
              "pin_argument_forwarding_source": "https://github.com/intel/pinplay-tools/blob/main/PinPoints/Test/sde-run.pinpoints.single-threaded.sh",
              "test_path": TEST_PATH, "new_full_experiment_executed": False,
              "scope": "Four unchanged Week 3 archival checks, including one source/three-fit exact replay. No E4 study execution."}
    try:
        assert sha(args.emulator_archive) == ARCHIVE_SHA256
        command = json.loads(args.command_file.read_text())
        assert command[:7] == [str(args.emulator.resolve()), "-p", "-follow_execv", "-p", "0", "-icl", "--"]
        assert Path(command[7]).resolve() == Path(sys.executable).resolve()
        assert Path(command[8]).resolve() == Path(__file__).resolve()
        assert command[9:] == sys.argv[1:]
        record.update(emulator_archive_sha256=sha(args.emulator_archive),
                      emulator_executable_sha256=sha(args.emulator),
                      emulator_version=args.version_file.read_text(), emulator_help_sha256=sha(args.help_file), command=command,
                      emulator_help_scope="Retained SDE Pin-tool help is informational; it is not an inventory of generic Pin launcher knobs.",
                      test_sha256=sha(root / TEST_PATH), script_sha256=sha(__file__))
        write(out / "source_provenance.json", record)
        cpu_controls = {n: os.environ.get(n) for n in
                        ("NPY_DISABLE_CPU_FEATURES", "OPENBLAS_CORETYPE", "MKL_ENABLE_INSTRUCTIONS")}
        assert all(v is None for v in cpu_controls.values()), cpu_controls
        assert os.environ.get("OPENBLAS_NUM_THREADS") == os.environ.get("OMP_NUM_THREADS") == "1"
        assert platform.python_version() == "3.12.14"
        import numpy as np
        # The original test imports this unmodified module, whose documented
        # import sets its three thread variables to 1. Capture that actual state.
        from ou_irregular import runner
        from ou_irregular.config import load_config
        from ou_irregular.experiments.provenance import environment
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            np.show_runtime()
        features = dict(np._core._multiarray_umath.__cpu_features__)
        runtime = {"pid": os.getpid(), "environment": environment(), "platform": platform.platform(),
                   "numpy_runtime": buffer.getvalue(), "numpy_build_config": str(np.__config__.CONFIG),
                   "cpu_features": features,
                   "cpu_dispatch": list(np._core._multiarray_umath.__cpu_dispatch__),
                   "cpu_control_environment": cpu_controls,
                   "execution_mode": record["execution_mode"]}
        write(out / "same_process_runtime.json", runtime)
        (out / "numpy_runtime.txt").write_text(buffer.getvalue())
        # Operational smoke replaces an invalid help-text check. It exercises
        # the exact auxiliary exec that failed before any historical test ran.
        # Numerical tests and runtime capture stay in this instrumented parent.
        probe_command = ["git", "rev-parse", "HEAD"]
        probe = subprocess.run(probe_command, text=True, capture_output=True, check=False)
        smoke = {"parent_pid": os.getpid(), "git_command": probe_command,
                 "git_returncode": probe.returncode, "git_stdout": probe.stdout,
                 "git_stderr": probe.stderr, "effective_AVX512F": features.get("AVX512F"),
                 "effective_AVX512_ICL": features.get("AVX512_ICL"),
                 "scope": "Auxiliary Git exec succeeds and parent NumPy exposes required emulated dispatch. This is a runtime smoke, not historical numerical replay certification."}
        write(out / "runtime_smoke.json", smoke)
        record["runtime_smoke"] = smoke
        assert probe.returncode == 0, smoke
        assert features.get("AVX512F") is True and features.get("AVX512_ICL") is True, features
        record["candidate_commit"] = probe.stdout.strip()
        record["git_status_before"] = git("status", "--porcelain", "--untracked-files=all")
        assert record["candidate_commit"] == FREEZE
        assert not record["git_status_before"]
        assert git("rev-parse", "--show-toplevel") == str(root)
        saved = json.loads((root / "artifacts/week3_validation_rerun/run_metadata.json").read_text())
        config = load_config(root / "artifacts/week3_validation_rerun/config.yaml")
        differences = runner.provenance_differences(saved, runner._provenance(config))
        record["historical_provenance_differences"] = differences
        assert not differences, differences
        record["effective_AVX512F_and_ICL_verified"] = True
        write(out / "source_provenance.json", record)
        # Reuse exactly the existing read-only outcome observer, rather than
        # replacing any test behavior or inferring outcomes from console text.
        from week8_release_regression import IDENTITY_PLUGIN
        observer = types.ModuleType("week8_sde_identity_observer")
        exec(IDENTITY_PLUGIN, observer.__dict__)
        os.environ["OU_RELEASE_COLLECTION_OUTPUT"] = str(out / "week3-isolation-nodeids.json")
        os.environ["OU_RELEASE_REPORT_OUTPUT"] = str(out / "week3-isolation-reports.json")
        import pytest
        scratch = out.with_name(out.name + "-pytest")
        assert not scratch.exists() and not scratch.is_symlink(), "pytest scratch must be fresh"
        pytest_args = ["-q", "-ra", TEST_PATH,
                       "--basetemp=" + str(scratch),
                       "--junitxml=" + str(out / "week3-isolation-junit.xml")]
        record["pytest_arguments"] = pytest_args
        record["pytest_pid"] = os.getpid()
        with (out / "pytest-week3-isolation.txt").open("w") as stream, contextlib.redirect_stdout(stream), contextlib.redirect_stderr(stream):
            code = int(pytest.main(pytest_args, plugins=[observer]))
        print((out / "pytest-week3-isolation.txt").read_text(), flush=True)
        record["pytest_returncode"] = code
        nodeids = json.loads((out / "week3-isolation-nodeids.json").read_text())
        reports = json.loads((out / "week3-isolation-reports.json").read_text())
        calls = Counter(r["nodeid"] for r in reports["reports"] if r["when"] == "call" and r["outcome"] == "passed")
        suites = list(ET.parse(out / "week3-isolation-junit.xml").getroot().iter("testsuite"))
        counts = {k: sum(int(s.attrib.get(k, 0)) for s in suites)
                  for k in ("tests", "errors", "failures", "skipped")}
        record["junit_counts"] = counts
        assert len(nodeids) == 4 and set(nodeids) == EXPECTED
        assert calls == Counter({n: 1 for n in EXPECTED})
        assert reports["exitstatus"] == code == 0
        assert counts == {"tests": 4, "errors": 0, "failures": 0, "skipped": 0}
        assert all(r["outcome"] == "passed" for r in reports["reports"])
        record["status"] = "passed"
    except BaseException as exc:
        record["exception"] = repr(exc)
        (out / "traceback.txt").write_text(traceback.format_exc())
        print(traceback.format_exc(), file=sys.stderr, flush=True)
    finally:
        try:
            record["git_status_after"] = git("status", "--porcelain", "--untracked-files=all")
        except BaseException as exc:
            # A diagnostic transport failure must not erase its own record.
            record["git_status_after_error"] = repr(exc)
            record["git_status_after"] = "unavailable"
        record["clean_before_and_after"] = not record.get("git_status_before", "missing") and not record["git_status_after"]
        if not record["clean_before_and_after"]:
            record["status"] = "failed"
        record["finished_utc"] = datetime.now(timezone.utc).isoformat()
        write(out / "source_provenance.json", record)
        write(out / "evidence_sha256.json", {str(p.relative_to(out)): {"sha256": sha(p), "bytes": p.stat().st_size}
              for p in sorted(out.rglob("*")) if p.is_file() and p.name != "evidence_sha256.json"})
    print(json.dumps({"status": record["status"], "junit_counts": record.get("junit_counts"),
                      "evidence": str(out), "execution_mode": record["execution_mode"]}, indent=2), flush=True)
    return 0 if record["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
