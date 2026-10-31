"""Run the unchanged release suite in five disjoint numerical-runtime groups.

Run from a clean, full-history scientific checkout:
    OU_LEGACY_PYTHON=/path/to/legacy/bin/python \
      python scripts/week8_release_regression.py /fresh/external/evidence

The caller installs the two recorded environments. This script installs nothing,
changes no source/test/archive bytes, and never starts a full E4 experiment.
Original logs, JUnit, complete node identities and failures are retained before
returning nonzero. Test fixtures include their prescribed small validation runs.
"""
import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import traceback
import xml.etree.ElementTree as ET


FREEZE = "0ff57919f5d8345f52d8689732ff94448ca3ca1d"
EXPECTED = {"collected": 2061, "passed": 2055, "skipped": 6}
HISTORICAL = ["tests/week7/test_plots_integrity.py",
              "tests/week7/test_plots_reference.py",
              "tests/week7/test_selection.py"]
GROUP_PATHS = {
    "core": ["tests", "--ignore=tests/week7", "--ignore=tests/week8",
             "--ignore=tests/week4/test_week3_isolation.py"],
    "week7": ["tests/week7", *["--ignore=" + name for name in HISTORICAL]],
    "week8": ["tests/week8"],
    "historical": HISTORICAL,
    "week3-isolation": ["tests/week4/test_week3_isolation.py"],
}

# A read-only pytest observer. Full nodeids avoid ambiguity in class/name or
# parameter strings; JUnit is independently checked for all execution counts.
IDENTITY_PLUGIN = '''"""Record test identities/outcomes without changing pytest behavior."""
import json, os
from pathlib import Path
_reports = []
def pytest_collection_finish(session):
    Path(os.environ["OU_RELEASE_COLLECTION_OUTPUT"]).write_text(
        json.dumps([item.nodeid for item in session.items], indent=2) + "\\n")
def pytest_runtest_logreport(report):
    _reports.append({"nodeid": report.nodeid, "when": report.when,
                     "outcome": report.outcome})
def pytest_sessionfinish(session, exitstatus):
    Path(os.environ["OU_RELEASE_REPORT_OUTPUT"]).write_text(
        json.dumps({"exitstatus": int(exitstatus), "reports": _reports},
                   indent=2) + "\\n")
'''

# Preimport the original Week3 runner exactly as the successful runner gate did.
# Its historical import sets MKL_NUM_THREADS=1. Restore that operational side
# effect before tests; scientific functions and assertions remain untouched.
# The real file and main guard keep multiprocessing spawn from reentering pytest.
PYTEST_ENTRY = '''"""Spawn-safe release test entry; no scientific/test monkeypatches."""
import json, os
from pathlib import Path
import sys
sys.path.insert(0, str(Path.cwd()))
def main():
    import ou_irregular.runner
    import pytest
    os.environ.pop("MKL_NUM_THREADS", None)
    os.environ.pop("VECLIB_MAXIMUM_THREADS", None)
    from ou_irregular.experiments.provenance import environment
    expected = json.loads(Path(os.environ["OU_RELEASE_RUNTIME"]).read_text())
    if environment() != expected:
        raise RuntimeError("test-entry runtime differs from recorded current runtime")
    return pytest.main(sys.argv[1:])
if __name__ == "__main__":
    raise SystemExit(main())
'''

RUNTIME_CAPTURE = '''import contextlib, io, json, os, platform, sys
import numpy as np
from ou_irregular.experiments.provenance import environment
buffer = io.StringIO()
with contextlib.redirect_stdout(buffer):
    np.show_runtime()
print(json.dumps({"environment": environment(), "executable": sys.executable,
    "platform": platform.platform(), "numpy_runtime": buffer.getvalue(),
    "cpu_features": np._core._multiarray_umath.__cpu_features__,
    "cpu_dispatch": list(np._core._multiarray_umath.__cpu_dispatch__),
    "cpu_control_environment": {name: os.environ.get(name) for name in
        ["NPY_DISABLE_CPU_FEATURES", "OPENBLAS_CORETYPE", "MKL_ENABLE_INSTRUCTIONS"]},
    "numpy_build_config": str(np.__config__.CONFIG)}, indent=2))
'''


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def read_json(path):
    return json.loads(path.read_text())


def git(root, *args):
    return subprocess.check_output(["git", *args], cwd=root, text=True).strip()


def run_logged(command, path, *, root, env):
    with path.open("w") as stream:
        result = subprocess.run(command, cwd=root, env=env, stdout=stream,
                                stderr=subprocess.STDOUT, check=False)
    return result.returncode


def capture_runtime(interpreter, label, root, out, env):
    result = subprocess.run([interpreter, "-c", RUNTIME_CAPTURE], cwd=root,
                            env=env, text=True, capture_output=True, check=False)
    (out / (label + "-runtime-stdout.txt")).write_text(result.stdout)
    (out / (label + "-runtime-stderr.txt")).write_text(result.stderr)
    write_json(out / (label + "-runtime-exit.json"), {"returncode": result.returncode})
    if result.returncode:
        raise RuntimeError(label + " numerical runtime capture failed; see retained logs")
    runtime = json.loads(result.stdout)
    write_json(out / (label + "-runtime.json"), runtime)
    code = run_logged([interpreter, "-m", "pip", "freeze"],
                      out / (label + "-dependencies.txt"), root=root, env=env)
    if code:
        raise RuntimeError(label + " dependency capture failed")
    return runtime["environment"]


def runtime_preflight(current, historical, root):
    checks = {}
    for name in ["artifacts/week5_e1/run_metadata.json",
                 "artifacts/week6_e2/run_metadata.json",
                 "artifacts/week7_recovery_evidence/floor_audit_current_capture/archive/run_metadata.json"]:
        saved = read_json(root / name)["provenance"]["environment"]
        checks[name] = {"expected_environment": saved, "matches_historical": saved == historical}
    name = "artifacts/week8_development/audit/ci/native_fixture/native_archive/run_metadata.json"
    saved = read_json(root / name)["provenance"]["environment"]
    checks[name] = {"expected_environment": saved, "matches_current": saved == current}
    return checks


def reconcile(out, exitcodes):
    collected = read_json(out / "collection-nodeids.json")
    errors, groups, all_nodeids, skips = [], {}, [], []
    if len(collected) != len(set(collected)):
        errors.append("full collection contains duplicate nodeids")
    if len(collected) != EXPECTED["collected"]:
        errors.append("full collection differs from required 2061 cases")
    for label in GROUP_PATHS:
        result = {"returncode": exitcodes.get(label), "paths": GROUP_PATHS[label]}
        groups[label] = result
        try:
            nodeids = read_json(out / (label + "-nodeids.json"))
            reports = read_json(out / (label + "-reports.json"))
            result["collected"] = len(nodeids)
            result["unique_collected"] = len(set(nodeids))
            all_nodeids.extend(nodeids)
            tree = ET.parse(out / (label + "-junit.xml"))
            suites = list(tree.getroot().iter("testsuite"))
            counts = {key: sum(int(s.attrib.get(key, 0)) for s in suites)
                      for key in ["tests", "errors", "failures", "skipped"]}
            counts["passed"] = counts["tests"] - counts["errors"] - counts["failures"] - counts["skipped"]
            result["junit"] = counts
            cases = list(tree.getroot().iter("testcase"))
            result["junit_case_count"] = len(cases)
            result["junit_identities"] = [t.attrib.get("classname", "") + "." + t.attrib.get("name", "") for t in cases]
            skips.extend({"group": label, "test": t.attrib.get("classname", "") + "." + t.attrib.get("name", ""),
                          "reason": t.find("skipped").attrib.get("message", "")}
                         for t in cases if t.find("skipped") is not None)
            executed = set(r["nodeid"] for r in reports["reports"])
            calls = Counter(r["nodeid"] for r in reports["reports"] if r["when"] == "call")
            setup_skips = {r["nodeid"] for r in reports["reports"] if r["when"] == "setup" and r["outcome"] == "skipped"}
            result["executed_nodeids"] = sorted(executed)
            if (executed != set(nodeids) or len(nodeids) != len(set(nodeids))
                    or counts["tests"] != len(nodeids) or len(cases) != len(nodeids)
                    or any(n != 1 for n in calls.values())
                    or set(calls) | setup_skips != set(nodeids)):
                errors.append(label + ": execution/collection/JUnit identity accounting differs")
            if counts["errors"] or counts["failures"] or exitcodes.get(label) != 0:
                errors.append(label + ": failed or errored tests; original logs retained")
        except Exception as exc:
            result["accounting_error"] = repr(exc)
            errors.append(label + ": incomplete or unreadable execution evidence")
    identities_match = Counter(all_nodeids) == Counter(collected) and len(set(all_nodeids)) == len(collected)
    if not identities_match:
        errors.append("five groups do not execute every full-collection node exactly once")
    totals = {key: sum(g.get("junit", {}).get(key, 0) for g in groups.values())
              for key in ["tests", "errors", "failures", "skipped", "passed"]}
    if totals != {"tests": 2061, "passed": 2055, "skipped": 6, "errors": 0, "failures": 0}:
        errors.append("final counts differ from 2061 accounted / 2055 passed / 6 skipped / zero failures")
    if len(skips) != 6 or not all("floor" in s["reason"] for s in skips):
        errors.append("skips are not exactly the six existing design-inapplicable floor cases")
    return {"groups": groups, "totals": totals, "skips": skips,
            "collected_test_count": len(collected),
            "all_collected_tests_accounted_for": len(collected) if identities_match else 0,
            "collection_execution_identity_match": identities_match,
            "unique_execution_identity_count": len(set(all_nodeids)), "errors": errors}


def execute(root, out, summary):
    env = os.environ.copy()
    env.pop("MKL_NUM_THREADS", None)
    env.pop("VECLIB_MAXIMUM_THREADS", None)
    if env.get("OPENBLAS_NUM_THREADS") != "1" or env.get("OMP_NUM_THREADS") != "1":
        raise RuntimeError("set OPENBLAS_NUM_THREADS=1 and OMP_NUM_THREADS=1 before starting Python")
    legacy = env.get("OU_LEGACY_PYTHON")
    if not legacy or not Path(legacy).is_file():
        raise RuntimeError("OU_LEGACY_PYTHON must name the installed historical interpreter")
    env["PYTHONPATH"] = os.pathsep.join([str(root), str(out)])
    env["OU_RELEASE_RUNTIME"] = str(out / "runtime.json")
    for label, variable in [("audit", "OU_E4_AUDIT_FIXTURE_ROOT"),
                            ("plots", "OU_E4_PLOTS_FIXTURE_ROOT"),
                            ("runner", "OU_E4_RUNNER_FIXTURE_ROOT")]:
        env[variable] = str(out / (label + "_fixture"))
    (out / "audit_fixture").mkdir()
    (out / "pytest_entry.py").write_text(PYTEST_ENTRY)
    (out / "week8_release_identity_plugin.py").write_text(IDENTITY_PLUGIN)
    cpu = subprocess.run(["lscpu"], cwd=root, text=True, capture_output=True, check=False)
    (out / "cpu.txt").write_text(cpu.stdout)
    (out / "cpu-stderr.txt").write_text(cpu.stderr)
    if cpu.returncode:
        raise RuntimeError("lscpu failed; release requires the recorded Linux/x86 runtime")
    summary["current_runtime"] = capture_runtime(sys.executable, "current", root, out, env)
    summary["historical_runtime"] = capture_runtime(legacy, "legacy", root, out, env)
    write_json(out / "runtime.json", summary["current_runtime"])
    checks = runtime_preflight(summary["current_runtime"], summary["historical_runtime"], root)
    write_json(out / "runtime-preflight.json", checks)
    if not all(all(v for k, v in item.items() if k.startswith("matches_")) for item in checks.values()):
        raise RuntimeError("recorded current/historical runtime preflight mismatch")
    collection_env = env.copy()
    collection_env["OU_RELEASE_COLLECTION_OUTPUT"] = str(out / "collection-nodeids.json")
    collection_env["OU_RELEASE_REPORT_OUTPUT"] = str(out / "collection-reports.json")
    collection_command = [sys.executable, "-m", "pytest", "--collect-only", "-q", "tests", "-p", "week8_release_identity_plugin"]
    summary["collection_command"] = collection_command
    code = run_logged(collection_command, out / "collection.txt", root=root, env=collection_env)
    summary["collection_returncode"] = code
    if code:
        raise RuntimeError("full test collection failed; original output retained")
    commands, active, codes = {}, [], {}
    try:
        for label, paths in GROUP_PATHS.items():
            # Match the prior gate's separate pytest temporary roots. Preserve
            # them alongside the evidence, without uploading thousands of
            # disposable corruption-test copies as publication artifacts.
            scratch = out.parent / (out.name + "-pytest-" + label)
            if scratch.exists() or scratch.is_symlink():
                raise RuntimeError("pytest scratch path already exists: " + str(scratch))
            launcher = ([legacy, "-m", "pytest"] if label == "historical" else
                        [sys.executable, "-m", "pytest"] if label == "week3-isolation" else
                        [sys.executable, str(out / "pytest_entry.py")])
            command = [*launcher, "-q", "-ra", *paths, "-p", "week8_release_identity_plugin",
                       "--basetemp=" + str(scratch),
                       "--junitxml=" + str(out / (label + "-junit.xml"))]
            commands[label] = command
            group_env = env.copy()
            group_env["OU_RELEASE_COLLECTION_OUTPUT"] = str(out / (label + "-nodeids.json"))
            group_env["OU_RELEASE_REPORT_OUTPUT"] = str(out / (label + "-reports.json"))
            stream = (out / ("pytest-" + label + ".txt")).open("w")
            try:
                process = subprocess.Popen(command, cwd=root, env=group_env, stdout=stream, stderr=subprocess.STDOUT)
            except BaseException:
                stream.close()
                raise
            active.append((label, process, stream))
        write_json(out / "test-commands.json", commands)
        for label, process, stream in active:
            codes[label] = process.wait()
            stream.close()
            print(label, codes[label], flush=True)
    finally:
        # Even an orchestration interruption preserves all launched child output.
        for label, process, stream in active:
            if process.poll() is None:
                process.terminate()
                try:
                    process.wait(timeout=30)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait()
            codes.setdefault(label, process.returncode)
            stream.close()
        write_json(out / "test-commands.json", commands)
        write_json(out / "test-exit.json", codes)
    summary.update(reconcile(out, codes))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path, help="new evidence directory outside the repository")
    args = parser.parse_args(argv)
    root = Path(git(Path.cwd(), "rev-parse", "--show-toplevel")).resolve()
    out = args.output.expanduser().absolute()
    if out.exists() or out.is_symlink():
        parser.error("output already exists; prior evidence will never be overwritten")
    if out.resolve().is_relative_to(root) or root.is_relative_to(out.resolve()):
        parser.error("output must be outside and not an ancestor of the checkout")
    out.mkdir(parents=True)
    started = time.monotonic()
    summary = {"started_utc": datetime.now(timezone.utc).isoformat(),
               "status": "failed", "scientific_freeze": FREEZE,
               "candidate_commit": git(root, "rev-parse", "HEAD"),
               "repository": str(root), "evidence": str(out),
               "group_collection_paths": GROUP_PATHS, "required_counts": EXPECTED,
               "strict_week3_note": "Exact Week3 assertions are unchanged and mandatory in an isolated process. CPU dispatch sensitivity remains a limitation; no skip, tolerance relaxation or CPU-feature override is introduced here.",
               "full_E4_started": False}
    try:
        before = git(root, "status", "--porcelain", "--untracked-files=all")
        (out / "git-status-before.txt").write_text(before + "\n")
        if before:
            raise RuntimeError("release regression requires a clean checkout")
        with (out / "scientific-source-diff.txt").open("w") as stream:
            subprocess.run(["git", "diff", "--exit-code", FREEZE, "HEAD", "--", "tests", "ou_irregular"], cwd=root, check=True,
                           stdout=stream, stderr=subprocess.STDOUT)
        execute(root, out, summary)
        after = git(root, "status", "--porcelain", "--untracked-files=all")
        (out / "git-status-after.txt").write_text(after + "\n")
        summary["checkout_remained_clean"] = after == ""
        if after:
            summary.setdefault("errors", []).append("checkout changed during regression")
        if not summary.get("errors"):
            summary["status"] = "passed"
    except BaseException as exc:
        summary["orchestration_exception"] = repr(exc)
        (out / "orchestration-traceback.txt").write_text(traceback.format_exc())
        summary["status"] = "failed"
    finally:
        summary["finished_utc"] = datetime.now(timezone.utc).isoformat()
        summary["elapsed_seconds"] = time.monotonic() - started
        write_json(out / "regression_summary.json", summary)
        records = {str(p.relative_to(out)): {"sha256": hashlib.sha256(p.read_bytes()).hexdigest(), "bytes": p.stat().st_size}
                   for p in sorted(out.rglob("*")) if p.is_file() and p.name != "evidence_sha256.json"}
        write_json(out / "evidence_sha256.json", records)
    print(json.dumps({"status": summary["status"], "summary": str(out / "regression_summary.json"),
                      "totals": summary.get("totals"), "errors": summary.get("errors", [])}, indent=2))
    return 0 if summary["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
