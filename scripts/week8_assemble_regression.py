"""Assemble a transparent multihost regression gate from retained evidence.

The native 41215391840 invocation remains failed and is copied in full. Its
four successful groups are combined with a separately executed, unchanged
four-case historical group. No test runs here, no original evidence changes,
and the resulting gate is explicitly not described as one successful run.
"""
import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import shutil
import types
import xml.etree.ElementTree as ET
import zipfile


FREEZE = "0ff57919f5d8345f52d8689732ff94448ca3ca1d"
NATIVE_RUN = 41215391840
NATIVE_ARTIFACT = 13048162246
NATIVE_CONTROL = "e54322849cb7c6d4993c46f278154cae57f44cdd"
NATIVE_ZIP_SHA = "6ba2b9a939b815771540cc99ab11e304f36438033013dc0e94b3078e9b275ea9"
HELPER_SHA = "4138b209030475b15ec5ae1bd8a986f9411caca1f1cbf1ec4b793e38174822ff"
NATIVE_WORKFLOW_SHA = "872a447cc7294abcb0497b87747358a9680756c75027c4716b0602dbc0d42174"
ORIGINAL_STRICT_CONTROL = "d5dee8da7012f868ebc106160bea3a9dd39c6dd2"
ORIGINAL_STRICT_SCRIPT_SHA = "04e829635253a9beb7c60c86e5b47066f1273b6b78fc73e3d69b7d6b967a58d5"
STRICT_TEST_SHA = "0628f727144568684205a5742f1746b19c3437bc1f75c4b6b0ba69720a32a0ea"
SDE_ARCHIVE_SHA = "94e97d623fec54385686e1e7ba65ebc9941748c05ee451423948334892bf2b50"
GOOD_NATIVE_GROUPS = {"core": (653, 653, 0), "week7": (775, 769, 6),
                      "week8": (492, 492, 0), "historical": (137, 137, 0)}
STRICT_PATH = "tests/week4/test_week3_isolation.py"
STRICT_FILES = ["week3-isolation-nodeids.json", "week3-isolation-reports.json",
                "week3-isolation-junit.xml", "pytest-week3-isolation.txt"]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read(path):
    return json.loads(path.read_text())


def write(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def inventory(root):
    require(root.is_dir() and not root.is_symlink(), "evidence root must be an ordinary directory")
    result = {}
    for p in sorted(root.rglob("*")):
        require(not p.is_symlink(), "evidence contains symbolic link: " + str(p))
        if p.is_file():
            result[str(p.relative_to(root))] = {"bytes": p.stat().st_size, "sha256": digest(p)}
        else:
            require(p.is_dir(), "evidence contains a nonregular object")
    return result


def verified_manifest(root):
    actual = inventory(root)
    require("evidence_sha256.json" in actual, "missing source evidence manifest")
    declared = read(root / "evidence_sha256.json")
    require(type(declared) is dict, "evidence manifest must be an object")
    for name in declared:
        p = PurePosixPath(name)
        require(not p.is_absolute() and ".." not in p.parts and str(p) == name,
                "noncanonical manifest path")
    require({n: r for n, r in actual.items() if n != "evidence_sha256.json"} == declared,
            "source manifest must cover every original file exactly")
    return actual


def verify_transport(path, metadata, extracted):
    if path is None:
        return {"checked_in_assembly": False,
                "scope": "No ZIP supplied to assembly. Caller transport verification remains separate; assembly checks complete extracted manifests."}
    require(path.is_file() and not path.is_symlink(), "transport must be an ordinary ZIP file")
    require(metadata["digest"] == "sha256:" + digest(path), "transport differs from API artifact digest")
    if "size_in_bytes" in metadata:
        require(path.stat().st_size == metadata["size_in_bytes"], "transport size differs from API artifact metadata")
    members = {}
    with zipfile.ZipFile(path) as archive:
        for entry in archive.infolist():
            name = PurePosixPath(entry.filename)
            require(not name.is_absolute() and ".." not in name.parts and str(name) == entry.filename.rstrip("/"),
                    "noncanonical ZIP member")
            if entry.is_dir():
                continue
            require(str(name) not in members, "duplicate ZIP member")
            data = archive.read(entry)
            members[str(name)] = {"bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}
    require(members == extracted, "ZIP members do not equal the safely extracted evidence")
    return {"checked_in_assembly": True, "sha256": digest(path),
            "bytes": path.stat().st_size, "regular_member_count": len(members)}


def artifact_metadata(path, expected_control):
    value = read(path)
    require(type(value) is dict and int(value["id"]) > 0, "invalid GitHub artifact metadata")
    require(value.get("expired") is False, "artifact metadata says expired")
    require(re.fullmatch(r"sha256:[0-9a-f]{64}", value.get("digest", "")) is not None,
            "artifact metadata lacks a SHA-256 transport digest")
    require(value["workflow_run"]["head_sha"] == expected_control, "artifact workflow head differs from declared control")
    return value


def load_helper(path):
    require(digest(path) == HELPER_SHA, "accounting helper differs from original executed control bytes")
    module = types.ModuleType("week8_original_regression_accounting")
    module.__file__ = str(path)
    exec(compile(path.read_text(), str(path), "exec"), module.__dict__)
    return module


def strict_validation(root, outer, runtime, expected_native_nodes, journal):
    emulated = "Intel SDE" in journal.get("execution_mode", "")
    require(read((outer if emulated else root) / "exit.json")["returncode"] == 0, "strict invocation did not pass")
    require(not (root / "hardware_unavailable.json").exists(), "strict artifact contains hardware-unavailable disposition")
    nodes = read(root / "week3-isolation-nodeids.json")
    require(len(nodes) == len(set(nodes)) == 4 and set(nodes) == set(expected_native_nodes),
            "strict invocation must collect exactly the original four identities")
    require(all(n.startswith(STRICT_PATH + "::") for n in nodes), "strict invocation includes another test path")
    reports = read(root / "week3-isolation-reports.json")
    require(reports["exitstatus"] == 0, "strict observer records a nonzero exit")
    require(all(r["outcome"] == "passed" for r in reports["reports"]), "strict invocation contains skip/failure/error")
    require(Counter((r["nodeid"], r["when"]) for r in reports["reports"]) ==
            Counter((n, stage) for n in nodes for stage in ["setup", "call", "teardown"]),
            "strict observer does not account for each test stage exactly once")
    xml = ET.parse(root / "week3-isolation-junit.xml")
    suites = list(xml.getroot().iter("testsuite"))
    counts = {key: sum(int(s.attrib.get(key, 0)) for s in suites)
              for key in ["tests", "errors", "failures", "skipped"]}
    require(counts == {"tests": 4, "errors": 0, "failures": 0, "skipped": 0}, "strict JUnit does not record four passes")
    require(len(list(xml.getroot().iter("testcase"))) == 4, "strict JUnit case identity count differs")
    if emulated:
        command = journal["command"]
        selection = journal["pytest_arguments"]
        require(command == read(outer / "command.json") and len(command) >= 6
                and command[1:4] == ["-no-follow-child", "-icl", "--"],
                "emulator command differs from retained parent-only instrumentation invocation")
        require(command.count("--help-file") == 1 and command.index("--help-file") + 1 < len(command),
                "emulator invocation lacks its retained driver-help file")
        help_path = outer / Path(command[command.index("--help-file") + 1]).name
        require(help_path.is_file() and digest(help_path) == journal["emulator_help_sha256"]
                and "-no-follow-child" in help_path.read_text(),
                "retained SDK driver help does not bind its documented no-follow-child option")
        parsed_path = outer / "internal-pin-command.txt"
        parsed = parsed_path.read_text()
        require(read(outer / "internal-pin-command-exit.json")["returncode"] == 0
                and bool(parsed.strip()) and command[4] in parsed and Path(command[5]).name in parsed
                and "-follow_execv" not in parsed and "-follow-execv" not in parsed,
                "SDK echo did not record the actual target command with exec-follow disabled")
        smoke_path = root / "runtime_smoke.json"
        smoke = read(smoke_path)
        smoke_binding = read(root / "evidence_sha256.json")["runtime_smoke.json"]
        require(smoke == journal["runtime_smoke"]
                and smoke_binding == {"sha256": digest(smoke_path), "bytes": smoke_path.stat().st_size},
                "actual runtime smoke is not bound to its source-provenance record and manifest")
        require(smoke["parent_pid"] == runtime["pid"] == journal["pytest_pid"]
                and smoke["git_command"] == ["git", "rev-parse", "HEAD"]
                and smoke["git_returncode"] == 0 and smoke["git_stdout"].strip() == FREEZE
                and smoke["effective_AVX512F"] is True and smoke["effective_AVX512_ICL"] is True,
                "same-parent runtime smoke did not establish auxiliary Git execution and effective required features")
        require(journal["pytest_pid"] == runtime["pid"] and journal["pytest_returncode"] == 0
                and journal["effective_AVX512F_and_ICL_verified"] is True
                and runtime["cpu_features"].get("AVX512_ICL") is True
                and runtime["execution_mode"] == journal["execution_mode"],
                "effective same-process emulator runtime is not bound to successful pytest execution")
        require(journal["test_sha256"] == STRICT_TEST_SHA and journal["emulator_archive_sha256"] == SDE_ARCHIVE_SHA
                and re.fullmatch(r"[0-9a-f]{64}", journal.get("emulator_executable_sha256", "")) is not None
                and bool(journal.get("emulator_version", "").strip())
                and (outer / "emulator-version.txt").read_text() == journal["emulator_version"],
                "emulator/test identity is incomplete or changed")
        require(journal["junit_counts"] == counts and not journal["historical_provenance_differences"],
                "emulator source-provenance record differs from successful original historical assertions")
    else:
        command = read(root / "command.json")
        selection = command
        require(type(command) is list and "-m" in command and "pytest" in command,
                "strict native command is not the original pytest invocation")
    require(type(selection) is list and selection.count(STRICT_PATH) == 1
            and not any(x.startswith(("-k", "--deselect", "--ignore")) for x in selection),
            "strict invocation changes original four-case selection")
    require(runtime["cpu_features"].get("AVX512F") is True, "strict runtime does not expose required AVX-512 capability")
    controls = runtime.get("cpu_control_environment", {})
    require(all(controls.get(n) is None for n in ["NPY_DISABLE_CPU_FEATURES", "OPENBLAS_CORETYPE", "MKL_ENABLE_INSTRUCTIONS"]),
            "strict runtime introduced an unreviewed CPU-dispatch override")
    return command


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--native", type=Path, required=True)
    p.add_argument("--strict", type=Path, required=True)
    p.add_argument("--native-metadata", type=Path, required=True)
    p.add_argument("--strict-metadata", type=Path, required=True)
    p.add_argument("--strict-control-commit", required=True)
    p.add_argument("--destination", type=Path, required=True, help="completion evidence directory; new child paths must not exist")
    p.add_argument("--native-zip", type=Path)
    p.add_argument("--strict-zip", type=Path)
    p.add_argument("--regression-script", type=Path, default=Path(__file__).with_name("week8_release_regression.py"))
    p.add_argument("--native-workflow", type=Path, default=Path(__file__).parents[1] / ".github/workflows/week8-regression.yml")
    p.add_argument("--strict-script", type=Path, help="defaults to native or SDE launcher according to the artifact")
    p.add_argument("--strict-workflow", type=Path, help="defaults to the native or SDE workflow according to the artifact")
    args = p.parse_args(argv)
    require(re.fullmatch(r"[0-9a-f]{40}", args.strict_control_commit) is not None, "strict control must be a full commit SHA")
    native, strict, destination = args.native.resolve(), args.strict.resolve(), args.destination.resolve()
    require(native != strict, "native and strict evidence must be distinct")
    for source in [native, strict]:
        require(not destination.is_relative_to(source) and not source.is_relative_to(destination),
                "destination and original inputs must not overlap")
    native_meta = artifact_metadata(args.native_metadata, NATIVE_CONTROL)
    strict_meta = artifact_metadata(args.strict_metadata, args.strict_control_commit)
    require(int(native_meta["id"]) == NATIVE_ARTIFACT and int(native_meta["workflow_run"]["id"]) == NATIVE_RUN,
            "native artifact/run identity differs")
    require(native_meta["name"] == "week8-publication-regression-" + str(NATIVE_RUN)
            and native_meta["digest"] == "sha256:" + NATIVE_ZIP_SHA, "native artifact name/digest differs")
    require(int(strict_meta["id"]) != NATIVE_ARTIFACT and strict_meta["name"].startswith(("week8-compatible-replay-", "week3-strict-sde-replay-"))
            and strict_meta["name"].endswith("-" + str(strict_meta["workflow_run"]["id"])), "strict artifact/run identity differs")
    strict_result = strict / "replay" if (strict / "replay").is_dir() else strict
    native_files = verified_manifest(native)
    verified_manifest(strict_result)
    strict_files = inventory(strict)
    if strict_result != strict:
        require(args.strict_zip is not None, "an outer artifact without a root manifest requires --strict-zip for full transport binding")
    require(len(native_files) == 271, "native invocation does not retain all 271 original files")
    native_transport = verify_transport(args.native_zip, native_meta, native_files)
    strict_transport = verify_transport(args.strict_zip, strict_meta, strict_files)
    helper = load_helper(args.regression_script)
    require((native / "week8_release_identity_plugin.py").read_bytes() == helper.IDENTITY_PLUGIN.encode(),
            "native identity-observer bytes differ from the original helper")
    require(digest(args.native_workflow) == NATIVE_WORKFLOW_SHA, "native workflow differs from original executed control")
    original = read(native / "regression_summary.json")
    codes = read(native / "test-exit.json")
    native_check = helper.reconcile(native, codes)
    require(original["status"] == "failed" and original["candidate_commit"] == original["scientific_freeze"] == FREEZE,
            "native summary lost original failed scientific-freeze identity")
    require(original["checkout_remained_clean"] is True and original["collection_execution_identity_match"] is True,
            "native invocation did not preserve a clean complete collection")
    require(original["totals"] == native_check["totals"] == {"tests": 2061, "passed": 2054, "skipped": 6, "failures": 1, "errors": 0},
            "native raw evidence differs from retained 2054/6/1 failure")
    for key, value in native_check.items():
        require(original[key] == value, "native summary does not reconcile with original raw evidence: " + key)
    for label, (tests, passed, skipped) in GOOD_NATIVE_GROUPS.items():
        require(codes[label] == 0 and native_check["groups"][label]["junit"] ==
                {"tests": tests, "passed": passed, "skipped": skipped, "failures": 0, "errors": 0},
                "native group is not fully successful: " + label)
    require(codes["week3-isolation"] == 1 and native_check["groups"]["week3-isolation"]["junit"] ==
            {"tests": 4, "passed": 3, "skipped": 0, "failures": 1, "errors": 0}, "native strict group failure accounting changed")
    for name in ["git-status-before.txt", "git-status-after.txt", "scientific-source-diff.txt"]:
        require(not (native / name).read_text().strip(), "native scientific checkout was modified")
    native_runtime = read(native / "current-runtime.json")
    runtime_name = "same_process_runtime.json" if (strict_result / "same_process_runtime.json").is_file() else "runtime.json"
    strict_runtime = read(strict_result / runtime_name)
    journal = strict_result / "source_provenance.json"
    if journal.exists():
        source_provenance = read(journal)
        if "candidate_commit" in source_provenance:
            require(source_provenance["candidate_commit"] == source_provenance["scientific_freeze"] == FREEZE
                    and source_provenance["git_status_before"] == source_provenance["git_status_after"] == ""
                    and source_provenance["clean_before_and_after"] is True
                    and source_provenance["test_path"] == STRICT_PATH and source_provenance["status"] == "passed",
                    "strict source/cleanliness journal does not bind the original unchanged source")
        else:
            require(source_provenance["scientific_commit"] == FREEZE and source_provenance["checkout_clean_before"] is True
                    and source_provenance["checkout_clean_after"] is True and source_provenance["executed_test_path"] == STRICT_PATH,
                    "strict source/cleanliness journal does not bind the original unchanged source")
        require(source_provenance.get("control_commit", args.strict_control_commit) == args.strict_control_commit,
                "strict source journal control differs from artifact control")
    else:
        require(args.strict_control_commit == ORIGINAL_STRICT_CONTROL,
                "a revised strict control must retain source_provenance.json")
        source_provenance = {"scientific_commit": FREEZE, "control_commit": args.strict_control_commit,
                             "basis": "Original pinned workflow checks out the scientific freeze; original strict launcher asserts HEAD and cleanliness before/after tests. Its runtime capture did not separately save checkout SHA."}
    emulated = "Intel SDE" in source_provenance.get("execution_mode", "")
    strict_script = args.strict_script or Path(__file__).with_name("week8_sde_replay.py" if emulated else "week8_strict_replay.py")
    strict_workflow = args.strict_workflow or Path(__file__).parents[1] / (".github/workflows/week8-sde-replay.yml" if emulated else ".github/workflows/week8-compatible-replay.yml")
    if not journal.exists():
        require(digest(strict_script) == ORIGINAL_STRICT_SCRIPT_SHA, "original strict script bytes differ")
    script_hash = source_provenance.get("script_sha256", source_provenance.get("strict_script_sha256"))
    if script_hash is not None:
        require(script_hash == digest(strict_script), "strict script bytes differ from source journal")
    if not emulated:
        require((strict_result / "week8_release_identity_plugin.py").read_bytes() == helper.IDENTITY_PLUGIN.encode(),
                "strict identity-observer bytes differ from original helper")
    expected_environment = json.loads(json.dumps(original["current_runtime"]))
    require(native_runtime["environment"] == expected_environment, "native current runtime changed")
    if emulated:
        # This is the captured state after importing the unchanged historical
        # runner, which itself sets MKL_NUM_THREADS=1. Record this narrow actual
        # import side effect; all other software/thread fields must still agree.
        require(expected_environment["thread_environment"]["MKL_NUM_THREADS"] is None,
                "native current runtime has an unexpected MKL setting")
        expected_environment["thread_environment"]["MKL_NUM_THREADS"] = "1"
    require(strict_runtime["environment"] == expected_environment, "strict software/thread environment differs beyond the declared original-runner import")
    command = strict_validation(strict_result, strict, strict_runtime, read(native / "week3-isolation-nodeids.json"), source_provenance)
    attempts = destination / "regression_attempts"
    native_target = attempts / str(NATIVE_RUN)
    strict_target = attempts / ("strict_artifact_" + str(strict_meta["id"]))
    final = destination / "final_regression"
    report_path = destination / "regression_assembly.json"
    require(all(not x.exists() and not x.is_symlink() for x in [native_target, strict_target, final, report_path]),
            "assembly target already exists; nothing may overwrite prior evidence")
    destination.mkdir(parents=True, exist_ok=True)
    attempts.mkdir(exist_ok=True)
    shutil.copytree(native, native_target)
    shutil.copytree(strict, strict_target)
    require(inventory(native_target) == native_files and inventory(strict_target) == strict_files, "preserved attempt copy differs")
    shutil.copytree(native, final)
    for name in STRICT_FILES:
        source = strict_result / name
        if name == "pytest-week3-isolation.txt" and not source.is_file():
            source = strict / "week3-isolation.txt"
        shutil.copyfile(source, final / name)
    codes["week3-isolation"] = 0
    write(final / "test-exit.json", codes)
    commands = read(native / "test-commands.json")
    commands["week3-isolation"] = command
    write(final / "test-commands.json", commands)
    shutil.copyfile(strict_result / runtime_name, final / "week3-compatible-runtime.json")
    shutil.copyfile(strict / "dependencies.txt", final / "week3-compatible-dependencies.txt")
    shutil.copyfile(strict / ("host-lscpu.txt" if emulated else "cpu.txt"), final / "week3-compatible-cpu.txt")
    shutil.copyfile(args.regression_script, final / "regression_script_executed.py")
    shutil.copyfile(args.native_workflow, final / "regression_workflow.yml")
    shutil.copyfile(strict_script, final / "strict_script_executed.py")
    shutil.copyfile(strict_workflow, final / "strict_workflow.yml")
    shutil.copyfile(args.native_metadata, final / "native_artifact_metadata.json")
    shutil.copyfile(args.strict_metadata, final / "strict_artifact_metadata.json")
    accounting = helper.reconcile(final, codes)
    require(accounting["errors"] == [] and accounting["totals"] == {"tests": 2061, "passed": 2055, "skipped": 6, "failures": 0, "errors": 0},
            "composite raw evidence does not account for every required test exactly once")
    summary = {k: v for k, v in original.items() if k not in native_check and k not in
               {"status", "started_utc", "finished_utc", "elapsed_seconds", "repository", "evidence", "strict_week3_note"}}
    group_origins = {}
    for label in helper.GROUP_PATHS:
        compatible = label == "week3-isolation"
        base = strict_target if compatible else native_target
        group_origins[label] = {"artifact_id": strict_meta["id"] if compatible else native_meta["id"],
                                "workflow_run_id": strict_meta["workflow_run"]["id"] if compatible else NATIVE_RUN,
                                "control_commit": args.strict_control_commit if compatible else NATIVE_CONTROL,
                                "scientific_commit": FREEZE,
                                "attempt_evidence": str(base.relative_to(destination)),
                                "runtime_capture": str((base / strict_result.relative_to(strict) / runtime_name if compatible else
                                                        base / ("legacy-runtime.json" if label == "historical" else "current-runtime.json")).relative_to(destination)),
                                "execution_mode": source_provenance.get("execution_mode", "native hosted CPU") if compatible else "native hosted CPU"}
    summary.update(accounting)
    summary.update(status="passed", gate_kind="composite_five_disjoint_groups_across_separate_CPU_runtimes",
                   single_successful_invocation=False, assembled_utc=datetime.now(timezone.utc).isoformat(),
                   candidate_commit=FREEZE, scientific_freeze=FREEZE, checkout_remained_clean=True,
                   evidence=str(final), original_failed_invocation_retained=True,
                   native_original_totals=original["totals"], group_provenance=group_origins,
                   native_execution_started_utc=original["started_utc"], native_execution_finished_utc=original["finished_utc"],
                   strict_source_provenance=source_provenance,
                   strict_execution_mode=source_provenance.get("execution_mode", "native hosted CPU"),
                   strict_effective_environment=strict_runtime["environment"],
                   strict_runtime_smoke_sha256=digest(strict_result / "runtime_smoke.json") if emulated else None,
                   strict_internal_pin_command_sha256=digest(strict / "internal-pin-command.txt") if emulated else None,
                   strict_runtime_capture_scope="Same pytest process after unchanged historical-runner import; MKL_NUM_THREADS is 1 from that module's original import side effect." if emulated else "Fresh interpreter runtime before unchanged historical pytest subprocess.",
                   assembly_scope="Four complete successful native groups (2057 tests:2051 passes and six existing skips) plus the entire new four-pass historical group. The original native invocation remains failed (2054 passes/six skips/one failure). No test, seed, hash assertion or tolerance changed; this record is not one successful invocation.",
                   strict_week3_note="Original exact historical assertions passed in a separately captured compatible CPU runtime. Native AVX2-only failure and every original invocation file remain separately retained; no universal bitwise-replay claim is made.")
    write(final / "regression_summary.json", summary)
    # This new manifest covers the assembled gate including provenance additions.
    write(final / "evidence_sha256.json", {n: r for n, r in inventory(final).items() if n != "evidence_sha256.json"})
    verified_manifest(final)
    require(inventory(native) == native_files and inventory(strict) == strict_files, "original inputs changed during assembly")
    require(inventory(native_target) == native_files and inventory(strict_target) == strict_files, "preserved attempts changed during assembly")
    record = {"passed": True, "assembled_utc": summary["assembled_utc"], "scientific_freeze": FREEZE,
              "native_failed_attempt": {"artifact_metadata": native_meta, "transport": native_transport,
                                        "retained_directory": str(native_target.relative_to(destination)), "files": native_files},
              "strict_successful_attempt": {"artifact_metadata": strict_meta, "transport": strict_transport,
                                            "retained_directory": str(strict_target.relative_to(destination)), "files": strict_files},
              "group_provenance": group_origins, "composite_totals": accounting["totals"],
              "complete_unique_node_identity_match": accounting["collection_execution_identity_match"],
              "strict_runtime_smoke_sha256": digest(strict_result / "runtime_smoke.json") if emulated else None,
              "strict_internal_pin_command_sha256": digest(strict / "internal-pin-command.txt") if emulated else None,
              "original_attempts_unchanged": True, "no_new_test_execution_in_assembly": True,
              "single_successful_invocation": False, "assembly_script_sha256": digest(Path(__file__))}
    write(report_path, record)
    print(json.dumps({"passed": True, "final_regression": str(final), "assembly_record": str(report_path),
                      "totals": accounting["totals"], "single_successful_invocation": False}, indent=2))


if __name__ == "__main__":
    main()
