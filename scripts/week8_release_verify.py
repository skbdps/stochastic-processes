"""Read-only final release verification from a clean full-history checkout.

Checks completed saved E4, every original artifact byte, full regression
evidence, all preexisting Git blobs, current docs and their relative links.
No new experiment, fit, numerical audit or figure render is performed.
"""
import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import subprocess
import sys
from urllib.parse import unquote
import zipfile

FREEZE = "0ff57919f5d8345f52d8689732ff94448ca3ca1d"
ZIP_SHA = "e9a780206925c4085ff4bb47cf988c7ccec66ce40db038f3fdef82703fe52513"
DOCS = ["README.md", "artifacts/README.md", "docs/week8_completion.md", "docs/week8_progress.md",
        "docs/week8_results_interpretation.md", "docs/week8_reproduction.md", "docs/week8_publication_recovery.md",
        "docs/week8_ci.md"]


def sha(data):
    return hashlib.sha256(data).hexdigest()


def git(*args):
    return subprocess.check_output(["git", *args], text=True).strip()


def tree(ref):
    return {path: identity for line in git("ls-tree", "-r", ref).splitlines()
            for identity, path in [line.split("\t", 1)]}


def evidence_path(base, relative):
    """Accept only an ordinary canonical path beneath this evidence root."""
    name = PurePosixPath(relative)
    assert not name.is_absolute() and ".." not in name.parts and str(name) == relative
    target = base / str(name)
    assert target.resolve().is_relative_to(base.resolve()) and target.exists(), relative
    assert not target.is_symlink(), relative
    return target


def verify_composite(root, evidence, regression, summary, accounting):
    """Read original attempts again; never replace a failed invocation record."""
    from scripts import week8_assemble_regression as assembly
    from scripts.week8_release_regression import reconcile, GROUP_PATHS

    record_path = evidence / "regression_assembly.json"
    record = json.loads(record_path.read_text())
    assert record["passed"] is True and record["scientific_freeze"] == FREEZE
    assert record["single_successful_invocation"] is False
    assert record["original_attempts_unchanged"] is True and record["no_new_test_execution_in_assembly"] is True
    assert record["complete_unique_node_identity_match"] is True
    assert record["composite_totals"] == accounting["totals"] == summary["totals"]
    assert record["assembly_script_sha256"] == assembly.digest(root / "scripts/week8_assemble_regression.py")
    assert summary["gate_kind"] == "composite_five_disjoint_groups_across_separate_CPU_runtimes"
    assert summary["single_successful_invocation"] is False
    assert summary["original_failed_invocation_retained"] is True
    assert summary["scientific_freeze"] == summary["candidate_commit"] == FREEZE
    assert summary["group_provenance"] == record["group_provenance"]
    assert set(summary["group_provenance"]) == set(GROUP_PATHS)

    native_record = record["native_failed_attempt"]
    strict_record = record["strict_successful_attempt"]
    native = evidence_path(evidence, native_record["retained_directory"])
    strict = evidence_path(evidence, strict_record["retained_directory"])
    assert native != strict and native.is_dir() and strict.is_dir()
    assert not native.is_relative_to(strict) and not strict.is_relative_to(native)
    native_files = assembly.verified_manifest(native)
    strict_files = assembly.inventory(strict)
    assert len(native_files) == 271 and native_files == native_record["files"]
    assert strict_files == strict_record["files"]
    strict_result = strict / "replay" if (strict / "replay").is_dir() else strict
    assembly.verified_manifest(strict_result)
    native_meta = assembly.artifact_metadata(regression / "native_artifact_metadata.json", assembly.NATIVE_CONTROL)
    strict_head = strict_record["artifact_metadata"]["workflow_run"]["head_sha"]
    strict_meta = assembly.artifact_metadata(regression / "strict_artifact_metadata.json", strict_head)
    assert native_meta == native_record["artifact_metadata"] and strict_meta == strict_record["artifact_metadata"]
    assert int(native_meta["id"]) == assembly.NATIVE_ARTIFACT
    assert int(native_meta["workflow_run"]["id"]) == assembly.NATIVE_RUN
    assert native_meta["digest"] == "sha256:" + assembly.NATIVE_ZIP_SHA
    for attempt_record in (native_record, strict_record):
        transport = attempt_record["transport"]
        metadata = attempt_record["artifact_metadata"]
        assert transport["checked_in_assembly"] is True
        assert metadata["digest"] == "sha256:" + transport["sha256"]
        assert transport["regular_member_count"] == len(attempt_record["files"])
        if "size_in_bytes" in metadata:
            assert transport["bytes"] == metadata["size_in_bytes"]

    original = json.loads((native / "regression_summary.json").read_text())
    native_codes = json.loads((native / "test-exit.json").read_text())
    original_accounting = reconcile(native, native_codes)
    failed_totals = {"tests": 2061, "passed": 2054, "skipped": 6, "failures": 1, "errors": 0}
    assert original["status"] == "failed" and original["totals"] == failed_totals
    assert summary["native_original_totals"] == failed_totals
    assert original["candidate_commit"] == original["scientific_freeze"] == FREEZE
    assert original["checkout_remained_clean"] is True
    for key, value in original_accounting.items():
        assert original[key] == value, key
    assert original_accounting["totals"] == failed_totals
    for label, (tests, passed, skipped) in assembly.GOOD_NATIVE_GROUPS.items():
        assert native_codes[label] == 0
        assert original_accounting["groups"][label]["junit"] == {
            "tests": tests, "passed": passed, "skipped": skipped, "errors": 0, "failures": 0}
    assert native_codes["week3-isolation"] == 1
    assert original_accounting["groups"]["week3-isolation"]["junit"] == {
        "tests": 4, "passed": 3, "skipped": 0, "errors": 0, "failures": 1}
    for name in ("git-status-before.txt", "git-status-after.txt", "scientific-source-diff.txt"):
        assert not (native / name).read_text().strip(), name

    # The native data, collection and four good groups remain exact copies.
    replacements = set(assembly.STRICT_FILES) | {
        "test-exit.json", "test-commands.json", "regression_summary.json", "evidence_sha256.json"}
    for name, identity in native_files.items():
        if name not in replacements:
            target = evidence_path(regression, name)
            assert target.stat().st_size == identity["bytes"] and assembly.digest(target) == identity["sha256"], name
    commands = json.loads((regression / "test-commands.json").read_text())
    native_commands = json.loads((native / "test-commands.json").read_text())
    final_codes = json.loads((regression / "test-exit.json").read_text())
    assert final_codes == {label: 0 for label in GROUP_PATHS}
    for label in assembly.GOOD_NATIVE_GROUPS:
        assert commands[label] == native_commands[label], label
    for name in assembly.STRICT_FILES:
        source = strict_result / name
        if name == "pytest-week3-isolation.txt" and not source.is_file():
            source = strict / "week3-isolation.txt"
        assert source.read_bytes() == (regression / name).read_bytes(), name

    journal = json.loads((strict_result / "source_provenance.json").read_text())
    assert journal == summary["strict_source_provenance"]
    assert journal["status"] == "passed"
    assert journal["candidate_commit"] == journal["scientific_freeze"] == FREEZE
    assert journal["git_status_before"] == journal["git_status_after"] == ""
    assert journal["clean_before_and_after"] is True and journal["test_path"] == assembly.STRICT_PATH
    assert journal["control_commit"] == strict_head
    runtime_name = "same_process_runtime.json" if (strict_result / "same_process_runtime.json").is_file() else "runtime.json"
    runtime = json.loads((strict_result / runtime_name).read_text())
    emulated = "Intel SDE" in journal["execution_mode"]
    assert summary["strict_execution_mode"] == journal["execution_mode"]
    assert summary["strict_effective_environment"] == runtime["environment"]
    reference_environment = json.loads(json.dumps(original["current_runtime"]))
    if emulated:
        assert reference_environment["thread_environment"]["MKL_NUM_THREADS"] is None
        reference_environment["thread_environment"]["MKL_NUM_THREADS"] = "1"
    assert runtime["environment"] == reference_environment
    strict_command = assembly.strict_validation(strict_result, strict, runtime,
                                               json.loads((native / "week3-isolation-nodeids.json").read_text()), journal)
    if emulated:
        smoke_hash = assembly.digest(strict_result / "runtime_smoke.json")
        assert summary["strict_runtime_smoke_sha256"] == record["strict_runtime_smoke_sha256"] == smoke_hash
        parsed_pin_hash = assembly.digest(strict / "internal-pin-command.txt")
        assert summary["strict_internal_pin_command_sha256"] == record["strict_internal_pin_command_sha256"] == parsed_pin_hash
    assert commands["week3-isolation"] == strict_command
    assert journal["script_sha256"] == assembly.digest(regression / "strict_script_executed.py")
    assert (strict_result / runtime_name).read_bytes() == (regression / "week3-compatible-runtime.json").read_bytes()
    assert (strict / "dependencies.txt").read_bytes() == (regression / "week3-compatible-dependencies.txt").read_bytes()
    assert (strict / ("host-lscpu.txt" if emulated else "cpu.txt")).read_bytes() == (regression / "week3-compatible-cpu.txt").read_bytes()

    # Bind retained launcher/workflow bytes to the actual control commits,
    # rather than relying only on filenames or a later equivalent launcher.
    controls = [("regression_script_executed.py", assembly.NATIVE_CONTROL, "scripts/week8_release_regression.py"),
                ("regression_workflow.yml", assembly.NATIVE_CONTROL, ".github/workflows/week8-regression.yml"),
                ("strict_script_executed.py", strict_head, "scripts/week8_sde_replay.py" if emulated else "scripts/week8_strict_replay.py"),
                ("strict_workflow.yml", strict_head, ".github/workflows/week8-sde-replay.yml" if emulated else ".github/workflows/week8-compatible-replay.yml")]
    for retained_name, commit, source in controls:
        assert (regression / retained_name).read_bytes() == subprocess.check_output(["git", "show", commit + ":" + source]), retained_name
    assert assembly.digest(regression / "regression_script_executed.py") == assembly.HELPER_SHA
    assert assembly.digest(regression / "regression_workflow.yml") == assembly.NATIVE_WORKFLOW_SHA

    runtime_links = {}
    for label, origin in summary["group_provenance"].items():
        compatible = label == "week3-isolation"
        base = strict if compatible else native
        metadata = strict_meta if compatible else native_meta
        assert origin["artifact_id"] == metadata["id"] and origin["workflow_run_id"] == metadata["workflow_run"]["id"]
        assert origin["control_commit"] == metadata["workflow_run"]["head_sha"] and origin["scientific_commit"] == FREEZE
        assert evidence_path(evidence, origin["attempt_evidence"]) == base
        actual_runtime = evidence_path(evidence, origin["runtime_capture"])
        expected_runtime = strict_result / runtime_name if compatible else base / (
            "legacy-runtime.json" if label == "historical" else "current-runtime.json")
        assert actual_runtime == expected_runtime
        assert origin["execution_mode"] == (journal["execution_mode"] if compatible else "native hosted CPU")
        capture = json.loads(actual_runtime.read_text())
        expected_env = runtime["environment"] if compatible else original[
            "historical_runtime" if label == "historical" else "current_runtime"]
        assert capture["environment"] == expected_env
        runtime_links[label] = origin["runtime_capture"]
    return {"verified": True, "single_successful_invocation": False,
            "native_failed_attempt_totals": failed_totals,
            "preserved_input_attempt_files": {"native": len(native_files), "strict": len(strict_files)},
            "strict_execution_mode": journal["execution_mode"], "per_group_runtime_captures": runtime_links,
            "strict_same_parent_runtime_smoke_sha256": summary["strict_runtime_smoke_sha256"],
            "strict_parsed_internal_pin_command_sha256": summary["strict_internal_pin_command_sha256"],
            "assembly_record_sha256": sha(record_path.read_bytes()),
            "all_2061_raw_test_identities_reconciled": accounting["collection_execution_identity_match"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = Path.cwd().resolve()
    assert git("rev-parse", "--show-toplevel") == str(root)
    assert not git("status", "--porcelain", "--untracked-files=all")
    assert not args.output.exists() and not args.output.resolve().is_relative_to(root)
    sys.path.insert(0, str(root))
    from ou_irregular.week8 import runner, plots
    from ou_irregular.experiments.provenance import utc_now, environment
    run = root / "artifacts/week8_e4"
    evidence = root / "artifacts/week8_evidence/completion"
    before = {str(p.relative_to(run)): sha(p.read_bytes()) for p in run.rglob("*") if p.is_file()}
    checked = runner.verify_run(run, require_completed=True)
    old, new = tree(FREEZE), tree("HEAD")
    changed = {n for n, identity in old.items() if new.get(n) != identity}
    assert changed == {"README.md"}, changed
    assert not git("diff", FREEZE, "HEAD", "--", "tests", "ou_irregular")
    original_zip = evidence / "original_artifact.zip"
    assert original_zip.stat().st_size == 35558419 and sha(original_zip.read_bytes()) == ZIP_SHA
    declared = json.loads((evidence / "original_member_hashes.json").read_text())
    members = {}
    with zipfile.ZipFile(original_zip) as archive:
        for entry in archive.infolist():
            if entry.is_dir():
                continue
            name = PurePosixPath(entry.filename)
            assert not name.is_absolute() and ".." not in name.parts and str(name) not in members
            data = archive.read(entry)
            members[str(name)] = {"sha256": sha(data), "bytes": len(data)}
            if name.parts[0] == "week8_e4":
                target = run / str(PurePosixPath(*name.parts[1:]))
                if str(name) == "week8_e4/run_metadata.json":
                    assert data == (evidence / "original_run_metadata.json").read_bytes()
                    original_metadata = json.loads(data)
                    sealed = checked["metadata"]
                    assert {k for k in original_metadata if original_metadata[k] != sealed[k]} == {"status"}
                    assert set(sealed) - set(original_metadata) == {"completed_utc", "visual_review_sha256", "artifact_sha256"}
                else:
                    assert target.read_bytes() == data, str(name)
            elif name.parts[0] == "week8_e4-figure-renderer":
                target = root / "artifacts" / str(name)
                assert target.read_bytes() == data, str(name)
            else:
                assert name.parts[0] == "evidence"
                target = root / "artifacts/week8_evidence/full_run_41200456847" / str(PurePosixPath(*name.parts[1:]))
                assert target.read_bytes() == data, str(name)
    assert len(members) == 124 and members == declared
    completion = json.loads((evidence / "completion_verification.json").read_text())
    assert completion["passed"] and completion["independent_full_reaudit_passed"]
    independent = json.loads((evidence / "independent_reaudit.json").read_text())
    assert independent["status"] == "passed" and independent["full_completion_audit"] and independent["replay_performed"]
    for key in ("planned_paths", "planned_fit_records", "source_paths_verified", "source_array_hashes_verified",
                "schedule_trace_hashes_verified", "replayed_fits", "replayed_source_paths", "pipeline_errors",
                "retained_failure_records", "point_valid_records", "interval_valid_records", "full_completion_audit",
                "raw_input_sha256", "run_identity_sha256", "source_dependency_files_verified", "tables_verified"):
        assert independent[key] == checked["audit"][key], key
    assert completion["all_36_replot_files_byte_identical"] and completion["sealed_archive_sha256"] == before
    reproduced = {n: sha((evidence / "external_replot/figures" / n).read_bytes()) for n in checked["figures"]["figures"]}
    assert len(reproduced) == 36 and reproduced == checked["figures"]["figures"]
    plots.verify_figures(run, runner.saved_design(run), extra_files=runner.WEEK_FILES, destination=evidence / "external_replot")
    regression = evidence / "final_regression"
    test_summary = json.loads((regression / "regression_summary.json").read_text())
    assert test_summary["status"] == "passed" and test_summary["checkout_remained_clean"]
    assert test_summary["totals"] == {"tests": 2061, "passed": 2055, "skipped": 6, "failures": 0, "errors": 0}
    assert test_summary["unique_execution_identity_count"] == test_summary["all_collected_tests_accounted_for"] == 2061
    assert test_summary["candidate_commit"] == FREEZE
    from scripts.week8_assemble_regression import verified_manifest
    complete_regression_files = verified_manifest(regression)
    manifest = json.loads((regression / "evidence_sha256.json").read_text())
    assert set(complete_regression_files) == set(manifest) | {"evidence_sha256.json"}
    assert (regression / "regression_script_executed.py").read_bytes() == (root / "scripts/week8_release_regression.py").read_bytes()
    from scripts.week8_release_regression import reconcile
    raw_reconciliation = reconcile(regression, json.loads((regression / "test-exit.json").read_text()))
    assert not raw_reconciliation["errors"] and raw_reconciliation["totals"] == test_summary["totals"]
    assert raw_reconciliation["all_collected_tests_accounted_for"] == 2061
    for key, value in raw_reconciliation.items():
        assert test_summary[key] == value, key
    for name, record in manifest.items():
        p = regression / name
        assert p.is_file() and p.stat().st_size == record["bytes"] and sha(p.read_bytes()) == record["sha256"], name
    composite = verify_composite(root, evidence, regression, test_summary, raw_reconciliation)
    reviewed = json.loads((run / "visual_review.json").read_text())
    assert len(reviewed["reviewed_pngs"]) == 18 and "AI" in reviewed["review_kind"]
    assert reviewed["figure_sha256"] == {n: sha((run / "figures" / n).read_bytes()) for n in reviewed["reviewed_pngs"]}
    # Earlier archives are preserved by exact Git object identity; their saved
    # completion states and denominators remain separate from E4.
    earlier = {}
    for name in ("week5_e1", "week6_e2", "week7_followup", "week7_floor", "week8_validation"):
        m = json.loads((root / "artifacts" / name / "run_metadata.json").read_text())
        assert m["status"] == "completed", name
        earlier[name] = "completed; prior Git bytes preserved"
    links = []
    for name in DOCS:
        text = (root / name).read_text()
        assert "PENDING_" not in text, name
        for href in re.findall(r"\[[^\]\n]*\]\(([^)]+)\)", text):
            href = href.strip().split("#", 1)[0]
            if not href or "://" in href or href.startswith("mailto:"):
                continue
            target = ((root / name).parent / unquote(href)).resolve()
            assert target.exists(), (name, href)
            links.append({"document": name, "target": href})
    after = {str(p.relative_to(run)): sha(p.read_bytes()) for p in run.rglob("*") if p.is_file()}
    assert before == after and not git("status", "--porcelain", "--untracked-files=all")
    result = {"passed": True, "checked_utc": utc_now(), "candidate_commit": git("rev-parse", "HEAD"),
              "scientific_freeze": FREEZE, "environment": environment(), "saved_completed_E4_verified": True,
              "original_member_files_verified": 124, "original_zip_sha256": ZIP_SHA,
              "all_prior_git_objects_except_intentionally_updated_readme_identical": True,
              "prior_git_files": len(old), "preserved_git_objects": {n: v for n, v in old.items() if n != "README.md"},
              "scientific_source_and_tests_unchanged": True, "earlier_archives": earlier,
              "regression_evidence_files_verified": len(manifest), "test_counts": test_summary["totals"],
              "composite_regression_provenance": composite,
              "all_36_retained_external_figure_bytes_identical": True, "all_18_actual_AI_reviews_bound": True,
              "documentation_relative_links_checked": links, "archive_unchanged_by_verification": True,
              "clean_checkout": True, "new_full_experiment_executed": False}
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k not in {"preserved_git_objects", "documentation_relative_links_checked"}}, indent=2))


if __name__ == "__main__":
    main()
