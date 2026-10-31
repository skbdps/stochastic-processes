"""Publish the retained successful scientific stages of recovery run 41214670935.

No source reconstruction, fitting, or rendering is performed. The failed
packaging path is replaced by exact saved-byte checks and native Git staging.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import shutil
import stat
import subprocess
import zipfile

FREEZE = "0ff57919f5d8345f52d8689732ff94448ca3ca1d"
ZIP_SHA = "e9a780206925c4085ff4bb47cf988c7ccec66ce40db038f3fdef82703fe52513"


def digest(data):
    return hashlib.sha256(data).hexdigest()


def snap(root):
    return {str(p.relative_to(root)): digest(p.read_bytes()) for p in sorted(root.rglob("*")) if p.is_file()}


def write(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def extract(path, target):
    target.mkdir(exist_ok=False)
    seen = set()
    with zipfile.ZipFile(path) as archive:
        for info in archive.infolist():
            name = PurePosixPath(info.filename)
            assert not name.is_absolute() and ".." not in name.parts and not stat.S_ISLNK(info.external_attr >> 16)
            if info.is_dir():
                continue
            assert str(name) not in seen
            seen.add(str(name))
            p = target / str(name)
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_bytes(archive.read(info))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--transport", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--control", type=Path, required=True)
    args = parser.parse_args()
    from ou_irregular.week8 import runner, plots
    from ou_irregular.experiments.provenance import utc_now, environment
    import numpy as np
    assert subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip() == FREEZE
    assert not subprocess.check_output(["git", "status", "--porcelain"], text=True).strip()
    metadata = json.loads(args.metadata.read_text())
    assert metadata["workflow_run"]["id"] == 41214670935
    assert metadata["workflow_run"]["head_sha"] == "2ee68378cc2c651784301fc59640ac554d00d8cb"
    assert digest((args.control / "release_inputs/recovery_script_executed.py").read_bytes()) == "88ae29b7b2568a2ecb40f419867d4d762045e36b3b0909e30d58af33c08265e9"
    assert metadata["name"] == "original-e4-recovery-41214670935" and not metadata["expired"]
    assert metadata["digest"] == "sha256:" + digest(args.transport.read_bytes())
    extract(args.transport, args.work)
    matches = list(args.work.rglob("original_member_hashes.json"))
    assert len(matches) == 1
    completion = matches[0].parent
    base = completion.parent
    original = base / "original"
    run = base / "week8_e4"
    renderer = base / "week8_e4-figure-renderer"
    original_zip = completion / "original_artifact.zip"
    assert original_zip.stat().st_size == 35558419 and digest(original_zip.read_bytes()) == ZIP_SHA
    fresh_original = args.work / "fresh_original_readback"
    extract(original_zip, fresh_original)
    fresh = snap(fresh_original)
    assert len(fresh) == 124 and fresh == snap(original)
    records = json.loads(matches[0].read_text())
    assert records == {n: {"sha256": h, "bytes": (fresh_original / n).stat().st_size} for n, h in fresh.items()}
    before = snap(run)
    raw = snap(original / "week8_e4")
    assert set(before) - set(raw) == {"visual_review.json"}
    assert {n for n in raw if raw[n] != before[n]} == {"run_metadata.json"}
    assert snap(renderer) == snap(original / "week8_e4-figure-renderer")
    checked = runner.verify_run(run, require_completed=True)
    a = checked["audit"]
    independent = json.loads((completion / "independent_reaudit.json").read_text())
    assert independent["status"] == "passed" and independent["full_completion_audit"] and independent["replay_performed"]
    # Repeat report must agree on all scientific fields; only lifecycle-bound
    # root input hashes may differ between separate audit invocations.
    for key in ("planned_paths", "planned_fit_records", "source_paths_verified", "source_array_hashes_verified",
                "schedule_trace_hashes_verified", "replayed_fits", "replayed_source_paths", "pipeline_errors",
                "retained_failure_records", "point_valid_records", "interval_valid_records", "full_completion_audit",
                "raw_input_sha256", "run_identity_sha256", "source_dependency_files_verified", "tables_verified"):
        assert independent[key] == a[key], key
    assert a["planned_paths"] == 14000 and a["planned_fit_records"] == 42000
    review = json.loads((run / "visual_review.json").read_text())
    assert review == json.loads((args.control / "release_inputs/week8_visual_review.json").read_text())
    runner.seal_run(run, review)
    assert snap(run) == before
    # Recovery rendered successfully before its path-lookup packaging failure.
    # Check its actual output, without another render or numerical audit.
    external = base / "external_replot/figures"
    wanted = checked["figures"]["figures"]
    got = {name: digest((external / name).read_bytes()) for name in wanted}
    assert len(got) == 36 and got == wanted
    plots.verify_figures(run, runner.saved_design(run), extra_files=runner.WEEK_FILES, destination=base / "external_replot")
    assert snap(run) == before and snap(original) == fresh
    write(completion / "completion_verification.json", {"passed": True, "checked_utc": utc_now(),
          "scientific_freeze": FREEZE, "original_workflow_run": "41200456847", "original_zip_sha256": ZIP_SHA,
          "numerical_reaudit_workflow_run": "41214670935", "publication_workflow_run": os.environ["GITHUB_RUN_ID"],
          "source_paths": 14000, "fit_records": 42000, "endpoint_replays": 42,
          "saved_verification_passed": True, "independent_full_reaudit_passed": True,
          "guarded_seal_passed": True, "idempotent_seal_read_only": True,
          "original_extraction_unchanged": True, "only_seal_metadata_and_review_added": True,
          "completed_archive_unchanged_by_replot": True, "all_36_replot_files_byte_identical": True,
          "figure_sha256": got, "sealed_archive_sha256": before,
          "new_full_experiment_executed": False,
          "packaging_disposition": "Recovery rendered externally, then used the wrong lookup path (external/name versus external/figures/name). This saved-byte continuation verifies all36 correct paths. Ignored artifact paths require explicit forced staging. Neither issue changes scientific source or outcomes.",
          "reaudit_runtime_capture_limitation": "Recovery host captured package/thread environment and NumPy CPU features/dispatch, but not full lscpu/NumPy build output; original full invocation and final-validation runtime records are separate evidence."})
    write(completion / "replot_verification.json", {"files": got, "all_byte_identical": True,
          "source_archive_unchanged": True, "external_only": True, "generated_by_recovery_run": "41214670935",
          "checked_by_publication_run": os.environ["GITHUB_RUN_ID"]})
    shutil.copytree(base / "external_replot", completion / "external_replot")
    shutil.copytree(base / "external_replot-figure-renderer", completion / "external_replot-figure-renderer")
    shutil.copy2(args.metadata, completion / "recovery_artifact_metadata.json")
    shutil.copy2(original / "week8_e4/run_metadata.json", completion / "original_run_metadata.json")
    for name in ("original_publication_failure.json", "recovery_packaging_failure.json", "original_workflow_reference.json", "recovery_script_executed.py"):
        shutil.copy2(args.control / "release_inputs" / name, completion / name)
    for name in ("week8_recover_original.py", "week8_publish_recovered.py"):
        shutil.copy2(args.control / "scripts" / name, completion / name)
    for name in ("week8-recover-original.yml", "week8-publish-recovered.yml"):
        shutil.copy2(args.control / ".github/workflows" / name, completion / name)
    logs = list(args.work.rglob("recovery.txt"))
    assert len(logs) == 1
    shutil.copy2(logs[0], completion / "recovery.txt")
    write(completion / "publication_runtime.json", {"environment": environment(), "cpu_features": np._core._multiarray_umath.__cpu_features__,
          "cpu_dispatch": np._core._multiarray_umath.__cpu_dispatch__, "workflow_run_id": os.environ["GITHUB_RUN_ID"],
          "control_commit": os.environ["GITHUB_SHA"]})
    root = Path.cwd()
    targets = [(run, root / "artifacts/week8_e4"), (renderer, root / "artifacts/week8_e4-figure-renderer"),
               (original / "evidence", root / "artifacts/week8_evidence/full_run_41200456847"),
               (completion, root / "artifacts/week8_evidence/completion")]
    inventory = {}
    for source, target in targets:
        assert not target.exists(), target
        shutil.copytree(source, target)
        assert snap(source) == snap(target)
        inventory.update({str(target.relative_to(root) / n): h for n, h in snap(target).items()})
    runner.verify_run(root / "artifacts/week8_e4", require_completed=True)
    write(args.work / "publication_inventory.json", inventory)
    print("Verified original bytes, retained full re-audit, seal, all36 external figures; ready for native Git staging", flush=True)


if __name__ == "__main__":
    main()
