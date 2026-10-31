"""Recover the fixed original E4 Actions artifact; never execute a new study.

Run from the clean, full-history scientific freeze using its pinned runtime.
All original bytes remain in the ZIP and in a separate untouched extraction.
The mutable copy changes only through the existing guarded visual seal.
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


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def snapshot(root):
    return {str(p.relative_to(root)): sha(p) for p in sorted(root.rglob("*")) if p.is_file()}


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--zip", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--review", type=Path, required=True)
    args = parser.parse_args()
    from ou_irregular.week8 import runner
    from ou_irregular.experiments.provenance import environment, utc_now
    import numpy as np

    assert subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip() == FREEZE
    assert not subprocess.check_output(["git", "status", "--porcelain"], text=True).strip()
    assert args.zip.stat().st_size == 35558419 and sha(args.zip) == ZIP_SHA
    args.work.mkdir(exist_ok=False)
    original = args.work / "original"
    original.mkdir()
    members = {}
    with zipfile.ZipFile(args.zip) as archive:
        for info in archive.infolist():
            name = PurePosixPath(info.filename)
            assert not name.is_absolute() and ".." not in name.parts
            assert name.parts[0] in {"week8_e4", "week8_e4-figure-renderer", "evidence"}
            assert not stat.S_ISLNK(info.external_attr >> 16)
            if info.is_dir():
                continue
            assert str(name) not in members
            data = archive.read(info)
            target = original / str(name)
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(data)
            members[str(name)] = {"bytes": len(data), "sha256": sha(target)}
    assert len(members) == 124
    untouched = snapshot(original)
    completion = args.work / "completion"
    completion.mkdir()
    write(completion / "original_member_hashes.json", members)
    shutil.copy2(args.zip, completion / "original_artifact.zip")
    write(completion / "runtime.json", {"environment": environment(), "cpu_features": np._core._multiarray_umath.__cpu_features__,
          "cpu_dispatch": np._core._multiarray_umath.__cpu_dispatch__, "workflow_run_id": os.environ.get("GITHUB_RUN_ID"),
          "control_commit": os.environ.get("CONTROL_COMMIT"), "started_utc": utc_now(), "scientific_freeze": FREEZE})
    run = args.work / "week8_e4"
    shutil.copytree(original / "week8_e4", run)
    shutil.copytree(original / "week8_e4-figure-renderer", args.work / "week8_e4-figure-renderer")
    # Saved verifier comes first: it neither samples nor fits.
    checked = runner.verify_run(run)
    assert checked["metadata"]["status"] == "audited_pending_visual_review"
    a = checked["audit"]
    assert a["planned_paths"] == a["source_paths_verified"] == 14000
    assert a["planned_fit_records"] == a["point_valid_records"] == a["interval_valid_records"] == 42000
    assert a["source_array_hashes_verified"] == 28000 and a["schedule_trace_hashes_verified"] == 42000
    assert a["replayed_fits"] == 42 and a["replayed_source_paths"] == 14
    assert a["pipeline_errors"] == a["retained_failure_records"] == 0 and a["full_completion_audit"]
    write(completion / "saved_verification.json", {"passed": True, "checked_utc": utc_now(), "audit_sha256": sha(run / "audit.json"),
          "original_zip_sha256": ZIP_SHA, "member_count": len(members), "status": checked["metadata"]["status"]})
    before = snapshot(run)
    print("Saved original verification PASSED; starting read-only independent full numerical re-audit", flush=True)
    repeat = runner.audit_study(run, replay=True)
    write(completion / "independent_reaudit.json", repeat)
    assert snapshot(run) == before
    for key in ("planned_paths", "planned_fit_records", "source_paths_verified", "source_array_hashes_verified",
                "schedule_trace_hashes_verified", "replayed_fits", "replayed_source_paths", "pipeline_errors",
                "retained_failure_records", "point_valid_records", "interval_valid_records", "full_completion_audit"):
        assert repeat[key] == a[key], key
    review = json.loads(args.review.read_text())
    assert review["figure_sha256"] == {name: sha(run / "figures" / name) for name in checked["figures"]["figures"] if name.endswith(".png")}
    assert review["figure_manifest_sha256"] == sha(run / "figures/figure_manifest.json")
    sealed = runner.seal_run(run, review)
    assert sealed["status"] == "completed"
    runner.verify_run(run, require_completed=True)
    after_seal = snapshot(run)
    assert set(after_seal) - set(before) == {"visual_review.json"}
    assert {n for n in before if before[n] != after_seal[n]} == {"run_metadata.json"}
    # Repeated sealing must be read-only.
    runner.seal_run(run, review)
    assert snapshot(run) == after_seal
    external = args.work / "external_replot"
    print("Guarded seal PASSED; reproducing all saved-table figures externally", flush=True)
    runner.replot(run, external)
    hashes = checked["figures"]["figures"]
    assert len(hashes) == 36
    reproduced = {name: sha(external / "figures" / name) for name in hashes}
    assert reproduced == hashes
    assert snapshot(run) == after_seal and snapshot(original) == untouched
    write(completion / "completion_verification.json", {"passed": True, "checked_utc": utc_now(),
          "scientific_freeze": FREEZE, "original_workflow_run": "41200456847", "original_zip_sha256": ZIP_SHA,
          "source_paths": 14000, "fit_records": 42000, "endpoint_replays": 42,
          "saved_verification_passed": True, "independent_full_reaudit_passed": True,
          "guarded_seal_passed": True, "idempotent_seal_read_only": True,
          "original_extraction_unchanged": True, "only_seal_metadata_and_review_added": True,
          "completed_archive_unchanged_by_replot": True, "all_36_replot_files_byte_identical": True,
          "figure_sha256": reproduced, "sealed_archive_sha256": after_seal,
          "new_full_experiment_executed": False})
    write(completion / "replot_verification.json", {"files": reproduced, "all_byte_identical": True,
          "source_archive_unchanged": True, "external_only": True})
    root = Path.cwd()
    targets = [(run, root / "artifacts/week8_e4"),
               (args.work / "week8_e4-figure-renderer", root / "artifacts/week8_e4-figure-renderer"),
               (original / "evidence", root / "artifacts/week8_evidence/full_run_41200456847"),
               (completion, root / "artifacts/week8_evidence/completion")]
    for source, target in targets:
        assert not target.exists(), target
        shutil.copytree(source, target)
        assert snapshot(source) == snapshot(target)
    shutil.copy2(original / "week8_e4/run_metadata.json", root / "artifacts/week8_evidence/completion/original_run_metadata.json")
    runner.verify_run(root / "artifacts/week8_e4", require_completed=True)
    print("Recovery, full numerical re-audit, actual-review seal and 36-file reproduction PASSED", flush=True)


if __name__ == "__main__":
    main()
