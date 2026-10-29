"""Restore a documented local freeze commit object for a fresh-clone audit.

The GitHub connector can publish the same tree with different commit metadata.
This helper verifies the archived object and the published tree before importing
the missing commit object. It changes no branch, working file or run artifact.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess


ROOT = Path(__file__).resolve().parents[1]


def restore(mapping_path, run_dir, *, root=ROOT):
    root = Path(root)
    mapping = json.loads(Path(mapping_path).read_text())
    metadata = json.loads((Path(run_dir) / "run_metadata.json").read_text())
    expected = mapping["local_recorded_commit"]
    if metadata["provenance"]["code_commit"] != expected:
        raise ValueError("mapping does not identify the recorded run commit")
    relative = Path(mapping["original_commit_object_file"])
    if relative.is_absolute() or ".." in relative.parts:
        raise ValueError("archived commit path must stay within the repository")
    raw = (root / relative).read_bytes()
    if hashlib.sha256(raw).hexdigest() != mapping["original_commit_object_sha256"]:
        raise ValueError("archived commit file digest mismatch")
    actual = hashlib.sha1(b"commit " + str(len(raw)).encode() + b"\0" + raw).hexdigest()
    if actual != expected:
        raise ValueError("archived Git commit identity mismatch")
    tree = mapping["identical_tree"]
    if raw.splitlines()[0] != f"tree {tree}".encode():
        raise ValueError("archived commit does not contain the declared tree")
    published_tree = subprocess.check_output(
        ["git", "rev-parse", mapping["published_freeze_commit"] + "^{tree}"],
        cwd=root, text=True).strip()
    if published_tree != tree:
        raise ValueError("published and recorded freeze trees differ")
    imported = subprocess.check_output(["git", "hash-object", "-w", "-t", "commit", "--stdin"],
                                       input=raw, cwd=root).decode().strip()
    if imported != expected:
        raise ValueError("Git object import returned a different identity")
    return {"recorded_commit": imported, "published_freeze_commit": mapping["published_freeze_commit"],
            "identical_tree": tree, "refs_changed": False, "run_artifacts_changed": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mapping", type=Path, required=True)
    parser.add_argument("--run-dir", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(restore(args.mapping, args.run_dir), sort_keys=True))


if __name__ == "__main__":
    main()
