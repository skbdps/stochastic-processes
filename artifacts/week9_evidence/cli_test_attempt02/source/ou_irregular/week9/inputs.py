"""Pinned, read-only Week 9 source access without historical runtime replay.

The catalog is itself bound to this module and each entry to a blob in the
published baseline. This verifies saved bytes and their declared accounting;
it does not repeat the original numerical audits or certify a new CPU runtime.
"""
from dataclasses import dataclass
import gzip
import hashlib
from io import BytesIO
import json
import math
from pathlib import Path, PurePosixPath
from types import MappingProxyType
import subprocess

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
BASELINE_COMMIT = "b6edb85d1e7a21ed19fb54b5c71213ed04703309"
CONFIG_SHA256 = "ba4bc13e2734de8b5fa30a05b31cd735e861e9ee1eaa5807a3b234d2c6e8c3a9"
DEFAULT_CONFIG = "configs/week9_publication.json"
STUDIES = {
    "e1": "artifacts/week5_e1", "e2": "artifacts/week6_e2",
    "e3_selected": "artifacts/week7_selected",
    "e3_followup": "artifacts/week7_followup",
    "e3_floor": "artifacts/week7_floor", "e4": "artifacts/week8_e4",
}
METHODS = ("exact", "pfml", "euler")
PARAMETERS = ("theta", "mu", "sigma")


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _sha(data):
    return hashlib.sha256(data).hexdigest()


def _canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _relative(name):
    _require(isinstance(name, str) and bool(name), "source path must be a nonempty string")
    p = PurePosixPath(name)
    _require(bool(p.parts) and not p.is_absolute() and ".." not in p.parts and "\\" not in name,
             "source path must stay within the repository")
    _require(str(p) == name and ".git" not in p.parts, "source path must be canonical")
    return p


def _path(root, name):
    p = _relative(name)
    current = root
    for part in p.parts:
        current = current / part
        _require(not current.is_symlink(), "source symlink is forbidden: " + name)
    _require(current.is_file(), "missing source file: " + name)
    return current


def _file_names(root, relative):
    """Also reject directory symlinks, which rglob must never follow."""
    directory = root / _relative(relative)
    _require(directory.is_dir() and not directory.is_symlink(), "missing/linked source directory: " + relative)
    result = set()
    for p in directory.rglob("*"):
        _require(not p.is_symlink(), "source symlink is forbidden: " + str(p))
        if p.is_file():
            result.add(p.relative_to(root).as_posix())
    return result


def _git_tree(root, prefixes):
    try:
        result = subprocess.run(["git", "ls-tree", "-rz", BASELINE_COMMIT, "--", *prefixes],
                                cwd=root, capture_output=True, check=True)
    except (OSError, subprocess.CalledProcessError) as exc:
        raise ValueError("full-history checkout with the pinned baseline commit is required") from exc
    tree = {}
    for item in result.stdout.split(b"\0"):
        if item:
            header, name = item.split(b"\t", 1)
            mode, kind, blob = header.decode().split()
            _require(kind == "blob" and mode in ("100644", "100755"), "baseline source is not a regular Git blob")
            tree[name.decode()] = (mode, blob)
    return tree


@dataclass(frozen=True)
class SourceBundle:
    root: Path
    config: dict
    config_path: Path
    _entries: dict
    _config_canonical: str

    def _check_config(self):
        _require(_canonical(self.config) == self._config_canonical, "in-memory source catalog was modified")
        _require(not self.config_path.is_symlink() and self.config_path.is_file(), "source catalog is missing or linked")
        _require(_sha(self.config_path.read_bytes()) == CONFIG_SHA256, "pinned source catalog digest mismatch")

    def study_dir(self, study):
        self._check_config()
        _require(study in STUDIES, "unknown study: " + str(study))
        return self.root / self.config["studies"][study]["path"]

    def _read_file(self, relative):
        _require(relative in self._entries, "file is outside the pinned source catalog: " + str(relative))
        data = _path(self.root, relative).read_bytes()
        entry = self._entries[relative]
        _require(len(data) == entry["size_bytes"] and _sha(data) == entry["sha256"],
                 "source digest mismatch: " + relative)
        blob = hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()
        _require(blob == entry["git_blob"], "source Git blob mismatch: " + relative)
        return data

    def read_bytes(self, relative):
        """Read and hash the exact bytes returned, avoiding a hash/read gap."""
        self._check_config()
        return self._read_file(relative)

    def json(self, relative):
        return json.loads(self.read_bytes(relative))

    def table(self, study, filename):
        self.study_dir(study)
        _relative(filename)
        _require("/" not in filename and filename.endswith(".csv"), "table must be a direct CSV in its study")
        relative = self.config["studies"][study]["path"] + "/" + filename
        return pd.read_csv(BytesIO(self.read_bytes(relative)), float_precision="round_trip")

    def snapshot(self):
        """Verify all cataloged source bytes and reject untracked archive additions."""
        self._check_config()
        for group in (*self.config["studies"].values(), *self.config["additional_sources"].values()):
            _require(_file_names(self.root, group["path"]) == set(group["files"]),
                     "source archive membership changed: " + group["path"])
        result = {name: _sha(self._read_file(name)) for name in sorted(self._entries)}
        self._check_config()
        return result

    def f1_records(self):
        """Read every original E1 fit; never regenerate an observation or fit."""
        study = self.config["studies"]["e1"]
        design = self.json(study["path"] + "/resolved_design.json")
        replications = design["replications"]
        cells = {c["cell_id"]: c for c in design["cells"]}
        rows, keys = [], set()
        for cell_id, cell in sorted(cells.items()):
            prefix = study["path"] + "/cells/" + cell_id
            manifest = self.json(prefix + "/cell_manifest.json")
            _require(manifest["cell"] == cell and manifest["status"] == "complete"
                     and manifest["completed_replications"] == replications,
                     "F1 cell manifest is incomplete or has changed identity")
            start = 0
            for part in manifest["partitions"]:
                _require(part["start"] == start and start < part["stop"] <= replications,
                         "F1 partition ranges are not a complete contiguous sequence")
                name = part["filename"]
                _relative(name)
                _require("/" not in name, "F1 partition filename must be direct")
                relative = prefix + "/" + name
                raw = self.read_bytes(relative)
                _require(_sha(raw) == part["sha256"], "F1 partition differs from saved cell manifest")
                records = [json.loads(line) for line in gzip.decompress(raw).splitlines() if line.strip()]
                _require(len(records) == part["rows"] == (part["stop"]-start)*len(METHODS),
                         "F1 partition row accounting mismatch")
                for row in records:
                    rep, method = row.get("rep"), row.get("estimator")
                    _require(type(rep) is int and start <= rep < part["stop"] and method in METHODS,
                             "F1 replication/method key is invalid")
                    key = (cell_id, rep, method)
                    _require(row.get("cell_id") == cell_id and key not in keys,
                             "F1 cell key mismatch or duplicate source fit")
                    _require(all(row.get(k) == value for k, value in cell.items()), "F1 design metadata mismatch")
                    _require(row.get("source_seed") == design["seed"], "F1 source seed mismatch")
                    _require(type(row.get("point_valid")) is bool and type(row.get("interval_valid")) is bool,
                             "F1 fit eligibility must be explicit Boolean values")
                    _require(not row.get("pipeline_error") and row.get("source_valid") is True,
                             "F1 source/pipeline failure is unexpected in this complete archive")
                    if row["point_valid"]:
                        _require(all(isinstance(row.get(p), (int, float)) and math.isfinite(row[p]) for p in PARAMETERS),
                                 "F1 valid fit has a nonfinite point estimate")
                    keys.add(key)
                    rows.append({**row, "_source_path": relative, "_source_sha256": part["sha256"]})
                start = part["stop"]
            _require(start == replications, "F1 cell is missing planned replications")
        expected = {(c, r, m) for c in cells for r in range(replications) for m in METHODS}
        _require(keys == expected and len(rows) == 42000, "F1 full planned key set is incomplete")
        return rows


def _check_bound_mapping(bundle, prefix, mapping, snapshots):
    _require(isinstance(mapping, dict), "saved binding map must be an object")
    for name, value in mapping.items():
        _relative(name)
        relative = prefix + "/" + name
        expected = value["sha256"] if isinstance(value, dict) else value
        _require(snapshots.get(relative) == expected, "saved manifest binding mismatch: " + relative)


def _verify_study(bundle, study, snapshots):
    spec = bundle.config["studies"][study]
    prefix, expected = spec["path"], spec["expected"]
    metadata = bundle.json(prefix + "/run_metadata.json")
    audit = bundle.json(prefix + "/audit.json")
    _require(audit.get("status") == "passed", "source audit has not passed: " + study)
    _require(metadata["audit_sha256"] == snapshots[prefix + "/audit.json"], "saved audit digest mismatch")
    for name, sha in metadata["provenance"]["source_files"].items():
        _require(snapshots.get(name) == sha, "original scientific dependency changed: " + name)
    if study == "e3_selected":
        _require(metadata["status"] == "completed_diagnostic_evidence", "selected diagnostics are incomplete")
        _require(audit["planned_sources"] == audit["retained_sources"] == audit["evaluated_sources"] == 53,
                 "selected source accounting changed")
        _require(audit["original_fit_records"] == audit["verified_fit_replays"] == 159 and not audit["failures"],
                 "selected diagnostic replay accounting changed")
        _require(not audit["headline_results_replaced"], "selected diagnostics must not replace headline results")
        _check_bound_mapping(bundle, prefix + "/sources", metadata["source_files"], snapshots)
        _require(metadata["source_files"] == audit["source_files"], "selected source maps disagree")
        return
    _require(metadata["status"] == "completed" and metadata["pipeline_errors"] == audit["pipeline_errors"] == 0,
             "source archive is not completed without pipeline errors: " + study)
    _require(metadata["visual_review_sha256"] == snapshots[prefix + "/visual_review.json"], "saved visual review digest mismatch")
    _check_bound_mapping(bundle, prefix, metadata["artifact_sha256"], snapshots)
    _check_bound_mapping(bundle, prefix, audit["audited_artifact_sha256"], snapshots)
    for key in ("planned_paths", "planned_fit_records", "point_valid_records", "interval_valid_records"):
        wanted = expected[{"planned_paths": "source_realizations", "planned_fit_records": "fit_records"}.get(key, key)]
        _require(audit[key] == wanted, "saved audit accounting mismatch: " + study + "/" + key)
    if study == "e3_floor":
        _require(audit["planned_independent_coupling_groups"] == expected["independent_coupling_groups"] == 8000,
                 "floor study independent-group accounting mismatch")
    design = bundle.json(prefix + "/resolved_design.json")
    _require(design == metadata["provenance"]["design"], "source resolved design/provenance mismatch")
    _require(tuple(design["estimators"]) == METHODS, "source estimator identities changed")
    cells = {c["cell_id"]: c for c in design["cells"]}
    _require(len(cells) == len(design["cells"]), "duplicate source design cell")
    _require(len(cells)*design["replications"] == expected["source_realizations"], "source path accounting changed")
    manifest = bundle.json(prefix + "/table_manifest.json")
    _check_bound_mapping(bundle, prefix, manifest, snapshots)
    summary = bundle.table(study, "parameter_summary.csv")
    actual = list(summary[["cell_id", "estimator", "parameter"]].itertuples(index=False, name=None))
    planned = {(c, m, p) for c in cells for m in METHODS for p in PARAMETERS}
    _require(len(actual) == len(planned) and set(actual) == planned, "source parameter summary is missing or duplicates planned keys")
    for row in summary.itertuples():
        cell = cells[row.cell_id]
        _require(all(getattr(row, k) == v for k, v in cell.items() if k in summary.columns)
                 and row.truth == cell["true_" + row.parameter], "source summary design coordinates changed")
        _require(row.n_planned == design["replications"] and 0 <= row.n_interval_valid <= row.n_point_valid <= row.n_planned,
                 "source summary eligibility accounting changed")
        _require((row.parameter == "mu") == (row.practical_classification == "not_applicable_zero_truth"),
                 "mu zero-truth classification scope changed")
    for metric, wanted in (("n_point_valid", expected["point_valid_records"]), ("n_interval_valid", expected["interval_valid_records"])):
        _require(int(summary.loc[summary.parameter == "theta", metric].sum()) == wanted,
                 "source valid-record totals changed")
    figures = bundle.json(prefix + "/figures/figure_manifest.json")
    _check_bound_mapping(bundle, prefix + "/figures", figures["figures"], snapshots)


def load_sources(config_path=None, *, root=None):
    """Verify the pinned source catalog, original Git membership and saved seals."""
    root = Path(root).resolve() if root is not None else ROOT
    config_path = Path(config_path) if config_path is not None else root / DEFAULT_CONFIG
    if not config_path.is_absolute():
        config_path = root / config_path
    _require(not config_path.is_symlink() and config_path.is_file(), "source catalog is missing or linked")
    raw = config_path.read_bytes()
    _require(_sha(raw) == CONFIG_SHA256, "pinned source catalog digest mismatch")
    config = json.loads(raw)
    _require(config["baseline_commit"] == BASELINE_COMMIT and config["schema_version"] == 1,
             "source baseline identity changed")
    _require(set(config["studies"]) == set(STUDIES), "source study membership changed")
    entries = dict(config["reference_files"])
    prefixes = list(entries)
    for study, spec in config["studies"].items():
        _require(spec["path"] == STUDIES[study], "source study path changed")
    for spec in (*config["studies"].values(), *config["additional_sources"].values()):
        _relative(spec["path"])
        prefixes.append(spec["path"])
        for name, entry in spec["files"].items():
            _require(name.startswith(spec["path"] + "/") and name not in entries,
                     "source file catalog has an unsafe or duplicate path")
            entries[name] = entry
    tree = _git_tree(root, prefixes)
    _require(set(tree) == set(entries), "source catalog differs from baseline Git tree membership")
    for name, entry in entries.items():
        _relative(name)
        _require(tree[name] == (entry["git_mode"], entry["git_blob"]), "source entry differs from baseline Git blob: " + name)
    bundle = SourceBundle(root, config, config_path,
                          MappingProxyType({k: MappingProxyType(v) for k, v in entries.items()}),
                          _canonical(config))
    snapshots = bundle.snapshot()
    for study in STUDIES:
        _verify_study(bundle, study, snapshots)
    selected = config["additional_sources"]["selected_figures"]["path"]
    figure_manifest = bundle.json(selected + "/figure_manifest.json")
    _check_bound_mapping(bundle, selected, figure_manifest["figures"], snapshots)
    return bundle


def guard_output(bundle, path):
    """Require a fresh directory destination disjoint from every protected input."""
    bundle._check_config()
    target = Path(path).expanduser()
    if not target.is_absolute():
        target = bundle.root / target
    _require(".." not in target.parts, "output path must not contain parent traversal")
    _require(not target.exists() and not target.is_symlink(), "output destination must be fresh")
    target = target.resolve()
    _require(target != bundle.root and target != bundle.root / ".git"
             and bundle.root / ".git" not in target.parents, "output destination overlaps repository metadata")
    protected = [bundle.root / spec["path"] for spec in
                 (*bundle.config["studies"].values(), *bundle.config["additional_sources"].values())]
    protected.extend(bundle.root / name for name in bundle.config["reference_files"])
    protected.append(bundle.config_path)
    for source in protected:
        source = source.resolve()
        _require(target != source and source not in target.parents and target not in source.parents,
                 "output destination overlaps protected source evidence")
    return target
