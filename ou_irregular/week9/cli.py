"""Build, verify and seal a derivative publication package from saved evidence.

The historical scientific verifiers and numerical services are deliberately
not invoked. A separate pinned data contract verifies their archived evidence.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import gzip
import hashlib
from importlib.metadata import distributions, version
import json
import os
from pathlib import Path
import platform
import re
import subprocess
from urllib.parse import urlsplit

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .inputs import load_sources, guard_output
from .plots import STEMS, build_figures, figure_bytes
from .tables import build_tables, build_claims


ROOT = Path(__file__).resolve().parents[2]
CONTRACT = "week9_verified_saved_evidence_publication_v1"
NOTES = ("week9_protocol.md", "week9_benchmarks.md", "week9_figure_plan.md",
         "week9_literature_comparison.md", "week9_reproduction.md")
EXCLUDED = {"manifest.json", "package_state.json", "visual_review.json"}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def plain(value):
    if isinstance(value, dict):
        return {str(k): plain(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [plain(v) for v in value]
    if isinstance(value, np.generic):
        return plain(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def json_bytes(value):
    return (json.dumps(plain(value), sort_keys=True, separators=(",", ":"),
                       ensure_ascii=False, allow_nan=False)+"\n").encode()


def sha(payload):
    return hashlib.sha256(payload).hexdigest()


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024*1024), b""):
            h.update(block)
    return h.hexdigest()


def write(path, payload):
    """Every build target is fresh; an existing different byte is never replaced."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        require(path.read_bytes() == payload, "refusing to replace different publication evidence: "+str(path))
        return
    with path.open("xb") as stream:
        stream.write(payload)
        stream.flush()
        os.fsync(stream.fileno())


def read_json(path):
    return json.loads(Path(path).read_text())


def csv_bytes(frame):
    return frame.to_csv(index=False, float_format="%.17g", lineterminator="\n").encode()


def table_payload(name, frame):
    payload = csv_bytes(frame)
    if name == "f1_ecdf.csv":
        return name+".gz", gzip.compress(payload, compresslevel=6, mtime=0)
    return name, payload


def renderer_identity(root):
    paths = [*sorted((root/"ou_irregular/week9").glob("*.py")),
             root/"configs/week9_publication.json", root/"requirements-week9.txt", root/"requirements-week9-lock.txt",
             *[root/"docs"/name for name in NOTES]]
    require(all(p.is_file() and not p.is_symlink() for p in paths), "renderer source/notes missing or linked")
    return {p.relative_to(root).as_posix(): digest(p) for p in paths}


def runtime():
    return {"python": platform.python_version(), "implementation": platform.python_implementation(),
            "system": platform.system(), "machine": platform.machine(),
            "packages": {name: version(name) for name in ("numpy", "scipy", "pandas", "matplotlib")},
            "package_closure": dict(sorted((d.metadata["Name"], d.version) for d in distributions())),
            "thread_environment": {name: os.environ.get(name) for name in
                                   ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")},
            "scope": "saved-data parsing, aggregation and new figure rendering; not historical numerical replay"}


def committed_renderer(root, identity, commit=None):
    """Bind the freeze to committed code without putting moving HEAD in figures."""
    if commit is None:
        commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip()
    require(re.fullmatch(r"[0-9a-f]{40}", commit) is not None, "invalid renderer commit")
    for relative, wanted in identity.items():
        result = subprocess.run(["git", "show", f"{commit}:{relative}"], cwd=root, capture_output=True)
        require(result.returncode == 0 and sha(result.stdout) == wanted,
                "renderer must be committed before finalization: "+relative)
    return commit


def publication_note(root, name, baseline_commit):
    """Keep new sibling-note links local; make old evidence links permanent."""
    source = root/"docs"/name
    text = source.read_text()
    def replace(match):
        target = match.group(1)
        parsed = urlsplit(target)
        if parsed.scheme or parsed.netloc or not parsed.path:
            return match.group(0)
        local = (source.parent/parsed.path).resolve()
        require(local.is_relative_to(root.resolve()) and local.is_file(), "note has an unresolved local link: "+target)
        if local.parent == root/"docs" and local.name in NOTES:
            rewritten = local.name
        else:
            relative = local.relative_to(root).as_posix()
            rewritten = f"https://github.com/skbdps/stochastic-processes/blob/{baseline_commit}/{relative}"
        if parsed.fragment:
            rewritten += "#"+parsed.fragment
        return "]("+rewritten+")"
    return re.sub(r"\]\(([^)]+)\)", replace, text).encode()


def captions(new_figures, copied):
    texts = {
        "F1_distributions_main": "F1. Empirical cumulative distributions of all eligible saved E1 estimates at nominal CV 0, 1 and 3. All three parameters are estimated jointly; truth is (1,0,0.5), q=0.5, n=1000 and R=2000 per cell. Dotted vertical lines denote truth. Counts show point-valid/planned fits; no confidence bands are claimed. Regular-grid exact/PFML estimates coincide, so their ECDFs overlap; line styles distinguish coincident methods. Selection is post-result presentation. All seven CV settings are supplied in supporting figures and f1_ecdf.csv.gz.",
        "S_F1_distributions_cv0_to_1": "Supporting F1, first page. Complete saved E1 ECDFs at nominal CV 0, 0.25, 0.5 and 1, with the same definitions as F1. Every eligible estimate is retained, including extremes.",
        "S_F1_distributions_cv1p5_to_3": "Supporting F1, second page. Complete saved E1 ECDFs at nominal CV 1.5, 2 and 3. These input Gamma gap laws have the registered floor; realized spans are not held fixed.",
        "F8_matched_mean": "F8. Exploratory E4 matched-marginal MMPP/IID comparisons for mu across r=1,5,20 and all methods. Ratios are MMPP/IID RMSE, differences are MMPP minus IID operational coverage in percentage points. Bars are pointwise 95% Monte Carlo intervals using independent-cell formulas; RMSE-ratio limits use the log-delta method. n=1000, q=0.5 and R=2000 per independent source cell. The contrast includes realized-span consequences; PFML/Euler Wald covariance is model-based. Displaying mu prominently is post-result presentation; all parameters remain supplied.",
        "S_calibration_e4_mean": "Supporting E4 calibration. Empirical SD and RMS reported SE for mu on the identical interval-valid subsets, with the same y-axis limits. These saved descriptive quantities do not have a derived Monte Carlo interval here. All methods, r values and both matched arms are retained. At r=20 both arms use the declared 5e-7 floor; MMPP labels denote floored MMPP-derived gaps. They do not establish a long-run covariance theorem.",
    }
    for n in (250, 1000, 4000):
        texts[f"F7_classification_n{n}"] = f"F7, n={n}. Complete registered E2 pointwise classifications over the 4 by 5 nominal q/CV grid for each method and theta/sigma. W means within the broad practical criteria; B:b, B:c and B:bc denote supported bias, coverage or both breaches; U remains unresolved. Classification uses unrounded bias Monte Carlo intervals and operational-coverage Wilson bounds. No interpolation, continuous boundary, monotonicity or simultaneous guarantee is asserted. The full 360 eligible E2 rows are in T1_e2.csv; mu's zero-truth relative-bias criterion is inapplicable."
        texts[f"S_calibration_e2_n{n}"] = f"Supporting E2 calibration, n={n}. Complete-grid empirical SD divided by RMS reported SE on the same interval-valid subset. Cell labels show the ratio; the common color scale is symmetric in log2 and saturates outside [0.5,2], without clipping the numerical labels. This is descriptive calibration, not a robust variance formula or a coverage explanation by itself. All counts and unrounded values are in calibration.csv."
    lines = ["# Week 9 publication captions", "", "## Newly rendered figures", ""]
    for stem in new_figures:
        require(stem in texts, "new figure lacks a caption")
        lines.extend([f"### {stem}", "", texts[stem], ""])
    lines.extend(["## Byte-identical source figures", "",
                  "These figures retain their original titles, definitions, uncertainty, source manifests and historical review identities. Their bytes are copied unchanged. They are supporting evidence and the existing F2–F6 figure families; each original review is copied beside the source figures. The source study's sealed report defines its original inspection scope.", ""])
    for item in copied:
        if item["path"].endswith(".png"):
            lines.append(f"- `{item['path']}` ← `{item['source_path']}` ({item['study']}); original review: `{item.get('original_review_path', 'see source manifest')}`.")
    return "\n".join(lines)+"\n"


def readable_claims(claims):
    lines = ["# Week 9 claims ledger", "", "Generated from the pinned saved evidence. Numerical wording follows the source row selections and unrounded values in claims.json; selected comparisons are pointwise.", ""]
    for claim in claims:
        lines.extend([f"## {claim['claim_id']} [{claim['status']}]", "", claim["text"], ""])
        for field in ("qualifications", "uncertainty", "comparison", "independent_unit"):
            if claim.get(field):
                lines.extend([f"{field.replace('_', ' ').capitalize()}: {claim[field]}", ""])
    return "\n".join(lines)


def claims_csv_bytes(claims):
    rows = [{k: json.dumps(plain(v), sort_keys=True, ensure_ascii=False, separators=(",", ":"))
             if isinstance(v, (dict, list, tuple)) else v for k,v in c.items()} for c in claims]
    columns = sorted({key for row in rows for key in row})
    return csv_bytes(pd.DataFrame(rows, columns=columns))


def readable_t1(tables):
    lines = ["# T1. Practical classifications within the investigated designs", "",
             "The registered criteria concern theta/sigma: operational coverage below 80% or absolute relative bias above 10%, assessed with their Monte Carlo intervals. W/B/U are pointwise classifications; W does not certify nominal 95% coverage. No continuous or maximum-safe-CV boundary is inferred. Mu's zero-truth relative-bias criterion is inapplicable.", "",
             "## Complete E2 group counts", "", "Each row covers all 20 q/CV designs at the stated sample size. Full values, uncertainty and source bindings are in T1_e2.csv (360 rows). The independent E1 inventory remains in T1_e1.csv (42 rows).", "",
             "| n | Method | Parameter | Within | Breach | Unresolved | Bias only | Coverage only | Both |",
             "| ---: | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |"]
    for row in tables["T1_counts.csv"].itertuples():
        lines.append(f"| {int(row.n)} | {row.estimator} | {row.parameter} | {row.within_criteria} | {row.supported_breach} | {row.unresolved} | {row.bias_only} | {row.coverage_only} | {row.bias_and_coverage} |")
    lines.extend(["", "## All six unresolved E1/E2 classifications", "",
                  "Intervals shown here are 95% Monte Carlo intervals for relative bias and 95% Wilson intervals for operational coverage, in percent. Classifications use unrounded values. Overlapping E1/E2 designs retain their independent study identities.", "",
                  "| Study | n | q | CV | Method | Parameter | Relative bias [MC] (%) | Coverage [Wilson] (%) |",
                  "| --- | ---: | ---: | ---: | --- | --- | --- | --- |"])
    for row in tables["unresolved.csv"].itertuples():
        bias = f"{100*row.relative_bias:.3f} [{100*row.relative_bias_mc_lower:.3f}, {100*row.relative_bias_mc_upper:.3f}]"
        coverage = f"{100*row.coverage_operational:.2f} [{100*row.operational_wilson_lower:.3f}, {100*row.operational_wilson_upper:.3f}]"
        lines.append(f"| {row.study} | {int(row.n)} | {row.theta_mean_gap:g} | {row.cv:g} | {row.estimator} | {row.parameter} | {bias} | {coverage} |")
    lines.extend(["", "Validity domain: stationary noise-free OU; truth (1,0,0.5); exogenous regular/Gamma observations; the registered n/q/CV grid and floor; conditional fits and natural-scale full-Hessian Wald intervals. Nominal span is (n-1)q in relaxation-time units. Different cells are independent and do not share their realized span. PFML/Euler covariance remains model-based. New real-data operating rules would require separate validation.", ""])
    return "\n".join(lines).encode()


def copy_source_figures(bundle, output=None):
    """Return the exact source-copy inventory; write only when output is supplied."""
    copied = []
    for study in bundle.config["studies"]:
        source = bundle.study_dir(study)
        if not (source/"figures").is_dir():
            continue
        review_path = source/"visual_review.json"
        review_rel = review_path.relative_to(bundle.root).as_posix()
        require(review_path.is_file(), "copied source figures have no original review")
        review_out = f"figures/source/{study}/original_visual_review.json"
        manifest_path = source/"figures/figure_manifest.json"
        manifest_rel = manifest_path.relative_to(bundle.root).as_posix()
        manifest_out = f"figures/source/{study}/original_figure_manifest.json"
        if output is not None:
            write(output/review_out, review_path.read_bytes())
            write(output/manifest_out, bundle.read_bytes(manifest_rel))
        for path in sorted((source/"figures").glob("*")):
            if path.suffix not in (".png", ".pdf"):
                continue
            target = f"figures/source/{study}/{path.name}"
            if output is not None:
                write(output/target, path.read_bytes())
            copied.append({"path": target, "sha256": digest(path), "kind": "byte_identical_source_copy",
                           "study": study, "source_path": path.relative_to(bundle.root).as_posix(),
                           "original_review_path": review_rel, "original_review_sha256": digest(review_path),
                           "review_record_in_package": review_out,
                           "original_figure_manifest_path": manifest_rel,
                           "original_figure_manifest_sha256": digest(manifest_path),
                           "figure_manifest_in_package": manifest_out,
                           "review_claim": "original saved review transferred only for identical figure bytes"})
    auxiliary = bundle.config.get("additional_sources", {}).get("selected_figures")
    if auxiliary:
        path = auxiliary["path"] if isinstance(auxiliary, dict) else auxiliary
        source = bundle.root/path
        review_rel = "artifacts/week7_evidence/selected_visual_review.json"
        review_out = "figures/source/e3_selected/original_visual_review.json"
        review_payload = bundle.read_bytes(review_rel)
        manifest_rel = (source/"figure_manifest.json").relative_to(bundle.root).as_posix()
        manifest_payload = bundle.read_bytes(manifest_rel)
        manifest_out = "figures/source/e3_selected/original_figure_manifest.json"
        if output is not None:
            write(output/review_out, review_payload)
            write(output/manifest_out, manifest_payload)
        for image in sorted(source.glob("*")):
            if image.suffix not in (".png", ".pdf"):
                continue
            target = "figures/source/e3_selected/"+image.name
            if output is not None:
                write(output/target, image.read_bytes())
            copied.append({"path": target, "sha256": digest(image), "kind": "byte_identical_source_copy",
                           "study": "e3_selected", "source_path": image.relative_to(bundle.root).as_posix(),
                           "original_review_path": review_rel, "original_review_sha256": sha(review_payload),
                           "review_record_in_package": review_out,
                           "original_figure_manifest_path": manifest_rel,
                           "original_figure_manifest_sha256": sha(manifest_payload),
                           "figure_manifest_in_package": manifest_out,
                           "review_claim": "original saved review transferred only for identical figure bytes"})
    return copied


def package_readme():
    return ("# Week 9 publication package\n\n"
            "This derivative package reads the completed E1–E4 studies. It contains publication figures, T1, full plotting data, claims, captions and benchmark/literature notes.\n\n"
            "- `package_state.json` records the lifecycle and binds the final visual review.\n"
            "- `manifest.json` binds every immutable package file.\n"
            "- `figure_manifest.json` distinguishes new layouts from byte-identical source copies.\n"
            "- `tables/T1.md` summarizes the complete classifications and all six unresolved cases.\n"
            "- `tables/T1_e2.csv`, `tables/T1_e1.csv` and `tables/unresolved.csv` retain all eligible classifications.\n"
            "- `claims.json` carries numerical values, source selectors, uncertainty and qualifications.\n"
            "- `source_catalog.json` pins baseline Git objects and all consumed saved bytes.\n"
            "- `notes/week9_reproduction.md` gives exact reproduction and review commands.\n\n"
            "The derivative verification does not extend the original numerical audit/replay scope. μ relative-to-truth metrics remain undefined. Reported comparisons are pointwise; classifications describe discrete investigated designs.\n").encode()


def build(config_path, output, *, root=ROOT):
    bundle = load_sources(config_path, root=root)
    destination = guard_output(bundle, output)
    before = bundle.snapshot()
    source_identity = renderer_identity(bundle.root)
    tables = build_tables(bundle)
    claims = build_claims(bundle, tables)
    destination.mkdir(parents=True)
    try:
        # Preserve the pinned catalog's exact bytes, including its serialization.
        # Its file digest, not just its parsed values, is part of the contract.
        catalog_bytes = bundle.config_path.read_bytes()
        write(destination/"source_catalog.json", catalog_bytes)
        write(destination/"source_snapshot.json", json_bytes(before))
        table_manifest = {}
        for name, frame in tables.items():
            filename, payload = table_payload(name, frame)
            path = "tables/"+filename
            write(destination/path, payload)
            table_manifest[name] = {"path": path, "rows": len(frame), "sha256": sha(payload),
                                    "columns": list(frame.columns)}
        write(destination/"tables/T1.md", readable_t1(tables))
        write(destination/"claims.json", json_bytes(claims))
        write(destination/"claims.csv", claims_csv_bytes(claims))
        write(destination/"claims.md", readable_claims(claims).encode())
        figures = build_figures(tables)
        new = []
        try:
            for stem, fig in figures.items():
                for extension in ("png", "pdf"):
                    path = f"figures/new/{stem}.{extension}"
                    payload = figure_bytes(fig, extension, stem)
                    write(destination/path, payload)
                    new.append({"path": path, "sha256": sha(payload), "kind": "new_week9_layout",
                                "caption_id": stem, "review_required": True})
        finally:
            for fig in figures.values():
                plt.close(fig)
        copied = copy_source_figures(bundle, destination)
        figure_manifest = {"schema_version": 1, "contract": CONTRACT, "figures": new+copied,
                           "new_pngs": [f["path"] for f in new if f["path"].endswith(".png")],
                           "copied_pngs": [f["path"] for f in copied if f["path"].endswith(".png")],
                           "table_bindings": table_manifest,
                           "selection": "post-result publication presentation; complete saved evidence retained"}
        write(destination/"figure_manifest.json", json_bytes(figure_manifest))
        write(destination/"captions.md", captions(figures, copied).encode())
        for note in NOTES:
            write(destination/"notes"/note, publication_note(bundle.root, note, bundle.config["baseline_commit"]))
        after = bundle.snapshot()
        require(before == after, "source bytes changed while building publication outputs")
        require(source_identity == renderer_identity(bundle.root), "renderer/notes changed during build")
        record = {"schema_version": 1, "contract": CONTRACT, "baseline_commit": bundle.config["baseline_commit"],
                  "source_catalog_sha256": sha(catalog_bytes), "source_snapshot_sha256": sha(json_bytes(before)),
                  "source_files_verified": len(before), "source_preservation_passed": True,
                  "renderer_files": source_identity, "runtime": runtime(), "tables": table_manifest,
                  "claims": len(claims), "new_pngs": len(figure_manifest["new_pngs"]),
                  "copied_pngs": len(figure_manifest["copied_pngs"]),
                  "new_simulated_sources": 0, "new_model_fits": 0,
                  "verification_scope": "pinned baseline bytes, saved contracts and derivative arithmetic; no new numerical source/fit replay"}
        write(destination/"build_record.json", json_bytes(record))
        write(destination/"README.md", package_readme())
        inventory = {p.relative_to(destination).as_posix(): digest(p) for p in sorted(destination.rglob("*")) if p.is_file()}
        write(destination/"manifest.json", json_bytes({"schema_version": 1, "contract": CONTRACT, "files": inventory}))
        write(destination/"package_state.json", json_bytes({"schema_version": 1, "status": "pending_visual_review",
            "manifest_sha256": digest(destination/"manifest.json"),
            "figure_manifest_sha256": digest(destination/"figure_manifest.json")}))
        return verify(destination, root=bundle.root)
    except Exception as exc:
        write(destination/"build_failure.json", json_bytes({"status": "failed_preserved", "error_type": type(exc).__name__,
                                                           "error": str(exc), "contract": CONTRACT}))
        raise


def verify(run_dir, *, root=ROOT, require_completed=False):
    output = Path(run_dir).resolve()
    require(output.is_dir() and not Path(run_dir).is_symlink(), "publication directory missing or linked")
    state, manifest = read_json(output/"package_state.json"), read_json(output/"manifest.json")
    require(state["manifest_sha256"] == digest(output/"manifest.json"), "package manifest changed")
    require(manifest.get("contract") == CONTRACT, "unsupported publication contract")
    expected = manifest["files"]
    require(all(not p.is_symlink() for p in output.rglob("*")), "publication contains a symlink")
    actual = {p.relative_to(output).as_posix() for p in output.rglob("*") if p.is_file()}
    allowed = set(expected)|{"manifest.json", "package_state.json"}
    if state.get("status") == "completed":
        allowed.add("visual_review.json")
    require(actual == allowed, "publication file membership changed")
    for relative, wanted in expected.items():
        path = output/relative
        require(not Path(relative).is_absolute() and ".." not in Path(relative).parts
                and not path.is_symlink() and path.resolve().is_relative_to(output), "unsafe package member")
        require(digest(path) == wanted, "publication file changed: "+relative)
    bundle = load_sources(output/"source_catalog.json", root=root)
    record = read_json(output/"build_record.json")
    require(record["source_snapshot_sha256"] == sha(json_bytes(bundle.snapshot()))
            == digest(output/"source_snapshot.json"), "publication source snapshot changed")
    require(record["renderer_files"] == renderer_identity(bundle.root), "publication renderer/notes changed")
    require(record["source_catalog_sha256"] == digest(output/"source_catalog.json"), "source catalog binding changed")
    require(record["new_simulated_sources"] == record["new_model_fits"] == 0, "unexpected scientific execution claim")
    # Rebuild numerical publication views from their pinned sources. This is
    # derivative arithmetic, not regeneration of an OU path or fitted interval.
    tables = build_tables(bundle)
    expected_tables = {}
    for name, frame in tables.items():
        filename, payload = table_payload(name, frame)
        require((output/"tables"/filename).read_bytes() == payload, "derived table arithmetic differs: "+name)
        expected_tables[name] = {"path": "tables/"+filename, "rows": len(frame), "sha256": sha(payload), "columns": list(frame.columns)}
    require(record.get("tables") == expected_tables, "build table accounting/bindings differ")
    require((output/"tables/T1.md").read_bytes() == readable_t1(tables), "readable T1 differs from unrounded tables")
    for note in NOTES:
        require((output/"notes"/note).read_bytes() == publication_note(bundle.root, note, bundle.config["baseline_commit"]),
                "publication note/link bindings changed: "+note)
    claims = build_claims(bundle, tables)
    require(record.get("claims") == len(claims), "build claim accounting differs")
    require((output/"claims.json").read_bytes() == json_bytes(claims), "claim arithmetic/source bindings differ")
    require((output/"claims.csv").read_bytes() == claims_csv_bytes(claims), "claim CSV differs from bound claims")
    require((output/"claims.md").read_bytes() == readable_claims(claims).encode(), "readable claims differ from bound claims")
    figures = read_json(output/"figure_manifest.json")
    require(state["figure_manifest_sha256"] == digest(output/"figure_manifest.json"), "figure manifest changed")
    copied = copy_source_figures(bundle)
    new = [{"path": f"figures/new/{stem}.{ext}", "sha256": digest(output/f"figures/new/{stem}.{ext}"),
            "kind": "new_week9_layout", "caption_id": stem, "review_required": True}
           for stem in STEMS for ext in ("png", "pdf")]
    require(figures.get("schema_version") == 1 and figures.get("contract") == CONTRACT
            and figures.get("figures") == new+copied, "publication figure role inventory differs")
    require(figures.get("new_pngs") == [f["path"] for f in new if f["path"].endswith(".png")]
            and figures.get("copied_pngs") == [f["path"] for f in copied if f["path"].endswith(".png")],
            "publication PNG review inventory differs")
    require(figures.get("table_bindings") == expected_tables, "figure table bindings differ")
    require(record.get("new_pngs") == len(figures["new_pngs"])
            and record.get("copied_pngs") == len(figures["copied_pngs"]), "build figure accounting differs")
    require(record.get("source_files_verified") == len(bundle.snapshot())
            and record.get("source_preservation_passed") is True
            and record.get("baseline_commit") == bundle.config["baseline_commit"], "build source accounting differs")
    fixed = {"source_catalog.json", "source_snapshot.json", "claims.json", "claims.csv", "claims.md",
             "figure_manifest.json", "captions.md", "build_record.json", "README.md", "tables/T1.md"}
    wanted_members = (fixed | {f["path"] for f in new+copied}
                      | {f["review_record_in_package"] for f in copied}
                      | {f["figure_manifest_in_package"] for f in copied}
                      | {v["path"] for v in expected_tables.values()} | {"notes/"+n for n in NOTES})
    require(set(expected) == wanted_members, "publication expected content membership differs")
    require((output/"README.md").read_bytes() == package_readme(), "publication index differs from its contract")
    for figure in figures["figures"]:
        require(digest(output/figure["path"]) == figure["sha256"], "figure binding changed")
        if figure["kind"] == "byte_identical_source_copy":
            require(digest(bundle.root/figure["source_path"]) == figure["sha256"], "source figure copy differs")
            require(digest(bundle.root/figure["original_review_path"])
                    == digest(output/figure["review_record_in_package"])
                    == figure["original_review_sha256"], "source review transfer changed")
            require(digest(bundle.root/figure["original_figure_manifest_path"])
                    == digest(output/figure["figure_manifest_in_package"])
                    == figure["original_figure_manifest_sha256"], "source figure manifest transfer changed")
    require((output/"captions.md").read_bytes() == captions(STEMS, copied).encode(), "captions differ from their recorded figure definitions")
    require(state.get("status") in ("pending_visual_review", "completed"), "unsupported package lifecycle")
    if state["status"] == "completed":
        require(state["visual_review_sha256"] == digest(output/"visual_review.json"), "visual review changed")
        validate_review(output, read_json(output/"visual_review.json"))
        require(isinstance(state.get("renderer_commit"), str), "completed renderer commit binding is missing")
        committed_renderer(bundle.root, record["renderer_files"], state.get("renderer_commit"))
    require(not require_completed or state["status"] == "completed", "publication is awaiting actual visual review")
    return {"contract": CONTRACT, "status": state["status"], "verified": True,
            "source_files": len(bundle.snapshot()), "immutable_files": len(expected),
            "tables": {k: len(v) for k,v in tables.items()}, "claims": record["claims"],
            "new_pngs": record["new_pngs"], "copied_pngs": record["copied_pngs"],
            "manifest_sha256": digest(output/"manifest.json"), "new_simulated_sources": 0, "new_model_fits": 0}


def validate_review(output, review):
    figures = read_json(output/"figure_manifest.json")
    require(review.get("figure_manifest_sha256") == digest(output/"figure_manifest.json"), "review is for different figures")
    require(isinstance(review.get("reviewer"), str) and review["reviewer"].strip(), "reviewer identity missing")
    require(review.get("review_type") == "actual_ai_raster_and_pdf_review", "review type must identify actual AI inspection")
    try:
        when = datetime.fromisoformat(review["reviewed_utc"].replace("Z", "+00:00"))
        require(when.tzinfo is not None and when <= datetime.now(timezone.utc), "invalid/future review timestamp")
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("actual review timestamp missing or invalid") from exc
    wanted = {f["path"]: f["sha256"] for f in figures["figures"]
              if f["kind"] == "new_week9_layout" and f["path"].endswith(".png")}
    require(review.get("reviewed_png_sha256") == wanted, "review omits or changes a new PNG")
    pdfs = {f["path"] for f in figures["figures"]
            if f["kind"] == "new_week9_layout" and f["path"].endswith(".pdf")}
    require(set(review.get("reviewed_pdf_pages", {})) == pdfs, "review omits a new PDF layout")
    for relative, value in review["reviewed_pdf_pages"].items():
        require(isinstance(value, dict) and value.get("page") == 1
                and re.fullmatch(r"[0-9a-f]{64}", value.get("rendered_page_sha256", "")) is not None
                and value.get("source_pdf_sha256") == digest(output/relative)
                and isinstance(value.get("render_command"), str) and value["render_command"]
                and isinstance(value.get("renderer"), str) and value["renderer"],
                "PDF page review binding is incomplete")
    require(review.get("all_panels_legible") is True and review.get("blocking_issues") == [], "visual review has blocking issues")
    require(isinstance(review.get("notes"), list) and review["notes"], "actual inspection notes are missing")


def finalize(run_dir, review_file, *, root=ROOT):
    output = Path(run_dir).resolve()
    result = verify(output, root=root)
    review_bytes = Path(review_file).read_bytes()
    review = json.loads(review_bytes)
    validate_review(output, review)
    if result["status"] == "completed":
        require((output/"visual_review.json").read_bytes() == review_bytes, "completed review cannot be replaced")
        return result
    renderer_commit = committed_renderer(Path(root), read_json(output/"build_record.json")["renderer_files"])
    write(output/"visual_review.json", review_bytes)
    state = read_json(output/"package_state.json")
    state.update(status="completed", visual_review_sha256=sha(review_bytes), completed_utc=datetime.now(timezone.utc).isoformat(),
                 renderer_commit=renderer_commit)
    # Only this derivative lifecycle changes. The immutable manifest excludes
    # state/review, and the existing pending state was verified immediately above.
    temporary = output/".package_state.tmp"
    with temporary.open("xb") as stream:
        stream.write(json_bytes(state))
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, output/"package_state.json")
    return verify(output, root=root, require_completed=True)


def compare(left, right, *, root=ROOT):
    verify(left, root=root)
    verify(right, root=root)
    a, b = Path(left), Path(right)
    ma, mb = read_json(a/"manifest.json"), read_json(b/"manifest.json")
    require(ma == mb, "external reproduction manifest differs")
    for relative, wanted in ma["files"].items():
        require(digest(a/relative) == digest(b/relative) == wanted, "external reproduction bytes differ: "+relative)
    return {"byte_identical": True, "immutable_files": len(ma["files"]),
            "manifest_sha256": digest(a/"manifest.json"),
            "scope": "all immutable derivative files, including PNG/PDF/table/claim bytes; lifecycle/review are separate"}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    command = sub.add_parser("build")
    command.add_argument("--sources", default=str(ROOT/"configs/week9_publication.json"))
    command.add_argument("--output", required=True)
    command = sub.add_parser("verify")
    command.add_argument("--run-dir", required=True)
    command.add_argument("--require-completed", action="store_true")
    command = sub.add_parser("finalize")
    command.add_argument("--run-dir", required=True)
    command.add_argument("--review-file", required=True)
    command = sub.add_parser("compare")
    command.add_argument("--left", required=True)
    command.add_argument("--right", required=True)
    args = parser.parse_args(argv)
    if args.command == "build":
        result = build(args.sources, args.output)
    elif args.command == "verify":
        result = verify(args.run_dir, require_completed=args.require_completed)
    elif args.command == "finalize":
        result = finalize(args.run_dir, args.review_file)
    else:
        result = compare(args.left, args.right)
    print(json_bytes(result).decode().strip())


if __name__ == "__main__":
    main()
