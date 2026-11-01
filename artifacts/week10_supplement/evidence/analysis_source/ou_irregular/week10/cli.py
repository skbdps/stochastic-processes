"""Explicit stages for the Week 10 explanatory supplement."""
import argparse
import json
from pathlib import Path

from ..experiments.provenance import file_digest, utc_now
from ..experiments.storage import atomic_json
from .design import full_design, validation_design
from .execution import audit, freeze, run_raw


def inventory(root):
    root = Path(root)
    result = {}
    for path in sorted(root.rglob("*")):
        if path.is_symlink():
            raise ValueError("supplement may not contain symbolic links")
        relative = str(path.relative_to(root))
        excluded = {"completion_seal.json", "horizon_full/.execution.lock",
                    "horizon_validation/.execution.lock"}
        if path.is_file() and relative not in excluded:
            result[relative] = file_digest(path)
    return result


def seal(package):
    """Bind completed evidence only after actual figure review and validation."""
    package = Path(package)
    if (package/"completion_seal.json").exists():
        raise FileExistsError("completed supplement is immutable")
    review = json.loads((package/"evidence/visual_review.json").read_text())
    validation = json.loads((package/"evidence/final_validation.json").read_text())
    figures = json.loads((package/"figures/figure_manifest.json").read_text())
    if review.get("status") != "passed" or validation.get("status") != "passed":
        raise ValueError("review and validation must actually pass before sealing")
    required_files = {str(p.relative_to(package)) for name in ("horizon_full", "span", "analysis", "figures", "interpretation")
                      for p in (package/name).rglob("*") if p.is_file() and p.name != ".execution.lock"}
    if not all((package/name).is_dir() for name in ("horizon_full", "span", "analysis", "figures", "interpretation")):
        raise ValueError("complete scientific evidence is required before sealing")
    bound = validation.get("validated_files", {})
    if not required_files.issubset(bound):
        raise ValueError("complete scientific evidence must be bound by final validation")
    for name, digest in bound.items():
        if Path(name).is_absolute() or ".." in Path(name).parts or file_digest(package/name) != digest:
            raise ValueError("stale or unsafe final validation binding")
    tests = validation.get("tests", {})
    if tests.get("passed", 0) < 1 or tests.get("failures") != 0 or tests.get("errors") != 0:
        raise ValueError("passing validation tests are required")
    audit_record = json.loads((package/"horizon_full/audit.json").read_text())
    expected = {"status": "passed", "sources_reconstructed": 14000, "fit_records": 42000,
                "endpoint_fit_replays": 42, "pipeline_errors": 0,
                "all_realized_spans_equal_horizon": True}
    if any(audit_record.get(k) != v for k, v in expected.items()):
        raise ValueError("complete full-study audit is required")
    from .execution import readback
    from .reporting import verify_tables
    from .span import verify_span_analysis
    from .interpretation import verify as verify_interpretation
    metadata, _ = readback(full_design(), package/"horizon_full")
    if audit_record["execution_identity_sha256"] != metadata["provenance"]["identity_sha256"]:
        raise ValueError("audit binds another execution")
    verify_tables(full_design(), package/"horizon_full", package/"analysis",
                  Path(__file__).resolve().parents[2]/"artifacts/week8_e4", replay_mechanism=True)
    verify_span_analysis(Path(__file__).resolve().parents[2]/"artifacts/week8_e4", package/"span")
    verify_interpretation(package)
    if review.get("figure_manifest_sha256") != file_digest(package/"figures/figure_manifest.json"):
        raise ValueError("review binds different figures")
    if (not figures["files"] or {Path(p).suffix for p in figures["files"]} != {".png", ".pdf"}
            or {Path(p).stem for p in figures["files"] if p.endswith(".png")}
            != {Path(p).stem for p in figures["files"] if p.endswith(".pdf")}):
        raise ValueError("complete PNG and PDF figure pairs are required")
    if set(review.get("inspected_files", {})) != set(figures["files"]):
        raise ValueError("every PNG and rendered PDF requires a recorded inspection")
    for name, digest in figures["files"].items():
        if file_digest(package/"figures"/name) != digest or review["inspected_files"][name] != digest:
            raise ValueError("inspected figure bytes changed")
    result = {"schema_version": 1, "status": "completed", "sealed_utc": utc_now(),
              "scope": "Week 10 explanatory supplement; original E4 remains unchanged",
              "files": inventory(package)}
    atomic_json(package/"completion_seal.json", result, immutable=True)
    return {"status": "completed", "sealed_files": len(result["files"]),
            "seal_sha256": file_digest(package/"completion_seal.json")}


def verify_seal(package):
    package = Path(package)
    seal = json.loads((package/"completion_seal.json").read_text())
    if seal.get("status") != "completed" or inventory(package) != seal["files"]:
        raise ValueError("supplement inventory changed after completion")
    return {"status": "passed", "sealed_files": len(seal["files"]),
            "seal_sha256": file_digest(package/"completion_seal.json")}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("freeze", "run", "audit", "span", "analyze", "verify", "seal"))
    parser.add_argument("--run-dir", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--e4-dir", type=Path, default=Path("artifacts/week8_e4"))
    parser.add_argument("--package", type=Path, default=Path("artifacts/week10_supplement"))
    parser.add_argument("--freeze-file", type=Path)
    parser.add_argument("--validation", action="store_true")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    design = validation_design() if args.validation else full_design()
    if args.stage == "freeze":
        if args.output is None:
            parser.error("freeze requires --output")
        result = freeze(design, extra_files=("requirements-week10-lock.txt",))
        atomic_json(args.output, result, immutable=True)
        print(json.dumps({"status": "frozen", "identity_sha256": result["identity_sha256"]}))
        return
    if args.stage in {"run", "audit", "analyze"} and args.run_dir is None:
        parser.error(args.stage + " requires --run-dir")
    if args.stage == "run":
        frozen = json.loads(args.freeze_file.read_text()) if args.freeze_file else None
        result = run_raw(design, args.run_dir, workers=args.workers, resume=args.resume,
                         frozen_provenance=frozen)
        result = {key: result[key] for key in ("status", "completed_cells", "pipeline_errors")}
    elif args.stage == "audit":
        result = audit(design, args.run_dir)
    elif args.stage == "span":
        from .span import write_span_analysis
        if args.output is None:
            parser.error("span requires --output")
        result = write_span_analysis(args.e4_dir, args.output)
    elif args.stage == "analyze":
        from .reporting import build
        if args.output is None:
            parser.error("analyze requires --output")
        result = build(design, args.run_dir, args.e4_dir, args.output)
    elif args.stage == "seal":
        result = seal(args.package)
    else:
        from .reporting import verify_tables
        from .span import verify_span_analysis
        from .interpretation import verify as verify_interpretation
        result = {"seal": verify_seal(args.package),
                  "tables": verify_tables(design, args.run_dir or args.package/"horizon_full",
                                          args.package/"analysis", args.e4_dir),
                  "span": verify_span_analysis(args.e4_dir, args.package/"span"),
                  "interpretation": verify_interpretation(args.package)}
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
