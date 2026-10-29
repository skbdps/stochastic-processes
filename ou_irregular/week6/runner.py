"""Week 6 E2 entry point; headline settings are the frozen registration.

Set OPENBLAS_NUM_THREADS=1 and OMP_NUM_THREADS=1 before starting Python for
parallel runs, and retain those settings for resume/replay. The numerical
environment, including explicit thread variables, is part of the run identity.
"""
import argparse
import json
from pathlib import Path

from ..experiments.audit import AUDIT_FILES, audit_run
from ..experiments.execution import run_raw
from ..experiments.provenance import ROOT
from ..experiments.storage import atomic_json
from ..experiments.workflow import analyze_and_audit, load_saved_design, seal_completed_run
from .config import load_config
from .plots import render_figures


WEEK_FILES = (*AUDIT_FILES, "ou_irregular/experiments/workflow.py", "ou_irregular/week6/__init__.py",
              "ou_irregular/week6/config.py", "ou_irregular/week6/reporting.py",
              "ou_irregular/week6/plots.py", "ou_irregular/week6/runner.py",
              "docs/week5_mathematics.md", "docs/week6_mathematics.md")


def _saved_design(path):
    design = load_saved_design(path)
    registered = load_config(ROOT / "configs/week6_e2.yaml")
    expected = registered if design.scope == "headline" else registered.validation(design.replications)
    if design.to_dict() != expected.to_dict():
        raise ValueError("this entry point accepts only the frozen E2 or its separate validation stream")
    return design


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    run = commands.add_parser("run", help="execute E2 or a separate development validation")
    run.add_argument("--config", type=Path, default=ROOT / "configs/week6_e2.yaml")
    run.add_argument("--output", type=Path, required=True)
    run.add_argument("--workers", type=int, default=1)
    run.add_argument("--partition-size", type=int, default=250)
    run.add_argument("--resume", action="store_true")
    run.add_argument("--validation", action="store_true")
    run.add_argument("--replications", type=int, help="validation only; never changes headline R")
    run.add_argument("--quiet", action="store_true")
    audit = commands.add_parser("audit", help="read-only complete arithmetic and fixed-path replay audit")
    audit.add_argument("--run-dir", type=Path, required=True)
    audit.add_argument("--output", type=Path, help="optional separate audit-report file")
    replot = commands.add_parser("replot", help="regenerate figures exclusively from saved tables")
    replot.add_argument("--run-dir", type=Path, required=True)
    replot.add_argument("--output", type=Path, help="optional separate figure directory")
    seal = commands.add_parser("finalize", help="seal audited evidence after recorded visual inspection")
    seal.add_argument("--run-dir", type=Path, required=True)
    seal.add_argument("--review-file", type=Path, required=True)
    args = parser.parse_args(argv)
    if args.command == "run":
        design = load_config(args.config)
        if args.validation:
            design = design.validation(4 if args.replications is None else args.replications)
        elif args.replications is not None:
            parser.error("--replications is allowed only with --validation; headline R is frozen at 2000")
        config_path = str(args.config.resolve().relative_to(ROOT))
        run_raw(design, args.output, extra_files=(*WEEK_FILES, config_path), workers=args.workers,
                partition_size=args.partition_size, resume=args.resume,
                headline_design_digest=design.digest if design.scope == "headline" else None,
                progress=not args.quiet)
        result = analyze_and_audit(args.output, design, render_figures)
        print(json.dumps({"status": "audited_pending_visual_review", "scope": design.scope,
                          "paths": result["planned_paths"], "fits": result["planned_fit_records"],
                          "output": str(args.output)}, sort_keys=True))
    else:
        design = _saved_design(args.run_dir)
        if args.command == "audit":
            result = audit_run(args.run_dir, design)
            if args.output is not None:
                if args.output.resolve().is_relative_to(args.run_dir.resolve()):
                    parser.error("write a repeat audit report outside the immutable run directory")
                atomic_json(args.output, result, immutable=True)
            print(json.dumps({"status": result["status"], "paths": result["planned_paths"],
                              "fits": result["planned_fit_records"], "replayed_paths": result["replayed_source_paths"]}))
        elif args.command == "replot":
            result = render_figures(args.run_dir, design, output_dir=args.output)
            print(json.dumps({"figures": len(result["figures"]), "design_sha256": design.digest}))
        else:
            review = json.loads(args.review_file.read_text())
            result = seal_completed_run(args.run_dir, design, review)
            print(json.dumps({"status": result["status"], "output": str(args.run_dir),
                              "artifact_files": len(result["artifact_sha256"])}))


if __name__ == "__main__":
    main()
