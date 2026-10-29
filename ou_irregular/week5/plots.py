"""First-pass E1 figures regenerated exclusively from saved complete tables."""
from io import BytesIO
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..experiments.provenance import file_digest, verify
from ..experiments.storage import atomic_bytes, atomic_json


STYLE = {"font.family": "DejaVu Sans", "font.size": 10, "axes.titlesize": 11,
         "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": "white",
         "axes.grid": True, "grid.alpha": .18, "axes.axisbelow": True,
         "pdf.fonttype": 42, "ps.fonttype": 42}
METHODS = (("exact", "Exact", "#0072B2", "o"), ("pfml", "PFML", "#D55E00", "s"),
           ("euler", "Euler", "#009E73", "^"))
PARAMETERS = (("theta", r"$\theta$"), ("mu", r"$\mu$"), ("sigma", r"$\sigma$"))
Z = 1.959963984540054


def _series(ax, group, metric, lo=None, hi=None, se=None):
    for name, label, color, marker in METHODS:
        selected = group[group.estimator == name].sort_values("cv")
        x, y = selected.cv.to_numpy(), selected[metric].to_numpy()
        ax.plot(x, y, label=label, color=color, marker=marker, markersize=4.8, linewidth=1.7)
        if se is not None:
            lower, upper = y-Z*selected[se].to_numpy(), y+Z*selected[se].to_numpy()
        else:
            lower, upper = selected[lo].to_numpy(), selected[hi].to_numpy()
        mask = np.isfinite(y) & np.isfinite(lower) & np.isfinite(upper)
        ax.errorbar(x[mask], y[mask], yerr=np.maximum(0, np.vstack((y[mask]-lower[mask], upper[mask]-y[mask]))),
                    fmt="none", color=color, capsize=2.4, linewidth=.9)
    cvs = sorted(group.cv.unique())
    ax.set_xticks(cvs, [f"{c:g}" for c in cvs], rotation=30, ha="right")
    ax.set_xlabel("Nominal gap CV")


def _finish(fig, axes, title, design, footnote, output, stem, *, top=.80):
    label = "REGISTERED E1" if design.scope == "headline" else "VALIDATION ONLY"
    cell = design.cells[0]
    fig.suptitle(title, fontsize=15, fontweight="semibold", y=.98)
    fig.text(.5, .925, f"{label}  |  q = {cell.q:g}  |  n = {cell.n:,} observations  |  R = {design.replications:,} per CV",
             ha="center", fontsize=10)
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles[:3], labels[:3], loc="upper center", bbox_to_anchor=(.5, .90),
               ncol=3, frameon=False)
    fig.subplots_adjust(left=.075, right=.985, bottom=.25 if fig.get_figheight() < 7 else .20,
                        top=top, wspace=.28, hspace=.40)
    fig.text(.075, .035, footnote, fontsize=8.5, color="#333333", va="bottom")
    for extension in ("png", "pdf"):
        stream = BytesIO()
        metadata = {"Title": title, "Creator": "OU irregular-sampling research",
                    "CreationDate": None, "ModDate": None} if extension == "pdf" else {"Software": "OU irregular-sampling research"}
        fig.savefig(stream, format=extension, dpi=170, metadata=metadata)
        atomic_bytes(output / f"{stem}.{extension}", stream.getvalue(), immutable=True)
    plt.close(fig)


def render_figures(run_dir, design, *, output_dir=None):
    run_dir = Path(run_dir)
    metadata = json.loads((run_dir / "run_metadata.json").read_text())
    verify(metadata["provenance"], expected_design=design)
    manifest = json.loads((run_dir / "table_manifest.json").read_text())
    for filename in ("parameter_summary.csv", "regular_comparisons.csv"):
        if file_digest(run_dir / filename) != manifest[filename]["sha256"]:
            raise ValueError("figure input table digest mismatch")
    summary = pd.read_csv(run_dir / "parameter_summary.csv", float_precision="round_trip")
    ratios = pd.read_csv(run_dir / "regular_comparisons.csv", float_precision="round_trip")
    if len({(c.n, c.q) for c in design.cells}) != 1:
        raise ValueError("E1 curves require a single n and q; use the separate full-grid plotter")
    if len(summary) != len(design.cells)*9 or len(ratios) != len(design.cells)*9:
        raise ValueError("incomplete figure inputs")
    output = Path(output_dir) if output_dir is not None else run_dir / "figures"
    output.mkdir(parents=True, exist_ok=True)
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(2, 3, figsize=(13.5, 8.2), squeeze=False)
        for col, (p, symbol) in enumerate(PARAMETERS):
            group = summary[summary.parameter == p]
            for row, scope in enumerate(("operational", "conditional")):
                ax = axes[row, col]
                _series(ax, group, f"coverage_{scope}", f"{scope}_wilson_lower", f"{scope}_wilson_upper")
                ax.axhline(.95, color="#555555", linestyle=":", linewidth=1)
                if p != "mu" and scope == "operational":
                    ax.axhline(.80, color="#999999", linestyle="--", linewidth=.9)
                ax.set_ylim(-.025, 1.025)
                ax.set_title(f"{symbol} | " + ("all planned paths" if scope == "operational" else "valid intervals"))
                if col == 0: ax.set_ylabel("Coverage")
        _finish(fig, axes, "F2  Coverage versus gap variability", design,
                "Bars: pointwise 95% Wilson Monte Carlo intervals. Dotted: nominal 95%; dashed: 80% practical coverage criterion.\n"
                "Operational coverage counts unavailable intervals as noncoverage. PFML/Euler intervals use model-based inverse Hessians.",
                output, "F2_coverage", top=.825)

        fig, axes = plt.subplots(1, 2, figsize=(11.2, 5.6), squeeze=False)
        group = summary[summary.parameter == "theta"]  # availability is joint across parameters
        for ax, kind in zip(axes.flat, ("point", "interval")):
            _series(ax, group, f"{kind}_availability", f"{kind}_availability_wilson_lower", f"{kind}_availability_wilson_upper")
            ax.set_title("Valid point fits" if kind == "point" else "Available Wald intervals")
            ax.set_ylim(-.025, 1.025)
            ax.set_ylabel("Fraction of planned paths")
        _finish(fig, axes, "F2 companion  Numerical availability", design,
                "Bars: pointwise 95% Wilson Monte Carlo intervals. Every planned path stays in the denominator.\n"
                "Validity is joint across the three parameters; the saved status table retains each exclusion reason.", output, "F2_availability")

        fig, axes = plt.subplots(1, 2, figsize=(11.2, 5.6), squeeze=False)
        for ax, (p, symbol) in zip(axes.flat, (PARAMETERS[0], PARAMETERS[2])):
            _series(ax, summary[summary.parameter == p], "relative_bias", se="relative_mcse_bias")
            ax.set_title(symbol)
            for value in (-.10, .10): ax.axhline(value, color="#999999", linestyle="--", linewidth=.9)
            ax.axhline(0, color="#555555", linestyle=":", linewidth=.8)
            ax.set_ylabel("Bias / absolute truth")
        _finish(fig, axes, "F3  Relative bias of theta and sigma", design,
                "Bars: approximate pointwise 95% Monte Carlo intervals, conditional on valid point fits. Dashed: +/-10% criteria.\n"
                "Lines join investigated CV values as visual guides; they do not estimate a universal smooth boundary.", output, "F3_relative_bias")

        fig, axes = plt.subplots(1, 2, figsize=(11.2, 5.6), squeeze=False)
        group = summary[summary.parameter == "mu"]
        for ax, metric, se, label in zip(axes.flat, ("bias", "rmse"), ("mcse_bias", "mcse_rmse"),
                                        ("Bias (original units)", "RMSE (original units)")):
            _series(ax, group, metric, se=se)
            ax.set_title(label)
            ax.set_ylabel("Original state units")
            ax.ticklabel_format(axis="y", style="sci", scilimits=(-3, 3))
            if metric == "bias": ax.axhline(0, color="#555555", linestyle=":", linewidth=.8)
        _finish(fig, axes, "F3 companion  Mean estimation in absolute units", design,
                "True mu is zero: truth-relative metrics are undefined. Bars: approximate pointwise 95% Monte Carlo intervals.\n"
                "Both metrics condition on valid point fits; the full table reports their planned and valid counts.", output, "F3_mu_absolute")

        fig, axes = plt.subplots(1, 3, figsize=(13.5, 5.6), squeeze=False)
        for ax, (p, symbol) in zip(axes.flat, PARAMETERS):
            group = ratios[ratios.parameter == p]
            _series(ax, group, "rmse_ratio", "rmse_ratio_mc_lower", "rmse_ratio_mc_upper")
            ax.set_yscale("log")
            zeros = group[group.rmse_ratio == 0]
            if len(zeros):
                labels = [f"{r.estimator}: CV={r.cv:g}" for r in zeros.itertuples()]
                ax.text(.02, .02, "Zero ratios (off log axis):\n" + "\n".join(labels),
                        transform=ax.transAxes, fontsize=7, va="bottom")
            ax.set_title(symbol)
            ax.axhline(1, color="#555555", linestyle=":", linewidth=1)
            ax.set_ylabel("RMSE / same-estimator regular RMSE")
        _finish(fig, axes, "F4  RMSE ratios to regular-design references", design,
                "Bars: pointwise log-delta Monte Carlo intervals; vertical scale is logarithmic. CV=0 is an identity, without independent-sample error.\n"
                "References match n, truth and nominal mean gap, not realized span. Cells are independent; contrasts share their regular baseline.",
                output, "F4_rmse_ratios")
    artifacts = {p.name: file_digest(p) for p in sorted(output.iterdir()) if p.suffix in (".png", ".pdf")}
    record = {"design_sha256": design.digest, "input_table_sha256": {
        name: manifest[name]["sha256"] for name in ("parameter_summary.csv", "regular_comparisons.csv")}, "figures": artifacts}
    atomic_json(output / "figure_manifest.json", record, immutable=True)
    return record
