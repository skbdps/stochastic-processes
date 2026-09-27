"""Headless exploratory figures for the Week 3 simulation harness.

These figures describe successful fits and separately expose fit failures.
They deliberately do not plot coverage: Week 3 supplies no confidence-interval
estimator, and missing coverage is not evidence of zero or nominal coverage.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd

from .metrics import add_twin_comparisons


_COLORS = {"exact": "#2166ac", "pfml": "#d6604d", "euler": "#1b9e77"}
_LABELS = {"exact": "Exact MLE", "pfml": "Mean-gap PFML", "euler": "Euler"}
_LINESTYLES = ("-", "--", ":", "-.")


def _draw_curves(ax, data: pd.DataFrame, metric: str, sample_sizes: list) -> None:
    """Use color for estimator and line style for n, preserving missing values."""
    for (estimator, n), group in data.groupby(["estimator", "n"], sort=True):
        group = group.sort_values("cv")
        ax.plot(
            group["cv"], group[metric], marker="o", markersize=4,
            color=_COLORS.get(estimator, "#666666"),
            linestyle=_LINESTYLES[sample_sizes.index(n) % len(_LINESTYLES)],
            linewidth=1.6,
        )
    ax.set_xlabel("Configured gap CV")
    ax.grid(alpha=0.22)
    ax.spines[["top", "right"]].set_visible(False)


def _legend(estimators: list, sample_sizes: list) -> list:
    handles = [Line2D([0], [0], color=_COLORS.get(estimator, "#666666"), lw=2,
                      label=_LABELS.get(estimator, estimator)) for estimator in estimators]
    if len(sample_sizes) > 1:
        handles.extend(Line2D([0], [0], color="#555555",
                             linestyle=_LINESTYLES[index % len(_LINESTYLES)],
                             label=f"n={int(n)}") for index, n in enumerate(sample_sizes))
    return handles


def _finish(fig, handles: list, title: str, path: Path) -> None:
    fig.suptitle(title, fontsize=13, y=0.995)
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.957),
               ncol=min(6, len(handles)), frameon=False, fontsize=9)
    fig.text(0.5, 0.012,
             "Week 3 smoke / exploratory results. Bias and RMSE conditional on success. "
             "No CI coverage computed.", ha="center", fontsize=8, color="#555555")
    fig.tight_layout(rect=(0.02, 0.05, 0.99, 0.89))
    fig.savefig(path, dpi=170, bbox_inches="tight")
    plt.close(fig)


def plot_summary(summary: pd.DataFrame, output_dir: str | Path) -> list[Path]:
    """Write bias, equidistant-twin RMSE-ratio, and fit-failure PNGs.

    Rows are parameters and columns are theta times nominal mean gap. With the
    Week 3 three-coarseness smoke design this gives a readable 3-by-3 figure.
    Theta/sigma bias is relative; mu bias is absolute because true mu is zero.
    Controls use nominal mean gaps; flooring can alter expected realized gaps.
    """
    if summary.empty:
        raise ValueError("cannot plot an empty simulation summary")
    required = {"cell_id", "estimator", "parameter", "n", "theta_mean_gap", "cv",
                "bias", "relative_bias", "failure_rate"}
    missing = required.difference(summary.columns)
    if missing:
        raise ValueError(f"summary missing required plot columns: {sorted(missing)}")
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    data = add_twin_comparisons(summary) if "rmse_ratio" not in summary else summary.copy()
    coarseness = sorted(data["theta_mean_gap"].unique())
    sample_sizes = sorted(data["n"].unique())
    estimators = sorted(data["estimator"].unique())
    if len(coarseness) > 3 or len(sample_sizes) > 3:
        raise ValueError("plot_summary supports up to three coarseness values and sample sizes; subset the design first")
    # A curve must not silently connect different truths/designs sharing CV.
    if data.duplicated(["theta_mean_gap", "cv", "n", "estimator", "parameter"]).any():
        raise ValueError("plot_summary requires one design per theta_mean_gap/CV/n/estimator/parameter")
    handles = _legend(estimators, sample_sizes)
    paths = []

    fig, axes = plt.subplots(3, len(coarseness), squeeze=False,
                             figsize=(max(8.5, 4.0 * len(coarseness)), 9.1))
    for row, parameter in enumerate(("theta", "mu", "sigma")):
        metric = "bias" if parameter == "mu" else "relative_bias"
        for column, gap in enumerate(coarseness):
            ax = axes[row, column]
            subset = data.loc[data["parameter"].eq(parameter) & data["theta_mean_gap"].eq(gap)]
            _draw_curves(ax, subset, metric, sample_sizes)
            ax.axhline(0, color="#555555", linewidth=0.8, alpha=0.5)
            ax.set_ylabel(f"{parameter}: {'absolute' if parameter == 'mu' else 'relative'} bias")
            if row == 0:
                ax.set_title(f"theta × nominal mean gap = {gap:g}", fontsize=10)
    path = output_dir / "week3_bias_vs_cv.png"
    _finish(fig, handles, "Parameter bias across sampling irregularity", path)
    paths.append(path)

    fig, axes = plt.subplots(3, len(coarseness), squeeze=False,
                             figsize=(max(8.5, 4.0 * len(coarseness)), 9.1))
    for row, parameter in enumerate(("theta", "mu", "sigma")):
        for column, gap in enumerate(coarseness):
            ax = axes[row, column]
            subset = data.loc[data["parameter"].eq(parameter) & data["theta_mean_gap"].eq(gap)]
            _draw_curves(ax, subset, "rmse_ratio", sample_sizes)
            ax.axhline(1, color="#555555", linewidth=0.8, alpha=0.5)
            ax.set_ylabel(f"{parameter}: RMSE / equidistant RMSE")
            if row == 0:
                ax.set_title(f"theta × nominal mean gap = {gap:g}", fontsize=10)
    path = output_dir / "week3_rmse_ratio_vs_cv.png"
    _finish(fig, handles, "RMSE relative to the same nominal-design equidistant twin", path)
    paths.append(path)

    # summarize_results has the same fit-level counts on each parameter row;
    # select theta once, rather than triple-counting the same planned fits.
    failures = data.loc[data["parameter"].eq("theta")]
    fig, axes = plt.subplots(1, len(coarseness), squeeze=False,
                             figsize=(max(8.5, 4.0 * len(coarseness)), 4.1))
    for column, gap in enumerate(coarseness):
        ax = axes[0, column]
        subset = failures.loc[failures["theta_mean_gap"].eq(gap)]
        _draw_curves(ax, subset, "failure_rate", sample_sizes)
        ax.set_ylim(-0.025, max(0.1, float(subset["failure_rate"].max()) * 1.15))
        ax.set_ylabel("Failed / all planned fits")
        ax.set_title(f"theta × nominal mean gap = {gap:g}", fontsize=10)
    path = output_dir / "week3_fit_failure_rate.png"
    _finish(fig, handles, "Fit failures retained in the Monte Carlo denominator", path)
    paths.append(path)
    return paths
