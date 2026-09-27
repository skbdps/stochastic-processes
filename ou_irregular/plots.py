"""Headless exploratory figures for the Week 3 simulation harness.

These figures describe point-valid fits and separately expose every exclusion.
Approximate Monte Carlo error bars describe uncertainty across simulated paths,
not parameter confidence intervals for an individual path.
They deliberately do not plot coverage: Week 3 supplies no confidence-interval
estimator, and missing coverage is not evidence of zero or nominal coverage.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

from .metrics import VALIDITY_POLICY, add_twin_comparisons


_COLORS = {"exact": "#2166ac", "pfml": "#d6604d", "euler": "#1b9e77"}
_LABELS = {"exact": "Exact MLE", "pfml": "Mean-gap PFML", "euler": "Euler"}
_LINESTYLES = ("-", "--", ":", "-.")


def _draw_curves(ax, data: pd.DataFrame, metric: str, sample_sizes: list,
                 mcse_column: str | None = None) -> None:
    """Preserve missing values; draw +/-1.96 MCSE only where it is estimable.

    Normal error bars are approximate Monte Carlo intervals for the plotted
    summary, not confidence intervals for a path estimate. They are descriptive
    at 20 replications and are not used as a Week 4 coverage calculation.
    """
    for (estimator, n), group in data.groupby(["estimator", "n"], sort=True):
        group = group.sort_values("cv")
        ax.plot(
            group["cv"], group[metric], marker="o", markersize=4,
            color=_COLORS.get(estimator, "#666666"),
            linestyle=_LINESTYLES[sample_sizes.index(n) % len(_LINESTYLES)],
            linewidth=1.6,
        )
        if mcse_column is not None and mcse_column in group:
            selected = np.isfinite(group[metric]) & np.isfinite(group[mcse_column])
            errors = group.loc[selected]
            if not errors.empty:
                ax.errorbar(errors["cv"], errors[metric],
                            yerr=1.96 * errors[mcse_column], fmt="none", capsize=3,
                            color=_COLORS.get(estimator, "#666666"), alpha=0.65,
                            elinewidth=1.0)
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


def _finish(fig, handles: list, title: str, path: Path, *, classified=False,
            uncertainty=False) -> None:
    fig.suptitle(title, fontsize=13, y=0.995)
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.957),
               ncol=min(6, len(handles)), frameon=False, fontsize=9)
    policy = ("Point metrics conditional on classified valid interior fits." if classified
              else "Historical policy: success and finite estimates.")
    uncertainty_note = (" Approx. bars: +/-1.96 MCSE (simulation uncertainty; exploratory small sample)."
                        if uncertainty else "")
    fig.text(0.5, 0.012,
             f"Week 3 exploratory results. {policy} No path CI coverage computed.\n"
             + uncertainty_note.strip(), ha="center", fontsize=8, color="#555555")
    fig.tight_layout(rect=(0.02, 0.05, 0.99, 0.89))
    fig.savefig(path, dpi=170, bbox_inches="tight")
    plt.close(fig)


def plot_summary(summary: pd.DataFrame, output_dir: str | Path) -> list[Path]:
    """Write bias/MCSE, same-estimator control ratios and exclusion-rate PNGs.

    Classified v2 summaries also get direct RMSE with +/-1.96 MCSE bars. V1
    records retain the original three-file surface and historical validity rule.
    Different CV cells are independent; no paired ratio uncertainty is plotted.

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
    policies = set(data["validity_policy"]) if "validity_policy" in data else set()
    if len(policies) > 1:
        raise ValueError("plot_summary cannot mix historical and classified validity policies")
    classified = policies == {VALIDITY_POLICY}
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
        mcse = "mcse_bias" if parameter == "mu" else "relative_mcse_bias"
        for column, gap in enumerate(coarseness):
            ax = axes[row, column]
            subset = data.loc[data["parameter"].eq(parameter) & data["theta_mean_gap"].eq(gap)]
            _draw_curves(ax, subset, metric, sample_sizes, mcse)
            ax.axhline(0, color="#555555", linewidth=0.8, alpha=0.5)
            ax.set_ylabel(f"{parameter}: {'absolute' if parameter == 'mu' else 'relative'} bias")
            if row == 0:
                ax.set_title(f"theta × nominal mean gap = {gap:g}", fontsize=10)
    path = output_dir / "week3_bias_vs_cv.png"
    _finish(fig, handles, "Parameter bias with approximate Monte Carlo uncertainty", path,
            classified=classified, uncertainty=True)
    paths.append(path)

    fig, axes = plt.subplots(3, len(coarseness), squeeze=False,
                             figsize=(max(8.5, 4.0 * len(coarseness)), 9.1))
    for row, parameter in enumerate(("theta", "mu", "sigma")):
        for column, gap in enumerate(coarseness):
            ax = axes[row, column]
            subset = data.loc[data["parameter"].eq(parameter) & data["theta_mean_gap"].eq(gap)]
            _draw_curves(ax, subset, "rmse_ratio", sample_sizes)
            ax.axhline(1, color="#555555", linewidth=0.8, alpha=0.5)
            # The shared title names the same-estimator CV=0 denominator;
            # repeating it on every axis would overlap adjacent panel labels.
            ax.set_ylabel(f"{parameter}: RMSE ratio")
            if row == 0:
                ax.set_title(f"theta × nominal mean gap = {gap:g}", fontsize=10)
    path = output_dir / "week3_rmse_ratio_vs_cv.png"
    _finish(fig, handles, "RMSE ratio: same estimator, same nominal design, regular-grid control", path,
            classified=classified)
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
        ax.set_ylabel("Excluded / all planned fits" if classified else "Failed / all planned fits")
        ax.set_title(f"theta × nominal mean gap = {gap:g}", fontsize=10)
    path = output_dir / "week3_fit_failure_rate.png"
    _finish(fig, handles, "All excluded fits retained in the Monte Carlo denominator", path,
            classified=classified)
    paths.append(path)
    if classified:
        fig, axes = plt.subplots(3, len(coarseness), squeeze=False,
                                 figsize=(max(8.5, 4.0 * len(coarseness)), 9.1))
        for row, parameter in enumerate(("theta", "mu", "sigma")):
            metric = "rmse" if parameter == "mu" else "relative_rmse"
            mcse = "mcse_rmse" if parameter == "mu" else "relative_mcse_rmse"
            for column, gap in enumerate(coarseness):
                ax = axes[row, column]
                subset = data.loc[data["parameter"].eq(parameter) & data["theta_mean_gap"].eq(gap)]
                _draw_curves(ax, subset, metric, sample_sizes, mcse)
                ax.axhline(0, color="#555555", linewidth=0.8, alpha=0.5)
                ax.set_ylabel(f"{parameter}: {'absolute' if parameter == 'mu' else 'relative'} RMSE")
                if row == 0:
                    ax.set_title(f"theta × nominal mean gap = {gap:g}", fontsize=10)
        path = output_dir / "week3_rmse_vs_cv.png"
        _finish(fig, handles, "RMSE with approximate Monte Carlo uncertainty", path,
                classified=True, uncertainty=True)
        paths.append(path)
    return paths
