"""Figures for the Week 4 checkpoint; no dependency on Week 3 plot schemas."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


class DiagnosticFigureWriter:
    def write(self, summary, bootstrap_summary, output_dir):
        output = Path(output_dir)
        output.mkdir(parents=True, exist_ok=True)
        paths = []
        plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
        data = summary.set_index("parameter").loc[["theta", "mu", "sigma"]]
        fig, ax = plt.subplots(figsize=(8, 5))
        x = np.arange(3)
        for shift, scope, color, label in ((-.055, "conditional", "#2166ac", "Among valid intervals"),
                                            (.055, "operational", "#d6604d", "All planned paths")):
            value = data[f"coverage_{scope}"].to_numpy()
            low, high = data[f"{scope}_wilson_lower"].to_numpy(), data[f"{scope}_wilson_upper"].to_numpy()
            ax.errorbar(x + shift, value, yerr=[value - low, high - value], fmt="o", capsize=5,
                        color=color, label=label)
        ax.axhspan(.94, .96, color="grey", alpha=.12, label="94–96% diagnostic reference")
        ax.axhline(.95, color="grey", lw=1, ls="--")
        ax.set(xticks=x, xticklabels=["theta", "mu", "sigma"], ylabel="Estimated coverage",
               title="Exact OU regular checkpoint · 95% natural-scale Wald intervals")
        lower = np.nanmin(data[["conditional_wilson_lower", "operational_wilson_lower"]].to_numpy())
        upper = np.nanmax(data[["conditional_wilson_upper", "operational_wilson_upper"]].to_numpy())
        ax.set_ylim(max(0, min(.94, lower) - .025), min(1, max(.96, upper) + .025))
        ax.grid(axis="y", alpha=.2)
        ax.legend(loc="lower left", frameon=False)
        n, valid = int(data.n_planned.iloc[0]), int(data.n_interval_valid.iloc[0])
        fig.text(.5, .025, f"{n} planned paths; {valid} valid intervals. Bars: 95% Wilson Monte Carlo intervals.\n"
                 "The reference band triggers investigation; it is not a target to force.", ha="center", fontsize=9)
        fig.tight_layout(rect=(0, .1, 1, 1))
        path = output / "week4_checkpoint_coverage.png"
        fig.savefig(path, dpi=160)
        plt.close(fig)
        paths.append(path)

        cells = bootstrap_summary[["cell_id", "cv"]].drop_duplicates().sort_values("cv")
        fig, axes = plt.subplots(len(cells), 3, figsize=(12, 3.8 * len(cells)), squeeze=False)
        names = ["exact", "pfml", "euler"]
        for row, cell in enumerate(cells.itertuples(index=False)):
            cell_rows = bootstrap_summary.loc[bootstrap_summary.cell_id == cell.cell_id].set_index("estimator").reindex(names)
            for col, parameter in enumerate(("theta", "mu", "sigma")):
                ax = axes[row, col]
                for shift, prefix, color, label in ((-.10, "wald", "#2166ac", "Wald"),
                                                    (.10, "bootstrap", "#d6604d", "Bootstrap percentile")):
                    lower = cell_rows[f"{prefix}_{parameter}_lower"].to_numpy(dtype=float)
                    upper = cell_rows[f"{prefix}_{parameter}_upper"].to_numpy(dtype=float)
                    center = (lower + upper) / 2
                    ax.errorbar(x + shift, center, yerr=(upper - lower) / 2, fmt="o", capsize=4, color=color, label=label)
                ax.scatter(x, cell_rows[parameter], color="black", marker="x", s=25, label="Point estimate", zorder=4)
                ax.set(xticks=x, xticklabels=["Exact", "PFML", "Euler"], ylabel=parameter,
                       title=f"CV={cell.cv:g} · {parameter}")
                ax.grid(axis="y", alpha=.2)
        handles, labels = axes[0, 0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(.5, .965), ncol=3, frameon=False)
        fig.suptitle(f"Fixed-path interval diagnostics · {len(cells)} prespecified cells", y=.999)
        draws = int(bootstrap_summary.planned_draws.iloc[0])
        fig.text(.5, .02, f"{draws} conditional bootstrap draws per method: timestamps and observed initial state fixed.\n"
                 "PFML/Euler bootstrap is diagnostic under misspecification. These plots do not measure bootstrap coverage.",
                 ha="center", fontsize=9)
        fig.tight_layout(rect=(0, .09, 1, .92))
        path = output / "week4_bootstrap_spotchecks.png"
        fig.savefig(path, dpi=160)
        plt.close(fig)
        paths.append(path)
        return paths
