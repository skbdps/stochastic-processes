"""Publication layouts built only from verified, saved experiment records.

These figures have a new renderer identity. They never change a historical
figure, numerical interval, eligibility decision, or practical classification.
"""
from io import BytesIO

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
import numpy as np


METHODS = (("exact", "Exact", "#0072B2", "o"),
           ("pfml", "PFML", "#D55E00", "s"),
           ("euler", "Euler", "#009E73", "^"))
PARAMETERS = (("theta", "θ"), ("mu", "μ"), ("sigma", "σ"))
STYLE = {**dict(matplotlib.rcParamsDefault), "backend": "Agg",
         "font.family": "DejaVu Sans", "font.size": 10,
         "axes.titlesize": 11, "axes.spines.top": False,
         "axes.spines.right": False, "axes.axisbelow": True,
         "figure.facecolor": "white", "pdf.fonttype": 42,
         "ps.fonttype": 42, "svg.hashsalt": "ou-week9-publication-v1"}
CLASSES = {"within_criteria": 0, "supported_breach": 1, "unresolved": 2}
CLASS_COLORS = ("#C9E4EE", "#EFB6A1", "#FFF0B0")
STEMS = ("F1_distributions_main", "S_F1_distributions_cv0_to_1", "S_F1_distributions_cv1p5_to_3",
         "F7_classification_n250", "S_calibration_e2_n250", "F7_classification_n1000", "S_calibration_e2_n1000",
         "F7_classification_n4000", "S_calibration_e2_n4000", "F8_matched_mean", "S_calibration_e4_mean")


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _legend(fig, *, y=.91, curves=False):
    handles = [Line2D([], [], color=color, marker=marker, label=label,
                      linewidth=1.6, markersize=5, linestyle=("-", "--", "-.")[i] if curves else "none")
               for i, (_, label, color, marker) in enumerate(METHODS)]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(.5, y),
               ncol=3, frameon=False)


def ecdf_figure(frame, cvs, title):
    """All retained eligible estimates in the declared CV selection, no KDE."""
    fig, axes = plt.subplots(len(cvs), 3, figsize=(12.6, 2.7*len(cvs)+2.1),
                             squeeze=False)
    fig.suptitle(title, fontsize=15, y=.98)
    fig.text(.5, .945, "Registered E1 | q = 0.5 | n = 1,000 | R = 2,000 per CV",
             ha="center", fontsize=10)
    _legend(fig, y=.925, curves=True)
    truths = {"theta": 1., "mu": 0., "sigma": .5}
    for row, cv in enumerate(cvs):
        for col, (parameter, symbol) in enumerate(PARAMETERS):
            ax = axes[row, col]
            group = frame[(frame.cv == cv) & (frame.parameter == parameter)]
            for method_index, (method, label, color, _) in enumerate(METHODS):
                values = group[group.estimator == method].sort_values(
                    ["estimate", "rep"], kind="stable")
                _require(len(values) > 0, "ECDF selected group is empty")
                _require(values.estimate.notna().all(), "ECDF estimate is missing")
                _require(values.n_point_valid.nunique() == 1
                         and int(values.n_point_valid.iloc[0]) == len(values),
                         "ECDF excludes eligible estimates")
                _require(np.allclose(values.ecdf, np.arange(1, len(values)+1)/len(values),
                                     rtol=0, atol=1e-15), "ECDF ranks changed")
                artist, = ax.step(values.estimate, values.ecdf, where="post",
                                  color=color, linewidth=1.5, label=label,
                                  linestyle=("-", "--", "-.")[method_index])
                artist.set_gid(f"ecdf:{cv:g}:{parameter}:{method}")
            ax.axvline(truths[parameter], color="#444444", linestyle=":", linewidth=1)
            ax.set_title(f"Nominal CV = {cv:g} | {symbol}")
            ax.set_xlabel(f"{symbol} estimate" + (" (original state units)" if parameter == "mu" else ""))
            if col == 0:
                ax.set_ylabel("Empirical cumulative probability")
            ax.set_ylim(0, 1.02)
            ax.grid(alpha=.18)
            counts = group.groupby("estimator")[["n_point_valid", "n_planned"]].first()
            _require(len(counts) == 3, "ECDF method missing")
            note = "; ".join(f"{m}: {int(r.n_point_valid)}/{int(r.n_planned)}"
                             for m, r in counts.iterrows())
            ax.text(.02, .035, note, transform=ax.transAxes, fontsize=7.4,
                    bbox={"facecolor": "white", "edgecolor": "none", "alpha": .85})
    fig.subplots_adjust(left=.075, right=.985, bottom=.12, top=.865,
                        hspace=.55, wspace=.30)
    fig.text(.075, .035,
             "Dotted lines: simulation truth. Curves use every point-valid saved estimate; counts are valid/planned.\n"
             "Descriptive distributions without confidence bands. Display selection is post-result presentation; all seven CVs are supplied.\n"
             "Regular-grid exact/PFML estimates coincide and their curves overlap. Line styles distinguish coincident methods.",
             fontsize=8.5, va="bottom")
    return fig


def classification_figure(frame, n):
    """Discrete complete-grid maps; no contour or interpolation convention."""
    fig, axes = plt.subplots(2, 3, figsize=(12.6, 8.3), squeeze=False)
    fig.suptitle(f"F7  Practical classifications | n = {n:,}", fontsize=15, y=.975)
    fig.text(.5, .929, "Registered E2 | R = 2,000 per cell | θ and σ only | pointwise Monte Carlo classification",
             ha="center", fontsize=10)
    cvs, qs = [0., .5, 1., 2., 3.], [.1, .5, 1., 2.]
    for row, (parameter, symbol) in enumerate((("theta", "θ"), ("sigma", "σ"))):
        for col, (method, label, _, _) in enumerate(METHODS):
            ax = axes[row, col]
            group = frame[(frame.n == n) & (frame.parameter == parameter)
                          & (frame.estimator == method)]
            _require(len(group) == 20, "classification map is not the complete 4 by 5 grid")
            image = np.empty((4, 5))
            for y, q in enumerate(qs):
                for x, cv in enumerate(cvs):
                    selected = group[(group.theta_mean_gap == q) & (group.cv == cv)]
                    _require(len(selected) == 1, "classification grid key is ambiguous")
                    value = selected.iloc[0]
                    code = value.practical_classification
                    _require(code in CLASSES, "unexpected practical classification")
                    image[y, x] = CLASSES[code]
                    annotation = {"within_criteria": "W", "unresolved": "U"}.get(code)
                    if annotation is None:
                        driver = str(value.breach_driver)
                        _require(driver in ("bias", "coverage", "both", "bias_and_coverage"), "breach driver is missing")
                        annotation = "B:"+{"bias": "b", "coverage": "c", "both": "bc", "bias_and_coverage": "bc"}[driver]
                    text = ax.text(x, y, annotation, ha="center", va="center",
                                   fontsize=10, color="#182027")
                    text.set_gid(f"class:{n}:{method}:{parameter}:{q:g}:{cv:g}")
            ax.imshow(image, origin="lower", vmin=-.5, vmax=2.5,
                      cmap=ListedColormap(CLASS_COLORS), interpolation="nearest", aspect="auto")
            ax.set_xticks(range(5), [f"{v:g}" for v in cvs])
            ax.set_yticks(range(4), [f"{v:g}" for v in qs])
            ax.set_xticks(np.arange(-.5, 5, 1), minor=True)
            ax.set_yticks(np.arange(-.5, 4, 1), minor=True)
            ax.grid(which="minor", color="white", linewidth=2)
            ax.tick_params(which="minor", bottom=False, left=False)
            ax.set_xlabel("Nominal gap CV")
            if col == 0:
                ax.set_ylabel("q = θ × nominal mean gap")
            ax.set_title(f"{label} | {symbol}")
    fig.subplots_adjust(left=.075, right=.985, bottom=.23, top=.855,
                        wspace=.26, hspace=.54)
    fig.text(.075, .055,
             "W: within criteria; B:b: bias breach; B:c: coverage breach; B:bc: both; U: unresolved.\n"
             "Criteria: operational coverage ≥80% and absolute relative bias ≤10%, assessed with the registered Monte Carlo intervals.\n"
             "These broad criteria do not certify 95% coverage. Cells are categorical design points; no continuous boundary is estimated.\n"
             "Natural-scale Wald intervals use full observed Hessians; PFML/Euler covariance is model-based. Flooring applies at CV ≥1.5.",
             fontsize=8.5, va="bottom")
    return fig


def matched_figure(frame):
    fig, axes = plt.subplots(1, 2, figsize=(12.6, 6.7))
    fig.suptitle("F8  Mean estimation under clustered observation times", fontsize=15, y=.975)
    fig.text(.5, .922, "Exploratory E4 | n = 1,000 | q = 0.5 | R = 2,000 independent paths per cell",
             ha="center", fontsize=10)
    _legend(fig, y=.905)
    metrics = (("rmse_ratio", "MMPP / IID RMSE", 1., 1.),
               ("operational_coverage_difference", "MMPP − IID coverage (percentage points)", 0., 100.))
    for ax, (metric, ylabel, null, scale) in zip(axes, metrics):
        for i, (method, label, color, marker) in enumerate(METHODS):
            group = frame[(frame.parameter == "mu") & (frame.estimator == method)].sort_values("ratio")
            _require(group.ratio.tolist() == [1, 5, 20], "matched μ contrasts are incomplete")
            values = group[metric].to_numpy()*scale
            lower, upper = (group[metric+suffix].to_numpy()*scale
                            for suffix in ("_mc_lower", "_mc_upper"))
            _require(np.all(lower <= values) and np.all(values <= upper), "invalid contrast interval")
            series = ax.errorbar(np.arange(3)+(i-1)*.14, values,
                                 yerr=np.vstack((values-lower, upper-values)),
                                 color=color, marker=marker, linestyle="none", capsize=4,
                                 label=label)
            series.lines[0].set_gid(f"matched:mu:{metric}:{method}")
        ax.axhline(null, color="#555555", linestyle=":", linewidth=1)
        ax.set_xticks(range(3), ["r = 1", "r = 5", "r = 20"])
        ax.set_xlim(-.4, 2.4)
        ax.set_xlabel("Intensity ratio of the matched gap-law pair")
        ax.set_ylabel(ylabel)
        ax.grid(axis="y", alpha=.2)
    fig.subplots_adjust(left=.085, right=.98, bottom=.28, top=.80, wspace=.33)
    fig.text(.075, .065,
             "Bars: pointwise 95% Monte Carlo intervals, using independent-cell formulas; positive RMSE-ratio limits use the log-delta method.\n"
             "MMPP and IID share the entire one-gap marginal and declared floor rule. Their contrast includes realized-span consequences.\n"
             "PFML μ shows increased RMSE and reduced coverage at r=5 and r=20; the comparison does not establish a mean-bias shift.\n"
             "PFML/Euler Wald covariance remains model-based. All parameters and matched contrasts are retained in the supplement.",
             fontsize=8.5, va="bottom")
    return fig


def calibration_grid(frame, n):
    fig, axes = plt.subplots(3, 3, figsize=(12.6, 10.8), squeeze=False)
    fig.suptitle(f"S  Empirical variability / reported uncertainty | E2 n = {n:,}", fontsize=14, y=.978)
    fig.text(.5, .936, "Descriptive ratio: empirical SD / RMS reported SE, using the same interval-valid subset",
             ha="center", fontsize=10)
    cvs, qs = [0., .5, 1., 2., 3.], [.1, .5, 1., 2.]
    for row, (parameter, symbol) in enumerate(PARAMETERS):
        for col, (method, label, _, _) in enumerate(METHODS):
            ax = axes[row, col]
            group = frame[(frame.study == "e2") & (frame.n == n)
                          & (frame.parameter == parameter) & (frame.estimator == method)]
            _require(len(group) == 20, "calibration grid incomplete")
            values = np.empty((4, 5))
            for y, q in enumerate(qs):
                for x, cv in enumerate(cvs):
                    selected = group[(group.theta_mean_gap == q) & (group.cv == cv)]
                    _require(len(selected) == 1, "calibration key ambiguous")
                    values[y, x] = selected.sd_over_rms_se.iloc[0]
            _require(np.isfinite(values).all() and (values > 0).all(), "calibration ratio unavailable")
            # A fixed symmetric log2 range treats half and double equally.
            ax.imshow(np.log2(values), origin="lower", vmin=-1, vmax=1,
                      cmap="RdBu_r", interpolation="nearest", aspect="auto")
            for y in range(4):
                for x in range(5):
                    text = ax.text(x, y, f"{values[y,x]:.2f}", ha="center", va="center",
                                   color="white" if abs(np.log2(values[y,x])) > .65 else "#182027", fontsize=8.5)
                    text.set_gid(f"calibration:{n}:{method}:{parameter}:{qs[y]:g}:{cvs[x]:g}")
            ax.set_xticks(range(5), [f"{v:g}" for v in cvs])
            ax.set_yticks(range(4), [f"{v:g}" for v in qs])
            ax.set_xlabel("Nominal gap CV")
            if col == 0:
                ax.set_ylabel("q")
            ax.set_title(f"{label} | {symbol}")
    fig.subplots_adjust(left=.07, right=.985, bottom=.14, top=.875, hspace=.55, wspace=.28)
    fig.text(.07, .035,
             "Cell text gives the full ratio rounded to two decimals. Blue: <1; white: 1; red: >1; color saturates at 0.5 and 2.\n"
             "This is a descriptive calibration diagnostic, without a derived ratio confidence interval or long-run covariance benchmark.\n"
             "It does not explain coverage by itself: bias, distribution shape and interval construction also matter. Counts are in calibration.csv.",
             fontsize=8.5, va="bottom")
    return fig


def calibration_e4(frame):
    fig, axes = plt.subplots(1, 2, figsize=(12.6, 6.6))
    fig.suptitle("S  Mean-estimator variability and reported uncertainty | E4", fontsize=14, y=.975)
    fig.text(.5, .925, "Exploratory E4 | μ in original state units | interval-valid subsets (2,000 per cell/method)",
             ha="center", fontsize=10)
    _legend(fig, y=.905)
    group = frame[(frame.study == "e4") & (frame.parameter == "mu")
                  & (frame.mechanism.isin(["mmpp", "iid"]))]
    order = [(r, m) for r in (1, 5, 20) for m in ("mmpp", "iid")]
    for ax, metric, ylabel in zip(axes, ("empirical_sd_interval_valid", "rms_reported_se"),
                                  ("Empirical SD", "RMS reported SE")):
        for i, (method, label, color, marker) in enumerate(METHODS):
            values = []
            for ratio, mechanism in order:
                selected = group[(group.estimator == method) & (group.ratio == ratio)
                                 & (group.mechanism == mechanism)]
                _require(len(selected) == 1, "E4 calibration group incomplete")
                values.append(selected[metric].iloc[0])
            line, = ax.plot(np.arange(6)+(i-1)*.12, values, marker=marker,
                            linestyle="none", color=color, label=label)
            line.set_gid(f"calibration:e4:mu:{metric}:{method}")
        ax.set_xticks(range(6), [f"{m.upper()}\nr={r}" for r,m in order], fontsize=8.5)
        ax.set_ylabel(ylabel)
        ax.grid(axis="y", alpha=.2)
    # Equal y limits permit a direct visual comparison of the two spreads.
    lower = min(ax.get_ylim()[0] for ax in axes)
    upper = max(ax.get_ylim()[1] for ax in axes)
    for ax in axes:
        ax.set_ylim(lower, upper)
    fig.subplots_adjust(left=.08, right=.98, bottom=.25, top=.80, wspace=.28)
    fig.text(.075, .06,
             "Points are saved descriptive summaries, without Monte Carlo intervals for these SD/SE quantities. Both panels use identical y limits.\n"
             "Empirical SD and RMS reported SE use the same interval-valid subset. All parameters/designs are retained in calibration.csv.\n"
             "The working-model uncertainty does not constitute a robust dependence-aware covariance calculation.",
             fontsize=8.5, va="bottom")
    return fig


def build_figures(tables):
    """Return new publication layouts; caller owns serialization and closing."""
    figures = {}
    with plt.rc_context(STYLE):
        frame = tables["f1_ecdf.csv"]
        figures["F1_distributions_main"] = ecdf_figure(frame, (0., 1., 3.), "F1  Estimator distributions at three E1 designs")
        figures["S_F1_distributions_cv0_to_1"] = ecdf_figure(frame, (0., .25, .5, 1.), "S  Complete E1 distributions | nominal CV 0 to 1")
        figures["S_F1_distributions_cv1p5_to_3"] = ecdf_figure(frame, (1.5, 2., 3.), "S  Complete E1 distributions | nominal CV 1.5 to 3")
        for n in (250, 1000, 4000):
            figures[f"F7_classification_n{n}"] = classification_figure(tables["T1_e2.csv"], n)
            figures[f"S_calibration_e2_n{n}"] = calibration_grid(tables["calibration.csv"], n)
        figures["F8_matched_mean"] = matched_figure(tables["e4_matched.csv"])
        figures["S_calibration_e4_mean"] = calibration_e4(tables["calibration.csv"])
    _require(tuple(figures) == STEMS, "publication figure set changed")
    return figures


def figure_bytes(fig, extension, title):
    _require(extension in ("png", "pdf"), "unsupported publication figure format")
    metadata = {"Title": title, "Creator": "OU Week 9 saved-evidence publication renderer",
                "CreationDate": None, "ModDate": None} if extension == "pdf" else {
                    "Software": "OU Week 9 saved-evidence publication renderer"}
    stream = BytesIO()
    with plt.rc_context(STYLE):
        fig.savefig(stream, format=extension, dpi=170, metadata=metadata)
    return stream.getvalue()
