"""Deterministic Week 7 figures from bound saved evidence, without new fits.

Pure builders keep plotted values inspectable. Archive adapters separately
enforce audit/input identity and immutable rendering; renderer provenance is
not retroactively inserted into an earlier raw-run dependency manifest.
"""
from io import BytesIO
from dataclasses import dataclass, replace
from importlib.metadata import version
from pathlib import Path
import gzip
import json
import math
import re
import subprocess
import textwrap
import zlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, MaxNLocator
import numpy as np
import pandas as pd

from ..experiments.design import Cell
from ..experiments import audit as ordinary_audit
from ..experiments.provenance import ROOT, CORE_FILES, canonical, capture, digest_bytes, file_digest, verify
from ..experiments.storage import atomic_bytes, atomic_json, run_lock
from . import execution, floor_audit, path_diagnostics, reporting
from .execution import require_design
from .path_diagnostics import _ManifestGuard, _check_incident_block, _reject_archive_ancestors


FILES = tuple(sorted(set(floor_audit.FILES) | {
    "ou_irregular/week7/plots.py", "docs/week7_plots.md"}))
CONTRACT = "week7_saved_evidence_png_pdf_v1"
FOLLOWUP_STEMS = ("W7_followup_coverage", "W7_followup_availability")
FLOOR_STEMS = ("W7_floor_paired_bias", "W7_floor_paired_coverage", "W7_floor_rmse_ratio", "W7_floor_spacing")
SELECTED_CASES = ((.5, 465, "W7_selected_E2_cv0p5_rep000465"),
                  (0., 1424, "W7_selected_E2_cv0_rep001424"))
METHODS = (("exact", "Exact", "#0072B2", "o"), ("pfml", "PFML", "#D55E00", "s"),
           ("euler", "Euler", "#009E73", "^"))
PARAMETERS = (("theta", r"$\theta$"), ("mu", r"$\mu$"), ("sigma", r"$\sigma$"))
FLOORS = ((1e-8, r"$\epsilon=10^{-8}$", "#0072B2", "o"),
          (1e-6, r"$\epsilon=10^{-6}$", "#333333", "s"),
          (1e-4, r"$\epsilon=10^{-4}$", "#D55E00", "^"))
STYLE = {"font.family": "DejaVu Sans", "font.size": 9.5, "axes.titlesize": 11,
         "axes.spines.top": False, "axes.spines.right": False,
         "figure.facecolor": "white", "axes.grid": True, "grid.alpha": .18,
         "axes.axisbelow": True, "pdf.fonttype": 42, "ps.fonttype": 42}
# Every effective setting starts from Matplotlib's versioned defaults. Ambient
# notebook rc mutations (TeX, transparency, fonts, DPI) cannot change evidence.
RENDER_SETTINGS = {**dict(matplotlib.rcParamsDefault), **STYLE, "backend": "Agg"}
Z = 1.959963984540054
CV_AXIS_LABEL = "Nominal CV (equally spaced design cells)"
SPACING_FIELDS = ("realized_mean_gap", "realized_cv", "realized_span", "clip_fraction", "min_gap")


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _numbers(values):
    result = np.asarray([np.nan if value is None else value for value in values], dtype=float)
    _require(not np.isinf(result).any(), "infinite plotted value is not a missing-value marker")
    return result


def _keys(table, names, expected, label):
    _require(set(names) <= set(table), "missing figure key columns: " + label)
    keys = list(table[list(names)].itertuples(index=False, name=None))
    _require(len(keys) == len(expected) and len(set(keys)) == len(keys) and set(keys) == expected,
             "incomplete or duplicate figure key set: " + label)


def _count(value, maximum, label):
    _require(isinstance(value, (int, float, np.integer, np.floating))
             and not isinstance(value, (bool, np.bool_)) and np.isfinite(value)
             and int(value) == value and 0 <= value <= maximum, "invalid count: " + label)
    return int(value)


def _coordinates(table, design, label):
    _require({"cell_id", "n", "theta_mean_gap", "cv"} <= set(table), "missing table coordinate fields: " + label)
    cells = {cell.identity: cell for cell in design.cells}
    for row in table.itertuples(index=False):
        _require(row.cell_id in cells, "unknown table cell: " + label)
        cell = cells[row.cell_id]
        actual = (row.n, row.theta_mean_gap, row.cv)
        _require(all(isinstance(value, (int, float, np.integer, np.floating)) and not isinstance(value, (bool, np.bool_))
                     and np.isfinite(value) for value in actual)
                 and actual == (cell.n, cell.q, cell.cv), "table canonical coordinates changed: " + label)


def _spacing_inputs(tables, design):
    spacing = tables["spacing_summary.csv"].copy()
    columns = {"cell_id", "n", "theta_mean_gap", "cv", "n_planned",
               *(field + "_" + suffix for field in SPACING_FIELDS for suffix in ("n_available", "mean", "mcse", "min", "max"))}
    _require(set(spacing) == columns, "spacing table header fields changed")
    _keys(spacing, ("cell_id",), {(cell.identity,) for cell in design.cells}, "spacing summary")
    _coordinates(spacing, design, "spacing")
    for row in spacing.to_dict("records"):
        _require(_count(row["n_planned"], design.replications, "spacing planned") == design.replications,
                 "spacing source count is not R")
        for field in SPACING_FIELDS:
            _require(_count(row[field+"_n_available"], design.replications, "spacing available") == design.replications,
                     "each plotted spacing mean must use all R saved sources")
            mean, low, high, se = [row[field+"_"+suffix] for suffix in ("mean", "min", "max", "mcse")]
            _require(np.isfinite([mean, low, high]).all() and low <= mean <= high,
                     "saved spacing extrema/mean are unavailable or inconsistent")
            _require((pd.isna(se) if design.replications < 2 else np.isfinite(se) and se >= 0),
                     "spacing Monte Carlo uncertainty has inconsistent availability")
    mapping = {cell.identity: cell for cell in design.cells}
    spacing["floor_fraction"] = spacing.cell_id.map(lambda cid: mapping[cid].floor_fraction)
    return spacing


def _standard_inputs(tables, design):
    """Validate complete saved summary keys before any renderer-state creation."""
    _require(set(tables) == set(floor_audit.TABLE_KEYS), "the six standard summary tables are required")
    summary = _summary(tables, design)
    single = {(cell.identity, method) for cell in design.cells for method, *_ in METHODS}
    parameter = {(cid, method, param) for cid, method in single for param, _ in PARAMETERS}
    _keys(tables["regular_comparisons.csv"], ("cell_id", "estimator", "parameter"), parameter, "regular comparisons")
    pairs = (("exact", "pfml"), ("exact", "euler"), ("pfml", "euler"))
    paired = {(cell.identity, a, b, param) for cell in design.cells for a, b in pairs for param, _ in PARAMETERS}
    _keys(tables["paired_comparisons.csv"], ("cell_id", "estimator_a", "estimator_b", "parameter"), paired, "paired comparisons")
    _keys(tables["numerical_diagnostics.csv"], ("cell_id", "estimator"), single, "numerical diagnostics")
    for name, table in tables.items():
        _coordinates(table, design, name)
    _spacing_inputs(tables, design)
    outcomes = tables["outcome_counts.csv"]
    required = {"cell_id", "n", "theta_mean_gap", "cv", "estimator", "point_status", "interval_status",
                "count", "n_planned", "fraction", "wilson_lower", "wilson_upper"}
    _require(set(outcomes) == required and not outcomes.duplicated(list(floor_audit.TABLE_KEYS["outcome_counts.csv"])).any(),
             "outcome table header or unique status keys changed")
    _require(set(zip(outcomes.cell_id, outcomes.estimator)) == single, "outcome table omits a complete cell/method group")
    for row in outcomes.to_dict("records"):
        _require(_count(row["n_planned"], design.replications, "outcome planned") == design.replications
                 and _count(row["count"], design.replications, "outcome count") > 0,
                 "outcome table has invalid denominators")
        _require(all(type(row[key]) is str and bool(row[key]) for key in ("point_status", "interval_status")),
                 "outcome status must be a nonempty string")
    _require(outcomes.groupby(["cell_id", "estimator"])["count"].sum().eq(design.replications).all(),
             "outcome status counts do not sum to R in every group")
    joint = summary[summary.parameter == "theta"].set_index(["cell_id", "estimator"])
    numerical = tables["numerical_diagnostics.csv"]
    required_numerical = {"cell_id", "estimator", "n", "theta_mean_gap", "cv", "n_planned", "n_interval_valid",
                          *(field+"_"+suffix for field in ordinary_audit.NUMERICAL_FIELDS
                            for suffix in ("n_finite", "n_nonfinite", "max", "median", "p95"))}
    _require(set(numerical) == required_numerical, "numerical table header fields changed")
    for row in numerical.to_dict("records"):
        _require(_count(row["n_planned"], design.replications, "numerical planned") == design.replications
                 and _count(row["n_interval_valid"], design.replications, "numerical interval count")
                 == joint.loc[(row["cell_id"], row["estimator"]), "n_interval_valid"],
                 "numerical table denominators changed")
        for field in ordinary_audit.NUMERICAL_FIELDS:
            _require(_count(row[field+"_n_finite"], design.replications, "finite diagnostics")
                     + _count(row[field+"_n_nonfinite"], design.replications, "nonfinite diagnostics") == design.replications,
                     "numerical diagnostic counts do not sum to R")
    return summary


def _summary(tables, design):
    data = tables["parameter_summary.csv"]
    expected = {(c.identity, method, parameter) for c in design.cells for method, *_ in METHODS
                for parameter, _ in PARAMETERS}
    _keys(data, ("cell_id", "estimator", "parameter"), expected, "parameter summary")
    cells = {cell.identity: cell for cell in design.cells}
    for row in data.itertuples(index=False):
        cell = cells[row.cell_id]
        _require(row.n == cell.n and row.theta_mean_gap == cell.q and row.cv == cell.cv,
                 "summary grid coordinates disagree with the canonical design")
        values = [row.n_planned, row.n_point_valid, row.n_interval_valid, row.n_covered]
        _require(all(not isinstance(v, (bool, np.bool_)) and np.isfinite(v) and int(v) == v for v in values)
                 and row.n_planned == design.replications and 0 <= row.n_covered <= row.n_interval_valid
                 <= row.n_point_valid <= row.n_planned, "figure denominator/count mismatch")
    for _, group in data.groupby(["cell_id", "estimator"]):
        _require(group[["n_planned", "n_point_valid", "n_interval_valid"]].nunique().eq(1).all(),
                 "joint validity counts differ across parameter rows")
    return data


def _plot(ax, x, y, lower, upper, *, label, color, marker, gid, log_ratio=False):
    """Plot saved estimates and endpoints; NA is a visible gap, never zero."""
    x, y, lower, upper = np.asarray(x, float), _numbers(y), _numbers(lower), _numbers(upper)
    _require(x.shape == y.shape == lower.shape == upper.shape, "plot estimate/interval shapes differ")
    finite = np.isfinite(y)
    bounds = finite & np.isfinite(lower) & np.isfinite(upper)
    _require(np.all(lower[bounds] <= y[bounds] + 1e-12) and np.all(upper[bounds] >= y[bounds] - 1e-12),
             "saved Monte Carlo endpoints do not contain their estimate")
    zeros = finite & (y == 0) if log_ratio else np.zeros(len(y), dtype=bool)
    if log_ratio:
        _require(np.all(y[finite] >= 0), "negative RMSE ratio")
        bounds &= (lower > 0) & (upper > 0) & (y > 0)
    plotted = y.copy()
    plotted[zeros] = np.nan
    line, = ax.plot(x, plotted, color=color, marker=marker, markersize=4.6, linewidth=1.5, label=label)
    line.set_gid(gid)
    bars = ax.vlines(x[bounds], lower[bounds], upper[bounds], color=color, linewidth=.9)
    bars.set_gid(gid + ":ci")
    for xx in x[~finite]:
        ax.text(xx, .025, "NA", color=color, fontsize=7, ha="center", va="bottom", transform=ax.get_xaxis_transform())
    for xx in x[zeros]:
        ax.text(xx, .025, "0", color=color, fontsize=8, ha="center", va="bottom", transform=ax.get_xaxis_transform())
    return line


def _series(ax, data, cvs, metric, lower, upper, *, scale=1., log_ratio=False):
    for index, (method, label, color, marker) in enumerate(METHODS):
        selected = data[data.estimator == method].set_index("cv").loc[list(cvs)]
        x = np.arange(len(cvs)) + (index-1)*.055
        _plot(ax, x, scale*selected[metric].to_numpy(float), scale*selected[lower].to_numpy(float),
              scale*selected[upper].to_numpy(float), label=label, color=color, marker=marker,
              gid=method, log_ratio=log_ratio)
    ax.set_xticks(range(len(cvs)), [f"{c:g}" for c in cvs])
    ax.set_xlim(-.28, len(cvs)-.72)
    ax.set_xlabel(CV_AXIS_LABEL)
    if log_ratio:
        ax.set_yscale("log")
        # Near-unit ratios otherwise receive long scientific tick labels that
        # can push the axis title off-canvas. Keep the log transform, but choose
        # compact decimal precision from the visible positive range.
        def compact(value, _):
            lo, hi = ax.get_ylim()
            span = abs(hi-lo)/max(abs(lo), abs(hi), np.finfo(float).tiny)
            digits = min(12, max(4, int(math.ceil(-math.log10(max(span, 1e-12))))+3))
            return f"{value:.{digits}g}"
        ax.yaxis.set_major_formatter(FuncFormatter(compact))
        ax.yaxis.set_minor_formatter(FuncFormatter(compact))
        ax.tick_params(axis="y", which="both", labelsize=8)


def _finish(fig, axes, title, subtitle, note, *, bottom=.20, top=.84, legend=True):
    fig.suptitle(title, fontsize=15, fontweight="semibold", y=.98)
    fig.text(.5, .94, subtitle, ha="center", fontsize=10)
    if legend:
        handles, labels = np.asarray(axes, dtype=object).flat[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(.5, .915), ncol=3, frameon=False, fontsize=9)
    fig.subplots_adjust(left=.075, right=.985, bottom=bottom, top=top, wspace=.30, hspace=.47)
    fig.text(.075, .025, note, fontsize=8.2, color="#333333", va="bottom", linespacing=1.4)


def _subtitle(design, *, floor=False):
    label = "EXPLORATORY E3" if design.scope == "exploratory" else "VALIDATION ONLY"
    cell = design.cells[0]
    unit = "per physical floor variant" if floor else "per independent CV cell"
    return f"{label}  |  q={cell.q:g}, n={cell.n:,} observations  |  R={design.replications:,} {unit}"


def build_followup_figures(tables, design):
    """Two complete saved-table follow-up figures, with no duplicated counts."""
    require_design("followup", design)
    summary = _summary(tables, design)
    cvs = (0., .5)
    figures = {}
    with plt.rc_context(RENDER_SETTINGS):
        fig, axes = plt.subplots(2, 3, figsize=(14.2, 8.8), squeeze=False)
        for column, (parameter, symbol) in enumerate(PARAMETERS):
            for row, scope in enumerate(("operational", "conditional")):
                ax = axes[row, column]
                ax.set_gid(f"coverage:{scope}:{parameter}")
                group = summary[summary.parameter == parameter]
                _series(ax, group, cvs, f"coverage_{scope}", f"{scope}_wilson_lower", f"{scope}_wilson_upper", scale=100)
                ax.axhline(95, color="#555555", linestyle=":", linewidth=1)
                ax.set_ylim(-3, 103)
                ax.set_ylabel("Coverage (%)")
                ax.set_title(f"{symbol} | " + ("all R planned paths" if scope == "operational" else "valid intervals only"))
        _finish(fig, axes, "Week 7  Independent follow-up coverage", _subtitle(design),
                "Bars: pointwise 95% Wilson Monte Carlo intervals. Operational coverage counts unavailable intervals as noncoverage.\n"
                "Conditional coverage uses each estimator's interval-valid denominator; undefined values are gaps labelled NA. Dotted: nominal 95%.\n"
                "PFML/Euler Wald covariance is model-based, not robust. These exploratory paths are separate from the unchanged E1/E2 headlines.")
        figures[FOLLOWUP_STEMS[0]] = fig
        fig, axes = plt.subplots(2, 2, figsize=(12.6, 9.2), squeeze=False)
        joint = summary[summary.parameter == "theta"]
        for ax, kind in zip(axes[0], ("point", "interval")):
            ax.set_gid("availability:" + kind)
            _series(ax, joint, cvs, kind + "_availability", kind + "_availability_wilson_lower",
                    kind + "_availability_wilson_upper", scale=100)
            ax.set_ylim(-3, 103)
            ax.set_ylabel("Available / R (%)")
            ax.set_title("Valid point fits" if kind == "point" else "Available Wald intervals")
        for ax, category in zip(axes[1], ("point_invalid", "interval_unavailable_given_point")):
            ax.set_gid("failures:" + category)
            maximum = 0
            for index, (method, label, color, _) in enumerate(METHODS):
                group = joint[joint.estimator == method].set_index("cv").loc[list(cvs)]
                counts = (group.n_planned-group.n_point_valid if category == "point_invalid"
                          else group.n_point_valid-group.n_interval_valid).to_numpy(int)
                maximum = max(maximum, int(max(counts)))
                bars = ax.bar(np.arange(2)+(index-1)*.23, counts, width=.21, color=color, label=label)
                for bar, count in zip(bars, counts):
                    bar.set_gid(method)
                    ax.annotate(f"{count:,}", (bar.get_x()+bar.get_width()/2, count),
                                xytext=(0, 3), textcoords="offset points", ha="center", fontsize=8)
            ax.set_ylim(0, max(1, maximum)*1.22)
            ax.set_xticks(range(2), ["0", "0.5"])
            ax.set_xlabel(CV_AXIS_LABEL)
            ax.set_ylabel("Number of planned paths")
            ax.yaxis.set_major_locator(MaxNLocator(integer=True))
            ax.set_title("Point-invalid paths" if category == "point_invalid" else "Valid point, unavailable interval")
        _finish(fig, axes, "Week 7  Availability and retained failures", _subtitle(design),
                "Top: joint point/interval availability with pointwise 95% Wilson intervals. Bottom counts are disjoint categories, each counted once per method.\n"
                "Point-invalid + valid-point/interval-unavailable = all unavailable intervals. Counts are not tripled across the three parameter rows.\n"
                "Zero-height bars are explicitly labelled 0. Statistical failures and their finite candidates remain in the saved records.", bottom=.19)
        figures[FOLLOWUP_STEMS[1]] = fig
    return figures


def _epsilon_inputs(epsilon, design):
    _require(tuple(epsilon.columns) == floor_audit.EPSILON_COLUMNS, "all 71 ordered epsilon fields are required")
    expected = {(c.identity, method, parameter) for c in design.cells if c.floor_fraction != 1e-6
                for method, *_ in METHODS for parameter, _ in PARAMETERS}
    _keys(epsilon, ("cell_id", "estimator", "parameter"), expected, "paired epsilon table")
    for cell in design.cells:
        if cell.floor_fraction == 1e-6:
            continue
        group = epsilon[epsilon.cell_id == cell.identity]
        reference = replace(cell, floor_fraction=1e-6)
        _require(group.cv.eq(cell.cv).all() and group.floor_fraction.eq(cell.floor_fraction).all()
                 and group.reference_floor_fraction.eq(1e-6).all() and group.n_planned.eq(design.replications).all()
                 and group.reference_cell_id.eq(reference.identity).all() and group.seed_anchor_cell_id.eq(reference.identity).all()
                 and group.source_seed.eq(design.seed).all() and group.source_coupling_policy.eq(floor_audit.COUPLING).all()
                 and group.difference_direction.eq("epsilon_minus_registered_floor").all()
                 and group.ratio_direction.eq("epsilon_over_registered_floor").all(), "paired plot direction or design changed")
        counts = [column for column in group.columns if column.startswith("n_")]
        for row in group.to_dict("records"):
            _require(all(not isinstance(row[key], (bool, np.bool_)) and np.isfinite(row[key])
                         and int(row[key]) == row[key] and 0 <= row[key] <= design.replications for key in counts),
                     "paired figure counts must be integral, non-Boolean and within R")
            _require(row["n_interval_valid_a"] <= row["n_point_valid_a"]
                     and row["n_interval_valid_b"] <= row["n_point_valid_b"]
                     and row["n_common_point_valid"] <= min(row["n_point_valid_a"], row["n_point_valid_b"])
                     and row["n_common_point_valid"] >= max(0, row["n_point_valid_a"]+row["n_point_valid_b"]-design.replications)
                     and row["n_common_interval_valid"] <= min(row["n_interval_valid_a"], row["n_interval_valid_b"], row["n_common_point_valid"])
                     and row["n_common_interval_valid"] >= max(0, row["n_interval_valid_a"]+row["n_interval_valid_b"]-design.replications)
                     and row["n_covered_a"] <= row["n_interval_valid_a"] and row["n_covered_b"] <= row["n_interval_valid_b"]
                     and row["n_both_cover"] <= row["n_common_interval_valid"]
                     and row["n_both_cover"] + row["n_a_only_cover"] == row["n_covered_a"]
                     and row["n_both_cover"] + row["n_b_only_cover"] == row["n_covered_b"]
                     and sum(row[key] for key in ("n_both_cover", "n_a_only_cover", "n_b_only_cover", "n_neither_cover")) == design.replications,
                     "paired own/common/coverage count relationships changed")
        for _, method in group.groupby("estimator"):
            joint = ["n_planned", "n_point_valid_a", "n_point_valid_b", "n_interval_valid_a", "n_interval_valid_b",
                     "n_common_point_valid", "n_common_interval_valid"]
            _require(method[joint].nunique().eq(1).all(), "paired counts disagree across parameters")
    return epsilon


def _common_note(epsilon):
    lines = ["Common-valid S, in Exact/PFML/Euler order (all parameter rows share these counts):"]
    for floor in (1e-8, 1e-4):
        parts = []
        for cv in (0., 1.5, 2., 3.):
            group = epsilon[(epsilon.floor_fraction == floor) & (epsilon.cv == cv) & (epsilon.parameter == "theta")].set_index("estimator")
            parts.append(f"CV={cv:g}: " + "/".join(f"{int(group.loc[method, 'n_common_point_valid']):,}" for method, *_ in METHODS))
        lines.append(f"epsilon={floor:g}: " + " | ".join(parts))
    return "\n".join(lines)


def build_floor_figures(tables, epsilon, design):
    """Four full-grid sensitivity figures; CV=0 remains one shared control."""
    require_design("floor", design)
    _summary(tables, design)
    epsilon = _epsilon_inputs(epsilon, design)
    cvs = (0., 1.5, 2., 3.)
    figures = {}
    with plt.rc_context(RENDER_SETTINGS):
        for kind, stem, title in (("bias", FLOOR_STEMS[0], "Week 7  Paired floor effects on estimation bias"),
                                   ("coverage", FLOOR_STEMS[1], "Week 7  Paired floor effects on operational coverage"),
                                   ("ratio", FLOOR_STEMS[2], "Week 7  Paired RMSE ratios to the registered floor")):
            fig, axes = plt.subplots(2, 3, figsize=(14.6, 10.0), squeeze=False)
            for row, floor in enumerate((1e-8, 1e-4)):
                for column, (parameter, symbol) in enumerate(PARAMETERS):
                    ax = axes[row, column]
                    ax.set_gid(f"{kind}:{floor:g}:{parameter}")
                    group = epsilon[(epsilon.floor_fraction == floor) & (epsilon.parameter == parameter)]
                    if kind == "bias":
                        metric = "bias_difference" if parameter == "mu" else "relative_bias_difference"
                        scale = 1. if parameter == "mu" else 100.
                        ylabel = "Bias difference (state units)" if parameter == "mu" else "Relative-bias difference (pp)"
                    elif kind == "coverage":
                        metric, scale, ylabel = "operational_coverage_difference", 100., "Coverage difference (percentage points)"
                    else:
                        metric, scale, ylabel = "rmse_ratio", 1., "RMSE(epsilon) / RMSE(1e-6)"
                    _series(ax, group, cvs, metric, metric + "_mc_lower", metric + "_mc_upper", scale=scale, log_ratio=kind == "ratio")
                    ax.axhline(1 if kind == "ratio" else 0, color="#555555", linestyle=":", linewidth=1)
                    ax.set_ylabel(ylabel)
                    ax.set_title(f"{symbol} | epsilon={floor:g} vs 1e-6")
            if kind == "coverage":
                note = (f"All R={design.replications:,} planned paired paths enter each operational-coverage difference; unavailable intervals count as noncoverage.\n"
                        "Bars: approximate pointwise 95% paired normal Monte Carlo intervals, without clipping to +/-100 percentage points.\n"
                        "Positive means higher coverage under that epsilon than under epsilon=1e-6. CV=0 is an exact shared control, not an independent comparison.\n"
                        "Collapsed normal-MC intervals reflect zero observed variation, not proof of population equality.\n"
                        "PFML/Euler intervals use model-based working covariance. Undefined values are NA gaps; absent bars mean unavailable MC uncertainty.")
            else:
                note = ("Bars: approximate pointwise 95% paired Monte Carlo intervals; both metrics use the same common-valid point subset.\n"
                        + ("Theta/sigma bias differences are percentage points of truth-relative bias; mu uses absolute units. No practical-breach threshold applies to a contrast.\n"
                           if kind == "bias" else "Ratio bars use the saved asymmetric log-delta limits. The vertical scale is log; a 0 label at the baseline marks a zero ratio off that scale.\n")
                        + "CV=0 is one shared control. NA labels mark undefined estimates; absent bars mark undefined uncertainty.\n" + _common_note(epsilon))
            _finish(fig, axes, title, _subtitle(design, floor=True), note, bottom=.25 if kind != "coverage" else .20)
            figures[stem] = fig
        spacing = _spacing_inputs(tables, design)
        # The preserved spacing table does not contain floor_fraction. Join it
        # by canonical cell identity rather than inferring it from row order.
        fig, axes = plt.subplots(2, 2, figsize=(12.8, 8.8), squeeze=False)
        specifications = (("realized_mean_gap", "Mean realized gap", "Time units", False),
                          ("clip_fraction", "Mean fraction of clipped gaps", "Clipped gaps (%)", False),
                          ("min_gap", "Mean of each path's minimum gap", "Time units (symlog around zero)", True),
                          ("realized_span", "Mean realized observation span", "Time units", False))
        for ax, (field, title, ylabel, logarithmic) in zip(axes.flat, specifications):
            _require(spacing[field + "_n_available"].eq(design.replications).all(),
                     "each plotted spacing mean must use all R saved sources")
            ax.set_gid("spacing:" + field)
            scale = 100. if field == "clip_fraction" else 1.
            for index, (floor, label, color, marker) in enumerate(FLOORS):
                group = spacing[spacing.floor_fraction == floor].set_index("cv").loc[list(cvs)]
                y = scale*group[field + "_mean"].to_numpy(float)
                se = scale*group[field + "_mcse"].to_numpy(float)
                _plot(ax, np.arange(4)+(index-1)*.055, y, y-Z*se, y+Z*se,
                      label=label, color=color, marker=marker, gid=f"epsilon:{floor:g}")
            if logarithmic:
                # A normal MC interval for a positive mean can cross zero in
                # small validation samples. A signed scale retains that saved
                # interval instead of silently clipping its lower endpoint.
                ax.set_yscale("symlog", linthresh=1e-8)
            ax.set_xticks(range(4), [f"c={c:g}" for c in cvs])
            ax.set_xlim(-.28, 3.28)
            ax.set_title(title)
            ax.set_xlabel(CV_AXIS_LABEL)
            ax.set_ylabel(ylabel)
        _finish(fig, axes, "Week 7  Realized spacing under the three floor settings", _subtitle(design, floor=True),
                "Each summary uses R source paths, not three estimator copies. Bars: approximate pointwise 95% Monte Carlo intervals for each mean.\n"
                "Minimum-gap panel shows the mean of path minima, not a global minimum; its symlog axis retains normal-MC bounds crossing zero (linear within +/-1e-8).\n"
                "n counts observations; K=n-1=999 transitions. No renormalization or realized-span matching; realized CV remains in the complete saved table.\n"
                "All three epsilon variants share Gamma/initial/innovation draws within each CV and replication; the curves are paired. CV=0 is exactly identical.", bottom=.22)
        figures[FLOOR_STEMS[3]] = fig
    return figures


def _selected_case(payload):
    item = payload["selection"]
    _require(item.get("experiment") == "E2" and type(item.get("rep")) is int,
             "selected figures require their fixed E2 source identities")
    for cv, rep, stem in SELECTED_CASES:
        cell = Cell(n=250, q=2., cv=cv)
        if item["rep"] == rep and item.get("cell_id") == cell.identity:
            return cv, rep, stem
    raise ValueError("selected figure is outside the two fixed protocol cases")


def _ld(value, *, missing=False):
    if value is None and missing:
        return np.longdouble(np.nan)
    _require(type(value) is str, "saved extended scalar must retain its decimal string")
    result = np.longdouble(value)
    _require(np.isfinite(result), "nonfinite saved extended decimal")
    return result


def profile_series(payload):
    """Saved profile deltas, subtracting extended decimals before float plotting."""
    _selected_case(payload)
    diagnostic = payload.get("diagnostic", {})
    profile = diagnostic.get("exact_profile", {})
    if payload.get("status") != "evaluated" or profile.get("status") != "evaluated":
        return {"available": False, "reason": profile.get("reason", diagnostic.get("error", payload.get("status", "unavailable")))}
    best = profile["best_finite"]
    baseline = _ld(best["nll_extended"])
    x = _numbers(profile["grid_log_theta"])
    stored = [_ld(value, missing=True) for value in profile["grid_nll_extended"]]
    _require(len(stored) == len(x) and len(x) >= 3 and np.isfinite(x).all() and (np.diff(x) > 0).all(),
             "saved profile grid is incomplete or unordered")
    # Convert only completed differences. Binary64 absolute NLL values can
    # coincide even when the saved longdouble differences are informative.
    delta = np.asarray(np.asarray(stored, dtype=np.longdouble)-baseline, dtype=float)
    differences = profile["differences"]
    values = {name: float(_ld(differences[name + "_extended"])) for name in
              ("stored_profile_minus_best", "stored_full_minus_best", "brownian_limit_minus_best", "iid_limit_minus_best")}
    best_x = float(best["log_theta"])
    stored_x = float(profile["stored_profile"]["log_theta"])
    resolution = float(_ld(profile["resolution_extended"]))
    _require(np.isfinite([best_x, stored_x, resolution]).all() and resolution > 0,
             "saved profile coordinates or resolution are unavailable")
    lower = max(float(x[0]), min(best_x, stored_x)-.15)
    upper = min(float(x[-1]), max(best_x, stored_x)+.15)
    _require(lower < upper, "profile zoom interval is empty")
    return {"available": True, "x": x, "delta_nll": delta, "best_x": best_x,
            "stored_x": stored_x, "resolution": resolution, "differences": values,
            "zoom_limits": (lower, upper), "zoom_mask": (x >= lower) & (x <= upper),
            "bounds": profile["bounds"], "disposition": profile["disposition"]}


def _profile_panel(ax, values, *, zoom=False):
    mask = values["zoom_mask"] if zoom else np.ones(len(values["x"]), dtype=bool)
    line, = ax.plot(values["x"][mask], values["delta_nll"][mask], color="#0072B2", linewidth=1.4, label="Saved LD profile grid")
    line.set_gid("profile-grid")
    ax.scatter([values["best_x"]], [0.], marker="D", color="#0072B2", s=32, label="Best finite candidate", zorder=5).set_gid("profile-best")
    differences = values["differences"]
    for name, label, color, marker in (("stored_profile_minus_best", "Stored theta; nuisances profiled", "#E69F00", "o"),
                                      ("stored_full_minus_best", "Stored parameters; LD profile representation", "#D55E00", "x")):
        ax.scatter([values["stored_x"]], [differences[name]], marker=marker, color=color, s=38, zorder=6, label=label).set_gid(name)
    if not zoom:
        for name, label, color in (("brownian_limit_minus_best", "Brownian + free-drift limit", "#009E73"),
                                   ("iid_limit_minus_best", "iid limiting likelihood", "#CC79A7")):
            ax.axhline(differences[name], color=color, linestyle="--", linewidth=1, label=label).set_gid(name)
        ax.set_xlim(values["x"][0], values["x"][-1])
        ax.set_yscale("symlog", linthresh=values["resolution"])
    else:
        ax.set_xlim(*values["zoom_limits"])
    ax.axhline(values["resolution"], color="#777777", linestyle=":", linewidth=1, label="Saved diagnostic resolution").set_gid("profile-resolution")
    ax.axhline(0, color="#777777", linewidth=.6)
    ax.set_xlabel(r"Working $\log(\theta_u)$")
    ax.set_ylabel("Conditional NLL minus best finite NLL")
    ax.set_title("Stored / best neighbourhood; saved grid only" if zoom else "Complete fixed profile grid (signed log y)")
    ax.legend(fontsize=6.5, frameon=False, loc="best")


def _ladder_traces(diagnostic):
    result = []
    for method, label, color, marker in METHODS:
        fit = diagnostic["fits"][method]
        if fit.get("status") == "evaluated":
            result.append((method, label + " float64", color, marker, fit["float64_ladder"]))
    exact = diagnostic["fits"]["exact"]
    if exact.get("status") == "evaluated" and "independent_longdouble_ladder" in exact:
        result.append(("exact_longdouble", "Exact LD; same normalized data", "#CC79A7", "D", exact["independent_longdouble_ladder"]))
    for _, _, _, _, ladder in result:
        _require(len(ladder["entries"]) == 4 and len(ladder["adjacent"]) == 3,
                 "selected figure requires all four saved Hessian steps and three comparisons")
        _require([entry["relative_step"] for entry in ladder["entries"]] == [.001, .0005, .00025, .000125],
                 "saved Hessian step sequence changed")
        _require([(pair["coarse_relative_step"], pair["fine_relative_step"]) for pair in ladder["adjacent"]]
                 == [(.001, .0005), (.0005, .00025), (.00025, .000125)], "saved adjacent Hessian pair identities changed")
        _require(ladder.get("replacement_interval_constructed") is False, "diagnostic ladder must not replace an interval")
    return result


def _diagnostic_line(ax, x, values, *, gid, label, color, marker, linestyle="-"):
    values = _numbers(values)
    line, = ax.plot(x, values, color=color, marker=marker, markersize=4, linewidth=1.3, linestyle=linestyle, label=label)
    line.set_gid(gid)
    for xx in np.asarray(x)[~np.isfinite(values)]:
        ax.text(xx, .03, "NA", transform=ax.get_xaxis_transform(), fontsize=6.5, color=color, ha="center")
    return line


def build_selected_figure(payload, *, validation):
    """One fixed case from saved diagnostics; no profile or Hessian evaluation."""
    _require(type(validation) is bool and payload.get("headline_results_replaced") is False,
             "selected diagnostic scope/replacement flag is invalid")
    cv, rep, _ = _selected_case(payload)
    literal = payload.get("literal_validation_fixture") is True
    _require(not literal or validation, "literal inputs cannot be labelled as exploratory source replay")
    label = ("LITERAL VALIDATION FIXTURE — not an original E2 replay" if literal else
             "VALIDATION ONLY — no original E2 scientific conclusion" if validation else
             "EXPLORATORY E3 — saved diagnostics of the original E2 source")
    with plt.rc_context(RENDER_SETTINGS):
        fig, axes = plt.subplots(2, 3, figsize=(15.4, 11.8), squeeze=False)
        for ax, name in zip(axes.flat, ("profile:global", "profile:zoom", "hessian:eigenvalues",
                                        "hessian:condition", "hessian:score", "hessian:adjacent")):
            ax.set_gid(name)
        values = profile_series(payload)
        diagnostic = payload.get("diagnostic", {})
        if payload.get("status") != "evaluated":
            for ax in axes.flat:
                ax.set_axis_off()
                ax.text(.5, .5, "Saved diagnostics unavailable\n" + str(values["reason"]), ha="center", va="center", wrap=True, transform=ax.transAxes)
            note = "The saved failure remains unresolved. No source, fit, profile, Hessian or replacement interval was generated for this figure."
        else:
            _require(diagnostic.get("headline_results_replaced") is False, "selected diagnostics replaced headline results")
            if values["available"]:
                _profile_panel(axes[0, 0], values)
                _profile_panel(axes[0, 1], values, zoom=True)
            else:
                for ax in axes[0, :2]:
                    ax.text(.5, .5, "Profile unavailable\n" + str(values["reason"]), ha="center", va="center", transform=ax.transAxes)
            traces = _ladder_traces(diagnostic)
            eigen_ax = axes[0, 2]
            eigen_colors = ("#0072B2", "#D55E00", "#009E73")
            for method, label_text, _, marker, ladder in traces:
                if method not in {"exact", "exact_longdouble"}:
                    continue
                for index, color in enumerate(eigen_colors):
                    series = [None if row.get("eigenvalues") is None else row["eigenvalues"][index] for row in ladder["entries"]]
                    _diagnostic_line(eigen_ax, range(4), series, gid=f"{method}:lambda{index+1}",
                                     label=f"lambda{index+1} " + ("LD" if method.endswith("longdouble") else "float64"),
                                     color=color, marker=marker, linestyle="--" if method.endswith("longdouble") else "-")
            eigen_ax.set_yscale("symlog", linthresh=1e-6)
            eigen_ax.axhline(0, color="#777777", linewidth=.7)
            eigen_ax.set_title("Exact full-Hessian eigenvalues")
            eigen_ax.set_ylabel("Signed eigenvalue (symlog)")
            eigen_ax.legend(fontsize=7, ncol=2, frameon=False)
            for ax, field, title, threshold in ((axes[1, 0], "condition_number", "Full-Hessian condition number", 1e10),
                                                 (axes[1, 1], "standardized_score", "Standardized score at stored point", .01)):
                for method, label_text, color, marker, ladder in traces:
                    _diagnostic_line(ax, range(4), [row.get(field) for row in ladder["entries"]],
                                     gid=method, label=label_text, color=color, marker=marker)
                ax.set_yscale("symlog", linthresh=1. if field == "condition_number" else 1e-6)
                ax.axhline(threshold, color="#777777", linestyle=":", linewidth=1, label="Registered gate")
                ax.set_title(title)
                ax.legend(fontsize=6.6, frameon=False, loc="best")
            adjacent = axes[1, 2]
            for method, label_text, color, marker, ladder in traces:
                for field, style in (("whitened_spectral", "-"), ("relative_frobenius", "--")):
                    _diagnostic_line(adjacent, range(3), [row.get(field) for row in ladder["adjacent"]],
                                     gid=method + ":" + field, label=label_text if style == "-" else "_nolegend_",
                                     color=color, marker=marker, linestyle=style)
            adjacent.set_yscale("symlog", linthresh=1e-6)
            adjacent.axhline(.01, color="#777777", linestyle=":", linewidth=1)
            adjacent.axhline(.005, color="#999999", linestyle=":", linewidth=.7)
            adjacent.set_title("Adjacent-step discrepancy: solid W; dashed F")
            adjacent.set_ylabel("Whitened spectral / relative Frobenius")
            adjacent.set_xticks(range(3), ["h → h/2", "h/2 → h/4", "h/4 → h/8"])
            adjacent.legend(fontsize=6.2, frameon=False, loc="best")
            for ax in (eigen_ax, axes[1, 0], axes[1, 1]):
                ax.set_xticks(range(4), [".001", ".0005", ".00025", ".000125"], rotation=20)
                ax.set_xlabel("Relative finite-difference step")
            precision = diagnostic["precision"]
            representations = diagnostic["representation_differences"]
            precision_note = (f"Saved longdouble: {precision['significand_bits']} significant binary digits "
                              f"({precision['fraction_bits']} fraction bits), {precision['storage_bytes']} storage bytes; not arbitrary precision.")
            representation_note = ("Registered Hessian ladders use identical float64-normalized inputs; the LD profile normalizes original values/gaps separately.\n"
                                   f"Saved representation differences: max |delta value|={representations['maximum_absolute_value_difference']:.3g}; "
                                   f"max |delta gap|={representations['maximum_absolute_gap_difference']:.3g}.")
            if values["available"]:
                bounds = values["bounds"]
                disposition = textwrap.fill(values["disposition"].replace("_", " "), width=150)
                profile_note = (f"Saved resolution={values['resolution']:.3g}; lower/upper tail targets reached: "
                                f"{bounds['lower_tail_reached']}/{bounds['upper_tail_reached']}. {disposition}\n"
                                "Finite searches do not prove a global optimum. Negative saved objective differences are retained, not clipped.")
            else:
                profile_note = "No available exact profile; its original point-validity or diagnostic limitation remains explicit."
            status = diagnostic["fits"]["exact"]
            note = (f"Original exact point/interval status: {status['original_point_status']} / {status['original_interval_status']}. "
                    "No replacement fit, interval or covariance.\n" + precision_note + "\n" + representation_note + "\n" + profile_note + "\n"
                    "h/2 is the registered matrix candidate subject to all original gates; h/4 and h/8 are diagnostic only. W gate=.01; F gate=.005; score gate=.01.")
        _finish(fig, axes, f"Week 7  Fixed E2 case: CV={cv:g}, replication {rep}", label + "  |  q=2, n=250", note,
                bottom=.255, top=.875, legend=False)
    return fig


def figure_bytes(figure, stem):
    """Narrow deterministic save boundary for failure injection and readback."""
    result = {}
    with plt.rc_context(RENDER_SETTINGS):
        for extension in ("png", "pdf"):
            stream = BytesIO()
            metadata = ({"Title": stem, "Creator": "OU irregular-sampling research",
                         "CreationDate": None, "ModDate": None} if extension == "pdf"
                        else {"Software": "OU irregular-sampling research"})
            figure.savefig(stream, format=extension, dpi=170, metadata=metadata)
            result[f"{stem}.{extension}"] = stream.getvalue()
    return result


def _digest(value):
    return digest_bytes(canonical(value).encode())


def _json_file(path):
    value, raw = floor_audit._json_file(Path(path))
    _require(type(value) is dict and raw == (canonical(value) + "\n").encode(),
             "renderer input must be a canonical JSON object")
    return value, raw


def _ordinary_path(path, *, directory=True):
    supplied = Path(path).absolute()
    _require(not any(item.is_symlink() for item in (supplied, *supplied.parents)),
             "symbolic links are forbidden in renderer paths")
    if supplied.exists():
        _require(supplied.is_dir() if directory else supplied.is_file(), "renderer artifact type changed")
    return supplied.resolve()


def _frozen_provenance(provenance, design, *, required, exact=None, frozen):
    """Check old raw provenance as recorded; new rendering has its own record."""
    files = provenance.get("source_files")
    _require(type(files) is dict and set(required) <= set(files)
             and (exact is None or set(files) == set(exact)), "renderer dependency membership changed")
    _require(all(type(name) is str and not Path(name).is_absolute() and ".." not in Path(name).parts
                 and type(value) is str and re.fullmatch(r"[0-9a-f]{64}", value) for name, value in files.items()),
             "unsafe dependency name or hash")
    verify(provenance, expected_design=design)
    dirty = provenance.get("uncommitted_dependency_files")
    _require(type(dirty) is list and len(dirty) == len(set(dirty))
             and all(type(name) is str and name in files for name in dirty), "invalid recorded dirty dependency list")
    commit = provenance.get("code_commit")
    _require(type(commit) is str and re.fullmatch(r"[0-9a-f]{40}", commit), "invalid recorded Git commit")
    if frozen:
        _require(dirty == [], "full exploratory rendering requires committed frozen dependencies")
        for name, expected in files.items():
            recorded = subprocess.run(["git", "show", f"{commit}:{name}"], cwd=ROOT, capture_output=True)
            _require(recorded.returncode == 0 and digest_bytes(recorded.stdout) == expected,
                     "recorded frozen commit does not contain dependency: " + name)


def _file_map(directory):
    return {str(path.relative_to(directory)): file_digest(path) for path in sorted(directory.rglob("*")) if path.is_file()}


def _saved_endpoint_descriptors(output, design, binding):
    """Read only fixed saved endpoint rows; do not regenerate their arrays."""
    descriptors = []
    endpoints = sorted({0, design.replications-1})
    for cell in design.cells:
        selected = {}
        for part in binding["cells"][cell.identity]["partitions"]:
            wanted = {rep for rep in endpoints if part["start"] <= rep < part["stop"]}
            if not wanted:
                continue
            with gzip.open(output / "cells" / cell.identity / part["filename"], "rt", encoding="utf-8") as stream:
                for line in stream:
                    row = floor_audit.strict_json(line, "saved endpoint partition")
                    _require(type(row) is dict and type(row.get("rep")) is int, "saved endpoint replication type changed")
                    if row["rep"] in wanted:
                        key = (row["rep"], row.get("estimator"))
                        _require(row.get("cell_id") == cell.identity and key not in selected
                                 and key[1] in floor_audit.METHODS, "saved endpoint key changed")
                        selected[key] = row
        _require(set(selected) == {(rep, method) for rep in endpoints for method in floor_audit.METHODS},
                 "saved endpoint triplet is incomplete")
        for rep in endpoints:
            triplet = [selected[rep, method] for method in floor_audit.METHODS]
            hashes = {name: triplet[0].get(name) for name in ("times_sha256", "observations_sha256")}
            _require(all(type(value) is str and re.fullmatch(r"[0-9a-f]{64}", value) for value in hashes.values())
                     and all({name: row.get(name) for name in hashes} == hashes for row in triplet),
                     "saved endpoint source hashes are malformed or not shared")
            descriptors.append({"cell_id": cell.identity, "replication": rep, **hashes,
                                "estimators_checked": list(floor_audit.METHODS)})
    return descriptors


def _audit_tolerances(kind):
    common = {"arithmetic_rtol": 2e-10, "arithmetic_atol": 2e-12, "replay_rtol": 2e-10, "replay_atol": 2e-12,
              "covariance_rtol": 2e-9, "covariance_atol": "2e-12 * max(1, expected matrix Frobenius norm)"}
    if kind == "floor":
        return {**common, "source_arrays_and_diagnostics": "exact types, values and SHA-256 array bytes"}
    return {**common,
            "span_roundoff": "8*gamma_(K+8)*max(span,K*nominal_gap)+8 ulps; gamma_m=m*eps/(1-m*eps)",
            "minimum_gap_roundoff": "4 ulps at the realized-span scale"}


def _verify_run_audit(output, design, metadata, raw_binding):
    saved, _ = _json_file(output / "audit.json")
    artifacts = ordinary_audit.artifact_fingerprint(output)
    common = {
        "schema_version": 1, "status": "passed", "full_completion_audit": True, "replay_performed": True,
        "design_sha256": design.digest, "run_identity_sha256": metadata["provenance"]["identity_sha256"],
        "code_commit": metadata["provenance"]["code_commit"],
        "source_dependency_files_verified": len(metadata["provenance"]["source_files"]),
        "audited_artifact_sha256": artifacts,
        "audited_run_metadata_sha256": ordinary_audit.semantic_metadata_digest(metadata),
        "cells_verified": len(design.cells), "planned_paths": len(design.cells)*design.replications,
        "planned_fit_records": 3*len(design.cells)*design.replications, "pipeline_errors": 0,
        "replayed_source_paths": len(design.cells)*len({0, design.replications-1}),
        "replayed_fits": 3*len(design.cells)*len({0, design.replications-1}),
    }
    _require(all(key in saved and canonical(saved[key]) == canonical(value) for key, value in common.items()),
             "a full saved audit of these exact scientific inputs is required")
    _require(type(saved.get("replay_records")) is list and len(saved["replay_records"]) == common["replayed_source_paths"],
             "saved audit endpoint descriptors are incomplete")
    kind = metadata["execution_kind"]
    selection = ("replications 0 and R-1 in every physical floor variant, all estimators" if kind == "floor"
                 else "replications 0 and R-1 in every cell, all estimators")
    _require(saved.get("replay_selection") == selection
             and canonical(saved.get("tolerances")) == canonical(_audit_tolerances(kind)),
             "saved audit fixed replay selection or tolerance declaration changed")
    _require(canonical(saved["replay_records"]) == canonical(_saved_endpoint_descriptors(output, design, raw_binding)),
             "saved audit endpoint identities, estimator scope or array hashes changed")
    if metadata["execution_kind"] == "floor":
        required = {
            "audit_contract": "week7_independent_paired_floor_audit_v1",
            "scientific_overrides": ["paired_floor_source"], "source_policy": metadata["source_policy"],
            "raw_input_sha256": _digest(raw_binding),
            "planned_source_realizations": common["planned_paths"],
            "planned_independent_coupling_groups": 4*design.replications,
            "source_replay_performed": True, "source_paths_verified": common["planned_paths"],
            "source_array_hashes_verified": 2*common["planned_paths"],
        }
        _require(all(key in saved and canonical(saved[key]) == canonical(value) for key, value in required.items()),
                 "saved audit does not disclose the paired-floor source override and full source scope")
    else:
        _require("audit_contract" not in saved and "scientific_overrides" not in saved,
                 "follow-up figures require the unchanged registered-source audit contract")
    closed = metadata["status"] in {"audited_pending_visual_review", "completed"}
    if closed or "audit_sha256" in metadata:
        _require(metadata.get("audit_sha256") == file_digest(output / "audit.json"), "saved root audit seal changed")
    return saved, artifacts


def _verify_completion(output, metadata):
    if metadata["status"] != "completed":
        return
    files = _file_map(output)
    expected = {name: value for name, value in files.items() if name not in {"run_metadata.json", ".execution.lock"}}
    _require(canonical(metadata.get("artifact_sha256")) == canonical(expected), "completed artifact seal changed")
    review, _ = _json_file(output / "visual_review.json")
    figure, _ = _json_file(output / "figures/figure_manifest.json")
    pngs = {name for name in figure["figures"] if name.endswith(".png")}
    _require(metadata.get("visual_review_sha256") == file_digest(output / "visual_review.json")
             and review.get("all_panels_legible") is True and type(review.get("reviewed_pngs")) is list
             and len(review["reviewed_pngs"]) == len(set(review["reviewed_pngs"]))
             and set(review["reviewed_pngs"]) == pngs
             and review.get("figure_manifest_sha256") == file_digest(output / "figures/figure_manifest.json")
             and type(review.get("reviewer")) is str and bool(review["reviewer"].strip()),
             "completed visual review is not bound to the current figures")


@dataclass(frozen=True)
class SavedFigureInputs:
    kind: str
    source: Path
    destination: Path
    stems: tuple
    metadata: dict
    root_bytes: bytes
    binding: dict
    snapshot: dict
    frozen: bool
    tables: dict
    epsilon: object = None
    payloads: object = None


def _load_run_inputs(kind, output_dir, design):
    """Read a full saved audit and table/raw bindings, without numerical replay."""
    require_design(kind, design)
    output = _ordinary_path(output_dir)
    _require(output.is_dir(), "saved run directory is missing")
    execution._root_layout(output, design)
    metadata, root_bytes = _json_file(output / "run_metadata.json")
    expected = {
        "schema_version": 1, "execution_kind": kind, "execution_contract": execution.EXECUTION_CONTRACT,
        "project": "ou-week7-paired-floor" if kind == "floor" else "ou-irregular-experiments",
        "procedure": execution._procedure(kind), "source_policy": execution.source_policy(kind, design),
        "scientific_overrides": ["paired_floor_source"] if kind == "floor" else [],
        "require_frozen": design.scope == "exploratory", "pipeline_errors": 0,
        "completed_cells": {cell.identity: {"cell_id": cell.identity, "paths": design.replications,
                                           "rows": 3*design.replications, "pipeline_errors": 0} for cell in design.cells},
    }
    _require(all(key in metadata and canonical(metadata[key]) == canonical(value) for key, value in expected.items()),
             "saved run scientific identity or complete accounting changed")
    _require(metadata.get("status") in {"raw_complete_pending_audit", "implementation_review_required",
                                       "audited_pending_visual_review", "completed"}, "figures require complete audited raw evidence")
    required = set(CORE_FILES) | (set(floor_audit.FILES) if kind == "floor"
                                  else set(reporting.FILES) | set(ordinary_audit.AUDIT_FILES))
    _frozen_provenance(metadata["provenance"], design, required=required, frozen=expected["require_frozen"])
    raw = reporting.raw_input_binding(output, design, metadata)
    standard = reporting._verify_standard_manifest(output, design, required=True)
    epsilon_manifest = reporting._verify_epsilon_manifest(output, kind, design, raw, required=kind == "floor")
    audit, artifacts = _verify_run_audit(output, design, metadata, raw)
    tables = {name: floor_audit._table(output / name) for name in floor_audit.TABLE_KEYS}
    summary = _standard_inputs(tables, design)
    joint = summary[summary.parameter == "theta"]
    audit_counts = {
        "partitions_verified": sum(len(cell["partitions"]) for cell in raw["cells"].values()),
        "point_valid_records": int(joint.n_point_valid.sum()),
        "interval_valid_records": int(joint.n_interval_valid.sum()),
        "retained_failure_records": standard["failed_fits.csv.gz"]["rows"],
        "tables_verified": {**{name: standard[name]["rows"] for name in floor_audit.TABLE_KEYS},
                            **({"epsilon_comparisons.csv": 72} if kind == "floor" else {})},
    }
    _require(all(key in audit and canonical(audit[key]) == canonical(value) for key, value in audit_counts.items()),
             "saved audit table/count declarations differ from the bound evidence")
    epsilon = floor_audit._table(output / "epsilon_comparisons.csv") if kind == "floor" else None
    if kind == "floor":
        _epsilon_inputs(epsilon, design)
    _verify_completion(output, metadata)
    binding = {
        "kind": kind, "design": design.to_dict(), "design_sha256": design.digest,
        "run_identity_sha256": metadata["provenance"]["identity_sha256"],
        "semantic_metadata_sha256": ordinary_audit.semantic_metadata_digest(metadata),
        "audit_sha256": file_digest(output / "audit.json"), "saved_audit_status": audit["status"],
        "saved_audit_full_completion": True, "audited_artifact_sha256": artifacts,
        "raw_input_sha256": _digest(raw), "table_manifest_sha256": file_digest(output / "table_manifest.json"),
        "input_table_sha256": {name: standard[name]["sha256"] for name in reporting.STANDARD_FILES},
        "epsilon_inputs": ({"epsilon_comparisons.csv": epsilon_manifest["table"]["sha256"],
                            "epsilon_table_manifest.json": file_digest(output / "epsilon_table_manifest.json")}
                           if kind == "floor" else {}),
    }
    snapshot = {**artifacts, "audit.json": binding["audit_sha256"]}
    _require((output / "run_metadata.json").read_bytes() == root_bytes, "root changed during saved-input readback")
    return SavedFigureInputs(kind, output, output / "figures", FOLLOWUP_STEMS if kind == "followup" else FLOOR_STEMS,
                             metadata, root_bytes, binding, snapshot, expected["require_frozen"], tables, epsilon)


def _load_selected_inputs(selected_dir, output_dir):
    """Validate closed retained diagnostics; failures remain review-required."""
    source, output = _ordinary_path(selected_dir), _ordinary_path(output_dir)
    _require(source.is_dir() and source != output and not source.is_relative_to(output)
             and not output.is_relative_to(source), "selected figures must be outside and disjoint from the closed selected archive")
    metadata, root_bytes = _json_file(source / "run_metadata.json")
    _require(type(metadata.get("schema_version")) is int and metadata["schema_version"] == 1
             and metadata.get("project") == "ou-week7-path-diagnostics"
             and metadata.get("status") in {"completed_diagnostic_evidence", "diagnostic_review_required"}
             and type(metadata.get("require_frozen")) is bool, "selected evidence must be closed with an explicit freeze policy")
    selection, selection_bytes = _json_file(source / "selection.json")
    _require(type(selection.get("sources")) is list and bool(selection["sources"]), "closed selection is missing")
    keys = [path_diagnostics._source_key(item) for item in selection["sources"]]
    _require(len(keys) == len(set(keys)), "duplicate selected-source identity")
    descriptor, descriptor_bytes = _json_file(source / "diagnostic_design.json")
    expected_descriptor = path_diagnostics._diagnostic_design(selection, selection_bytes)
    _require(canonical(descriptor) == canonical(expected_descriptor.to_dict()), "selected diagnostic descriptor changed")
    frozen = metadata["require_frozen"]
    required = set(CORE_FILES) | set(path_diagnostics.FILES)
    if frozen:
        fixed = "artifacts/week7_evidence/selection.json"
        _require((ROOT / fixed).read_bytes() == selection_bytes, "full selected input differs from the fixed repository selection")
        required.add(fixed)
    _frozen_provenance(metadata["provenance"], expected_descriptor, required=required, frozen=frozen)
    path_diagnostics._check_output_layout(source, selection, require_complete=True)
    originals, payloads, all_hashes = {}, {}, {}
    for item in selection["sources"]:
        filename = path_diagnostics._filename(item)
        payload, _ = _json_file(source / "sources" / filename)
        floor_audit._nested_json(payload)
        originals[path_diagnostics._source_key(item)] = payload["original_records"]
        all_hashes[filename] = file_digest(source / "sources" / filename)
        if item["experiment"] == "E2" and any(item["cell_id"] == Cell(250, 2., cv).identity
                                              and item["rep"] == rep for cv, rep, _ in SELECTED_CASES):
            _, _, stem = _selected_case(payload)
            payloads[stem] = payload
    _require(set(payloads) == {stem for _, _, stem in SELECTED_CASES}, "one of the two fixed selected figure sources is missing")
    # This is the unchanged structural/hash/full-record readback helper, not
    # a source generator, profile calculation, Hessian calculation or fit.
    expected_audit = path_diagnostics._audit_selected_files(source, selection, originals, metadata)
    audit, _ = _json_file(source / "audit.json")
    _require(canonical(audit) == canonical(expected_audit)
             and metadata.get("audit_sha256") == file_digest(source / "audit.json"), "closed selected audit binding changed")
    _require(metadata["status"] == ("completed_diagnostic_evidence" if audit["status"] == "passed" else "diagnostic_review_required"),
             "selected failure status was not retained honestly")
    snapshot = {name: value for name, value in _file_map(source).items() if name != ".execution.lock"}
    binding = {
        "kind": "selected", "design": descriptor, "design_sha256": expected_descriptor.digest,
        "run_identity_sha256": metadata["provenance"]["identity_sha256"],
        "root_metadata_sha256": digest_bytes(root_bytes), "audit_sha256": file_digest(source / "audit.json"),
        "saved_audit_status": audit["status"], "saved_audit_failures": audit["failures"],
        "diagnostic_design_sha256": digest_bytes(descriptor_bytes), "selection_sha256": digest_bytes(selection_bytes),
        "source_payload_sha256": all_hashes,
        "selected_sources": {stem: {"selection": payload["selection"], "status": payload["status"],
                                    "payload_sha256": all_hashes[path_diagnostics._filename(payload["selection"])]}
                             for stem, payload in payloads.items()},
        "input_table_sha256": {}, "epsilon_inputs": {}, "input_artifact_sha256": snapshot,
    }
    _require((source / "run_metadata.json").read_bytes() == root_bytes, "selected root changed during readback")
    return SavedFigureInputs("selected", source, output, tuple(stem for _, _, stem in SELECTED_CASES),
                             metadata, root_bytes, binding, snapshot, frozen, {}, payloads=payloads)


def _settings_value(value):
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, (list, tuple)):
        return [_settings_value(item) for item in value]
    # The only nonprimitive default is the color Cycler; its values are plain
    # lists, so use its public value mapping rather than process-specific repr.
    if hasattr(value, "by_key"):
        return {name: _settings_value(item) for name, item in value.by_key().items()}
    raise ValueError("unrecorded effective Matplotlib setting type: " + type(value).__name__)


def renderer_environment():
    """Bind the complete normalized settings and font/compression implementation."""
    from matplotlib import ft2font
    font_root = Path(matplotlib.get_data_path()) / "fonts/ttf"
    return {
        "backend": "Agg", "matplotlib": matplotlib.__version__, "freetype": ft2font.__freetype_version__,
        "pillow": version("pillow"), "fonttools": version("fonttools"), "zlib": zlib.ZLIB_RUNTIME_VERSION,
        "effective_rc_settings": {key: _settings_value(value) for key, value in RENDER_SETTINGS.items()},
        "bundled_ttf_sha256": {path.name: file_digest(path) for path in sorted(font_root.glob("*.ttf"))},
        "png_dpi": 170, "pdf_creation_date": None, "pdf_modification_date": None,
    }


@dataclass(frozen=True)
class _RendererDesign:
    serialized: str

    def to_dict(self):
        return json.loads(self.serialized)

    @property
    def digest(self):
        return digest_bytes(self.serialized.encode())


def _descriptor(prepared):
    return {
        "schema_version": 1, "rendering_contract": CONTRACT, "kind": prepared.kind,
        "input_binding": prepared.binding, "require_frozen": prepared.frozen,
        "figure_filenames": sorted(f"{stem}.{extension}" for stem in prepared.stems for extension in ("png", "pdf")),
        "renderer_environment": renderer_environment(),
    }


def _dependency_set(extra_files):
    result = set(CORE_FILES) | set(FILES) | set(extra_files)
    _require(all(type(name) is str and not Path(name).is_absolute() and ".." not in Path(name).parts
                 for name in result), "unsafe renderer dependency path")
    return result


def _validate_renderer_state(state, descriptor, extra_files):
    _require(type(state) is dict and set(state) == {"schema_version", "rendering_contract", "descriptor",
                                                  "descriptor_sha256", "provenance"}
             and type(state.get("schema_version")) is int and state["schema_version"] == 1
             and state.get("rendering_contract") == CONTRACT
             and canonical(state.get("descriptor")) == canonical(descriptor)
             and state.get("descriptor_sha256") == _digest(descriptor), "renderer state/input identity changed")
    design = _RendererDesign(canonical(descriptor))
    dependencies = _dependency_set(extra_files)
    _frozen_provenance(state["provenance"], design, required=dependencies, exact=dependencies,
                       frozen=descriptor["require_frozen"])
    return state


def renderer_state_path(output_dir):
    """Immutable sibling to the run root or the external selected output."""
    output = _ordinary_path(output_dir)
    return output.with_name(output.name + "-figure-renderer") / "renderer_state.json"


def _state_layout(path):
    _ordinary_path(path.parent)
    _reject_archive_ancestors(path.parent)
    if path.parent.exists():
        _require(path.parent.is_dir() and {child.name for child in path.parent.iterdir()} == {"renderer_state.json"},
                 "renderer state directory is missing its descriptor or contains unexpected evidence")
        _require(path.is_file() and not path.is_symlink(), "renderer state must be a regular file")


def _figure_layout(prepared, *, complete=False):
    output = prepared.destination
    if not output.exists():
        _require(not complete, "declared figure directory is missing")
        return None
    _ordinary_path(output)
    expected = {f"{stem}.{extension}" for stem in prepared.stems for extension in ("png", "pdf")}
    allowed = expected | {"figure_manifest.json"} | ({".execution.lock"} if prepared.kind == "selected" else set())
    _require(all(child.name in allowed and child.is_file() and not child.is_symlink() for child in output.iterdir()),
             "figure directory contains a noncanonical name, nonregular file or symbolic link")
    path = output / "figure_manifest.json"
    if not path.exists():
        _require(not complete, "figure manifest is missing")
        return None
    manifest, _ = _json_file(path)
    _require(type(manifest.get("figures")) is dict and set(manifest["figures"]) == expected,
             "figure manifest does not declare the complete canonical PNG/PDF set")
    present = {child.name for child in output.iterdir()} - {".execution.lock"}
    _require(present == expected | {"figure_manifest.json"}, "declared figure is missing or an extra figure is present")
    for name, digest in manifest["figures"].items():
        _require(type(digest) is str and re.fullmatch(r"[0-9a-f]{64}", digest)
                 and file_digest(output / name) == digest, "declared figure bytes changed")
    return manifest


def _manifest(prepared, state, figures):
    return {
        "schema_version": 1, "rendering_contract": CONTRACT, "kind": prepared.kind,
        "design_sha256": prepared.binding["design_sha256"],
        "run_identity_sha256": prepared.binding["run_identity_sha256"],
        "input_binding": prepared.binding,
        # Compatibility with the historical figure sealer: this map contains
        # only the seven standard table-manifest entries, never epsilon files.
        "input_table_sha256": prepared.binding["input_table_sha256"],
        "epsilon_inputs": prepared.binding["epsilon_inputs"],
        "renderer": state, "renderer_state_sha256": _digest(state),
        "figures": figures, "new_sources_or_fits": False, "headline_results_replaced": False,
        "visual_review_performed": False,
    }


def _verify_manifest(prepared, descriptor, extra_files):
    manifest = _figure_layout(prepared, complete=True)
    state = manifest.get("renderer")
    _validate_renderer_state(state, descriptor, extra_files)
    expected = _manifest(prepared, state, manifest["figures"])
    _require(canonical(manifest) == canonical(expected), "figure manifest scope, input binding or renderer identity changed")
    return manifest


def _current_snapshot(prepared):
    if prepared.kind == "selected":
        path_diagnostics._check_output_layout(prepared.source, prepared.binding["design"]["selection"], require_complete=True)
        return {name: value for name, value in _file_map(prepared.source).items() if name != ".execution.lock"}
    execution._root_layout(prepared.source, require_design(prepared.kind, _saved_design(prepared.binding["design"])))
    return {**ordinary_audit.artifact_fingerprint(prepared.source),
            "audit.json": file_digest(prepared.source / "audit.json")}


def _saved_design(value):
    from ..experiments.design import RunDesign
    cells = tuple(Cell(n=row["n"], q=row["theta_mean_gap"], cv=row["cv"], theta=row["true_theta"],
                       mu=row["true_mu"], sigma=row["true_sigma"], floor_fraction=row["floor_fraction"],
                       floor_cv_threshold=row["floor_cv_threshold"]) for row in value["cells"])
    design = RunDesign(value["name"], value["seed"], value["replications"], cells, value["scope"])
    _require(canonical(design.to_dict()) == canonical(value), "saved renderer design changed")
    return design


def _check_inputs(prepared, guard=None, stage="readback"):
    if guard is not None:
        guard.check(stage + ":before_inputs")
    else:
        _require((prepared.source / "run_metadata.json").read_bytes() == prepared.root_bytes,
                 "root metadata changed during figure verification")
    _require(canonical(_current_snapshot(prepared)) == canonical(prepared.snapshot),
             "saved raw/table/audit/selected input bytes changed during rendering")
    if guard is not None:
        guard.check(stage + ":after_inputs")
    else:
        _require((prepared.source / "run_metadata.json").read_bytes() == prepared.root_bytes,
                 "root metadata changed during figure verification")


def _prepare_state(prepared, descriptor, extra_files):
    parent = prepared.destination if prepared.kind == "selected" else prepared.source
    path = renderer_state_path(parent)
    _state_layout(path)
    if path.exists():
        state, _ = _json_file(path)
        return path, _validate_renderer_state(state, descriptor, extra_files), False
    existing = {child.name for child in prepared.destination.iterdir()} if prepared.destination.exists() else set()
    _require(not (existing - {".execution.lock"}), "partial figures lack their original immutable renderer state")
    design = _RendererDesign(canonical(descriptor))
    provenance = capture(design, _dependency_set(extra_files), require_frozen=prepared.frozen)
    state = {"schema_version": 1, "rendering_contract": CONTRACT, "descriptor": descriptor,
             "descriptor_sha256": design.digest, "provenance": provenance}
    _validate_renderer_state(state, descriptor, extra_files)
    return path, state, True


def _state_unchanged(path, expected):
    _state_layout(path)
    _require(path.read_bytes() == expected, "immutable renderer state changed during rendering")


def _write_prepared(prepared, loader, extra_files):
    descriptor = _descriptor(prepared)
    existing = _figure_layout(prepared)
    if existing is not None:
        manifest = _verify_manifest(prepared, descriptor, extra_files)
        _check_inputs(prepared)
        return manifest
    _require(prepared.kind == "selected" or prepared.metadata["status"] != "completed",
             "completed archives permit verification only; missing figures cannot be recreated")
    lock_root = prepared.destination if prepared.kind == "selected" else prepared.source
    # Read-only preflight precedes directory creation, lock acquisition and any
    # rendering. Archived input roots themselves are never modified here.
    execution._output_path(lock_root)
    _check_incident_block(prepared.source)
    if prepared.kind == "selected":
        _require(not renderer_state_path(lock_root).parent.is_relative_to(prepared.source),
                 "renderer state must be disjoint from selected evidence")
    _prepare_state(prepared, descriptor, extra_files)
    _check_inputs(prepared)
    initial_guard = _ManifestGuard(prepared.source / "run_metadata.json", prepared.root_bytes)
    with run_lock(lock_root):
        initial_guard.check("after_figure_lock")
        prepared = loader()
        descriptor = _descriptor(prepared)
        existing = _figure_layout(prepared)
        if existing is not None:
            manifest = _verify_manifest(prepared, descriptor, extra_files)
            _check_inputs(prepared)
            return manifest
        path, state, create = _prepare_state(prepared, descriptor, extra_files)
        guard = _ManifestGuard(prepared.source / "run_metadata.json", prepared.root_bytes)
        figures = {}
        try:
            _check_inputs(prepared, guard, "before_renderer_state")
            if create:
                atomic_json(path, state, immutable=True)
            state_bytes = path.read_bytes()
            _require(state_bytes == (canonical(state) + "\n").encode(), "renderer state changed at creation")
            if prepared.kind == "selected":
                figures = {stem: build_selected_figure(prepared.payloads[stem], validation=not prepared.frozen)
                           for stem in prepared.stems}
            else:
                design = _saved_design(prepared.binding["design"])
                figures = (build_followup_figures(prepared.tables, design) if prepared.kind == "followup"
                           else build_floor_figures(prepared.tables, prepared.epsilon, design))
            _require(set(figures) == set(prepared.stems), "renderer builders omitted or added a planned figure")
            hashes = {}
            for stem in prepared.stems:
                _check_inputs(prepared, guard, "before_render:" + stem)
                _state_unchanged(path, state_bytes)
                payloads = figure_bytes(figures[stem], stem)
                _require(set(payloads) == {stem + ".png", stem + ".pdf"}
                         and all(type(value) is bytes and value for value in payloads.values()), "renderer save boundary returned an incomplete pair")
                _check_inputs(prepared, guard, "before_figure_commit:" + stem)
                _state_unchanged(path, state_bytes)
                _figure_layout(prepared)
                for name, payload in payloads.items():
                    atomic_bytes(prepared.destination / name, payload, immutable=True)
                    _require((prepared.destination / name).read_bytes() == payload, "figure changed during immutable publication")
                    hashes[name] = digest_bytes(payload)
                plt.close(figures[stem])
            _check_inputs(prepared, guard, "before_figure_manifest")
            _state_unchanged(path, state_bytes)
            _figure_layout(prepared)
            _validate_renderer_state(state, descriptor, extra_files)
            manifest = _manifest(prepared, state, hashes)
            atomic_json(prepared.destination / "figure_manifest.json", manifest, immutable=True)
            _verify_manifest(prepared, descriptor, extra_files)
            _check_inputs(prepared, guard, "after_figure_manifest")
            return manifest
        except BaseException:
            if not guard.drifted:
                guard.check("figure_render_failure")
            raise
        finally:
            for figure in figures.values():
                plt.close(figure)


def write_figures(kind, output_dir, design, *, extra_files=()):
    """Render/adopt run figures without changing raw, audit or lifecycle bytes."""
    extra_files = tuple(extra_files)
    loader = lambda: _load_run_inputs(kind, output_dir, design)
    return _write_prepared(loader(), loader, extra_files)


def verify_figures(kind, output_dir, design, *, extra_files=()):
    """Pure complete readback; the embedded renderer identity is self-contained."""
    prepared = _load_run_inputs(kind, output_dir, design)
    manifest = _verify_manifest(prepared, _descriptor(prepared), tuple(extra_files))
    _check_inputs(prepared)
    return manifest


def write_selected_figures(selected_dir, output_dir, *, extra_files=()):
    """Render the two fixed cases outside their immutable diagnostic archive."""
    extra_files = tuple(extra_files)
    loader = lambda: _load_selected_inputs(selected_dir, output_dir)
    return _write_prepared(loader(), loader, extra_files)


def verify_selected_figures(selected_dir, output_dir, *, extra_files=()):
    """Pure selected-figure readback, including retained review-required scope."""
    prepared = _load_selected_inputs(selected_dir, output_dir)
    manifest = _verify_manifest(prepared, _descriptor(prepared), tuple(extra_files))
    _check_inputs(prepared)
    return manifest
