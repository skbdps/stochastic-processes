"""Complete E2 heatmaps from saved, verified reporting tables.

The axes are the investigated discrete q/CV cells. No interpolation, smoothing,
or outcome-dependent grid selection is applied. Operational coverage uses all
planned paths; an unavailable interval is a noncovering outcome. Conditional
coverage is undefined when no interval is valid, and is masked rather than
silently replaced by the operational zero in that case.
"""
from io import BytesIO
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np

from ..experiments.provenance import file_digest
from ..experiments.storage import atomic_bytes, atomic_json
from .config import REGISTERED
from .reporting import load_tables, write_grid_reports


METHODS = (("exact", "Exact"), ("pfml", "PFML"), ("euler", "Euler"))
PARAMETERS = (("theta", r"$\theta$"), ("mu", r"$\mu$"), ("sigma", r"$\sigma$"))
CLASS_LABELS = {"within_criteria": "W", "supported_breach": "B", "unresolved": "U"}
REPORT_FILES = ("classification_table.csv", "classification_counts.csv",
                "failed_cells.csv", "report_manifest.json")
STYLE = {"font.family": "DejaVu Sans", "font.size": 10, "axes.titlesize": 11,
         "figure.facecolor": "white", "axes.grid": False, "pdf.fonttype": 42, "ps.fonttype": 42}
NA_COLOR = "#d8d8d8"

OPERATIONAL_NOTE = (
    "Operational coverage: all planned paths; unavailable intervals count as noncoverage. Colors show point estimates on a fixed 0–100% scale.\n"
    "W = within criteria; B = supported breach; U = unresolved. These pointwise labels combine operational-coverage Wilson bounds AND relative-bias MC intervals.\n"
    "A breach can therefore reflect bias even with high coverage. Mu has zero truth and no relative-bias classification. Complete uncertainty is in the saved tables.\n"
    "Intervals are nominal 95% natural-scale Wald intervals; PFML/Euler covariance is model-based. Discrete grid cells do not define a smooth or simultaneous boundary."
)
CONDITIONAL_NOTE = (
    "Conditional coverage: covering intervals / S valid intervals. S is shown in every cell; NA with gray fill means S = 0. A colored 0.00% means S > 0 and no coverage.\n"
    "Every planned q/CV cell is retained. Compare with operational coverage and availability; conditioning on valid intervals can select a different subset.\n"
    "Intervals are nominal 95% symmetric natural-scale Wald intervals; PFML/Euler use model-based inverse-Hessian covariance.\n"
    "Saved tables retain pointwise Wilson Monte Carlo intervals. Colors show estimates, without interpolation or simultaneous guarantees."
)
AVAILABILITY_NOTE = (
    "Availability: valid count / all planned paths. S is the valid count for the indicated row. Point and interval validity are joint across theta, mu and sigma.\n"
    "An interval can be unavailable despite a valid point fit. These categories overlap; do not add point-invalid and interval-unavailable counts.\n"
    "Numerical unavailability differs from statistical noncoverage. The complete failed-cell report retains every cell/method, including zero-failure rows.\n"
    "Colors use the same fixed 0–100% scale. Saved tables retain pointwise Wilson Monte Carlo intervals and detailed outcome counts."
)


def _layout(design):
    """A later experiment must not accidentally inherit an E1-sized plot."""
    ns = tuple(REGISTERED["n"])
    qs = tuple(REGISTERED["q"])
    cvs = tuple(REGISTERED["cv"])
    expected = {(n, q, cv) for n in ns for q in qs for cv in cvs}
    actual = [(c.n, c.q, c.cv) for c in design.cells]
    if len(actual) != len(expected) or set(actual) != expected:
        raise ValueError("Week 6 figures require the complete registered E2 n/q/CV grid")
    if any((c.theta, c.mu, c.sigma, c.floor_fraction, c.floor_cv_threshold)
           != (1.0, 0.0, 0.5, REGISTERED["floor_fraction"], REGISTERED["floor_cv_threshold"])
           for c in design.cells):
        raise ValueError("Week 6 figures require the E2 truth and floor settings")
    if (design.scope not in ("headline", "validation")
            or (design.scope == "headline" and (design.name, design.seed, design.replications)
                != ("week6_e2", REGISTERED["seed"], REGISTERED["replications"]))):
        raise ValueError("Week 6 figures require the registered E2 identity or an explicit validation scope")
    return ns, qs, cvs


def _matrix(summary, n, method, parameter, column, qs, cvs):
    """Select a complete matrix with explicit q rows and CV columns.

    Checking the Cartesian keys before reindexing prevents a pivot or reindex
    from disguising a missing group as an unavailable statistical quantity.
    """
    group = summary[(summary.n == n) & (summary.estimator == method)
                    & (summary.parameter == parameter)]
    keys = list(group[["theta_mean_gap", "cv"]].itertuples(index=False, name=None))
    expected = {(q, cv) for q in qs for cv in cvs}
    if len(keys) != len(expected) or set(keys) != expected:
        raise ValueError(f"incomplete heatmap matrix: n={n}, {method}, {parameter}")
    indexed = group.set_index(["theta_mean_gap", "cv"])[column]
    return np.asarray([[indexed.loc[(q, cv)] for cv in cvs] for q in qs])


def _counts(summary, n, method, parameter, qs, cvs):
    counts = [_matrix(summary, n, method, parameter, name, qs, cvs).astype(float)
              for name in ("n_planned", "n_point_valid", "n_interval_valid", "n_covered")]
    planned, point, interval, covered = counts
    if (any(not np.isfinite(a).all() or not (a == np.floor(a)).all() for a in counts)
            or np.any(planned <= 0) or np.any(point > planned) or np.any(interval > point)
            or np.any(covered > interval) or np.any(covered < 0) or np.any(interval < 0)):
        raise ValueError("heatmap validity/coverage counts are inconsistent")
    return planned, point, interval, covered


def _coverage_matrix(summary, n, method, parameter, qs, cvs, scope):
    planned, _, interval, covered = _counts(summary, n, method, parameter, qs, cvs)
    if scope not in ("operational", "conditional"):
        raise ValueError("unknown coverage scope")
    denominator = planned if scope == "operational" else interval
    expected = np.divide(covered, denominator, out=np.full_like(covered, np.nan), where=denominator > 0)
    values = _matrix(summary, n, method, parameter, f"coverage_{scope}", qs, cvs).astype(float)
    # A valid zero-coverage cell is distinct from a zero-denominator cell. Do
    # not use fillna(0), clipping, or a color scale computed from finite cells.
    if (np.isinf(values).any() or np.any(np.isfinite(values) & ((values < 0) | (values > 1)))
            or not np.allclose(values, expected, rtol=1e-12, atol=1e-14, equal_nan=True)):
        raise ValueError("heatmap coverage does not match its saved denominator")
    labels = _matrix(summary, n, method, parameter, "practical_classification", qs, cvs)
    annotations = np.empty(values.shape, dtype=object)
    for index in np.ndindex(values.shape):
        value = values[index]
        if scope == "conditional":
            annotations[index] = ("NA" if np.isnan(value) else f"{100*value:.2f}%") + f"\nS={int(interval[index]):,}"
        elif parameter == "mu":
            if labels[index] != "not_applicable_zero_truth":
                raise ValueError("mu must not receive a truth-relative practical classification")
            annotations[index] = f"{100*value:.2f}%"
        else:
            if labels[index] not in CLASS_LABELS:
                raise ValueError("unknown practical classification in heatmap")
            # The audited joint classification is copied, never inferred from
            # the rounded displayed coverage or from conditional coverage.
            annotations[index] = f"{100*value:.2f}%\n{CLASS_LABELS[labels[index]]}"
    return np.ma.masked_invalid(100*values), annotations


def _availability_matrix(summary, n, method, kind, qs, cvs):
    if kind not in ("point", "interval"):
        raise ValueError("unknown availability scope")
    # write_grid_reports verifies joint counts across all three parameters.
    # Each path is counted once here, using theta's copy of those joint counts.
    planned, point, interval, _ = _counts(summary, n, method, "theta", qs, cvs)
    valid = point if kind == "point" else interval
    values = _matrix(summary, n, method, "theta", f"{kind}_availability", qs, cvs).astype(float)
    if not np.isfinite(values).all() or not np.allclose(values, valid/planned, rtol=1e-12, atol=1e-14):
        raise ValueError("heatmap availability does not match its saved count")
    annotations = np.asarray([[f"{100*values[i, j]:.2f}%\nS={int(valid[i, j]):,}"
                               for j in range(len(cvs))] for i in range(len(qs))])
    return np.ma.array(100*values), annotations


def _draw_panel(ax, values, annotations, qs, cvs, title):
    cmap = matplotlib.colormaps["viridis"].copy()
    cmap.set_bad(NA_COLOR)
    image = ax.imshow(values, interpolation="nearest", origin="upper", aspect="auto",
                      cmap=cmap, norm=Normalize(vmin=0, vmax=100))
    ax.set_xticks(range(len(cvs)), [f"{value:g}" for value in cvs])
    ax.set_yticks(range(len(qs)), [f"{value:g}" for value in qs])
    ax.set_xlabel("Nominal gap CV")
    ax.set_ylabel(r"$q = \theta\,\overline{\Delta}_{\rm nominal}$")
    ax.set_title(title, pad=7)
    ax.set_xticks(np.arange(len(cvs)+1)-.5, minor=True)
    ax.set_yticks(np.arange(len(qs)+1)-.5, minor=True)
    ax.grid(which="minor", color="white", linewidth=.8)
    ax.tick_params(which="minor", bottom=False, left=False)
    for index in np.ndindex(values.shape):
        rgba = image.cmap.get_bad() if np.ma.is_masked(values[index]) else image.to_rgba(values[index])
        luminance = np.dot(rgba[:3], [.2126, .7152, .0722])
        ax.text(index[1], index[0], annotations[index], ha="center", va="center",
                fontsize=9.2, linespacing=1.35, color="#171717" if luminance > .52 else "white")
    return image


def _save(fig, output, stem, title):
    for extension in ("png", "pdf"):
        stream = BytesIO()
        metadata = ({"Title": title, "Creator": "OU irregular-sampling research",
                     "CreationDate": None, "ModDate": None} if extension == "pdf"
                    else {"Software": "OU irregular-sampling research"})
        fig.savefig(stream, format=extension, dpi=170, metadata=metadata)
        atomic_bytes(output / f"{stem}.{extension}", stream.getvalue(), immutable=True)


def _render_one(output, design, n, qs, cvs, stem, title, note, panels, *, availability=False):
    rows = 2 if availability else 3
    fig, axes = plt.subplots(rows, 3, figsize=(14.4, 8.8 if availability else 11.7), squeeze=False)
    try:
        for ax, (values, annotations, panel_title) in zip(axes.flat, panels):
            image = _draw_panel(ax, values, annotations, qs, cvs, panel_title)
        top, bottom = (.84, .23) if availability else (.88, .18)
        fig.subplots_adjust(left=.065, right=.895, bottom=bottom, top=top, wspace=.30, hspace=.43)
        color_ax = fig.add_axes([.929, bottom+.025, .016, top-bottom-.05])
        colorbar = fig.colorbar(image, cax=color_ax, ticks=[0, 20, 40, 60, 80, 100])
        colorbar.set_label("Availability (%)" if availability else "Coverage (%)", labelpad=10)
        label = "REGISTERED E2" if design.scope == "headline" else "VALIDATION ONLY"
        fig.suptitle(title, fontsize=16, fontweight="semibold", y=.975)
        fig.text(.5, .935 if not availability else .92,
                 f"{label}  |  n = {n:,} observations  |  R = {design.replications:,} per (q, CV) cell",
                 ha="center", fontsize=10.5)
        fig.text(.065, .035, note, va="bottom", fontsize=8.5, color="#333333", linespacing=1.6)
        _save(fig, output, stem, title)
    finally:
        plt.close(fig)


def render_figures(run_dir, design, *, output_dir=None):
    """Write nine PNG/PDF pairs and a manifest binding every derived report.

    Rendering is deterministic in the recorded environment. Existing outputs
    are accepted only if their bytes are identical; a different result requires
    a separate output directory, not an overwrite of scientific evidence.
    """
    ns, qs, cvs = _layout(design)
    root = Path(run_dir)
    tables, inputs = load_tables(root, design)
    summary = tables["parameter_summary.csv"]
    prepared = []
    # Validate every matrix before writing any figure. This includes all three
    # parameters even when their interval availability happens to be identical.
    for n in ns:
        for scope, prefix, title, note in (
            ("operational", "F5", "F5  Operational coverage and practical classifications", OPERATIONAL_NOTE),
            ("conditional", "F6", "F6  Coverage conditional on available intervals", CONDITIONAL_NOTE),
        ):
            panels = [(*_coverage_matrix(summary, n, method, parameter, qs, cvs, scope), f"{symbol} | {label}")
                      for parameter, symbol in PARAMETERS for method, label in METHODS]
            prepared.append((n, f"{prefix}_{scope}_n{n}", title, note, panels, False))
        panels = [(*_availability_matrix(summary, n, method, kind, qs, cvs), f"{title} | {label}")
                  for kind, title in (("point", "Valid point fits"), ("interval", "Valid intervals"))
                  for method, label in METHODS]
        prepared.append((n, f"F6_availability_n{n}", "F6 companion  Point-fit and interval availability",
                         AVAILABILITY_NOTE, panels, True))

    names = {f"{item[1]}.{extension}" for item in prepared for extension in ("png", "pdf")} | set(REPORT_FILES)
    output = Path(output_dir) if output_dir is not None else root / "figures"
    if output.exists() and {p.name for p in output.iterdir()} - names - {"figure_manifest.json"}:
        raise ValueError("unexpected artifact in Week 6 figure directory")
    output.mkdir(parents=True, exist_ok=True)
    report = write_grid_reports(root, design, output_dir=output)
    if report["input_table_sha256"] != inputs:
        raise ValueError("report inputs changed during figure preparation")
    with plt.rc_context(STYLE):
        for n, stem, title, note, panels, availability in prepared:
            _render_one(output, design, n, qs, cvs, stem, title, note, panels, availability=availability)
    current_manifest = json.loads((root / "table_manifest.json").read_text())
    if any(file_digest(root / name) != digest or current_manifest[name]["sha256"] != digest
           for name, digest in inputs.items()):
        raise ValueError("figure source tables changed during rendering")
    record = {
        "design_sha256": design.digest,
        "input_table_sha256": inputs,
        # The frozen workflow checks the exact figure-directory file set. Its
        # generic 'figures' map deliberately includes all four report artifacts.
        "figures": {name: file_digest(output / name) for name in sorted(names)},
        "grid": {"n": list(ns), "q_rows": list(qs), "cv_columns": list(cvs)},
        "coverage_color_limits_percent": [0, 100],
        "classification_labels": CLASS_LABELS,
        "classification_scope": "audited joint operational-coverage and relative-bias criteria; pointwise",
        "conditional_unavailable": "NA with gray fill when interval-valid S is zero",
    }
    atomic_json(output / "figure_manifest.json", record, immutable=True)
    return record
