"""E4 single-axis figures from immutable, fully audited saved tables.

No sampling, fitting, source reconstruction or summary recalculation occurs here.
The saved full audit is bound to every raw/table byte and its fixed endpoint
identities. Rendering has its own frozen provenance; it never rewrites an older
raw-run dependency manifest. Each figure is a separate PNG/PDF pair.
"""
from dataclasses import dataclass
from importlib.metadata import version
from io import BytesIO
import json
import math
from pathlib import Path
import re
import subprocess
import textwrap
import zlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..experiments.provenance import CORE_FILES, ROOT, canonical, capture, digest_bytes, file_digest, verify
from ..experiments.storage import atomic_bytes, atomic_json, run_lock
from . import audit, execution, reporting
from .design import parse_design, validate_design

FILES = tuple(sorted(set(audit.FILES) | {"ou_irregular/week8/plots.py", "docs/week8_plots.md"}))
CONTRACT = "week8_saved_evidence_single_axis_png_pdf_v1"
METHODS = (("exact", "Exact", "o"), ("pfml", "PFML", "s"), ("euler", "Euler", "^"))
PARAMETERS = (("theta", r"$\theta$"), ("mu", r"$\mu$"), ("sigma", r"$\sigma$"))
STEMS = tuple([f"W8_coverage_{p}" for p, _ in PARAMETERS] + ["W8_availability"] +
              [f"W8_matched_{metric}_{p}" for p, _ in PARAMETERS for metric in ("bias", "coverage", "rmse_ratio")] +
              ["W8_source_mean", "W8_source_adjacent_product", "W8_source_cv", "W8_source_span_sd", "W8_source_clipping"])
FILENAMES = frozenset(f"{stem}.{ext}" for stem in STEMS for ext in ("png", "pdf"))
# No palette or Matplotlib style is prescribed. Versioned defaults are retained;
# only typography, canvas geometry, output font encoding and backend are fixed.
SETTINGS = {**dict(matplotlib.rcParamsDefault), "backend": "Agg", "font.family": ["DejaVu Sans"],
            "font.size": 10., "axes.titlesize": 13., "pdf.fonttype": 42, "ps.fonttype": 42}
Z = 1.959963984540054
_CLOSED = {"audited_pending_visual_review", "completed"}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(value):
    return digest_bytes(canonical(value).encode())


def _json(path):
    return execution._json_file(path)[0]


def _row(frame, **keys):
    mask = np.ones(len(frame), dtype=bool)
    for name, value in keys.items():
        mask &= frame[name].eq(value).to_numpy()
    selected = frame.loc[mask]
    require(len(selected) == 1, "figure data do not identify exactly one planned row")
    return selected.iloc[0]


def _ordered(frame, cells, **keys):
    return pd.DataFrame([_row(frame, cell_id=c.identity, **keys) for c in cells])


def validate_tables(tables, design):
    """Pure builder boundary: require complete schemas/keys/identities, not a new audit."""
    validate_design(design)
    require(set(tables) == set(reporting.TABLE_KEYS), "figure builders require all eight uncompressed tables")
    cells = {c.identity: c.to_dict() for c in design.cells}
    for name, table in tables.items():
        require(tuple(table.columns) == reporting.TABLE_COLUMNS[name], "figure table schema changed: "+name)
        keys = list(reporting.TABLE_KEYS[name])
        require(not table.duplicated(keys).any() and table.cell_id.isin(cells).all(), "figure keys changed")
        for field in next(iter(cells.values())):
            wanted = table.cell_id.map(lambda cid: cells[cid][field])
            require(table[field].eq(wanted).all(), "figure scientific metadata changed: "+field)
    summary = reporting._summary_input(tables["parameter_summary.csv"], design)
    for kind in ("regular", "matched_iid"):
        name = kind+"_comparisons.csv"
        expected = reporting._planned_keys(name, design, {})
        actual = set(tables[name][list(reporting.TABLE_KEYS[name])].itertuples(index=False, name=None))
        require(actual == expected, "figure contrast keys incomplete")
    for name in ("source_summary.csv", "spacing_summary.csv"):
        require(set(tables[name].cell_id) == set(cells) and len(tables[name]) == len(cells), "source figure cells incomplete")
    return summary


def _axes(design, title, ylabel, *, matched=False):
    fig, ax = plt.subplots(figsize=(10.6, 7.0))
    fig.subplots_adjust(left=.11, right=.97, bottom=.31, top=.79)
    fig.suptitle(title, fontsize=15, y=.965)
    role = "VALIDATION ONLY — not E4 results" if design.scope == "validation" else "E4 exploratory experiment"
    fig.text(.5, .916, f"{role}  |  R = {design.replications:,} per cell  |  n = 1,000  |  θm = 0.5", ha="center", fontsize=10)
    if matched:
        ax.set_xticks(range(3), ["r = 1", "r = 5", "r = 20"])
        ax.set_xlabel("Matched marginal-gap pairs (independent MMPP and IID paths)")
    else:
        ax.set_xticks(range(7), ["MMPP\nr = 1", "MMPP\nr = 5", "MMPP-derived\nr = 20", "IID\nr = 1", "IID\nr = 5", "IID floored\nr = 20", "Regular\ngap m"])
        ax.set_xlabel("Seven distinct, independently simulated cells")
    ax.set_ylabel(ylabel)
    ax.grid(axis="y", alpha=.2)
    return fig, ax


def _note(fig, text):
    lines = textwrap.wrap(text, width=119)
    require(len(lines) <= 6, "figure caption exceeds its reserved space")
    fig.text(.07, .075, "\n".join(lines), ha="left", va="bottom", fontsize=9)


def _series(ax, x, y, lower, upper, *, label, marker, gid, ratio=False):
    x, y, lower, upper = (np.asarray(a, dtype=float) for a in (x, y, lower, upper))
    require(x.shape == y.shape == lower.shape == upper.shape, "figure coordinate lengths differ")
    require(np.isfinite(x).all() and not any(np.isinf(a).any() for a in (y, lower, upper)), "infinite figure values")
    finite = np.isfinite(y)
    interval = finite & np.isfinite(lower) & np.isfinite(upper)
    require(np.all(lower[interval] <= y[interval]) and np.all(y[interval] <= upper[interval]), "figure error-bar direction changed")
    require(not ratio or np.all(y[finite] >= 0), "RMSE ratios cannot be negative")
    shown = finite & (y > 0) if ratio else finite
    require(not ratio or np.all(lower[interval & shown] > 0), "positive log-ratio bounds must be positive")
    # One errorbar call per named series keeps the default color cycle stable.
    errors = np.vstack((np.where(interval, y-lower, 0.), np.where(interval, upper-y, 0.)))
    artist = ax.errorbar(x[shown], y[shown], yerr=errors[:, shown], fmt=marker, linestyle="none", capsize=3, label=label)
    artist.lines[0].set_gid(gid)
    for i in np.flatnonzero(~finite):
        ax.text(x[i], .035, "NA", ha="center", fontsize=8, transform=ax.get_xaxis_transform())
    if ratio:
        for i in np.flatnonzero(finite & ~shown):
            ax.text(x[i], .035, "0", ha="center", fontsize=8, transform=ax.get_xaxis_transform())
    for i in np.flatnonzero(shown & ~interval):
        ax.annotate("†", (x[i], y[i]), xytext=(4, 5), textcoords="offset points", fontsize=9)
    return artist


def _legend(ax, columns=3):
    ax.legend(loc="lower center", bbox_to_anchor=(.5, 1.015), ncol=columns, fontsize=9, frameon=False)


def build_figures(tables, design):
    """Build exactly 18 independent single-axis figures, exposing saved series by gid."""
    summary = validate_tables(tables, design)
    figures = {}
    try:
        with plt.rc_context(SETTINGS):
            for parameter, symbol in PARAMETERS:
                stem = f"W8_coverage_{parameter}"
                fig, ax = _axes(design, f"Operational coverage of {symbol}", "Coverage (%)")
                figures[stem] = fig
                for index, (method, label, marker) in enumerate(METHODS):
                    data = _ordered(summary, design.cells, estimator=method, parameter=parameter)
                    _series(ax, np.arange(7)+(index-1)*.20, 100*data.coverage_operational,
                            100*data.operational_wilson_lower, 100*data.operational_wilson_upper,
                            label=label, marker=marker, gid=method+"/operational_coverage")
                ax.axhline(95, linestyle="--", linewidth=1, label="Nominal 95%")
                ax.set_ylim(-2, 103); _legend(ax, 4)
                _note(fig, "Operational coverage divides covered intervals by all planned paths; unavailable intervals count as noncoverage. "
                      "Bars are pointwise 95% Wilson Monte Carlo intervals. PFML/Euler intervals are model-based. "
                      "Exact counts, conditional coverage and failures remain in the saved tables. No simultaneous-coverage claim is made.")
            fig, ax = _axes(design, "Point-estimate and interval availability", "Available records (%)")
            figures["W8_availability"] = fig
            for index, (method, label, marker) in enumerate(METHODS):
                data = _ordered(summary, design.cells, estimator=method, parameter="theta")
                for kind, shift, mark in (("point", -.055, marker), ("interval", .055, "x")):
                    _series(ax, np.arange(7)+(index-1)*.24+shift, 100*data[kind+"_availability"],
                            100*data[kind+"_availability_wilson_lower"], 100*data[kind+"_availability_wilson_upper"],
                            label=label+" "+kind, marker=mark, gid=method+"/"+kind+"_availability")
            ax.set_ylim(-2, 103); _legend(ax, 3)
            _note(fig, "Both denominators are all R planned paths. Availability is joint across the three parameters; optimizer success "
                  "alone does not establish point or interval validity. Bars are pointwise 95% Wilson Monte Carlo intervals. "
                  "The full ordered failure export retains every unavailable record, even when none occur.")
            contrasts = tables["matched_iid_comparisons.csv"]
            pairs = tuple(c for c in design.cells if c.mechanism == "mmpp")
            for parameter, symbol in PARAMETERS:
                for metric in ("bias", "coverage", "rmse_ratio"):
                    stem = f"W8_matched_{metric}_{parameter}"
                    field, low, high, factor, title, ylabel = {
                        "bias": ("excess_bias", "excess_bias_mc_lower", "excess_bias_mc_upper", 1., f"Matched-control difference in {symbol} bias", "Bias: MMPP minus IID (parameter units)"),
                        "coverage": ("operational_coverage_difference", "operational_coverage_difference_mc_lower", "operational_coverage_difference_mc_upper", 100., f"Matched-control difference in {symbol} coverage", "Coverage: MMPP minus IID (percentage points)"),
                        "rmse_ratio": ("rmse_ratio", "rmse_ratio_mc_lower", "rmse_ratio_mc_upper", 1., f"Matched-control {symbol} RMSE ratio", "RMSE: MMPP / IID (log scale)"),
                    }[metric]
                    fig, ax = _axes(design, title, ylabel, matched=True); figures[stem] = fig
                    for index, (method, label, marker) in enumerate(METHODS):
                        data = _ordered(contrasts, pairs, estimator=method, parameter=parameter)
                        _series(ax, np.arange(3)+(index-1)*.18, factor*data[field], factor*data[low], factor*data[high],
                                label=label, marker=marker, gid=method+"/"+field, ratio=metric == "rmse_ratio")
                    ax.axhline(1 if metric == "rmse_ratio" else 0, linestyle="--", linewidth=1)
                    if metric == "rmse_ratio": ax.set_yscale("log")
                    _legend(ax)
                    _note(fig, "The full one-gap marginal and floor transform match within each MMPP/IID pair; paths are independent. "
                          "Realized span is not held fixed. Bias/RMSE use each cell's own valid points; operational coverage uses all paths. "
                          "Bars show pointwise Monte Carlo uncertainty (normal differences; log-delta ratios). NA: undefined point; †: uncertainty unavailable; 0: zero ratio off log scale.")
            source = _ordered(tables["source_summary.csv"], design.cells)
            for metric, suffix, title, ylabel in (("mean", "mean_gap", "Mean observation gap", "Mean gap (model time)"),
                                                  ("adjacent_product", "adjacent_product", "Adjacent-gap cross-product", "Mean adjacent product (model time squared)"),
                                                  ("cv", "cv", "Average within-path gap CV", "Average within-path CV")):
                fig, ax = _axes(design, title, ylabel); figures["W8_source_"+metric] = fig
                for offset, kind, label, marker in ((-.12, "raw", "Raw gaps", "o"), (.12, "transformed", "After floor", "s")):
                    values = source[f"{kind}_{suffix}_mean"]
                    errors = Z*source[f"{kind}_{suffix}_mcse"]
                    _series(ax, np.arange(7)+offset, values, values-errors, values+errors,
                            label=label, marker=marker, gid=kind+"/"+suffix)
                if metric != "cv":
                    artist, = ax.plot(np.arange(7), source["expected_transformed_"+suffix], "x", linestyle="none", label="Analytic after-floor target")
                    artist.set_gid("analytic/transformed_"+suffix)
                _legend(ax)
                note = ("Average within-path CV is not the CV of the pooled marginal gap law; finite-path variability and dependence can affect it. "
                        if metric == "cv" else "Analytical targets describe the declared gap law, not new simulation draws. ")
                _note(fig, note+"Bars are 1.96 Monte Carlo standard errors across independent source paths, never across dependent gaps or repeated estimators. "
                      "Only r = 20 is floored; no renormalization or hidden-chain evolution through floor padding occurs. Near-overlap of raw and floored points is retained.")
            fig, ax = _axes(design, "Variation in realized observation span", "Sample standard deviation of span (model time)")
            figures["W8_source_span_sd"] = fig
            values = source.realized_span_mcse*np.sqrt(source.realized_span_n_available)
            artist, = ax.plot(np.arange(7), values, "o", linestyle="none", label="Across-path span SD")
            artist.set_gid("realized/span_sd"); _legend(ax, 1)
            _note(fig, "Sample SD uses independent source paths and is recovered as stored MCSE × √R (ddof = 1). No uncertainty interval for the SD is implied. "
                  "Matching a marginal gap law does not hold the realized observation span fixed. This figure separates span variability from a change in the nominal mean gap.")
            fig, ax = _axes(design, "Observed posthoc floor clipping", "Clipped gaps per million generated gaps")
            figures["W8_source_clipping"] = fig
            artist, = ax.plot(np.arange(7), 1e6*source.clip_fraction_mean, "o", linestyle="none", label="Observed clipped fraction")
            artist.set_gid("source/clipping"); _legend(ax, 1)
            ax.set_ylim(0, max(1., 1.15*float((1e6*source.clip_fraction_mean).max())))
            counts = ", ".join(str(int(x)) for x in source.clip_count_total)
            _note(fig, f"Clip counts in displayed cell order: {counts}; each cell contains {999*design.replications:,} generated gaps. "
                  "No binomial interval treats dependent MMPP gaps as independent. Only r = 20 has a positive floor (10⁻⁶m); zero observed clips is not evidence of zero clipping probability. "
                  "Full raw, transformed and stored-time diagnostics remain in source_summary.csv.")
        require(set(figures) == set(STEMS) and all(len(fig.axes) == 1 for fig in figures.values()), "figure plan incomplete or multiple panels introduced")
        return figures
    except BaseException:
        for fig in figures.values(): plt.close(fig)
        raise


def figure_bytes(figure, stem):
    result = {}
    with plt.rc_context(SETTINGS):
        for extension in ("png", "pdf"):
            stream = BytesIO()
            metadata = ({"Title": stem, "Creator": "OU irregular-sampling research", "CreationDate": None, "ModDate": None}
                        if extension == "pdf" else {"Software": "OU irregular-sampling research"})
            figure.savefig(stream, format=extension, dpi=170, metadata=metadata)
            result[stem+"."+extension] = stream.getvalue()
    return result


def _effective(value):
    if value is None or isinstance(value, (str, int, float, bool)): return value
    if isinstance(value, (list, tuple)): return [_effective(x) for x in value]
    if hasattr(value, "by_key"): return {k: _effective(v) for k, v in value.by_key().items()}
    raise ValueError("unrecorded renderer setting")


def renderer_environment():
    from matplotlib import ft2font
    fonts = Path(matplotlib.get_data_path())/"fonts/ttf"
    return {"matplotlib": matplotlib.__version__, "freetype": ft2font.__freetype_version__, "pillow": version("pillow"),
            "fonttools": version("fonttools"), "zlib": zlib.ZLIB_RUNTIME_VERSION, "png_dpi": 170,
            "settings": {k: _effective(v) for k, v in SETTINGS.items()},
            "font_sha256": {p.name: file_digest(p) for p in sorted(fonts.glob("*.ttf"))},
            "pdf_creation_date": None, "pdf_modification_date": None}


def _scientific_snapshot(output):
    return {str(p.relative_to(output)): file_digest(p) for p in sorted(output.rglob("*"))
            if p.is_file() and p.name != ".execution.lock" and not p.is_relative_to(output/"figures")
            and p.name != "visual_review.json"}


def _endpoints_and_counts(output, design, binding):
    endpoints, point, interval, partitions = [], 0, 0, 0
    coordinates = {"point_coordinate_rows": 0, "interval_coordinate_rows": 0, "coordinate_scalars_checked": 0}
    for cell in sorted(design.cells, key=lambda c: c.identity):
        rows, cb = audit.load_cell(output, cell, design, binding["run_identity"])
        binding["cells"][cell.identity] = cb
        partitions += len(cb["partitions"])
        for row in rows:
            point += row["point_valid"]; interval += row["interval_valid"]
            for kind, names in (("point", ("state_center", "state_scale", "time_scale")),
                                ("interval", ("state_offset", "state_scale", "time_scale"))):
                diagnostic = audit.strict_json(row[kind+"_diagnostics_json"])
                present = [name in diagnostic for name in names]
                require(not any(present) or all(present), "incomplete saved coordinate triple")
                if all(present):
                    coordinates[kind+"_coordinate_rows"] += 1; coordinates["coordinate_scalars_checked"] += 3
        for rep in (0, design.replications-1):
            triplet = [row for row in rows if row["rep"] == rep]
            require(len(triplet) == 3 and [r["estimator"] for r in triplet] == list(audit.METHODS), "saved endpoint incomplete")
            hashes = {key: triplet[0][key] for key in ("times_sha256", "observations_sha256")}
            require(all(all(row[k] == v for k, v in hashes.items()) for row in triplet), "saved endpoint hash mismatch")
            endpoints.append({"cell_id": cell.identity, "replication": rep, **hashes, "estimators_checked": list(audit.METHODS)})
    return endpoints, {"point_valid_records": point, "interval_valid_records": interval, "partitions_verified": partitions, **coordinates}


@dataclass(frozen=True)
class Inputs:
    output: Path
    design: object
    metadata: dict
    root_bytes: bytes
    tables: dict
    binding: dict
    snapshot: dict


def load_inputs(output, design):
    """Pure no-sampling readback of a full saved audit and all its bound evidence."""
    validate_design(design)
    output = execution._ordinary_path(output)
    execution._root_layout(output, design)
    root_bytes = (output/"run_metadata.json").read_bytes()
    before = _scientific_snapshot(output)
    preliminary = _json(output/"run_metadata.json")
    files = preliminary["provenance"]["source_files"]
    require(type(files) is dict and set(CORE_FILES) | set(audit.FILES) <= set(files), "raw run lacks the audited dependency closure")
    raw_extra = tuple(sorted(set(files) - set(CORE_FILES) - set(audit.FILES)))
    metadata, raw = audit._metadata(output, design, raw_extra)
    saved = _json(output/"audit.json")
    tables = {name: pd.read_csv(output/name, float_precision="round_trip") for name in reporting.TABLE_KEYS}
    validate_tables(tables, design)
    manifest = audit._table_manifest(output)
    endpoints, counts = _endpoints_and_counts(output, design, raw)
    binding = _json(output/"table_binding.json")
    expected_binding = {"schema_version": 1, "project": "ou-week8-e4-reporting", "reporting_contract": reporting.REPORTING_CONTRACT,
                        "run_identity": raw["run_identity"], "design_sha256": design.digest, "raw_inputs": raw, "raw_input_sha256": digest(raw)}
    require(canonical(binding) == canonical(expected_binding), "saved raw/table binding changed")
    paths = 7*design.replications
    artifacts = {k: v for k, v in before.items() if k not in {"run_metadata.json", "audit.json"}}
    expected = {"schema_version": 1, "audit_contract": audit.AUDIT_CONTRACT, "status": "passed", "full_completion_audit": True,
                "design_sha256": design.digest, "run_identity_sha256": raw["run_identity"], "code_commit": metadata["provenance"]["code_commit"],
                "source_dependency_files_verified": len(files), "audited_artifact_sha256": artifacts,
                "audited_run_metadata_sha256": audit.semantic_metadata_digest(metadata), "raw_input_sha256": digest(raw),
                "scientific_overrides": [], "source_policy": metadata["source_policy"], "cells_verified": 7,
                "planned_paths": paths, "planned_source_realizations": paths, "planned_independent_source_paths": paths,
                "planned_fit_records": 3*paths, "source_replay_performed": True, "source_paths_verified": paths,
                "source_array_hashes_verified": 2*paths, "schedule_trace_hashes_verified": 3*paths, "pipeline_errors": 0,
                "retained_failure_records": manifest["failed_fits.csv.gz"]["rows"], **counts,
                "tables_verified": {name: entry["rows"] for name, entry in manifest.items()}, "replay_performed": True,
                "replayed_source_paths": 14, "replayed_fits": 42, "replay_records": endpoints,
                "replay_selection": "replications 0 and R-1 in every E4 cell, all estimators",
                "tolerances": {"arithmetic_rtol": 2e-10, "arithmetic_atol": 2e-12, "replay_rtol": 2e-10, "replay_atol": 2e-12,
                  "source_hashes": "exact SHA-256 of both arrays, both gap vectors and full phase/switch trace",
                  "covariance_rtol": 2e-9, "covariance_atol": "2e-12 * max(1, expected matrix Frobenius norm)"}}
    require(all(k in saved and canonical(saved[k]) == canonical(v) for k, v in expected.items()), "full saved audit does not bind these exact inputs/counts/endpoints")
    require(metadata["status"] not in _CLOSED or metadata.get("audit_sha256") == file_digest(output/"audit.json"), "root audit seal changed")
    # Completed archives must retain their full seal; partial figures are allowed
    # only while the raw/audited run is still awaiting actual visual completion.
    if metadata["status"] == "completed": audit._closed_bindings(output, metadata, audit._snapshot(output, design))
    identity = {"design": design.to_dict(), "design_sha256": design.digest, "run_identity_sha256": raw["run_identity"],
                "raw_input_sha256": digest(raw), "audited_artifact_sha256": artifacts,
                "semantic_metadata_sha256": audit.semantic_metadata_digest(metadata), "audit_sha256": file_digest(output/"audit.json"),
                "input_table_sha256": {name: entry["sha256"] for name, entry in manifest.items()}}
    require((output/"run_metadata.json").read_bytes() == root_bytes and _scientific_snapshot(output) == before, "saved inputs changed during readback")
    return Inputs(output, design, metadata, root_bytes, tables, identity, before)


@dataclass(frozen=True)
class RendererDesign:
    descriptor: dict
    def to_dict(self): return self.descriptor
    @property
    def digest(self): return digest(self.descriptor)


def _dependencies(extra_files):
    return tuple(sorted(set(CORE_FILES) | set(FILES) | set(execution._dependencies(extra_files))))


def _descriptor(inputs):
    return {"schema_version": 1, "rendering_contract": CONTRACT, "input_binding": inputs.binding,
            "require_frozen": inputs.design.scope == "exploratory", "figure_filenames": sorted(FILENAMES),
            "renderer_environment": renderer_environment()}


def renderer_state_path(output):
    output = execution._ordinary_path(output)
    return output.with_name(output.name+"-figure-renderer")/"renderer_state.json"


def _state_layout(path):
    execution._ordinary_path(path.parent)
    execution._reject_archive_ancestors(path.parent)
    if path.parent.exists():
        require(path.parent.is_dir() and {p.name for p in path.parent.iterdir()} == {path.name}
                and path.is_file() and not path.is_symlink(), "renderer state layout changed")


def _validate_state(state, descriptor, extra_files):
    require(type(state) is dict and set(state) == {"schema_version", "rendering_contract", "descriptor", "descriptor_sha256", "provenance"}
            and type(state["schema_version"]) is int and state["schema_version"] == 1 and state["rendering_contract"] == CONTRACT
            and canonical(state["descriptor"]) == canonical(descriptor) and state["descriptor_sha256"] == digest(descriptor), "renderer state identity changed")
    provenance = state["provenance"]
    require(type(provenance) is dict and set(provenance["source_files"]) == set(_dependencies(extra_files)), "renderer dependencies changed")
    verify(provenance, expected_design=RendererDesign(descriptor))
    dirty = provenance.get("uncommitted_dependency_files")
    require(type(dirty) is list and len(dirty) == len(set(dirty)) and set(dirty) <= set(provenance["source_files"]), "renderer dirty list changed")
    commit = provenance.get("code_commit")
    require(type(commit) is str and re.fullmatch(r"[0-9a-f]{40}", commit), "renderer commit missing")
    if descriptor["require_frozen"]:
        require(not dirty, "full E4 renderer must be committed")
        for name, wanted in provenance["source_files"].items():
            old = subprocess.run(["git", "show", f"{commit}:{name}"], cwd=ROOT, capture_output=True)
            require(old.returncode == 0 and digest_bytes(old.stdout) == wanted, "renderer freeze changed: "+name)
    return state


def _figure_layout(output, complete=False):
    directory = output/"figures"
    if not directory.exists():
        require(not complete, "figure directory missing"); return None
    require(directory.is_dir() and not directory.is_symlink(), "figure directory invalid")
    require(all(p.name in FILENAMES | {"figure_manifest.json"} and p.is_file() and not p.is_symlink() for p in directory.iterdir()), "noncanonical figure artifact")
    path = directory/"figure_manifest.json"
    if not path.exists():
        require(not complete, "figure manifest missing"); return None
    manifest = _json(path)
    require(type(manifest.get("figures")) is dict and set(manifest["figures"]) == FILENAMES
            and {p.name for p in directory.iterdir()} == FILENAMES | {path.name}, "figure manifest/file set incomplete")
    require(all(file_digest(directory/name) == wanted for name, wanted in manifest["figures"].items()), "saved figure bytes changed")
    return manifest


def _check_inputs(inputs, guard=None):
    if guard: guard.check("figures:input_precheck")
    require((inputs.output/"run_metadata.json").read_bytes() == inputs.root_bytes
            and _scientific_snapshot(inputs.output) == inputs.snapshot, "raw/table/audit bytes changed during rendering")
    if guard: guard.check("figures:input_postcheck")


def _manifest(inputs, state, figures):
    return {"schema_version": 1, "rendering_contract": CONTRACT, "design_sha256": inputs.design.digest,
            "input_binding": inputs.binding, "renderer": state, "renderer_state_sha256": digest(state), "figures": figures,
            "new_sources_or_fits": False, "visual_review_performed": False, "axes_per_figure": 1}


def _verify(inputs, extra_files, target=None):
    manifest = _figure_layout(inputs.output if target is None else target, True)
    state = _validate_state(manifest.get("renderer"), _descriptor(inputs), extra_files)
    require(canonical(manifest) == canonical(_manifest(inputs, state, manifest["figures"])), "figure manifest binding changed")
    _check_inputs(inputs)
    return manifest


def _destination(inputs, destination):
    if destination is None: return inputs.output
    target = execution._output_path(destination)
    require(target != inputs.output and not target.is_relative_to(inputs.output)
            and not inputs.output.is_relative_to(target), "external replot destination overlaps source archive")
    if target.exists():
        require(target.is_dir() and all(p.name in {"figures", ".execution.lock"} for p in target.iterdir()),
                "external replot destination contains unrelated evidence")
    return target


def verify_figures(output, design, *, extra_files=(), destination=None):
    """Complete read-only verification; destination optionally identifies a replot."""
    extra_files = _dependencies(extra_files)
    inputs = load_inputs(output, design)
    return _verify(inputs, extra_files, _destination(inputs, destination))


def write_figures(output, design, *, extra_files=(), destination=None):
    """Publish/adopt exact figures after whole-candidate validation, never overwrite."""
    extra_files = _dependencies(extra_files)
    inputs = load_inputs(output, design)
    target = _destination(inputs, destination)
    if _figure_layout(target) is not None: return _verify(inputs, extra_files, target)
    require(destination is not None or inputs.metadata["status"] != "completed", "closed evidence cannot recreate figures")
    execution._output_path(target)
    state_path = renderer_state_path(target)
    _state_layout(state_path)
    require(state_path.exists() or not (target/"figures").exists() or not any((target/"figures").iterdir()),
            "partial figures lack original renderer state")
    initial = execution._ManifestGuard(inputs.output/"run_metadata.json", inputs.root_bytes)
    with run_lock(target):
        initial.check("figures:after_lock")
        inputs = load_inputs(inputs.output, design)
        _destination(inputs, destination)
        if _figure_layout(target) is not None: return _verify(inputs, extra_files, target)
        descriptor = _descriptor(inputs)
        _state_layout(state_path)
        state = (_json(state_path) if state_path.exists() else
                 {"schema_version": 1, "rendering_contract": CONTRACT, "descriptor": descriptor, "descriptor_sha256": digest(descriptor),
                  "provenance": capture(RendererDesign(descriptor), _dependencies(extra_files), require_frozen=descriptor["require_frozen"])})
        _validate_state(state, descriptor, extra_files)
        guard = execution._ManifestGuard(inputs.output/"run_metadata.json", inputs.root_bytes)
        figures = {}
        try:
            _check_inputs(inputs, guard)
            figures = build_figures(inputs.tables, design)
            candidates = {}
            for stem in STEMS:
                pair = figure_bytes(figures[stem], stem)
                require(set(pair) == {stem+".png", stem+".pdf"} and all(type(v) is bytes and v for v in pair.values()), "renderer pair incomplete")
                candidates.update(pair)
                plt.close(figures[stem])
            require(set(candidates) == FILENAMES, "figure plan changed")
            _check_inputs(inputs, guard); _destination(inputs, destination); _figure_layout(target)
            # Check ALL existing partial bytes before writing ANY missing output.
            for name, payload in candidates.items():
                path = target/"figures"/name
                require(not path.exists() or path.read_bytes() == payload, "different immutable partial figure: "+name)
            _validate_state(state, descriptor, extra_files)
            atomic_json(state_path, state, immutable=True)
            state_bytes = state_path.read_bytes()
            require(state_bytes == (canonical(state)+"\n").encode(), "renderer state changed at creation")
            for name, payload in candidates.items():
                guard.check("figures:before_write:"+name)
                require(state_path.read_bytes() == state_bytes, "renderer state changed")
                atomic_bytes(target/"figures"/name, payload, immutable=True)
            _check_inputs(inputs, guard)
            require(state_path.read_bytes() == state_bytes, "renderer state changed before final manifest")
            hashes = {name: digest_bytes(payload) for name, payload in candidates.items()}
            require(all(file_digest(target/"figures"/name) == value for name, value in hashes.items()), "figure bytes changed before manifest")
            manifest = _manifest(inputs, state, hashes)
            atomic_json(target/"figures/figure_manifest.json", manifest, immutable=True)
            result = _verify(load_inputs(inputs.output, design), extra_files, target)
            guard.check("figures:complete")
            return result
        except BaseException:
            if not guard.drifted: guard.check("figures:failure")
            raise
        finally:
            for fig in figures.values(): plt.close(fig)
