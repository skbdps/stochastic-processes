import json
import shutil
from dataclasses import replace

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from ou_irregular.experiments.provenance import file_digest
from ou_irregular.experiments.workflow import verify_figures
from ou_irregular.week6 import plots
from ou_irregular.week6.reporting import INPUTS, load_tables
from test_reporting import full_grid_tables


def copy_inputs(source, destination):
    for name in (*INPUTS, "run_metadata.json", "table_manifest.json"):
        shutil.copyfile(source/name, destination/name)


def replace_summary(root, table):
    table.to_csv(root/"parameter_summary.csv", index=False, float_format="%.17g")
    path = root/"table_manifest.json"
    manifest = json.loads(path.read_text())
    manifest["parameter_summary.csv"] = {"rows": len(table), "sha256": file_digest(root/"parameter_summary.csv")}
    path.write_text(json.dumps(manifest))


def test_full_grid_outputs_orientation_and_exact_regeneration(tmp_path, full_grid_tables, monkeypatch):
    root, design = full_grid_tables
    tables, inputs = load_tables(root, design)
    summary = tables["parameter_summary.csv"]
    captured = {}
    original_save = plots._save

    def capture(fig, output, stem, title):
        panels = [ax for ax in fig.axes if ax.images]
        captured[stem] = {
            "values": [ax.images[0].get_array().copy() for ax in panels],
            "limits": [(ax.images[0].norm.vmin, ax.images[0].norm.vmax) for ax in panels],
            "xticks": [[t.get_text() for t in ax.get_xticklabels()] for ax in panels],
            "yticks": [[t.get_text() for t in ax.get_yticklabels()] for ax in panels],
            "text": "\n".join(t.get_text() for t in fig.texts),
        }
        original_save(fig, output, stem, title)

    monkeypatch.setattr(plots, "_save", capture)
    report = plots.render_figures(root, design)
    stems = {f"{prefix}_n{n}" for n in (250, 1000, 4000)
             for prefix in ("F5_operational", "F6_conditional", "F6_availability")}
    expected = {f"{stem}.{extension}" for stem in stems for extension in ("png", "pdf")} | set(plots.REPORT_FILES)
    assert set(report["figures"]) == expected
    assert len(report["figures"]) == 22
    assert set(captured) == stems
    assert report["input_table_sha256"] == inputs
    assert verify_figures(root, design) == report
    assert set(p.name for p in (root/"figures").iterdir()) == expected | {"figure_manifest.json"}
    for stem, figure in captured.items():
        assert "VALIDATION ONLY" in figure["text"] and "R = 2 per" in figure["text"]
        assert "REGISTERED E2" not in figure["text"]
        assert all(value.shape == (4, 5) for value in figure["values"])
        assert all(limits == (0, 100) for limits in figure["limits"])
        assert all(ticks == ["0", "0.5", "1", "2", "3"] for ticks in figure["xticks"])
        assert all(ticks == ["0.1", "0.5", "1", "2"] for ticks in figure["yticks"])
        n = int(stem.rsplit("n", 1)[1])
        if "availability" in stem:
            assert len(figure["values"]) == 6
            specifications = [(kind+"_availability", "theta", method)
                              for kind in ("point", "interval") for method in ("exact", "pfml", "euler")]
        else:
            assert len(figure["values"]) == 9
            scope = "operational" if stem.startswith("F5") else "conditional"
            specifications = [("coverage_"+scope, parameter, method)
                              for parameter in ("theta", "mu", "sigma") for method in ("exact", "pfml", "euler")]
        # Check every displayed entry against its scientific key, including
        # opposite corners. This detects a transpose even when shape is 4x5.
        for values, (metric, parameter, method) in zip(figure["values"], specifications):
            group = summary[(summary.n == n) & (summary.parameter == parameter) & (summary.estimator == method)]
            indexed = {(row.theta_mean_gap, row.cv): getattr(row, metric) for row in group.itertuples()}
            expected_values = [[100*indexed[(q, cv)] for cv in (0, .5, 1, 2, 3)] for q in (.1, .5, 1, 2)]
            np.testing.assert_allclose(values.filled(np.nan), expected_values, rtol=0, atol=0, equal_nan=True)
    report_manifest = json.loads((root/"figures/report_manifest.json").read_text())
    assert report_manifest["input_table_sha256"] == inputs
    for name in plots.REPORT_FILES:
        assert report["figures"][name] == file_digest(root/"figures"/name)
    first = {p.name: p.read_bytes() for p in (root/"figures").iterdir()}
    assert all(b"/CreationDate" not in data and b"/ModDate" not in data
               for name, data in first.items() if name.endswith(".pdf"))
    assert plots.render_figures(root, design) == report
    assert {p.name: p.read_bytes() for p in (root/"figures").iterdir()} == first
    # The complete figure manifest must bind report bytes as well as PNG/PDFs.
    # Corrupt a separate copy, preserving the validation evidence used above.
    copy_inputs(root, tmp_path)
    shutil.copytree(root/"figures", tmp_path/"figures")
    assert len(verify_figures(tmp_path, design)["figures"]) == 22
    with (tmp_path/"figures/classification_counts.csv").open("a") as stream:
        stream.write("\n")
    with pytest.raises(ValueError, match="figure digest"):
        verify_figures(tmp_path, design)


def test_conditional_unavailable_is_gray_na_and_observed_zero_is_colored(full_grid_tables):
    root, design = full_grid_tables
    tables, _ = load_tables(root, design)
    summary = tables["parameter_summary.csv"].copy()
    ns, qs, cvs = plots._layout(design)
    base = (summary.n == ns[0]) & (summary.estimator == "exact") & (summary.parameter == "theta")
    unavailable = base & (summary.theta_mean_gap == qs[0]) & (summary.cv == cvs[0])
    zero = base & (summary.theta_mean_gap == qs[0]) & (summary.cv == cvs[1])
    # Deliberately constructed display fixture, not a scientific rerun.
    summary.loc[unavailable, ["n_interval_valid", "n_covered", "coverage_operational"]] = 0
    summary.loc[unavailable, "coverage_conditional"] = np.nan
    summary.loc[zero, ["n_point_valid", "n_interval_valid"]] = design.replications
    summary.loc[zero, ["n_covered", "coverage_operational", "coverage_conditional"]] = 0
    values, annotations = plots._coverage_matrix(summary, ns[0], "exact", "theta", qs, cvs, "conditional")
    assert values.mask[0, 0] and not values.mask[0, 1]
    assert annotations[0, 0] == "NA\nS=0"
    assert annotations[0, 1] == "0.00%\nS=2" and values[0, 1] == 0
    operational, operational_text = plots._coverage_matrix(summary, ns[0], "exact", "theta", qs, cvs, "operational")
    assert not np.ma.getmaskarray(operational)[0, 0]
    assert operational[0, 0] == 0 and operational_text[0, 0].startswith("0.00%")
    fig, ax = plt.subplots()
    try:
        image = plots._draw_panel(ax, values, annotations, qs, cvs, "Synthetic denominator test")
        colors = image.to_rgba(values)
        np.testing.assert_array_equal(colors[0, 0], image.cmap.get_bad())
        assert not np.array_equal(colors[0, 0], colors[0, 1])
        assert {"NA\nS=0", "0.00%\nS=2"} <= {t.get_text() for t in ax.texts}
    finally:
        plt.close(fig)


@pytest.mark.parametrize("invalid", ["nan_with_valid_intervals", "zero_without_intervals", "out_of_range"])
def test_conditional_denominator_mismatches_are_rejected(full_grid_tables, invalid):
    root, design = full_grid_tables
    tables, _ = load_tables(root, design)
    summary = tables["parameter_summary.csv"].copy()
    ns, qs, cvs = plots._layout(design)
    mask = ((summary.n == ns[0]) & (summary.estimator == "exact") & (summary.parameter == "theta")
            & (summary.theta_mean_gap == qs[0]) & (summary.cv == cvs[0]))
    if invalid == "zero_without_intervals":
        summary.loc[mask, ["n_interval_valid", "n_covered", "coverage_conditional"]] = 0
    elif invalid == "nan_with_valid_intervals":
        summary.loc[mask, "coverage_conditional"] = np.nan
    else:
        summary.loc[mask, "coverage_conditional"] = 1.1
    with pytest.raises(ValueError, match="coverage does not match"):
        plots._coverage_matrix(summary, ns[0], "exact", "theta", qs, cvs, "conditional")


@pytest.mark.parametrize("corruption", ["missing_group", "duplicate_with_same_length"])
def test_plot_inputs_require_full_keys_before_writing(tmp_path, full_grid_tables, corruption):
    root, design = full_grid_tables
    copy_inputs(root, tmp_path)
    summary = pd.read_csv(tmp_path/"parameter_summary.csv", float_precision="round_trip")
    if corruption == "missing_group":
        summary = summary.iloc[9:]
    else:
        summary = pd.concat([summary.iloc[1:], summary.iloc[[1]]], ignore_index=True)
    replace_summary(tmp_path, summary)
    with pytest.raises(ValueError, match="full planned key set"):
        plots.render_figures(tmp_path, design)
    assert not (tmp_path/"figures").exists()


def test_plotter_rejects_a_shrunken_design_before_io(tmp_path, full_grid_tables):
    _, design = full_grid_tables
    partial = replace(design, cells=design.cells[:-1])
    with pytest.raises(ValueError, match="complete registered E2"):
        plots.render_figures(tmp_path, partial)
    assert not list(tmp_path.iterdir())


def test_development_run_cannot_acquire_a_registered_label(full_grid_tables):
    _, design = full_grid_tables
    with pytest.raises(ValueError, match="registered E2 identity"):
        plots._layout(replace(design, scope="headline"))


def test_changed_source_table_bytes_are_rejected(tmp_path, full_grid_tables):
    root, design = full_grid_tables
    copy_inputs(root, tmp_path)
    with (tmp_path/"parameter_summary.csv").open("a") as stream:
        stream.write("\n")
    with pytest.raises(ValueError, match="input digest mismatch"):
        plots.render_figures(tmp_path, design)


def test_classification_is_copied_without_rethresholding_coverage(full_grid_tables):
    root, design = full_grid_tables
    tables, _ = load_tables(root, design)
    summary = tables["parameter_summary.csv"].copy()
    ns, qs, cvs = plots._layout(design)
    mask = ((summary.n == ns[0]) & (summary.estimator == "exact") & (summary.parameter == "theta")
            & (summary.theta_mean_gap == qs[0]) & (summary.cv == cvs[0]))
    summary.loc[mask, ["n_point_valid", "n_interval_valid", "n_covered"]] = design.replications
    summary.loc[mask, ["coverage_operational", "coverage_conditional"]] = 1
    summary.loc[mask, "practical_classification"] = "supported_breach"
    _, annotations = plots._coverage_matrix(summary, ns[0], "exact", "theta", qs, cvs, "operational")
    assert annotations[0, 0] == "100.00%\nB"
    _, mu_text = plots._coverage_matrix(summary, ns[0], "exact", "mu", qs, cvs, "operational")
    assert all(label.endswith("%") and "\n" not in label for label in mu_text.flat)
