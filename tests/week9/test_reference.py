"""Independent arithmetic and artist checks against original saved experiment bytes.

References use stdlib CSV/gzip/statistics and the written statistical definitions;
production classification, selectors, contrast and ECDF helpers are not used to
construct the expected values. No source simulation, fitting or replay is run.
"""
from collections import Counter, defaultdict
import csv
from decimal import Decimal
import gzip
import json
import math
from pathlib import Path
import statistics

import matplotlib.pyplot as plt
import numpy as np
import pytest

from ou_irregular.week9.inputs import load_sources
from ou_irregular.week9 import plots, tables

ROOT = Path(__file__).resolve().parents[2]
Z = 1.959963984540054
METHODS = ("exact", "pfml", "euler")


def saved_csv(study, name="parameter_summary.csv"):
    with (ROOT / "artifacts" / study / name).open(newline="") as stream:
        return list(csv.DictReader(stream))


def raw_records(study, cells=None):
    root = ROOT / "artifacts" / study / "cells"
    for cell_dir in sorted(root.iterdir()):
        if cells is not None and cell_dir.name not in cells:
            continue
        with (cell_dir / "cell_manifest.json").open() as stream:
            manifest = json.load(stream)
        for part in manifest["partitions"]:
            with gzip.open(cell_dir / part["filename"], "rt") as stream:
                for line in stream:
                    if line.strip():
                        yield json.loads(line)


@pytest.fixture(scope="module")
def publication():
    return tables.build_tables(load_sources())


@pytest.fixture(scope="module")
def e1_reference():
    values = defaultdict(list)
    truth = {"theta": 1., "mu": 0., "sigma": .5}
    for row in raw_records("week5_e1"):
        assert row["point_valid"] is True
        for parameter in truth:
            values[(row["cell_id"], row["estimator"], parameter)].append(
                (row[parameter], row["rep"], row["cv"], truth[parameter]))
    return {key: sorted(rows, key=lambda r: (r[0], r[1])) for key, rows in values.items()}


def independent_class(row):
    """Apply the registered inequalities to exact decimal CSV endpoints."""
    if row["parameter"] == "mu":
        return "not_applicable_zero_truth", "not_applicable"
    lo = Decimal(row["operational_wilson_lower"])
    hi = Decimal(row["operational_wilson_upper"])
    lower_bias = Decimal(row["relative_bias_mc_lower"])
    upper_bias = Decimal(row["relative_bias_mc_upper"])
    bad_coverage = hi < Decimal("0.8")
    bad_bias = lower_bias > Decimal("0.1") or upper_bias < Decimal("-0.1")
    if bad_coverage or bad_bias:
        return "supported_breach", ("bias_and_coverage" if bad_coverage and bad_bias else
                                    "coverage" if bad_coverage else "bias")
    if lo >= Decimal("0.8") and lower_bias >= Decimal("-0.1") and upper_bias <= Decimal("0.1"):
        return "within_criteria", "none"
    return "unresolved", "uncertain"


def wilson(successes, count):
    p = successes/count
    denominator = 1 + Z*Z/count
    center = (p + Z*Z/(2*count))/denominator
    radius = Z*math.sqrt(p*(1-p)/count + Z*Z/(4*count*count))/denominator
    return max(0., center-radius), min(1., center+radius)


def test_all_402_headline_classifications_and_drivers_match_independent_arithmetic(publication):
    expected_totals = {"e1": Counter(within_criteria=19, supported_breach=22, unresolved=1),
                       "e2": Counter(within_criteria=160, supported_breach=195, unresolved=5)}
    for study, directory in (("e1", "week5_e1"), ("e2", "week6_e2")):
        originals = [r for r in saved_csv(directory) if r["parameter"] != "mu"]
        view = publication["T1_"+study+".csv"].set_index(["cell_id", "estimator", "parameter"])
        assert len(view) == len(originals) == (42 if study == "e1" else 360)
        counts = Counter()
        for original in originals:
            expected_class, expected_driver = independent_class(original)
            key = original["cell_id"], original["estimator"], original["parameter"]
            result = view.loc[key]
            assert original["practical_classification"] == result.practical_classification == expected_class
            assert result.breach_driver == expected_driver
            successes, planned = int(original["n_covered"]), int(original["n_planned"])
            lower, upper = wilson(successes, planned)
            assert float(original["operational_wilson_lower"]) == pytest.approx(lower, abs=3e-15)
            assert float(original["operational_wilson_upper"]) == pytest.approx(upper, abs=3e-15)
            assert result.coverage_operational == successes/planned
            counts[expected_class] += 1
        assert counts == expected_totals[study]


def test_all_six_unresolved_identities_remain_distinct_across_studies(publication):
    result = publication["unresolved.csv"]
    actual = set(result[["study", "n", "theta_mean_gap", "cv", "estimator", "parameter"]]
                 .itertuples(index=False, name=None))
    assert actual == {
        ("e1", 1000, .5, .5, "pfml", "sigma"),
        ("e2", 250, .1, 0., "euler", "theta"),
        ("e2", 250, .1, .5, "euler", "theta"),
        ("e2", 250, .1, 2., "euler", "sigma"),
        ("e2", 250, .1, 3., "pfml", "theta"),
        ("e2", 1000, .5, .5, "pfml", "sigma"),
    }
    assert len(result) == 6
    overlap = result[(result.n == 1000) & (result.estimator == "pfml") & (result.parameter == "sigma")]
    assert overlap.cell_id.nunique() == 1 and overlap.study.nunique() == 2
    assert overlap.coverage_operational.nunique() == 2


def test_failed_fit_export_accounts_for_85_points_one_rejection_and_46_sources(publication):
    with gzip.open(ROOT / "artifacts/week6_e2/failed_fits.csv.gz", "rt", newline="") as stream:
        failed = list(csv.DictReader(stream))
    assert len(failed) == 86
    assert sum(row["point_valid"] == "False" for row in failed) == 85
    assert sum(row["point_valid"] == "True" and row["interval_valid"] == "False" for row in failed) == 1
    assert len({(row["cell_id"], int(row["rep"])) for row in failed}) == 46
    per_method = Counter((row["cell_id"], row["estimator"]) for row in failed)
    view = publication["availability.csv"]
    e2 = view[view.study == "e2"]
    assert int(e2.n_point_invalid.sum()) == 85
    assert int(e2.n_point_valid_without_interval.sum()) == 1
    assert int(e2.n_interval_unavailable.sum()) == 86
    for row in e2.itertuples():
        assert row.n_interval_unavailable == per_method[(row.cell_id, row.estimator)]
        assert row.n_interval_unavailable == row.n_point_invalid + row.n_point_valid_without_interval


def test_ecdf_samples_ranks_and_truth_match_all_original_e1_saved_fits(publication, e1_reference):
    frame = publication["f1_ecdf.csv"]
    assert len(frame) == 126000
    assert set(frame.cv) == {0, .25, .5, 1, 1.5, 2, 3}
    groups = {key: group for key, group in frame.groupby(["cell_id", "estimator", "parameter"], sort=False)}
    assert set(groups) == set(e1_reference)
    for key, expected in e1_reference.items():
        group = groups[key].sort_values("rank")
        assert group.estimate.tolist() == [r[0] for r in expected]
        assert group.rep.tolist() == [r[1] for r in expected]
        assert group.cv.tolist() == [r[2] for r in expected]
        assert group.truth.tolist() == [r[3] for r in expected]
        assert group["rank"].tolist() == list(range(1, 2001))
        assert group.ecdf.tolist() == [rank/2000 for rank in range(1, 2001)]
        assert group.n_planned.eq(2000).all() and group.n_point_valid.eq(2000).all()


def test_mu_metrics_remain_absolute_and_unclassified(publication):
    summary = publication["parameter_summary.csv"]
    mu = summary[summary.parameter == "mu"]
    assert mu.truth.eq(0).all()
    for column in summary.columns:
        if column.startswith("relative_"):
            assert mu[column].isna().all(), column
    assert mu.practical_classification.eq("not_applicable_zero_truth").all()
    assert mu.breach_driver.eq("not_applicable").all()
    assert not (publication["T1_e2.csv"].parameter == "mu").any()


def test_classification_artist_labels_and_colors_encode_every_original_cell(publication):
    original = saved_csv("week6_e2")
    for n in (250, 1000, 4000):
        fig = plots.classification_figure(publication["T1_e2.csv"], n)
        try:
            texts = {text.get_gid(): text.get_text() for ax in fig.axes for text in ax.texts}
            expected = {}
            for row in original:
                if int(row["n"]) != n or row["parameter"] == "mu":
                    continue
                label, driver = independent_class(row)
                code = {"within_criteria": "W", "unresolved": "U"}.get(label)
                if code is None:
                    code = "B:"+{"bias": "b", "coverage": "c", "bias_and_coverage": "bc"}[driver]
                key = f"class:{n}:{row['estimator']}:{row['parameter']}:{float(row['theta_mean_gap']):g}:{float(row['cv']):g}"
                expected[key] = code
            assert texts == expected
            assert len(fig.axes) == 6 and len(texts) == 120
            for ax in fig.axes:
                image = ax.images[0]
                assert image.get_interpolation() == "nearest"
                assert image.get_array().shape == (4, 5)
                # Text-code and image-class correspondence is independent of palette color choice.
                annotations = {(round(t.get_position()[0]), round(t.get_position()[1])): t.get_text() for t in ax.texts}
                for y in range(4):
                    for x in range(5):
                        assert image.get_array()[y, x] == (0 if annotations[x,y] == "W" else 2 if annotations[x,y] == "U" else 1)
        finally:
            plt.close(fig)


def test_ecdf_artist_samples_and_truth_lines_match_independent_raw_e1(publication, e1_reference):
    frame = publication["f1_ecdf.csv"]
    by_coordinates = {(rows[0][2], key[1], key[2]): rows for key, rows in e1_reference.items()}
    for cvs in ((0., .25, .5, 1.), (1.5, 2., 3.)):
        fig = plots.ecdf_figure(frame, cvs, "Independent reference check")
        try:
            for ax in fig.axes:
                series = [line for line in ax.lines if line.get_gid()]
                assert len(series) == 3
                assert {line.get_label() for line in series} == {"Exact", "PFML", "Euler"}
                for line in series:
                    _, cv, parameter, method = line.get_gid().split(":")
                    expected = by_coordinates[(float(cv), method, parameter)]
                    assert list(line.get_xdata()) == [r[0] for r in expected]
                    assert list(line.get_ydata()) == [rank/2000 for rank in range(1, 2001)]
                    assert line.get_drawstyle() == "steps-post"
                truth = {"theta": 1., "mu": 0., "sigma": .5}[series[0].get_gid().split(":")[2]]
                truth_lines = [line for line in ax.lines if not line.get_gid()]
                assert len(truth_lines) == 1 and list(truth_lines[0].get_xdata()) == [truth, truth]
        finally:
            plt.close(fig)


@pytest.fixture(scope="module")
def e4_pfml_reference():
    design = json.loads((ROOT / "artifacts/week8_e4/resolved_design.json").read_text())
    cells = {c["cell_id"]: (int(c["ratio"]), c["mechanism"]) for c in design["cells"]
             if c["ratio"] in (5, 20) and c["mechanism"] in ("mmpp", "iid")}
    raw = defaultdict(list)
    for row in raw_records("week8_e4", set(cells)):
        if row["estimator"] == "pfml":
            assert row["point_valid"] is True and row["interval_valid"] is True
            raw[cells[row["cell_id"]]].append((row["mu"], row["mu_se"], row["mu_lower"] <= 0 <= row["mu_upper"]))
    result = {}
    for key, rows in raw.items():
        assert len(rows) == 2000
        errors = [r[0] for r in rows]
        squares = [v*v for v in errors]
        rmse = math.sqrt(statistics.mean(squares))
        result[key] = {"rmse": rmse, "rmse_mcse": statistics.stdev(squares)/(2*math.sqrt(2000)*rmse),
                       "coverage": sum(r[2] for r in rows)/2000,
                       "sd": statistics.stdev(errors),
                       "rms_se": math.sqrt(statistics.mean(r[1]*r[1] for r in rows))}
    return result


def test_e4_pfml_mu_contrasts_and_calibration_recompute_from_saved_raw_fits(publication, e4_pfml_reference):
    matched = publication["e4_matched.csv"]
    calibration = publication["calibration.csv"]
    for ratio in (5, 20):
        a, b = (e4_pfml_reference[(ratio, mechanism)] for mechanism in ("mmpp", "iid"))
        rmse_ratio = a["rmse"]/b["rmse"]
        se_log = math.hypot(a["rmse_mcse"]/a["rmse"], b["rmse_mcse"]/b["rmse"])
        cov_difference = a["coverage"]-b["coverage"]
        cov_se = math.sqrt((a["coverage"]*(1-a["coverage"])+b["coverage"]*(1-b["coverage"]))/2000)
        expected = {"rmse_ratio": rmse_ratio,
                    "rmse_ratio_mc_lower": rmse_ratio*math.exp(-Z*se_log),
                    "rmse_ratio_mc_upper": rmse_ratio*math.exp(Z*se_log),
                    "operational_coverage_difference": cov_difference,
                    "operational_coverage_difference_mc_lower": cov_difference-Z*cov_se,
                    "operational_coverage_difference_mc_upper": cov_difference+Z*cov_se}
        row = matched[(matched.ratio == ratio) & (matched.estimator == "pfml") & (matched.parameter == "mu")].iloc[0]
        for field, value in expected.items():
            assert row[field] == pytest.approx(value, rel=8e-13, abs=2e-15), field
        assert row.rmse_ratio_mc_lower > 1 and row.operational_coverage_difference_mc_upper < 0
        for mechanism, reference in (("mmpp", a), ("iid", b)):
            row = calibration[(calibration.study == "e4") & (calibration.ratio == ratio)
                              & (calibration.mechanism == mechanism) & (calibration.estimator == "pfml")
                              & (calibration.parameter == "mu")].iloc[0]
            assert row.empirical_sd_interval_valid == pytest.approx(reference["sd"], rel=8e-13)
            assert row.rms_reported_se == pytest.approx(reference["rms_se"], rel=8e-13)
            assert row.sd_over_rms_se == pytest.approx(reference["sd"]/reference["rms_se"], rel=8e-13)
            assert row.n_interval_valid == 2000


def test_matched_plot_uses_source_values_and_asymmetric_interval_endpoints(publication):
    original = saved_csv("week8_e4", "matched_iid_comparisons.csv")
    fig = plots.matched_figure(publication["e4_matched.csv"])
    try:
        for axis, metric, scale in zip(fig.axes, ("rmse_ratio", "operational_coverage_difference"), (1., 100.)):
            containers = {c.lines[0].get_gid(): c for c in axis.containers}
            for index, method in enumerate(METHODS):
                rows = sorted((r for r in original if r["estimator"] == method and r["parameter"] == "mu"), key=lambda r: float(r["ratio"]))
                series = containers[f"matched:mu:{metric}:{method}"]
                values = np.array([float(r[metric])*scale for r in rows])
                np.testing.assert_array_equal(series.lines[0].get_ydata(), values)
                segments = series.lines[2][0].get_segments()
                for i, (segment, row) in enumerate(zip(segments, rows)):
                    assert segment[0,0] == pytest.approx(i+(index-1)*.14)
                    assert segment[1,0] == segment[0,0]
                    assert segment[0,1] == pytest.approx(float(row[metric+"_mc_lower"])*scale, abs=2e-15)
                    assert segment[1,1] == pytest.approx(float(row[metric+"_mc_upper"])*scale, abs=2e-15)
                if metric == "rmse_ratio" and method == "pfml":
                    for row in rows[1:]:
                        assert not math.isclose(float(row[metric])-float(row[metric+"_mc_lower"]),
                                                float(row[metric+"_mc_upper"])-float(row[metric]), rel_tol=1e-6)
    finally:
        plt.close(fig)


def test_calibration_artist_uses_rms_not_mean_se_and_identical_axis_scale(publication):
    original = saved_csv("week8_e4")
    fig = plots.calibration_e4(publication["calibration.csv"])
    try:
        order = [(ratio, mechanism) for ratio in (1, 5, 20) for mechanism in ("mmpp", "iid")]
        for axis, metric in zip(fig.axes, ("empirical_sd_interval_valid", "rms_reported_se")):
            lines = {line.get_gid(): line for line in axis.lines}
            for method in METHODS:
                expected = []
                for ratio, mechanism in order:
                    row, = [r for r in original if r["estimator"] == method and r["parameter"] == "mu"
                            and float(r["ratio"]) == ratio and r["mechanism"] == mechanism]
                    expected.append(float(row[metric]))
                assert list(lines[f"calibration:e4:mu:{metric}:{method}"].get_ydata()) == expected
        assert fig.axes[0].get_ylim() == fig.axes[1].get_ylim()
        pfml_mmpp = next(r for r in original if r["estimator"] == "pfml" and r["parameter"] == "mu"
                         and r["mechanism"] == "mmpp" and float(r["ratio"]) == 20)
        assert float(pfml_mmpp["mean_reported_se"]) != float(pfml_mmpp["rms_reported_se"])
    finally:
        plt.close(fig)
