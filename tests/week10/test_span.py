"""Independent arithmetic and preservation checks for span standardization."""
import gzip
import json

import numpy as np
import pandas as pd
import pytest

from ou_irregular.experiments.analysis import csv_bytes
from ou_irregular.experiments.provenance import file_digest
from ou_irregular.week10 import span
from ou_irregular.week10.span import (COMPACT_COLUMNS, ESTIMATORS, analyze_spans,
                                     stratum_indices, verify_span_analysis, write_span_analysis)


def records():
    rows = []
    # MMPP has three short and two long paths; IID has two short and three long.
    # Errors are identical within a stratum but overall mixtures differ.
    errors = {"mmpp": [1., 2., 3., 4., 6.], "iid": [1., 3., 3., 5., 7.]}
    spans = {"mmpp": [100., 150., 200., 600., 700.], "iid": [100., 200., 600., 700., 800.]}
    for arm in ("mmpp", "iid"):
        for rep, (span, error) in enumerate(zip(spans[arm], errors[arm])):
            for estimator in ESTIMATORS:
                row = {"cell_id": arm + "_cell", "rep": rep, "estimator": estimator,
                       "mechanism": arm, "ratio": 5., "n": 1000, "realized_span": span,
                       "point_valid": True, "interval_valid": True,
                       "point_status": "valid_interior", "interval_status": "valid",
                       "pipeline_error": "", "source_valid": True, "source_seed": 10,
                       "source_spawn_key": f"{arm}-{rep}", "times_sha256": f"time-{arm}-{rep}",
                       "observations_sha256": f"data-{arm}-{rep}"}
                for parameter, truth in (("theta", 1.), ("mu", 0.), ("sigma", .5)):
                    row.update({parameter: truth + error, "true_" + parameter: truth,
                                parameter + "_se": 1., parameter + "_lower": truth + error - 4,
                                parameter + "_upper": truth + error + 4})
                rows.append(row)
    return pd.DataFrame(rows, columns=COMPACT_COLUMNS)


def result(tables, name):
    frame = tables[name]
    return frame[(frame.estimator == "pfml") & (frame.parameter == "mu")].iloc[0]


def test_strata_boundary_convention_and_invalid_spans():
    assert stratum_indices([.1, 399.999, 400, 450, 475, 500, 525, 550, 600, 1000]).tolist() == [0, 0, 1, 2, 3, 4, 5, 6, 7, 7]
    for values in ([0], [-1], [np.inf], [np.nan]):
        with pytest.raises(ValueError):
            stratum_indices(values)
    with pytest.raises(ValueError):
        stratum_indices([10], edges=[0, 10, 10, np.inf])


def test_standardization_and_mc_variance_against_independent_hand_arithmetic():
    tables = analyze_spans(records(), edges=[0, 500, np.inf], minimum_count=2)
    row = result(tables, "span_standardized.csv")
    # Pooled stratum counts are 5 and 5, so the same 1/2 weights apply to both.
    assert row.mmpp_bias == pytest.approx(.5 * 2 + .5 * 5)
    assert row.iid_bias == pytest.approx(.5 * 2 + .5 * 5)
    assert row.bias_difference == pytest.approx(0)
    expected_mmpp_mse = .5 * (1 + 4 + 9) / 3 + .5 * (16 + 36) / 2
    expected_iid_mse = .5 * (1 + 9) / 2 + .5 * (9 + 25 + 49) / 3
    assert row.mmpp_mse == pytest.approx(expected_mmpp_mse)
    assert row.iid_mse == pytest.approx(expected_iid_mse)
    assert row.rmse_ratio == pytest.approx(np.sqrt(expected_mmpp_mse / expected_iid_mse))
    independent_mmpp_var = .25 * (np.var([1, 2, 3], ddof=1) / 3 + np.var([4, 6], ddof=1) / 2)
    independent_iid_var = .25 * (np.var([1, 3], ddof=1) / 2 + np.var([3, 5, 7], ddof=1) / 3)
    assert row.bias_difference_mcse ** 2 == pytest.approx(independent_mmpp_var + independent_iid_var)
    assert row.mmpp_coverage_operational == pytest.approx(.75)
    assert row.iid_coverage_operational == pytest.approx(2 / 3)
    assert row.mmpp_excluded_fraction == 0
    assert row.iid_excluded_fraction == 0
    raw = result(tables, "span_unadjusted.csv")
    assert raw.bias_difference == pytest.approx(3.2 - 3.8)


def test_failure_eligibility_and_outcome_independent_weights():
    baseline = records()
    modified = baseline.copy()
    invalid = (modified.mechanism == "mmpp") & (modified.rep == 0)
    modified.loc[invalid, ["point_valid", "interval_valid"]] = False
    modified.loc[invalid, "point_status"] = "invalid"
    modified.loc[invalid, "interval_status"] = "unavailable"
    for p in ("theta", "mu", "sigma"):
        modified.loc[invalid, [p, p + "_se", p + "_lower", p + "_upper"]] = np.nan
    a = analyze_spans(baseline, edges=[0, 500, np.inf], minimum_count=2)
    b = analyze_spans(modified, edges=[0, 500, np.inf], minimum_count=2)
    pd.testing.assert_frame_equal(a["span_support.csv"], b["span_support.csv"])
    row = result(b, "span_standardized.csv")
    assert row.mmpp_n_planned == 5
    assert row.mmpp_n_point_valid == 4
    assert row.mmpp_n_interval_valid == 4
    assert row.mmpp_bias == pytest.approx(.5 * 2.5 + .5 * 5)
    assert row.mmpp_coverage_operational == pytest.approx(.5 * (2 / 3) + .5 * .5)


def test_missing_common_support_retains_all_bins_and_excluded_mass():
    tables = analyze_spans(records(), edges=[0, 500, np.inf], minimum_count=3)
    row = result(tables, "span_standardized.csv")
    assert row.support_strata == 0
    assert row.mmpp_excluded_fraction == 1
    assert row.iid_excluded_fraction == 1
    assert np.isnan(row.rmse_ratio)
    assert len(tables["span_support.csv"]) == 2
    assert len(tables["span_strata.csv"]) == 2 * 2 * 3 * 3


@pytest.mark.parametrize("valid_points", [0, 1])
def test_valid_point_metrics_remain_undefined_or_uncertain_when_support_has_no_information(valid_points):
    data = records()
    missing = (data.mechanism == "mmpp") & (data.rep < 3 - valid_points)
    data.loc[missing, ["point_valid", "interval_valid"]] = False
    data.loc[missing, "point_status"] = "invalid"
    data.loc[missing, "interval_status"] = "unavailable"
    for p in ("theta", "mu", "sigma"):
        data.loc[missing, [p, p + "_se", p + "_lower", p + "_upper"]] = np.nan
    tables = analyze_spans(data, edges=[0, 500, np.inf], minimum_count=2)
    row = result(tables, "span_standardized.csv")
    assert row.support_strata == 2
    assert row.mmpp_n_planned == 5
    assert row.mmpp_n_point_valid == 2 + valid_points
    assert row.mmpp_coverage_operational == pytest.approx(.5 * (valid_points / 3) + .25)
    assert np.isnan(row.bias_difference_mcse)
    assert np.isnan(row.rmse_ratio_mc_lower)
    if valid_points == 0:
        assert np.isnan(row.mmpp_bias)
        assert np.isnan(row.mmpp_mse)
        assert np.isnan(row.rmse_ratio)
    else:
        assert np.isfinite(row.mmpp_bias)
        assert np.isfinite(row.rmse_ratio)


def test_empty_bins_reported_and_boundary_coverage_has_wilson_uncertainty():
    data = records()
    for p in ("theta", "mu", "sigma"):
        data[p + "_lower"] = data["true_" + p] + 10
        data[p + "_upper"] = data["true_" + p] + 11
    tables = analyze_spans(data, edges=[0, 500, 1000, np.inf], minimum_count=2)
    strata = tables["span_strata.csv"]
    nonempty = strata[strata.n_planned > 0]
    assert nonempty.coverage_operational.eq(0).all()
    assert nonempty.coverage_wilson_upper.gt(0).all()
    assert strata[strata.stratum == 2].n_planned.eq(0).all()
    standardized = result(tables, "span_standardized.csv")
    assert standardized.coverage_difference_mcse == 0
    assert np.isnan(standardized.coverage_difference_mc_lower)
    assert np.isnan(standardized.coverage_difference_mc_upper)


@pytest.mark.parametrize("kind", ["duplicate", "missing_estimator", "mismatched_source", "pipeline", "flags", "truth", "bounds"])
def test_rejects_corrupted_fit_or_pairing_records(kind):
    data = records()
    if kind == "duplicate":
        data = pd.concat([data, data.iloc[[0]]])
    elif kind == "missing_estimator":
        data = data.iloc[1:]
    elif kind == "mismatched_source":
        data.loc[0, "realized_span"] = 111
    elif kind == "pipeline":
        data.loc[0, "pipeline_error"] = "retained failure"
    elif kind == "flags":
        data["point_valid"] = data.point_valid.astype(object)
        data.loc[0, "point_valid"] = 1
    elif kind == "truth":
        data.loc[0, "true_mu"] = .1
    else:
        data.loc[0, "mu_lower"] = 100
    with pytest.raises(ValueError):
        analyze_spans(data, edges=[0, 500, np.inf], minimum_count=2)


def test_writer_refuses_overlap_or_existing_destination(tmp_path):
    root = tmp_path / "original"
    root.mkdir()
    with pytest.raises(ValueError):
        write_span_analysis(root, root / "new")
    with pytest.raises(ValueError):
        write_span_analysis(root, tmp_path)
    other = tmp_path / "existing"
    other.mkdir()
    with pytest.raises(FileExistsError):
        write_span_analysis(root, other)


def verification_fixture(tmp_path, monkeypatch):
    root, output = tmp_path / "source", tmp_path / "output"
    root.mkdir()
    output.mkdir()
    source = root / "source.json"
    source.write_text("original source identity")
    data = records()
    binding = {"contract": span.CONTRACT, "source_design": {"synthetic_fixture": True},
               "fit_records": len(data), "source_paths": 10,
               "source_files": {"source.json": {"sha256": file_digest(source), "bytes": source.stat().st_size}}}
    # The raw reader is tested by full archive validation; this fixture isolates
    # read-only export guards without simulating or loading a research archive.
    monkeypatch.setattr(span, "read_e4_records", lambda _: (data, binding))
    tables = analyze_spans(data)
    tables["span_fit_records.csv.gz"] = data
    outputs = {}
    for name, table in tables.items():
        payload = csv_bytes(table)
        if name.endswith(".gz"):
            payload = gzip.compress(payload, mtime=0)
        (output / name).write_bytes(payload)
        outputs[name] = {"sha256": file_digest(output / name), "rows": len(table), "bytes": len(payload)}
    manifest = {**binding, "span_edges": [0., 400., 450., 475., 500., 525., 550., 600., "infinity"],
                "stratum_convention": "left_closed_right_open", "minimum_count_each_arm": 30,
                "weighting": "pooled_counts_normalized_on_common_support",
                "uncertainty": "independent_arm_delta_method_conditional_on_observed_stratum_counts_and_weights",
                "zero_coverage_variance": "MCSE_zero_but_normal_confidence_bounds_undefined;_stratum_Wilson_bounds_retained",
                "scope": "post_result_descriptive_sensitivity_not_causal_adjustment",
                "original_unadjusted_arithmetic_checks": 216,
                "new_simulated_sources": 0, "new_fits": 0, "outputs": outputs}
    (output / "span_manifest.json").write_text(json.dumps(manifest))
    return root, output, manifest


def test_readonly_verifier_rederives_every_export(tmp_path, monkeypatch):
    root, output, _ = verification_fixture(tmp_path, monkeypatch)
    before = {p.name: p.read_bytes() for p in output.iterdir()}
    result = verify_span_analysis(root, output)
    assert result["rederived_exports"] == 6
    assert result["new_simulated_sources"] == result["new_fits"] == 0
    assert before == {p.name: p.read_bytes() for p in output.iterdir()}


@pytest.mark.parametrize("kind", ["source", "settings", "missing", "extra", "self_consistent_edit"])
def test_verifier_rejects_changed_input_settings_inventory_or_arithmetic(tmp_path, monkeypatch, kind):
    root, output, manifest = verification_fixture(tmp_path, monkeypatch)
    if kind == "source":
        (root / "source.json").write_text("changed source identity")
    elif kind == "settings":
        manifest["minimum_count_each_arm"] = 29
    elif kind == "missing":
        (output / "span_support.csv").unlink()
        del manifest["outputs"]["span_support.csv"]
    elif kind == "extra":
        (output / "unplanned.txt").write_text("extra output")
    else:
        target = output / "span_unadjusted.csv"
        table = pd.read_csv(target)
        table.loc[0, "bias_difference"] += 1
        target.write_bytes(csv_bytes(table))
        manifest["outputs"][target.name].update(sha256=file_digest(target), bytes=target.stat().st_size)
    (output / "span_manifest.json").write_text(json.dumps(manifest))
    with pytest.raises(ValueError):
        verify_span_analysis(root, output)
