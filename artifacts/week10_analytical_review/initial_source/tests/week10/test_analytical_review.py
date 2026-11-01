"""Independent checks of conditional-mixture and paired coverage arithmetic."""
import gzip
import json
import math

import numpy as np
import pandas as pd
import pytest

from ou_irregular.experiments.analysis import csv_bytes
from ou_irregular.experiments.provenance import file_digest
from ou_irregular.week10.analytical_review import (
    Z95, build, diagnostic_tables, known_width_coverage, paired_difference, verify,
)


def fixture_tables():
    paths, summaries = [], []
    truth = 2.5
    fixed_errors = np.array([0., 2., 3., -1.])
    fitted_errors = np.array([0., 1.5, 3., -2.])
    fitted_se = np.array([1., 1., 2., 2.])
    for mechanism, scale in (("mmpp", 1.), ("iid", .5)):
        cid = mechanism + "-5"
        for rep in range(4):
            paths.append({"cell_id": cid, "rep": rep, "mechanism": mechanism, "ratio": 5.,
                          "point_valid": True, "interval_valid": True,
                          "fixed_rho_mu": truth + fixed_errors[rep],
                          "pfml_mu": truth + fitted_errors[rep], "pfml_mu_se": fitted_se[rep],
                          "fixed_rho_conditional_variance": scale*(rep+1)**2,
                          "population_working_variance": 1.,
                          "times_sha256": cid + "-time-" + str(rep),
                          "observations_sha256": cid + "-values-" + str(rep)})
        summaries.append({"cell_id": cid, "estimator": "pfml", "parameter": "mu", "truth": truth,
                          "n_planned": 4, "n_point_valid": 4, "n_interval_valid": 4,
                          "n_covered": 4, "coverage_operational": 1.,
                          "bias": fitted_errors.mean(), "rmse": np.sqrt(np.mean(fitted_errors**2)),
                          "empirical_sd_point_valid": fitted_errors.std(ddof=1)})
    return pd.DataFrame(paths), pd.DataFrame(summaries)


def test_known_gaussian_width_has_exact_nominal_conditional_coverage():
    variance = np.array([.1, 1., 4., 100.])
    assert known_width_coverage(variance, variance) == pytest.approx(np.full(4, .95))
    # Normal central probability checked with the independent erf identity.
    assert float(known_width_coverage(4., 1.)) == pytest.approx(math.erf(Z95/(2*math.sqrt(2))))


@pytest.mark.parametrize("variance,working", [(0., 1.), (-1., 1.), (1., 0.),
                                                (np.nan, 1.), (1., np.inf)])
def test_nonpositive_or_nonfinite_known_widths_are_rejected(variance, working):
    with pytest.raises(ValueError):
        known_width_coverage(variance, working)


def test_mixture_coverage_averages_probabilities_before_comparing_to_marginal_normal():
    paths, summary = fixture_tables()
    row = diagnostic_tables(paths, summary)["cell_calibration.csv"].query("mechanism == 'mmpp'").iloc[0]
    probabilities = [math.erf(Z95 / (i*math.sqrt(2))) for i in (1, 2, 3, 4)]
    assert row.exact_mixture_known_width_coverage_estimate == pytest.approx(sum(probabilities)/4)
    assert row.mixture_coverage_schedule_mcse == pytest.approx(np.std(probabilities, ddof=1)/2)
    marginal = math.erf(Z95/math.sqrt(2*7.5))
    assert row.marginal_normal_rms_width_coverage_approximation == pytest.approx(marginal)
    assert abs(row.exact_mixture_known_width_coverage_estimate-marginal) > .05
    assert row.A_coverage_residual_conditional_mcse == pytest.approx(
        math.sqrt(sum(p*(1-p) for p in probabilities))/4)


def test_predicted_sd_and_ratio_uncertainty_uses_independent_schedule_arms():
    paths, summary = fixture_tables()
    tables = diagnostic_tables(paths, summary)
    row = tables["cell_calibration.csv"].query("mechanism == 'mmpp'").iloc[0]
    variance_se = np.std([1., 4., 9., 16.], ddof=1)/2
    assert row.predicted_variance == 7.5
    assert row.predicted_variance_mcse == pytest.approx(variance_se)
    assert row.predicted_sd_mcse == pytest.approx(variance_se/(2*math.sqrt(7.5)))
    assert row.predicted_sd_mc_lower == pytest.approx(math.sqrt(7.5)*math.exp(-Z95*variance_se/15))
    contrast = tables["matched_prediction_comparisons.csv"].iloc[0]
    logse = .5*math.sqrt(2)*(variance_se/7.5)
    assert contrast.predicted_sd_ratio == pytest.approx(math.sqrt(2))
    assert contrast.predicted_sd_ratio_log_mcse == pytest.approx(logse)
    assert contrast.predicted_sd_ratio_mc_upper == pytest.approx(math.sqrt(2)*math.exp(Z95*logse))
    assert contrast.observed_pfml_sd_ratio == contrast.observed_pfml_rmse_ratio == 1.


def test_shifted_truth_and_center_width_decomposition_keep_the_paired_indicators():
    paths, summary = fixture_tables()
    tables = diagnostic_tables(paths, summary)
    row = tables["cell_calibration.csv"].query("mechanism == 'mmpp'").iloc[0]
    assert row.fixed_rho_rmse == pytest.approx(math.sqrt(3.5))
    assert row.A_covered == row.B_covered == 2
    assert row.C_covered == 4
    paired = tables["paired_coverage.csv"].query("mechanism == 'mmpp'").set_index("contrast")
    assert paired.loc["B_minus_A", "coverage_difference"] == 0.
    assert paired.loc["B_minus_A", "paired_mcse"] == pytest.approx(math.sqrt(2/3)/2)
    assert paired.loc["C_minus_B", "coverage_difference"] == .5
    assert paired.loc["C_minus_B", "paired_mcse"] == pytest.approx(math.sqrt(1/3)/2)
    assert paired.loc["C_minus_A", "coverage_difference"] == (
        paired.loc["C_minus_B", "coverage_difference"] + paired.loc["B_minus_A", "coverage_difference"])
    assert paired.loc["B_minus_A", "first_only_covered"] == 1
    assert paired.loc["B_minus_A", "second_only_covered"] == 1


def test_identical_paired_indicators_have_zero_difference_uncertainty():
    result = paired_difference([0, 1, 0, 1], [0, 1, 0, 1])
    assert result["coverage_difference"] == result["paired_mcse"] == 0.
    # Treating these as two independent binomials would incorrectly give >0.
    with pytest.raises(ValueError):
        paired_difference([True], [False])
    with pytest.raises(ValueError):
        paired_difference([0., .5], [0., 1.])


@pytest.mark.parametrize("mutation", ["missing", "duplicate", "unavailable", "wrong_count", "wrong_truth"])
def test_review_rejects_changed_sources_denominators_or_truth(mutation):
    paths, summary = fixture_tables()
    if mutation == "missing":
        paths = paths.iloc[:-1]
    elif mutation == "duplicate":
        paths = pd.concat([paths, paths.iloc[:1]])
    elif mutation == "unavailable":
        paths.loc[0, "interval_valid"] = False
    elif mutation == "wrong_count":
        summary.loc[0, "n_covered"] = 3
    else:
        summary.loc[0, "truth"] = 1.5
    with pytest.raises(ValueError):
        diagnostic_tables(paths, summary)


def fake_package(tmp_path):
    package = tmp_path / "sealed"
    (package/"analysis").mkdir(parents=True)
    paths, summary = fixture_tables()
    (package/"analysis/schedule_predictions.csv.gz").write_bytes(gzip.compress(csv_bytes(paths), mtime=0))
    (package/"analysis/parameter_summary.csv").write_bytes(csv_bytes(summary))
    original = pd.DataFrame([{"mechanism": mechanism, "ratio": 5., "fixed_rho_mean_sd": scale,
                              "observed_sd": scale, "observed_rmse": scale+.1, "observed_bias": .1}
                             for mechanism, scale in (("mmpp", 2.), ("iid", 1.))])
    (package/"analysis/original_e4_mechanism.csv").write_bytes(csv_bytes(original))
    for cid, group in paths.groupby("cell_id"):
        part = package/"horizon_full/cells"/cid/"part-000000-000004.jsonl.gz"
        part.parent.mkdir(parents=True)
        records = []
        for row in group.to_dict(orient="records"):
            records.append({"cell_id": cid, "rep": row["rep"], "estimator": "pfml",
                            "point_valid": True, "interval_valid": True,
                            "mu": row["pfml_mu"], "mu_se": row["pfml_mu_se"], "true_mu": 2.5,
                            "mu_lower": row["pfml_mu"]-Z95*row["pfml_mu_se"],
                            "mu_upper": row["pfml_mu"]+Z95*row["pfml_mu_se"],
                            "times_sha256": row["times_sha256"], "observations_sha256": row["observations_sha256"]})
        part.write_bytes(gzip.compress(("\n".join(json.dumps(r) for r in records)+"\n").encode(), mtime=0))
    sealed_files = {str(p.relative_to(package)): file_digest(p) for p in package.rglob("*") if p.is_file()}
    (package/"completion_seal.json").write_text(json.dumps({"status": "completed", "files": sealed_files}))
    return package


def test_additive_review_binds_raw_intervals_and_preserves_the_entire_original_seal(tmp_path):
    package = fake_package(tmp_path)
    (package/"horizon_full/.execution.lock").touch()
    before = {str(p.relative_to(package)): file_digest(p) for p in package.rglob("*") if p.is_file()}
    output = tmp_path/"review"
    manifest = build(package, output)
    assert manifest["raw_readback"]["pfml_sources_checked"] == 8
    assert manifest["raw_readback"]["every_C_indicator_matches_raw_interval"] is True
    assert manifest["raw_readback"]["raw_C_covered_by_cell"] == {"iid-5": 4, "mmpp-5": 4}
    assert verify(package, output)["status"] == "passed"
    assert before == {str(p.relative_to(package)): file_digest(p) for p in package.rglob("*") if p.is_file()}
    with pytest.raises(FileExistsError):
        build(package, output)
    with pytest.raises(FileExistsError):
        build(package, package/"additional_review")


def test_review_rejects_mutated_sealed_input_and_mutated_output(tmp_path):
    package = fake_package(tmp_path)
    output = tmp_path/"review"
    build(package, output)
    reviewed = output/"paired_coverage.csv"
    original_bytes = reviewed.read_bytes()
    reviewed.write_bytes(original_bytes+b"changed\n")
    with pytest.raises(ValueError):
        verify(package, output)
    reviewed.write_bytes(original_bytes)
    with (package/"analysis/parameter_summary.csv").open("a") as stream:
        stream.write("changed\n")
    with pytest.raises(ValueError, match="sealed bytes"):
        verify(package, output)
