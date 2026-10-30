"""Author checks of the renderer's literal-input and displayed-domain contract.

No observation is generated and no estimator/auditor is run in this module.
Independent suites separately check numerical artist values and the archive
lifecycle. These focused checks cover unavailable evidence and display domains.
"""
from copy import deepcopy

import matplotlib.pyplot as plt
import numpy as np
import pytest

from ou_irregular.week7 import plots
from .test_plots_reference import (forbid_science, literal_floor_tables,
                                    literal_followup_tables, literal_selected_payload)


@pytest.fixture
def floor_input():
    return literal_floor_tables()


def test_canonical_filename_sets_are_eight_png_pdf_pairs():
    names = (*plots.FOLLOWUP_STEMS, *plots.FLOOR_STEMS, *(item[2] for item in plots.SELECTED_CASES))
    assert len(names) == len(set(names)) == 8
    assert len({stem+"."+suffix for stem in names for suffix in ("png", "pdf")}) == 16


@pytest.mark.parametrize("value", [True, 1.5, -1, 6])
def test_summary_rejects_nonintegral_boolean_and_out_of_range_counts(value):
    design, tables = literal_followup_tables()
    tables["parameter_summary.csv"]["n_point_valid"] = tables["parameter_summary.csv"]["n_point_valid"].astype(object)
    tables["parameter_summary.csv"].loc[0, "n_point_valid"] = value
    with pytest.raises(ValueError, match="count|denominator"):
        plots._summary(tables, design)


def test_summary_requires_the_whole_cell_method_parameter_cartesian_set():
    design, tables = literal_followup_tables()
    summary = tables["parameter_summary.csv"]
    tables["parameter_summary.csv"] = summary[summary.cell_id != design.cells[0].identity]
    with pytest.raises(ValueError, match="incomplete"):
        plots._summary(tables, design)


def test_epsilon_rejects_counts_that_violate_intersection_lower_bound(floor_input):
    design, _, epsilon = floor_input
    mask = (epsilon.cv == 1.5) & (epsilon.floor_fraction == 1e-8) & (epsilon.estimator == "exact")
    epsilon.loc[mask, "n_common_point_valid"] = 0
    epsilon.loc[mask, "n_common_interval_valid"] = 0
    with pytest.raises(ValueError, match="relationships"):
        plots._epsilon_inputs(epsilon, design)


@pytest.mark.parametrize("column,value", [
    ("reference_cell_id", "0"*64),
    ("source_seed", 20261009),
    ("ratio_direction", "registered_floor_over_epsilon"),
])
def test_epsilon_identity_and_direction_cannot_be_relabelled(floor_input, column, value):
    design, _, epsilon = floor_input
    epsilon.loc[0, column] = value
    with pytest.raises(ValueError, match="direction|design"):
        plots._epsilon_inputs(epsilon, design)


def test_spacing_caption_cannot_claim_r_when_a_field_uses_fewer_sources(floor_input):
    design, tables, epsilon = floor_input
    tables["spacing_summary.csv"].loc[0, "min_gap_n_available"] = design.replications-1
    try:
        with forbid_science(), pytest.raises(ValueError, match="all R"):
            plots.build_floor_figures(tables, epsilon, design)
    finally:
        plt.close("all")


def test_selected_literal_payload_cannot_acquire_original_e2_label():
    payload = literal_selected_payload(.5)
    with forbid_science(), pytest.raises(ValueError, match="literal"):
        plots.build_selected_figure(payload, validation=False)


def test_selected_case_identity_is_fixed():
    payload = literal_selected_payload(.5)
    payload["selection"]["rep"] += 1
    with forbid_science(), pytest.raises(ValueError, match="fixed"):
        plots.build_selected_figure(payload, validation=True)


def test_selected_unavailable_payload_retains_its_failure_and_six_empty_panels():
    payload = literal_selected_payload(.5)
    payload["status"] = "diagnostic_failure"
    payload["diagnostic"] = {"status": "failed", "error": "Literal retained arithmetic diagnostic failure"}
    with forbid_science():
        figure = plots.build_selected_figure(payload, validation=True)
    try:
        text = "\n".join(item.get_text() for item in figure.texts) + "\n" + "\n".join(
            item.get_text() for axis in figure.axes for item in axis.texts)
        assert "LITERAL VALIDATION FIXTURE" in text
        assert "unavailable" in text.lower()
        assert "Literal retained arithmetic diagnostic failure" in text
        assert "No source, fit, profile, Hessian or replacement interval" in text
        assert len(figure.axes) == 6
        assert all(not axis.lines and not axis.collections for axis in figure.axes)
    finally:
        plt.close(figure)


def test_selected_adjacent_pair_cannot_be_displayed_under_a_different_step():
    payload = literal_selected_payload(.5)
    payload["diagnostic"]["fits"]["exact"]["float64_ladder"]["adjacent"][0]["fine_relative_step"] = .00025
    try:
        with forbid_science(), pytest.raises(ValueError, match="pair identities"):
            plots.build_selected_figure(payload, validation=True)
    finally:
        plt.close("all")


def test_profile_requires_ordered_complete_saved_grid():
    payload = literal_selected_payload(.5)
    payload["diagnostic"]["exact_profile"]["grid_log_theta"] = [-1., 1., 0., 2.]
    with forbid_science(), pytest.raises(ValueError, match="unordered"):
        plots.profile_series(payload)


def test_rmse_tick_formatter_keeps_log_transform_and_compact_near_unit_labels(floor_input):
    design, tables, epsilon = floor_input
    mask = epsilon.parameter == "theta"
    epsilon.loc[mask, "rmse_ratio"] = 1.000025
    epsilon.loc[mask, "rmse_ratio_mc_lower"] = 1.00001
    epsilon.loc[mask, "rmse_ratio_mc_upper"] = 1.00004
    with forbid_science():
        figures = plots.build_floor_figures(tables, epsilon, design)
    try:
        figure = figures[plots.FLOOR_STEMS[2]]
        axis = next(axis for axis in figure.axes if axis.get_gid() == "ratio:1e-08:theta")
        assert axis.get_yscale() == "log"
        label = axis.yaxis.get_major_formatter()(1.000025, 0)
        assert "times" not in label and len(label) <= 10
        assert float(label) == pytest.approx(1.000025, abs=1e-10)
    finally:
        for figure in figures.values():
            plt.close(figure)


def test_profile_builder_does_not_mutate_saved_payload():
    payload = literal_selected_payload(0.)
    before = deepcopy(payload)
    with forbid_science():
        result = plots.profile_series(payload)
    assert payload == before
    assert np.any(result["delta_nll"] < 0)
