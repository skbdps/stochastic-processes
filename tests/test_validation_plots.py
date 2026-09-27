"""Check that exploratory figures show simulation uncertainty and exclusions."""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from ou_irregular.metrics import summarize_results
from ou_irregular import plots


def test_error_bars_are_1_96_times_mcse_not_path_intervals():
    data = pd.DataFrame({"estimator": ["exact"] * 3, "n": [250] * 3,
                         "cv": [0.0, 1.0, 2.0], "bias": [2.0, 4.0, np.nan],
                         "mcse_bias": [0.5, np.nan, np.nan]})
    fig, ax = plt.subplots()
    plots._draw_curves(ax, data, "bias", [250], "mcse_bias")
    assert len(ax.containers) == 1
    # Only CV=0 has both a summary and an estimable Monte Carlo SE.
    segments = ax.containers[0].lines[2][0].get_segments()
    assert len(segments) == 1
    np.testing.assert_allclose(segments[0], [[0, 2 - 1.96 * 0.5], [0, 2 + 1.96 * 0.5]])
    plt.close(fig)


def _all_invalid():
    return pd.DataFrame([
        {"cell_id": "all_invalid", "theta_mean_gap": 0.5, "cv": 0.0, "n": 250,
         "mean_gap": 0.5, "true_theta": 1.0, "true_mu": 0.0, "true_sigma": 0.5,
         "rep": rep, "estimator": "exact", "success": True, "optimizer_success": True,
         "fit_status": "boundary", "valid_for_point_summary": False,
         "reason": "known fixture", "solver": "fixture", "nll": 0.0,
         "theta": 100.0, "mu": 0.0, "sigma": 1.0}
        for rep in range(2)
    ])


def test_classified_plots_keep_all_invalid_cells_and_explain_scope(tmp_path, monkeypatch):
    captured = {}
    original = plots._finish

    def capture(fig, handles, title, path, **kwargs):
        captured[path.name] = {
            "title": title, "classified": kwargs.get("classified"),
            "uncertainty": kwargs.get("uncertainty", False),
            "ylabel": fig.axes[0].get_ylabel(),
            "curves": [line.get_ydata().copy() for line in fig.axes[0].lines],
        }
        original(fig, handles, title, path, **kwargs)

    monkeypatch.setattr(plots, "_finish", capture)
    summary = summarize_results(_all_invalid())
    paths = plots.plot_summary(summary, tmp_path)
    assert len(paths) == 4
    for path in paths:
        assert path.stat().st_size > 10_000
        assert path.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
        assert captured[path.name]["classified"]
    assert captured["week3_bias_vs_cv.png"]["uncertainty"]
    assert captured["week3_rmse_vs_cv.png"]["uncertainty"]
    assert not captured["week3_rmse_ratio_vs_cv.png"]["uncertainty"]
    assert "same estimator" in captured["week3_rmse_ratio_vs_cv.png"]["title"]
    failures = captured["week3_fit_failure_rate.png"]
    assert failures["ylabel"] == "Excluded / all planned fits"
    np.testing.assert_array_equal(failures["curves"][0], [1.0])
    assert np.isnan(captured["week3_bias_vs_cv.png"]["curves"][0]).all()


def test_plotting_rejects_mixed_historical_and_classified_policies(tmp_path):
    summary = summarize_results(_all_invalid())
    summary.loc[0, "validity_policy"] = "legacy_success_finite_v1"
    with pytest.raises(ValueError, match="cannot mix"):
        plots.plot_summary(summary, tmp_path)
