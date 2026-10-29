import numpy as np
import pandas as pd
import pytest

from ou_irregular.experiments.design import Cell, RunDesign, ESTIMATORS
from ou_irregular.experiments.summaries import ExperimentSummarizer, validate_records, spacing_summary


def fixture_records():
    design = RunDesign("test", 42, 3, (Cell(10, .5, 0),), "validation")
    records = []
    for rep in range(3):
        for name in ESTIMATORS:
            row = {**design.cells[0].to_dict(), "rep": rep, "estimator": name,
                   "point_valid": rep < 2, "interval_valid": rep == 0,
                   "interval_status": "valid" if rep == 0 else "unavailable",
                   "realized_span": 4.5}
            for p, truth in [("theta", 1), ("mu", 0), ("sigma", .5)]:
                row.update({p: truth + [0, 2, 100][rep], f"{p}_lower": truth-.1 if rep == 0 else np.nan,
                            f"{p}_upper": truth+.1 if rep == 0 else np.nan, f"{p}_se": .05 if rep == 0 else np.nan})
            records.append(row)
    return design, pd.DataFrame(records)


def test_hand_calculated_metrics_and_failure_denominators():
    design, data = fixture_records()
    summary = ExperimentSummarizer().summarize(data, design)
    theta = summary.query("estimator == 'exact' and parameter == 'theta'").iloc[0]
    assert theta.bias == 1
    assert theta.mcse_bias == 1
    assert theta.rmse == pytest.approx(np.sqrt(2))
    assert theta.mcse_rmse == pytest.approx(1/np.sqrt(2))
    assert theta.coverage_conditional == 1
    assert theta.coverage_operational == pytest.approx(1/3)
    assert theta.point_availability == pytest.approx(2/3)
    assert theta.interval_availability == pytest.approx(1/3)
    assert theta.conditional_wilson_lower < 1
    assert summary.query("parameter == 'mu'").relative_rmse.isna().all()
    assert spacing_summary(data).n_planned.tolist() == [3]


@pytest.mark.parametrize("failure", ["missing_group", "missing_rep", "duplicate", "metadata"])
def test_global_plan_integrity(failure):
    design, data = fixture_records()
    if failure == "missing_group":
        data = data[data.estimator != "pfml"]
    elif failure == "missing_rep":
        data = data.iloc[1:]
    elif failure == "duplicate":
        data = pd.concat([data, data.iloc[:1]])
    else:
        data.loc[0, "n"] = 42
    with pytest.raises(ValueError):
        validate_records(data, design)


def test_no_valid_fits_keep_operational_denominator():
    design, data = fixture_records()
    data["point_valid"] = False
    data["interval_valid"] = False
    summary = ExperimentSummarizer().summarize(data, design)
    assert summary.bias.isna().all() and summary.mcse_rmse.isna().all()
    assert summary.coverage_conditional.isna().all()
    assert summary.coverage_operational.eq(0).all()
    assert summary.operational_wilson_upper.gt(0).all()
