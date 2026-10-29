"""Complete-denominator Monte Carlo summaries, composed with frozen Week 4 CI policy."""
import itertools

import numpy as np
import pandas as pd

from ..week4.summaries import CoverageSummarizer, PARAMETERS, wilson_interval
from .design import ESTIMATORS


def validate_records(records, design):
    """Verify the whole Cartesian plan, including entirely absent groups."""
    data = pd.DataFrame(records)
    keys = ["cell_id", "rep", "estimator"]
    if not set(keys).issubset(data) or data.duplicated(keys).any():
        raise ValueError("missing keys or duplicate planned records")
    expected = set(itertools.product((c.identity for c in design.cells),
                                     range(design.replications), ESTIMATORS))
    if set(data[keys].itertuples(index=False, name=None)) != expected:
        raise ValueError("records do not equal the complete planned Cartesian key set")
    for cell in design.cells:
        group = data[data.cell_id == cell.identity]
        for name, value in cell.to_dict().items():
            if name not in group or not group[name].eq(value).all():
                raise ValueError(f"inconsistent scientific metadata: {name}")
    return data


def mean_mcse(values):
    """Independent-replication SE of a mean; undefined with fewer than two draws."""
    values = np.asarray(values, dtype=float)
    return float(values.std(ddof=1) / np.sqrt(len(values))) if len(values) > 1 else np.nan


def covered(group, parameter):
    truth = group[f"true_{parameter}"].to_numpy()
    return (group.interval_valid.to_numpy(dtype=bool)
            & (group[f"{parameter}_lower"].to_numpy() <= truth)
            & (truth <= group[f"{parameter}_upper"].to_numpy()))


class ExperimentSummarizer:
    """Extend, rather than modify, historical coverage and classification rules.

    Conditional point metrics use S valid fits. Delta-method RMSE uncertainty is
    sd(e²)/(2 sqrt(S) RMSE), since d sqrt(M)/dM=1/(2 sqrt(M)). This approximation
    is undefined at zero RMSE; neither that case nor S<2 receives a fabricated SE.
    Availability and operational coverage always retain all planned paths.
    """

    def summarize(self, records, design):
        data = validate_records(records, design)
        summary = CoverageSummarizer().summarize(data, planned_replications=design.replications)
        summary = summary.drop(columns="checkpoint_in_reference_band")
        groups = {(cid, name): group for (cid, name), group in data.groupby(["cell_id", "estimator"])}
        additions = []
        for row in summary.itertuples(index=False):
            group = groups[row.cell_id, row.estimator]
            p, truth = row.parameter, row.truth
            errors = group.loc[group.point_valid, p].to_numpy() - truth
            rmse_se = mean_mcse(errors ** 2) / (2 * row.rmse) if row.rmse > 0 else np.nan
            extra = {"mcse_rmse": rmse_se,
                     "relative_rmse": row.rmse / abs(truth) if truth else np.nan,
                     "relative_mcse_rmse": rmse_se / abs(truth) if truth else np.nan,
                     "empirical_sd_point_valid": float(errors.std(ddof=1)) if len(errors) > 1 else np.nan,
                     "bias_mc_lower": row.bias - 1.959963984540054 * row.mcse_bias,
                     "bias_mc_upper": row.bias + 1.959963984540054 * row.mcse_bias}
            # Empirical spread and reported SE are compared on the SAME eligible
            # interval-valid subset, not on two different selected populations.
            estimates = group.loc[group.interval_valid, p].to_numpy()
            se = group.loc[group.interval_valid, f"{p}_se"].to_numpy()
            extra.update(empirical_sd_interval_valid=float(estimates.std(ddof=1)) if len(estimates) > 1 else np.nan,
                         mean_reported_se=float(se.mean()) if len(se) else np.nan,
                         rms_reported_se=float(np.sqrt(np.mean(se ** 2))) if len(se) else np.nan)
            for kind, count in (("point", row.n_point_valid), ("interval", row.n_interval_valid)):
                proportion = count / row.n_planned
                lo, hi = wilson_interval(count, row.n_planned)
                extra.update({f"{kind}_availability_mcse": np.sqrt(proportion * (1-proportion) / row.n_planned),
                              f"{kind}_availability_wilson_lower": lo, f"{kind}_availability_wilson_upper": hi})
            extra.update({key: group[key].iloc[0] for key in
                          ("n", "theta_mean_gap", "cv", "mean_gap", "floor_fraction")})
            extra["nominal_dimensionless_span"] = (extra["n"] - 1) * extra["theta_mean_gap"]
            additions.append(extra)
        return pd.concat([summary.reset_index(drop=True), pd.DataFrame(additions)], axis=1)


def spacing_summary(records):
    """Each source schedule counted once, even though three estimators share it."""
    data = pd.DataFrame(records)
    source = data[data.estimator == "exact"]
    fields = ("realized_mean_gap", "realized_cv", "realized_span", "clip_fraction", "min_gap")
    rows = []
    for cell_id, group in source.groupby("cell_id"):
        row = {"cell_id": cell_id, "n_planned": len(group),
               **{k: group[k].iloc[0] for k in ("n", "theta_mean_gap", "cv")}}
        for name in fields:
            values = pd.to_numeric(group.get(name, pd.Series(dtype=float)), errors="coerce").dropna().to_numpy()
            row.update({f"{name}_n_available": len(values),
                        f"{name}_mean": float(values.mean()) if len(values) else np.nan,
                        f"{name}_mcse": mean_mcse(values),
                        f"{name}_min": float(values.min()) if len(values) else np.nan,
                        f"{name}_max": float(values.max()) if len(values) else np.nan})
        rows.append(row)
    return pd.DataFrame(rows)


def outcome_counts(records):
    data = pd.DataFrame(records)
    keys = ["cell_id", "n", "theta_mean_gap", "cv", "estimator", "point_status", "interval_status"]
    result = data.groupby(keys, dropna=False).size().rename("count").reset_index()
    totals = data.groupby(["cell_id", "estimator"]).size()
    result["n_planned"] = [totals[c, e] for c, e in zip(result.cell_id, result.estimator)]
    result["fraction"] = result["count"] / result.n_planned
    bounds = [wilson_interval(int(c), int(n)) for c, n in zip(result["count"], result.n_planned)]
    result["wilson_lower"], result["wilson_upper"] = zip(*bounds)
    return result
