"""Independent-design and paired-path Monte Carlo contrasts.

These are uncertainty estimates for simulation summaries, not new parameter
interval methods. Different CV cells are independent; equal replication labels
do not induce pairing. Estimators within a cell share their source arrays.
Point contrasts condition on each design's valid fits (independent comparison)
or on the common-valid intersection (paired comparison). Operational coverage
always encodes unavailable intervals as zero over all planned replications.
"""
from dataclasses import replace
from itertools import combinations

import numpy as np
import pandas as pd

from ..week4.summaries import PARAMETERS
from .design import ESTIMATORS
from .summaries import ExperimentSummarizer, covered, mean_mcse, validate_records


Z = 1.959963984540054


def _interval(row, name, value, se):
    row.update({name: value, f"{name}_mcse": se,
                f"{name}_mc_lower": value - Z * se,
                f"{name}_mc_upper": value + Z * se})


def _ratio(row, numerator, denominator, log_se):
    # A zero numerator has a defined ratio, but the square-root/log delta
    # derivatives are singular. Preserve the ratio and withhold uncertainty.
    ratio = numerator / denominator if np.isfinite(denominator) and denominator > 0 else np.nan
    row.update(rmse_ratio=ratio, rmse_ratio_mcse=ratio * log_se,
               rmse_ratio_log_mcse=log_se,
               rmse_ratio_mc_lower=ratio * np.exp(-Z * log_se),
               rmse_ratio_mc_upper=ratio * np.exp(Z * log_se),
               rmse_ratio_uncertainty="delta_method_log_ratio_95pct")


def _relative(row, prefix, truth):
    for suffix in ("", "_mcse", "_mc_lower", "_mc_upper"):
        row[f"relative_{prefix}{suffix}"] = row[f"{prefix}{suffix}"] / abs(truth) if truth else np.nan


def regular_reference_comparisons(records, design):
    """Compare each cell to its SAME-estimator CV=0 nominal-design control.

    For independent cells, Var(B_a-B_0)=Var(B_a)+Var(B_0). For positive
    RMSEs A,B, log(A/B) has delta variance Var(A)/A²+Var(B)/B². Exponentiating
    its approximate MC interval retains positive ratio endpoints. No matching
    of realized spans or pure attribution to gap CV is asserted.

    A CV=0 row compared with itself has an algebraically zero difference and
    unit ratio when defined, with zero MC uncertainty. It is never treated as
    two independent estimates. Missing metric values remain missing.
    """
    summary = ExperimentSummarizer().summarize(records, design)
    index = summary.set_index(["cell_id", "estimator", "parameter"])
    cell_map = {c.identity: c for c in design.cells}
    rows = []
    for a in summary.itertuples(index=False):
        cell = cell_map[a.cell_id]
        reference_id = replace(cell, cv=0).identity
        key = (reference_id, a.estimator, a.parameter)
        if key not in index.index:
            raise ValueError(f"missing exact nominal-design regular reference for {a.cell_id}")
        b = index.loc[key]
        self_reference = a.cell_id == reference_id
        row = {"cell_id": a.cell_id, "reference_cell_id": reference_id,
               "estimator": a.estimator, "parameter": a.parameter, "truth": a.truth,
               "n": cell.n, "theta_mean_gap": cell.q, "cv": cell.cv,
               "comparison_type": "self_reference" if self_reference else "independent_cells",
               "point_subset": "each_cell_own_valid_point_fits",
               "matching": "same_n_truth_nominal_mean_gap_and_floor_settings_not_realized_span"}
        for field in ("n_planned", "n_point_valid", "n_interval_valid", "n_covered", "bias", "rmse",
                      "mcse_bias", "mcse_rmse", "coverage_operational", "coverage_conditional"):
            row[f"cell_{field}"], row[f"reference_{field}"] = getattr(a, field), b[field]
        for name, field, uncertainty in (
                ("excess_bias", "bias", "mcse_bias"),
                ("operational_coverage_difference", "coverage_operational", "operational_mcse"),
                ("conditional_coverage_difference", "coverage_conditional", "conditional_mcse")):
            value = getattr(a, field) - b[field]
            se = (0.0 if np.isfinite(value) else np.nan) if self_reference else np.hypot(getattr(a, uncertainty), b[uncertainty])
            _interval(row, name, value, se)
        _relative(row, "excess_bias", a.truth)
        log_se = np.nan
        if a.rmse > 0 and b.rmse > 0:
            log_se = 0.0 if self_reference else np.hypot(a.mcse_rmse / a.rmse, b.mcse_rmse / b.rmse)
        _ratio(row, a.rmse, b.rmse, log_se)
        rows.append(row)
    return pd.DataFrame(rows)


def paired_estimator_comparisons(records, design):
    """A-minus-B comparisons on common paths, retaining paired covariance.

    Let e_A,e_B be errors on S common-valid paths, M_j=mean(e_j²), R_j=sqrt(M_j).
    Bias and MSE difference MCSEs are sd(e_A-e_B)/sqrt(S) and
    sd(e_A²-e_B²)/sqrt(S). RMSE-difference delta influence is
    e_A²/(2R_A)-e_B²/(2R_B). Log-RMSE-ratio influence is
    (e_A²/M_A-e_B²/M_B)/2. The sample SD of each influence divided by sqrt(S)
    includes its covariance automatically; constants from centering cancel.
    Delta derivatives are undefined when a required RMSE is zero. Sample Monte
    Carlo uncertainty is unavailable when S<2.

    Coverage comparison uses D_r=I_A,r-I_B,r in {-1,0,1} on all R planned paths.
    Its MCSE is sd(D)/sqrt(R), using sample variance. The discordant counts
    expose the pairing; adding independent binomial variances is incorrect.
    """
    data = validate_records(records, design)
    # Also enforce full point/interval semantic checks before Boolean masking.
    ExperimentSummarizer().summarize(data, design)
    rows = []
    for cell in design.cells:
        group = data[data.cell_id == cell.identity]
        methods = {name: group[group.estimator == name].set_index("rep").sort_index()
                   for name in ESTIMATORS}
        for name_a, name_b in combinations(ESTIMATORS, 2):
            a, b = methods[name_a], methods[name_b]
            common = a.point_valid.to_numpy() & b.point_valid.to_numpy()
            common_intervals = a.interval_valid.to_numpy() & b.interval_valid.to_numpy()
            count = int(common.sum())
            for p in PARAMETERS:
                truth = cell.to_dict()[f"true_{p}"]
                ea = a[p].to_numpy(dtype=float)[common] - truth
                eb = b[p].to_numpy(dtype=float)[common] - truth
                ma = float(np.mean(ea ** 2)) if count else np.nan
                mb = float(np.mean(eb ** 2)) if count else np.nan
                ra, rb = np.sqrt(ma), np.sqrt(mb)
                row = {"cell_id": cell.identity, "n": cell.n, "theta_mean_gap": cell.q, "cv": cell.cv,
                       "estimator_a": name_a, "estimator_b": name_b, "parameter": p, "truth": truth,
                       "comparison_type": "paired_shared_source_paths", "direction": "a_minus_b_and_a_over_b",
                       "n_planned": design.replications, "n_point_valid_a": int(a.point_valid.sum()),
                       "n_point_valid_b": int(b.point_valid.sum()), "n_common_point_valid": count,
                       "n_common_interval_valid": int(common_intervals.sum()),
                       "point_subset": "intersection_of_valid_point_fits",
                       "coverage_subset": "all_planned_paths_unavailable_counts_as_not_covering",
                       "bias_a_common": float(ea.mean()) if count else np.nan,
                       "bias_b_common": float(eb.mean()) if count else np.nan,
                       "rmse_a_common": ra, "rmse_b_common": rb,
                       "squared_error_covariance_common": float(np.cov(ea ** 2, eb ** 2, ddof=1)[0, 1]) if count > 1 else np.nan}
                for name, differences in (("bias_difference", ea - eb), ("mse_difference", ea ** 2 - eb ** 2)):
                    _interval(row, name, float(differences.mean()) if count else np.nan, mean_mcse(differences))
                _relative(row, "bias_difference", truth)
                rmse_se = mean_mcse(ea ** 2 / (2 * ra) - eb ** 2 / (2 * rb)) if ra > 0 and rb > 0 else np.nan
                _interval(row, "rmse_difference", ra - rb, rmse_se)
                log_se = mean_mcse((ea ** 2 / ma - eb ** 2 / mb) / 2) if ma > 0 and mb > 0 else np.nan
                _ratio(row, ra, rb, log_se)
                ca, cb = covered(a, p), covered(b, p)
                differences = ca.astype(float) - cb.astype(float)
                row.update(n_covered_a=int(ca.sum()), n_covered_b=int(cb.sum()),
                           n_both_cover=int((ca & cb).sum()), n_a_only_cover=int((ca & ~cb).sum()),
                           n_b_only_cover=int((~ca & cb).sum()), n_neither_cover=int((~ca & ~cb).sum()),
                           operational_coverage_a=float(ca.mean()), operational_coverage_b=float(cb.mean()))
                _interval(row, "operational_coverage_difference", float(differences.mean()), mean_mcse(differences))
                rows.append(row)
    return pd.DataFrame(rows)
