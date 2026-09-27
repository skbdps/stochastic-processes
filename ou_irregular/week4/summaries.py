"""Coverage accounting and predeclared practical criteria, outside inference.

Only experiment evaluation sees simulation truth. Inference services receive
observations and return intervals; these summaries count all planned outcomes.
"""
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.stats import norm


PARAMETERS = ("theta", "mu", "sigma")


def wilson_interval(successes, trials, level=0.95):
    """Binomial Wilson interval; no fabricated interval when denominator is zero."""
    for count in (successes, trials):
        if isinstance(count, (bool, np.bool_)) or not isinstance(count, (int, np.integer)) or count < 0:
            raise ValueError("binomial counts must be nonnegative integers")
    if isinstance(level, (bool, np.bool_)) or not np.isfinite(level) or not 0 < level < 1:
        raise ValueError("level must be strictly between zero and one")
    if not 0 <= successes <= trials:
        raise ValueError("binomial counts must satisfy 0 <= successes <= trials")
    if trials == 0:
        return np.nan, np.nan
    p, z = successes / trials, norm.ppf((1 + level) / 2)
    denominator = 1 + z * z / trials
    center = (p + z * z / (2 * trials)) / denominator
    half = z * np.sqrt(p * (1 - p) / trials + z * z / (4 * trials ** 2)) / denominator
    # Include the empirical proportion exactly, even at all-zero/all-one
    # endpoints where floating-point cancellation otherwise creates negative
    # plotting error bars of order 1e-19.
    return max(0., min(p, float(center - half))), min(1., max(p, float(center + half)))


@dataclass(frozen=True)
class ThresholdPolicy:
    coverage_threshold: float = 0.80
    relative_bias_threshold: float = 0.10

    def classify(self, coverage_lower, coverage_upper, bias_lower, bias_upper, truth):
        if truth == 0:
            return "not_applicable_zero_truth"
        coverage_breach = np.isfinite(coverage_upper) and coverage_upper < self.coverage_threshold
        bias_available = np.isfinite([bias_lower, bias_upper]).all()
        bias_breach = bias_available and (bias_lower > self.relative_bias_threshold
                                         or bias_upper < -self.relative_bias_threshold)
        if coverage_breach or bias_breach:
            return "supported_breach"
        if (np.isfinite(coverage_lower) and coverage_lower >= self.coverage_threshold
                and bias_available and bias_lower >= -self.relative_bias_threshold
                and bias_upper <= self.relative_bias_threshold):
            return "within_criteria"
        return "unresolved"


@dataclass(frozen=True)
class CoverageSummarizer:
    level: float = 0.95
    policy: ThresholdPolicy = ThresholdPolicy()

    def summarize(self, records, *, planned_replications):
        if (isinstance(planned_replications, (bool, np.bool_))
                or not isinstance(planned_replications, (int, np.integer)) or planned_replications < 1):
            raise ValueError("planned_replications must be a positive integer")
        data = pd.DataFrame(records)
        required = {"cell_id", "rep", "estimator", "point_valid", "interval_valid", "interval_status",
                    *PARAMETERS, *(f"true_{p}" for p in PARAMETERS),
                    *(f"{p}_{end}" for p in PARAMETERS for end in ("lower", "upper"))}
        if required.difference(data):
            raise ValueError(f"missing fields: {sorted(required.difference(data))}")
        if data.empty or data.duplicated(["cell_id", "rep", "estimator"]).any():
            raise ValueError("empty or duplicate planned records")
        if not data.rep.map(lambda x: isinstance(x, (int, np.integer))
                            and not isinstance(x, (bool, np.bool_)) and x >= 0).all():
            raise ValueError("replication IDs must be nonnegative integers")
        for column in ("point_valid", "interval_valid"):
            if not data[column].map(lambda x: isinstance(x, (bool, np.bool_))).all():
                raise ValueError(f"{column} must be Boolean")
        if (data.interval_valid & ~data.point_valid).any():
            raise ValueError("a valid interval requires a valid point fit")
        finite_parameters = np.isfinite(data[list(PARAMETERS)].to_numpy(dtype=float)).all(axis=1)
        if (data.point_valid & ~finite_parameters).any():
            raise ValueError("valid point estimates must be finite")
        if (data.point_valid & ((data.theta <= 0) | (data.sigma <= 0))).any():
            raise ValueError("valid theta and sigma estimates must be positive")
        for p in PARAMETERS:
            valid_bounds = (np.isfinite(data[f"{p}_lower"]) & np.isfinite(data[f"{p}_upper"])
                            & (data[f"{p}_lower"] <= data[f"{p}_upper"]))
            if (data.interval_valid & ~valid_bounds).any():
                raise ValueError("valid intervals require finite ordered endpoints")
        rows = []
        z = norm.ppf((1 + self.level) / 2)
        for (cell_id, estimator), group in data.groupby(["cell_id", "estimator"], sort=True):
            if set(group.rep) != set(range(planned_replications)) or len(group) != planned_replications:
                raise ValueError("all planned replication IDs must survive aggregation")
            total = len(group)
            point_mask = group.point_valid.to_numpy(dtype=bool)
            interval_mask = group.interval_valid.to_numpy(dtype=bool)
            n_point, n_interval = int(point_mask.sum()), int(interval_mask.sum())
            for p in PARAMETERS:
                if group[f"true_{p}"].nunique(dropna=False) != 1:
                    raise ValueError("truth must be constant within a scientific cell")
                truth = float(group[f"true_{p}"].iloc[0])
                if not np.isfinite(truth):
                    raise ValueError("truth must be finite")
                errors = group[p].to_numpy(dtype=float)[point_mask] - truth
                bias = float(errors.mean()) if n_point else np.nan
                rmse = float(np.sqrt(np.mean(errors ** 2))) if n_point else np.nan
                mcse = float(errors.std(ddof=1) / np.sqrt(n_point)) if n_point > 1 else np.nan
                relative_bias = bias / abs(truth) if truth else np.nan
                relative_mcse = mcse / abs(truth) if truth else np.nan
                covered = interval_mask & (group[f"{p}_lower"].to_numpy() <= truth) & (truth <= group[f"{p}_upper"].to_numpy())
                n_covered = int(covered.sum())
                conditional = n_covered / n_interval if n_interval else np.nan
                operational = n_covered / total
                cond_low, cond_high = wilson_interval(n_covered, n_interval, self.level)
                op_low, op_high = wilson_interval(n_covered, total, self.level)
                bias_low, bias_high = relative_bias - z * relative_mcse, relative_bias + z * relative_mcse
                rows.append({"cell_id": cell_id, "estimator": estimator, "parameter": p,
                             "truth": truth, "n_planned": total, "n_point_valid": n_point,
                             "n_interval_valid": n_interval, "n_covered": n_covered,
                             "point_availability": n_point / total, "interval_availability": n_interval / total,
                             "bias": bias, "rmse": rmse, "mcse_bias": mcse,
                             "relative_bias": relative_bias, "relative_mcse_bias": relative_mcse,
                             "relative_bias_mc_lower": bias_low, "relative_bias_mc_upper": bias_high,
                             "coverage_conditional": conditional, "coverage_operational": operational,
                             "conditional_mcse": np.sqrt(conditional * (1 - conditional) / n_interval) if n_interval else np.nan,
                             "operational_mcse": np.sqrt(operational * (1 - operational) / total),
                             "conditional_wilson_lower": cond_low, "conditional_wilson_upper": cond_high,
                             "operational_wilson_lower": op_low, "operational_wilson_upper": op_high,
                             "negative_lower_count": int((interval_mask & (group[f"{p}_lower"].to_numpy() < 0)).sum()),
                             "practical_classification": self.policy.classify(op_low, op_high, bias_low, bias_high, truth),
                             "checkpoint_in_reference_band": bool(n_interval and .94 <= conditional <= .96),
                             "point_metric_scope": "conditional_on_valid_point_fit",
                             "primary_coverage_scope": "all_planned_paths_unavailable_counts_as_not_covering"})
        return pd.DataFrame(rows)
