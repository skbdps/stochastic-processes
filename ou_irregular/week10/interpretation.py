"""Generate the supplement's numerical narrative with exact source bindings."""
from pathlib import Path
import json
import math

import pandas as pd

from ..experiments.provenance import canonical, file_digest
from ..experiments.storage import atomic_bytes, atomic_json


def render(package):
    package = Path(package)
    names = {"original": "span/span_unadjusted.csv", "standardized": "span/span_standardized.csv",
             "horizon": "analysis/matched_iid_comparisons.csv",
             "original_mechanism": "analysis/original_e4_mechanism.csv",
             "horizon_mechanism": "analysis/horizon_mechanism.csv",
             "parameters": "analysis/parameter_summary.csv", "counts": "analysis/source_summary.csv"}
    tables = {key: pd.read_csv(package/name, float_precision="round_trip") for key, name in names.items()}
    bindings = {key: {"path": name, "sha256": file_digest(package/name)} for key, name in names.items()}
    claims = {"schema_version": 1, "sources": bindings, "contrasts": [], "mechanism": []}
    rows = []
    labels = {"original": "Original fixed N", "standardized": "Span-standardized original",
              "horizon": "New fixed T"}
    for ratio in (5,20):
        for key in ("original", "standardized", "horizon"):
            a = tables[key].query("ratio == @ratio and estimator == 'pfml' and parameter == 'mu'").iloc[0]
            claims["contrasts"].append({"source": key, "row_key": {"ratio": ratio, "estimator": "pfml", "parameter": "mu"},
                                         "values": a.to_dict()})
            rows.append(f"| {ratio} | {labels[key]} | {a.rmse_ratio:.4f} [{a.rmse_ratio_mc_lower:.4f}, {a.rmse_ratio_mc_upper:.4f}] | "
                        f"{100*a.mmpp_coverage_operational:.2f}% / {100*a.iid_coverage_operational:.2f}% | "
                        f"{100*a.coverage_difference:.3f} [{100*a.coverage_difference_mc_lower:.3f}, {100*a.coverage_difference_mc_upper:.3f}] |")
    mechanism_rows = []
    for ratio in (5,20):
        for key, label, pred, observed, cov in (
            ("original_mechanism", "Original fixed N", "fixed_rho_mean_sd", "observed_rmse", "observed_coverage"),
            ("horizon_mechanism", "New fixed T", "predicted_fixed_rho_sd", "observed_pfml_rmse", "observed_pfml_coverage")):
            a = tables[key].query("ratio == @ratio and mechanism == 'mmpp'").iloc[0]
            b = tables[key].query("ratio == @ratio and mechanism == 'iid'").iloc[0]
            values={"ratio": ratio, "source": key, "predicted_sd_ratio": a[pred]/b[pred],
                    "observed_rmse_ratio": a[observed]/b[observed],
                    "approximate_mmpp_coverage": a.normal_working_coverage_approximation,
                    "observed_mmpp_coverage": a[cov]}
            claims["mechanism"].append(values)
            mechanism_rows.append(f"| {ratio} | {label} | {values['predicted_sd_ratio']:.4f} | {values['observed_rmse_ratio']:.4f} | "
                                  f"{100*values['approximate_mmpp_coverage']:.2f}% | {100*values['observed_mmpp_coverage']:.2f}% |")
    support = []
    for ratio in (1,5,20):
        a=tables["standardized"].query("ratio == @ratio and estimator == 'pfml' and parameter == 'mu'").iloc[0]
        support.append(f"| {ratio} | {100*a.mmpp_excluded_fraction:.2f}% | {100*a.iid_excluded_fraction:.2f}% | {a.mmpp_span_mean:.3f} | {a.iid_span_mean:.3f} |")
    null=tables["horizon"].query("ratio == 1 and estimator == 'pfml' and parameter == 'mu'").iloc[0]
    exact=tables["parameters"].query("estimator == 'exact'")
    counts=tables["counts"]
    countrows=[f"| {r.mechanism.upper()} | {int(r.ratio)} | {r.actual_n_mean:.3f} | {int(r.actual_n_min)}–{int(r.actual_n_max)} |" for r in counts.itertuples(index=False)]
    diagnostic=tables["horizon_mechanism"]
    max_z=max(abs(diagnostic.fixed_rho_mse_minus_predicted_variance/diagnostic.fixed_rho_mse_minus_predicted_variance_mcse))
    max_correction=max(diagnostic.fitted_minus_fixed_mu_rms/diagnostic.observed_pfml_sd)
    claims["context"]={"exact_coverage_min": exact.coverage_operational.min(),
                         "exact_coverage_max": exact.coverage_operational.max(),
                         "max_fixed_coefficient_mse_discrepancy_mcse_units": max_z,
                         "max_fitted_minus_fixed_rms_over_fitted_sd": max_correction,
                         "poisson_control": null.to_dict()}
    text=f"""# Week 10 supplement: results and interpretation

This is a post-result exploratory extension of E4, authorized after the Week 10
methods baseline. The original fixed-N study is preserved. Two additional
analyses address span variation: descriptive reweighting of the saved original
records, and an independent fixed-horizon experiment. A separate covariance
calculation explains the mean-estimation mechanism without fitting its
parameters to the observed error curves.

## 1. Equal-span results

The new experiment contains seven cells, 2,000 independent sources per cell,
14,000 sources and 42,000 estimator fits. Every observed span is exactly
T=499.5. All points and intervals were available; there were zero pipeline
errors. Observation counts vary, and the same three estimators and original
full-Hessian Wald inference are used.

| r | Analysis | PFML mean RMSE ratio, MMPP/IID [95% MC interval] | Coverage, MMPP / IID | Coverage difference, percentage points [95% MC interval] |
| --- | --- | --- | --- | --- |
{chr(10).join(rows)}

The fixed-horizon PFML mean RMSE and coverage contrasts at r=5 and r=20 have
pointwise MC intervals excluding equal performance. Thus variation in total
observed span is not necessary for the adverse PFML mean effect in these
tested designs. The r=5 point estimate is smaller than in original E4, but
these independent experiments differ in observation-count distribution and
endpoint convention. Their point estimates do not quantify a causal fraction
explained by span, and no formal attenuation claim follows from this table.

At r=1 the two parent processes both reduce to Poisson sampling. The new PFML
mean RMSE ratio is {null.rmse_ratio:.4f} [{null.rmse_ratio_mc_lower:.4f},
{null.rmse_ratio_mc_upper:.4f}], and its coverage difference is
{100*null.coverage_difference:.3f} pp [{100*null.coverage_difference_mc_lower:.3f},
{100*null.coverage_difference_mc_upper:.3f}]. Both intervals include no
difference. This is a useful negative control, not an equivalence test.

Across all seven cells and all three parameters, exact-likelihood coverage
ranges from {100*exact.coverage_operational.min():.2f}% to {100*exact.coverage_operational.max():.2f}%.
The complete parameter tables retain PFML/Euler theta and sigma failures and
Euler mean overcoverage. Those parameter-specific findings must not be
collapsed into a claim that every estimator fails because of clustering.

## 2. What the saved-span diagnostic establishes

Bins and the common-support rule were specified before generating the new
analysis tables, after seeing original E4. The same pooled span-bin weights
were used for both arms. The excluded fractions and weighted mean spans are:

| r | MMPP excluded | IID excluded | Weighted MMPP span | Weighted IID span |
| --- | --- | --- | --- | --- |
{chr(10).join(support)}

PFML mean differences persist under this coarse standardization. These
intervals condition on the observed support, bin counts and weights; they
omit uncertainty from estimating those choices. Within-bin span variation
and differences between schedules remain. Unequal tail exclusions change the
target population. The analysis is therefore corroborating descriptive
evidence, not a matched-span causal experiment or a decomposition of the
original effect into percentages due to different mechanisms.

## 3. A quantitative covariance mechanism

At an interior PFML solution, with K=N-1,

$$
\\widehat\\mu_P=\\overline X_{{0:K-1}}+
\\frac{{X_K-X_0}}{{K(1-\\widehat\\rho)}}.
$$

Thus interior observations have equal event-index weights. Bursts produce
many observations of closely related OU states; those extra observations
provide less independent information about the long-run mean than their
count suggests. Exogenous sampling does not systematically select high or
low OU states, so this is a precision and interval-calibration mechanism,
without a demonstrated mean-bias shift.

For a coefficient rho held fixed, the corresponding linear estimator has
weights -rho/[K(1-rho)] at the first observation, 1/K in the interior, and
1/[K(1-rho)] at the last. Its schedule-conditional variance is exactly

$$
V\\sum_{{i,j}} w_iw_j\\exp(-\\theta_0|t_i-t_j|),
\\qquad V=\\frac{{\\sigma_0^2}}{{2\\theta_0}}.
$$

In the original event-stationary fixed-count design, let B be the matrix
Laplace transform of the joint used-gap/ending-phase kernel at theta_0 and
phi=alpha B 1. The lag-h OU correlation is alpha B^h 1 for MMPP and phi^h
for IID gaps. The matched one-gap law makes lag one and the population
working variance agree; it does not equate higher-lag correlations. The
working AR(1) variance consequently misses additional serial covariance.
The floor-aware kernel and finite-count weighted sums are derived and tested
in the companion mathematical note.

| r | Design | Predicted mean SD ratio | Observed PFML RMSE ratio | Approximate MMPP coverage | Observed MMPP coverage |
| --- | --- | --- | --- | --- | --- |
{chr(10).join(mechanism_rows)}

The fixed-count prediction integrates the known generating law analytically.
For the new fixed-horizon study, the prediction averages exact conditional
variances over its saved schedules. The latter has MC uncertainty from those
schedules and uses the parent-law coefficient as a reference, without claiming
it is the exact finite-window pseudo-target. Both calculations are benchmarks
for a fixed-coefficient diagnostic; coverage additionally substitutes a
population working variance and a Gaussian approximation.

The diagnostic was also evaluated on every new OU path. Across the seven
cells, its observed MSE minus its average predicted variance is at most
{max_z:.2f} estimated MC standard errors from zero. The RMS difference between
the fitted PFML mean and the fixed-coefficient diagnostic is at most
{100*max_correction:.2f}% of the fitted mean's empirical SD. This directly
checks the approximation that links the covariance argument to the actual
jointly fitted estimator. It does not establish exact finite-sample coverage,
a general dependent-gap CLT, or a new robust confidence-interval method.

## 4. Observation-count and boundary scope

| Mechanism | r | Mean N | Observed N range |
| --- | --- | --- | --- |
{chr(10).join(countrows)}

Both endpoints are observed. The stochastic parent process starts at an event
with its Palm phase distribution; it is not initialized at an arbitrary
calendar-time origin. Arrivals strictly before T are retained and X(T) is
appended. Accordingly, the r=1 Poisson design has expected N=2+T/m=1001,
whereas the regular design has exactly 1000 observations. Event-Palm MMPP
initialization can also affect finite-window counts near the origin.

The terminal observation interval is not a fresh parent interarrival and is
not subjected to the parent-gap floor. Stopping and endpoint insertion can
change the finite-window retained-gap distribution even when the underlying
complete-gap laws match. Random N and its distribution remain features of
this contrast. These results establish persistence at a common observed
span under the declared boundary rule, not a universal dependence-only
effect across all possible observation designs.

## 5. Research conclusion and remaining manuscript work

The supplemental evidence strengthens the original E4 interpretation:
clustered exogenous schedules can reduce PFML mean precision and lead its
working Wald intervals to undercover, even when every path spans the same
time interval. Higher-lag OU covariance provides a quantitative explanation
of the measured size. The matched one-gap population target for the mean
remains correct, while working uncertainty need not reflect actual spread.

All uncertainty statements are pointwise and the supplement is exploratory.
No simulation seed, estimator, interval formula or replication count was
changed to obtain these outcomes. Week 11 still needs the integrated E1–E4
Results, practical guidance, broader related work and full discussion, followed
by academic review. Week 12 remains revision, final layout and submission
preparation. This supplement does not certify journal novelty or acceptance.
"""
    return text, claims


def write(package, output):
    output=Path(output)
    if output.exists():
        raise FileExistsError("interpretation output must be new")
    text, claims=render(package)
    output.mkdir(parents=True)
    atomic_bytes(output/"results.md",text.encode(),immutable=True)
    atomic_json(output/"claims.json",claims,immutable=True)
    atomic_json(output/"manifest.json",{"files":{p.name:file_digest(p) for p in sorted(output.iterdir())}},immutable=True)
    return claims


def verify(package):
    package=Path(package)
    output=package/"interpretation"
    text,claims=render(package)
    if (output/"results.md").read_bytes()!=text.encode():
        raise ValueError("generated interpretation differs from its source tables")
    if canonical(json.loads((output/"claims.json").read_text()))!=canonical(claims):
        raise ValueError("numerical claim bindings differ from their source tables")
    manifest=json.loads((output/"manifest.json").read_text())
    if set(p.name for p in output.iterdir())!=set(manifest["files"])|{"manifest.json"}:
        raise ValueError("interpretation inventory changed")
    for name,digest in manifest["files"].items():
        if file_digest(output/name)!=digest:
            raise ValueError("interpretation digest changed")
    return {"status":"passed","contrast_claims":len(claims["contrasts"]),
            "mechanism_claims":len(claims["mechanism"]),"context_rederived":True}
