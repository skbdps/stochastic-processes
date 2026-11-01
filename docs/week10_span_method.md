# Week 10: saved E4 span sensitivity

This is a post-result supplementary analysis of the original fixed-count E4
experiment. It preserves the original archive and its unconditional contrasts.
The outcome-independent analysis specification was frozen in the
[Week 10 supplement protocol](week10_supplement_protocol.md) before producing
the supplementary effect tables. It does not retroactively preregister E4.

## Question and strata

E4 fixes 1,000 observations and has nominal span 499.5, while realized span is
random. Its MMPP/IID comparison includes changes to the span distribution.
This analysis asks how the observed contrast changes after giving both arms
the same coarse span distribution. It complements the new fixed-horizon study.

Use intervals `[0,400)`, `[400,450)`, `[450,475)`, `[475,500)`, `[500,525)`,
`[525,550)`, `[550,600)`, and `[600,infinity)`, separately at ratios 1, 5 and 20.
All bins, including empty bins, appear in the output. All three estimators and
all three parameters are retained. The regular reference is retained in the
compact source export but does not enter the matched stochastic contrasts.

Let `n_ah` be the number of planned paths in arm `a` and stratum `h`. A stratum
enters standardization only if both arms have at least 30 planned paths. This
support rule uses schedule counts, never the estimated effects, fit validity
or coverage outcomes. On the selected support `S`, both arms use

$$
w_h=\frac{n_{Mh}+n_{Ih}}{\sum_{g\in S}(n_{Mg}+n_{Ig})}.
$$

The output records every support decision, shared weight, each arm's excluded
fraction, and actual within-bin span means and ranges. An empty support yields
undefined standardized estimates, with excluded fractions one.

## Estimands and uncertainty

Within each bin, bias and MSE use point-valid errors. RMSE is the square root
of MSE. Operational coverage counts every planned path; an unavailable
interval counts as noncoverage. Conditional coverage and availability are
also reported. No failure is silently removed from the coverage denominator.

Standardized bias and MSE are weighted means of the corresponding within-bin
valid-point means. Standardized RMSE is the square root of standardized MSE;
it is not the weighted mean of bin RMSEs. When point availability differs,
these are standardized conditional-on-valid-point metrics, with availability
shown explicitly. Operational coverage is the weighted mean of within-bin
operational proportions. The MMPP-minus-IID difference and MMPP/IID RMSE ratio
are reported alongside the unadjusted original-sample contrasts.

Monte Carlo uncertainty conditions on the observed stratum counts, selected
support and weights. For a within-bin sample mean `m_ah` with estimated
variance `v_ah`, the standardized mean has estimated variance

$$
\widehat{\operatorname{Var}}\!\left(\sum_h w_hm_{ah}\right)
=\sum_h w_h^2v_{ah}.
$$

For bias and MSE, `v_ah` is the unbiased sample variance of errors or squared
errors divided by the number of valid point fits. For operational coverage,
it is the binomial plug-in variance `p_ah(1-p_ah)/n_ah`. Independent-arm
variances add for differences. RMSE-ratio intervals use the log delta method,
with gradient `(+1/(2 MSE_M), -1/(2 MSE_I))`. All intervals are pointwise.

These intervals exclude uncertainty from estimating the weights and selecting
support. They are conditional descriptive intervals, not unconditional
repeated-study guarantees. The binwise coverage table also gives Wilson
intervals, including nonzero uncertainty for zero or complete observed
coverage. If an estimated coverage-difference variance is zero, its MCSE is
reported as zero but its normal confidence bounds remain undefined; no exact
certainty is inferred from an empirical boundary.

## Interpretation limits

Conditioning on total span changes the distribution of individual gaps and
their dependence. Even within a shared bin the arms may have different mean
spans and other schedule properties. Standardization neither holds span
exactly fixed nor identifies a causal direct effect of clustering. Binwise
results with small counts remain visible even when excluded from the pooled
contrast. Results cannot justify choosing narrower bins or new cut points
after observing the desired effect.

The fixed-horizon experiment answers a related but different design question:
it fixes the observation window and lets the observation count vary. Its
boundary and sampling conventions are specified separately in the supplement
protocol. Agreement across these analyses strengthens the interpretation;
disagreement is an informative result to report.

## Reproduction and evidence

`ou_irregular.week10.span.read_e4_records` verifies the full original seven-cell
design, all planned fit keys, original partition hashes and shared-source
identities. It reads 42,000 fits on 14,000 sources. The writer records source
file hashes, creates a separate immutable output directory, and verifies the
source hashes again afterward. It generates no new schedules, observations,
fits or intervals. It cannot overwrite an existing analysis or nest output
inside the original E4 archive.

`verify_span_analysis` reads the same immutable sources again, checks the
complete settings and output inventory, and rederives every exported table
and the compact fit export. It rejects changed source bytes, missing or extra
exports, edited support rules, and changed arithmetic even if an edited
output's manifest digest has also been changed. The verifier writes nothing.

The outputs are complete bin summaries, bin contrasts, support weights,
standardized contrasts, unadjusted contrasts and a compressed compact fit
export. The compact export retains all estimators, parameters, validity flags,
source keys and observation/time hashes; the original full raw records remain
authoritative. Synthetic tests independently check the weighted means and
variances, boundary conventions, failure denominators, missing support and
invalid-record rejection.
