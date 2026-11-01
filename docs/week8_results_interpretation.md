# Week 8 E4 results and mathematical interpretation

Prepared 31 October 2026 from the recovered original full E4 tables. Numerical
claims below trace to
[parameter_summary.csv](../artifacts/week8_e4/parameter_summary.csv),
[matched_iid_comparisons.csv](../artifacts/week8_e4/matched_iid_comparisons.csv),
[source_summary.csv](../artifacts/week8_e4/source_summary.csv) and
[spacing_summary.csv](../artifacts/week8_e4/spacing_summary.csv). These files
belong to the experiment executed at
`0ff57919f5d8345f52d8689732ff94448ca3ca1d`; the
[completion record](week8_completion.md) separates original execution from
later recovery, actual figure inspection and release checks.

## What was compared

Each of seven independent cells has R=2,000 source paths, n=1,000 observations,
truth (1,0,0.5) and nominal q=0.5. MMPP and IID cells at each r=1,5,20 share the
entire event-stationary one-gap law and the declared floor transformation.
The MMPP retains arrival phase between gaps; IID resets to an independent Palm
phase before every gap. Different cells are independent, including r=1,
where MMPP and IID are independent draws from the same exponential-gap law.
The regular control is a separate constant-gap design. Estimators within a
source path are paired.

All 42,000 points and intervals are available. Bias/RMSE therefore use 2,000
points per cell/method/parameter, and operational and conditional coverage
coincide. Unavailable intervals would count as noncoverage operationally;
their absence in this run does not remove that rule. Relative mu metrics are
undefined because its truth is zero.

## Complete parameter overview

The following percentages are relative bias for theta/sigma and operational
coverage. Mu bias is in its natural parameter units. Full MCSEs, Wilson
intervals, RMSEs, reported versus empirical standard errors and classifications
remain in the 63-row parameter table. Coverage MCSEs describe variation across
independent simulated paths, not uncertainty for a single fitted path.

| Design | Method | theta bias (%) | sigma bias (%) | mu bias | Coverage theta / mu / sigma (%) |
| --- | --- | ---: | ---: | ---: | --- |
| MMPP r=1 | exact | +1.333 | -0.027 | -0.000406 | 95.55 / 94.65 / 94.90 |
| MMPP r=1 | pfml | -17.807 | -9.972 | -0.000491 | 30.75 / 94.80 / 5.75 |
| MMPP r=1 | euler | -32.662 | -16.536 | -0.000371 | 0.05 / 97.60 / 0.00 |
| MMPP r=5 | exact | +1.429 | +0.099 | +0.000282 | 95.55 / 94.55 / 95.30 |
| MMPP r=5 | pfml | -25.928 | -14.587 | -0.000009 | 7.60 / 92.75 / 0.25 |
| MMPP r=5 | euler | -38.341 | -14.616 | +0.000397 | 0.00 / 98.55 / 0.00 |
| MMPP r=20 | exact | +1.100 | +0.020 | -0.000389 | 95.20 / 94.45 / 95.50 |
| MMPP r=20 | pfml | -45.636 | -26.977 | -0.000419 | 0.00 / 90.70 / 0.00 |
| MMPP r=20 | euler | -52.736 | -10.488 | -0.000284 | 0.00 / 99.20 / 0.35 |
| IID r=1 | exact | +1.238 | +0.083 | +0.000058 | 95.25 / 94.90 / 95.40 |
| IID r=1 | pfml | -17.982 | -9.933 | +0.000039 | 30.45 / 95.05 / 5.75 |
| IID r=1 | euler | -32.827 | -16.435 | +0.000246 | 0.00 / 97.70 / 0.00 |
| IID r=5 | exact | +1.732 | +0.115 | +0.000209 | 94.70 / 94.80 / 95.10 |
| IID r=5 | pfml | -25.813 | -14.576 | +0.000472 | 7.05 / 95.45 / 0.20 |
| IID r=5 | euler | -38.177 | -14.666 | +0.000403 | 0.00 / 98.40 / 0.00 |
| IID r=20 | exact | +1.265 | -0.034 | -0.000020 | 94.80 / 95.05 / 94.10 |
| IID r=20 | pfml | -45.895 | -27.024 | +0.000739 | 0.00 / 94.55 / 0.00 |
| IID r=20 | euler | -52.672 | -10.573 | -0.000108 | 0.00 / 99.00 / 0.25 |
| Regular | exact | +1.205 | +0.205 | +0.001124 | 95.50 / 94.25 / 94.80 |
| Regular | pfml | +1.205 | +0.205 | +0.001124 | 95.50 / 94.25 / 94.80 |
| Regular | euler | -20.684 | -20.532 | +0.001124 | 2.35 / 94.25 / 0.00 |

Exact theta/sigma meet the adopted practical criteria in every design: 14
within-criteria labels. The regular PFML fits coincide with exact and add two
within-criteria labels. Stochastic PFML and all Euler theta/sigma combinations
give 26 supported breaches; there are no unresolved E4 labels. All 21 mu
classifications remain inapplicable. These criteria use a coverage threshold
of 80% and a relative-bias threshold of 10% with prespecified MC uncertainty;
passing them does not assert precisely 95% coverage.

Exact relative theta bias ranges from about +1.10% to +1.73%; its sigma bias
is small. Exact coverage across all parameters/designs spans 94.10–95.55%,
with independent-path MCSEs around half a percentage point. This supports
the specified exact conditional-likelihood procedure within this design; it
does not establish universal robustness or eliminate finite-sample bias.
PFML/Euler theta/sigma intervals often have very low coverage of the true OU
parameters despite complete numerical availability. Those are retained
scientific findings, not pipeline exceptions.

## Dependence at matched full marginals

The next table reports **mu** using the saved independent-cell MMPP-minus-IID contrasts.
RMSE ratios are MMPP/IID. Coverage differences and MC limits are in percentage
points; the intervals are pointwise normal MC intervals for coverage and
log-delta MC intervals for positive RMSE ratios. Equal replication numbers
do not create cross-cell pairs. All components have 2,000 valid fits.

| r | Method | RMSE ratio [95% MC] | Coverage difference pp [95% MC] |
| ---: | --- | --- | --- |
| 1 | exact | 0.994528 [0.951207, 1.039822] | -0.25 [-1.629, +1.129] |
| 1 | pfml | 1.002571 [0.958997, 1.048125] | -0.25 [-1.610, +1.110] |
| 1 | euler | 1.000717 [0.957877, 1.045472] | -0.10 [-1.039, +0.839] |
| 5 | exact | 0.998990 [0.956918, 1.042912] | -0.25 [-1.642, +1.142] |
| 5 | pfml | 1.083874 [1.037545, 1.132272] | -2.70 [-4.158, -1.242] |
| 5 | euler | 0.987276 [0.945330, 1.031084] | +0.15 [-0.610, +0.910] |
| 20 | exact | 1.028231 [0.983061, 1.075476] | -0.60 [-1.982, +0.782] |
| 20 | pfml | 1.157863 [1.107265, 1.210774] | -3.85 [-5.466, -2.234] |
| 20 | euler | 1.007900 [0.962998, 1.054896] | +0.20 [-0.385, +0.785] |

PFML mu shows increased RMSE and reduced coverage at r=5 and r=20, with the
displayed pointwise MC intervals excluding their null values. At r=20 the
absolute RMSEs are 0.0353987 for MMPP and 0.0305724 for IID; coverage is
90.70% versus 94.55%. Its matched bias difference is −0.00115784 with MC
interval [−0.00320791,0.00089222], so this comparison does not establish a
mean-bias shift. At r=5 PFML mu coverage is 92.75% versus 95.45%.
Exact/Euler mu contrast intervals contain the null values in this design.

Every one of the 27 matched bias MC intervals contains zero. All theta/sigma
RMSE-ratio MC intervals contain one. A single r=20 exact-sigma coverage
contrast is +1.40 percentage points, with pointwise interval
[0.02457,2.77543]; this is an isolated contrast among many examined and is
not evidence of a general beneficial dependence effect. Other theta/sigma
coverage contrasts either include zero or have an empirical-zero normal
interval because both observed success counts are zero. A zero empirical
MCSE in those cases does not establish an exactly zero underlying coverage
probability or a certain null difference; the component Wilson intervals in
the parameter table retain uncertainty.

The matched full-marginal comparison assesses consequences of gap dependence
within the fixed design, including the distribution of realized span. It
does not condition on or hold each path's span fixed. The regular comparisons
match nominal design, and comparisons sharing their single regular reference
are not mutually independent. Neither those contrasts nor an r-curve alone
isolates dependence. There is no simultaneous or familywise inference claim.

## Population targets explain marginal effects, not finite-sample equality

The prospective [population-target derivation](week8_population_targets.md)
was recorded before fitted E4 observations. For actual gap Y, let
mbar=E[Y], L(s)=E[exp(−sY)], V=sigma0²/(2theta0) and
H(theta0)=E[(1−exp(−theta0Y))/Y]. Under the declared exogeneity and stationary
OU initialization, the expected working objectives give

$$
\theta_{P,*}=-\log L(\theta_0)/\bar m,\qquad
\sigma_{P,*}^2=2V\theta_{P,*},\qquad \mu_{P,*}=\mu_0,
$$

$$
\theta_{E,*}=(1-L(\theta_0))/\bar m,\qquad
\sigma_{E,*}^2=V\{2H(\theta_0)-(1-L(\theta_0))^2/\bar m\},\qquad
\mu_{E,*}=\mu_0.
$$

Euler uses each actual gap; dropping the H term would give the wrong random-gap
diffusion target. MMPP and matched IID share these targets because their full
one-gap marginals agree. The r-curve changes the marginal and hence the
targets: PFML theta targets are about 0.810930, 0.728444 and 0.531875, while
Euler theta targets are about 0.666667, 0.610526 and 0.467026. All summaries
continue to assess the true OU parameters; results were not recentered at
these pseudo-targets.

Matching expected objectives does not imply equal finite-sample variance,
mean-gap fluctuation, span variation or coverage. Serial dependence can change
covariance between score contributions. In particular, PFML's realized mean
gap remains part of each actual fit; the target calculation does not replace
it by the nominal mean. PFML/Euler full-Hessian Wald covariance is the
registered model-based working covariance, without a sandwich correction.
Shared targets do not make it robust to misspecification or dependence.
Probability-limit interpretations require appropriate ergodic laws, and no
general dependent-gap CLT or robust variance formula is established here.

## Source diagnostics, flooring and span

The source summaries take one record per path, avoiding triple counting
estimators and avoiding treatment of dependent gaps as independent samples.
Their MCSEs use independent paths. Raw gaps precede the floor; transformed
gaps precede cumulative summation; realized spacing is obtained by differencing
the stored timestamps used for OU simulation. Tiny cumulative roundoff can
distinguish transformed and stored-time minimum gaps.

| Design | Raw mean (MCSE) | Raw adjacent product (MCSE) | Realized span mean (MCSE) | Span SD | Clipped gaps |
| --- | --- | --- | --- | ---: | ---: |
| MMPP r=1 | 0.500032 (0.000348) | 0.250187 (0.000390) | 499.532 (0.348) | 15.546 | 0 |
| MMPP r=5 | 0.499722 (0.000676) | 0.312249 (0.000784) | 499.222 (0.675) | 30.192 | 0 |
| MMPP r=20 | 0.499642 (0.001144) | 0.520181 (0.002205) | 499.143 (1.143) | 51.122 | 7 |
| IID r=1 | 0.500832 (0.000370) | 0.251174 (0.000415) | 500.331 (0.370) | 16.536 | 0 |
| IID r=5 | 0.500759 (0.000455) | 0.251138 (0.000553) | 500.259 (0.454) | 20.317 | 0 |
| IID r=20 | 0.500297 (0.000749) | 0.251462 (0.001052) | 499.797 (0.748) | 33.445 | 10 |
| Regular | 0.500000 (0.000000) | 0.250000 (0.000000) | 499.500 (0.000) | 0.000 | 0 |

The raw adjacent-product expectations are 0.25 for IID, 0.3125 for MMPP r=5
and approximately 0.520970 for MMPP r=20. The two MMPP sample values are each
within 0.4 across-path MCSE of their respective expectations. The IID r=1
adjacent product is about 2.83 MCSE above 0.25; this descriptive full-study
fluctuation is retained and is not a reason to change seeds or institute a
new acceptance gate. MMPP r=20 has span SD about 51.12 versus IID
33.45; r=5 has about 30.19 versus 20.32. This illustrates why matching the
marginal and expected span does not fix realized-span variation. Mean
within-path sample CV is also distinct from pooled marginal CV; the source
CV figure labels that distinction.

Only r=20 is floored at 5e−7. Seven MMPP and ten IID raw gaps are clipped out
of 1,998,000 gaps per cell, versus theoretical clip probability about
3.88799e−6. The source retains exact counts and independent-path descriptive
uncertainty; dependent MMPP gaps are not assigned an independent-binomial
confidence interval. The defined transformed mean is about
0.500000000000972, without renormalization. No latent-clock padding is
inserted. Consequently the transformed MMPP-derived sequence is described
as floored MMPP-derived gaps rather than an unchanged-clock MMPP.

## Limits

E4 is exploratory, with one truth, q, observation count and finite ratio grid.
Its pointwise findings do not locate a universal operating boundary. Exact
conditional inference excludes the stationary initial density even though
simulation starts stationary. Natural-scale Wald bounds remain unclipped.
The numerical audit has stated shared dependencies and endpoint-only full-fit
replay. Raster inspection establishes legibility within its scope; it does
not establish independent correctness of the numerical model. Historical
hardware-sensitive byte replay and unattributed manifest incidents retain
their documented uncertainty. Weeks 5–7 outcomes and failure classifications
remain unchanged; this report does not certify personal mastery, novelty,
paper acceptance or later manuscript stages.
