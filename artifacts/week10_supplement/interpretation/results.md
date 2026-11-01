# Week 10 supplement: results and interpretation

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
| 5 | Original fixed N | 1.0839 [1.0375, 1.1323] | 92.75% / 95.45% | -2.700 [-4.158, -1.242] |
| 5 | Span-standardized original | 1.0813 [1.0326, 1.1323] | 92.61% / 95.30% | -2.690 [-4.252, -1.129] |
| 5 | New fixed T | 1.0560 [1.0101, 1.1039] | 92.65% / 94.50% | -1.850 [-3.369, -0.331] |
| 20 | Original fixed N | 1.1579 [1.1073, 1.2108] | 90.70% / 94.55% | -3.850 [-5.466, -2.234] |
| 20 | Span-standardized original | 1.1450 [1.0925, 1.2001] | 90.98% / 94.34% | -3.368 [-5.048, -1.687] |
| 20 | New fixed T | 1.1564 [1.1063, 1.2087] | 90.40% / 94.40% | -4.000 [-5.638, -2.362] |

The fixed-horizon PFML mean RMSE and coverage contrasts at r=5 and r=20 have
pointwise MC intervals excluding equal performance. Thus variation in total
observed span is not necessary for the adverse PFML mean effect in these
tested designs. The r=5 point estimate is smaller than in original E4, but
these independent experiments differ in observation-count distribution and
endpoint convention. Their point estimates do not quantify a causal fraction
explained by span, and no formal attenuation claim follows from this table.

At r=1 the two parent processes both reduce to Poisson sampling. The new PFML
mean RMSE ratio is 1.0184 [0.9736,
1.0653], and its coverage difference is
-0.750 pp [-2.071,
0.571]. Both intervals include no
difference. This is a useful negative control, not an equivalence test.

Across all seven cells and all three parameters, exact-likelihood coverage
ranges from 93.80% to 96.20%.
The complete parameter tables retain PFML/Euler theta and sigma failures and
Euler mean overcoverage. Those parameter-specific findings must not be
collapsed into a claim that every estimator fails because of clustering.

## 2. What the saved-span diagnostic establishes

Bins and the common-support rule were specified before generating the new
analysis tables, after seeing original E4. The same pooled span-bin weights
were used for both arms. The excluded fractions and weighted mean spans are:

| r | MMPP excluded | IID excluded | Weighted MMPP span | Weighted IID span |
| --- | --- | --- | --- | --- |
| 1 | 0.25% | 0.20% | 499.853 | 499.981 |
| 5 | 9.85% | 1.65% | 499.820 | 499.798 |
| 20 | 5.50% | 0.30% | 498.471 | 498.547 |

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
\widehat\mu_P=\overline X_{0:K-1}+
\frac{X_K-X_0}{K(1-\widehat\rho)}.
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
V\sum_{i,j} w_iw_j\exp(-\theta_0|t_i-t_j|),
\qquad V=\frac{\sigma_0^2}{2\theta_0}.
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
| 5 | Original fixed N | 1.0776 | 1.0839 | 93.11% | 92.75% |
| 5 | New fixed T | 1.0757 | 1.0560 | 93.11% | 92.65% |
| 20 | Original fixed N | 1.1624 | 1.1579 | 90.82% | 90.70% |
| 20 | New fixed T | 1.1602 | 1.1564 | 90.89% | 90.40% |

The fixed-count prediction integrates the known generating law analytically.
For the new fixed-horizon study, the prediction averages exact conditional
variances over its saved schedules. The latter has MC uncertainty from those
schedules and uses the parent-law coefficient as a reference, without claiming
it is the exact finite-window pseudo-target. Both calculations are benchmarks
for a fixed-coefficient diagnostic; coverage additionally substitutes a
population working variance and a Gaussian approximation.

The diagnostic was also evaluated on every new OU path. Across the seven
cells, its observed MSE minus its average predicted variance is at most
0.94 estimated MC standard errors from zero. The RMS difference between
the fitted PFML mean and the fixed-coefficient diagnostic is at most
0.77% of the fitted mean's empirical SD. This directly
checks the approximation that links the covariance argument to the actual
jointly fitted estimator. It does not establish exact finite-sample coverage,
a general dependent-gap CLT, or a new robust confidence-interval method.

## 4. Observation-count and boundary scope

| Mechanism | r | Mean N | Observed N range |
| --- | --- | --- | --- |
| MMPP | 1 | 1000.783 | 895–1099 |
| MMPP | 5 | 1003.398 | 815–1205 |
| MMPP | 20 | 1005.499 | 660–1341 |
| IID | 1 | 1000.456 | 902–1115 |
| IID | 5 | 1000.232 | 869–1141 |
| IID | 20 | 1002.335 | 792–1250 |
| REGULAR | 0 | 1000.000 | 1000–1000 |

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
