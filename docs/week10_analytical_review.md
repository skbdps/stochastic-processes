# Critical review of the Week 10 analytical argument

Dated 1 November 2026. The finite-design covariance identities survive critical
review. Their transfer to the fitted PFML mean and its reported intervals
needed substantive qualifications. This correction uses the saved experiment;
it changes no estimator, seed, replication, result, figure or original seal.
The current manuscript is [Draft 3](../manuscript/week10_draft3.md).

## 1. Main correction: stochastic closeness does not transfer moments

The fitted mean is a ratio,

$$
\widehat\mu_P=\overline X_{0:K-1}
 +\frac{X_K-X_0}{K(1-\widehat\rho)}.
$$

The earlier argument correctly derived its difference from a mean computed
at a fixed coefficient. However, smallness in probability, and small
differences on the saved paths, do not establish closeness in mean square.
The denominator can approach zero from the admissible side.

The [boundary proof](week10_pfml_boundary_moments.md) conditions on the first
K observations and varies the last. The slope is affine in that last response,
and the fitted mean generically has a nonzero simple pole where the slope
reaches one. The OU transition density is positive there. On the interior
side alone, integrating the absolute error to any power p>=1 diverges. For
example, with prefix (0,1,3) and last response 4.2-delta, the fit has
rho_hat=1-5 delta/14 and mu_hat=0.4+3.92/delta, while the residual sum of
squares approaches 0.56. It is not a degenerate-regression artifact.

This concerns the ideal real-arithmetic interior estimator. It does not mean
a finite-precision program literally has infinite moments. The numerical
point-valid rule nevertheless supplies no declared statistical separation
from rho=1. Interval rejection is a separate decision. On a regular mesh,
exact and Euler regression means also equal this ratio on the common interior
region; this is not a pathology unique to clustered PFML.

The consequence is precise: empirical bias, SD and RMSE remain descriptions
of the saved experiment, but the ordinary moment-based RMSE/bias MC intervals
cannot be presented as justified confidence statements about the ideal
estimator's population risk. The usual RMSE delta method needs finite second
and fourth moments. Those assumptions cannot be inferred from 2,000 moderate
realizations. The fixed-coefficient diagnostic has finite moments and its
covariance derivation remains valid. Coverage indicators are bounded, so the
original coverage contrasts retain their independent-replication, pointwise
MC interpretation.

No observations have been removed and no new slope cutoff is imposed.
A regularized estimator or alternative loss would require a separately
declared analysis. The defensible headline now emphasizes the observed
coverage deficit and uses fitted empirical spread as supporting description.

## 2. Correct coverage object: conditional Gaussian versus marginal mixture

For fixed coefficient and an exogenous schedule t, the diagnostic is exactly
Gaussian with mean mu_0 and variance v(t). Give it the known generating-law
working width z sqrt(w(t)). Its conditional coverage is exactly

$$
C(t)=2\Phi\{z\sqrt{w(t)/v(t)}\}-1.
$$

Averaging C(t) over schedules estimates the exact mixture coverage of this
known-width diagnostic. This is what the sealed fixed-horizon field named
`normal_working_coverage_approximation` computes. Its historical name is kept,
but calling this calculation itself a Gaussian approximation was inaccurate.
In contrast, applying the same formula to averaged variances approximates the
mixture and, when counts vary, also replaces random widths by their RMS.
Neither calculation is exact coverage of the fitted, studentized PFML interval.
The known width uses generating truth and is an explanatory benchmark, not a
proposed implementable interval.

The new [cell calibration table](../artifacts/week10_analytical_review/calibration/cell_calibration.csv)
reports both calculations under distinct names and adds schedule MC uncertainty.
At r=20, exact known-width mixture coverage is estimated as 95.0066% for IID
and 90.8908% for MMPP. These expected probabilities must not be confused with
coverage indicators observed on the particular saved OU paths.

## 3. Separate the fitted center from its estimated width

The exact regression-coordinate inverse observed-Hessian mean variance is

$$
\widehat v_{\mu,\mathrm{Hess}}
=\frac{\widehat v}{(1-\widehat\rho)^2}
 \left\{\frac1K+
 \frac{(\widehat\mu_P-\overline X_{0:K-1})^2}{S_{xx}}\right\},
\qquad \widehat v=\mathrm{SSE}/K.
$$

The full derivation includes joint estimation of slope and residual variance.
The deterministic tests also compare the formula with the actual normalized
numerical full-Hessian implementation. The original leading population
working-variance formula is correct, but it does not equal every fitted width.

Using the same saved source paths, define A as the fixed-coefficient center
with the known population width, B as the fitted center with that same known
width, and C as the actual fitted center and reported Hessian width. All
14,000 fitted centers, standard errors and C coverage indicators were checked
against the original raw partitions. All fits were available, so each cell's
denominator is 2,000 throughout.

| Fixed-T cell | A observed coverage | B observed coverage | C actual coverage | Width change C-B, pp [conservative pointwise 95% MC interval] |
| --- | ---: | ---: | ---: | ---: |
| IID, r=20 | 95.60% | 95.65% | 94.40% | -1.250 [-2.340, -0.148] |
| MMPP, r=20 | 90.85% | 90.80% | 90.40% | -0.400 [-1.630, 0.833] |

Across all seven cells, center-only changes are at most 0.05 percentage points.
The IID r=20 width-only change has paired MCSE 0.342 pp. The full C-A change
there is -1.20 pp, with paired MCSE 0.338 pp. Therefore, the small observed
center discrepancy does not establish interval calibration. The
[complete paired table](../artifacts/week10_analytical_review/calibration/paired_coverage.csv)
uses the per-source indicator difference, including discordant counts;
it does not treat A, B and C as independent samples. Its confidence limits
subtract simultaneous 97.5% two-sided Clopper–Pearson limits for the two
discordance probabilities. The union bound guarantees at least 95% coverage
for each paired difference without assuming the discordant counts independent.
For zero discordances in 2,000 paths the interval is approximately +/-0.2189 pp,
although the descriptive plug-in MCSE is zero. First-pass normal intervals
are retained under explicitly approximate column names and as preserved
initial output. These remain post-result, pointwise diagnostics: the internal
Bonferroni construction does not control multiplicity across different cells
or contrasts. No observed difference is not an equivalence result.

## 4. Compare diagnostic SD with empirical SD and report prediction uncertainty

| Design | r | Predicted fixed-coefficient SD ratio | Empirical fitted SD ratio | Empirical fitted RMSE ratio |
| --- | ---: | ---: | ---: | ---: |
| Original fixed N | 5 | 1.077596 | 1.084052 | 1.083874 |
| New fixed T | 5 | 1.075742 | 1.055928 | 1.055960 |
| Original fixed N | 20 | 1.162403 | 1.158120 | 1.157863 |
| New fixed T | 20 | 1.160223 | 1.156786 | 1.156378 |

All ratios compare MMPP with its corresponding IID control. The fixed-count
prediction is an analytic generating-law calculation. The fixed-horizon
predictions average conditional variances over the saved schedules, with
approximate pointwise schedule-MC intervals [1.074827, 1.076657] at r=5 and
[1.158359, 1.162090] at r=20. These intervals quantify sampling of schedules,
not fitted-estimator population-risk uncertainty.

The benchmark uses no OU responses but shares the experiment's schedules.
It is not an independent validation experiment. Comparing its point estimate
with an empirical RMSE interval is not a formal independent calibration test.
The original paired fixed-diagnostic squared-error-minus-conditional-variance
check does retain this dependence; all seven differences are within 0.945
estimated MCSEs of zero. The fitted-minus-fixed RMS discrepancy below 0.77%
of fitted empirical SD establishes agreement on saved paths only.
Sources: [fixed-N comparison](../artifacts/week10_analytical_review/calibration/original_sd_comparisons.csv)
and [fixed-T comparison](../artifacts/week10_analytical_review/calibration/matched_prediction_comparisons.csv).

## 5. Algebra, numerical scope and conclusions retained

The review checked the signed endpoint weights, conditional OU covariance,
innovations recursion, floor-aware matrix Laplace kernel, event-index
correlation, finite-count weighted sum, and population working information.
Independent 70-digit calculations at r=1,5,20 and N=2,7,1000 agree with the
adopted q=.5 finite-count predictions to relative error below 5e-16. A
subtractive finite-sum evaluation can lose precision for extreme off-grid
near-unit correlations; that observation does not affect this study. Any
extension there needs separate numerical validation. Signed endpoint weights
also preclude a universal monotonicity claim: the net variance increase is
verified for this calibration, not every conceivable dependent schedule.

The fixed-T actual PFML coverage deficits remain -1.85 pp at r=5
(pointwise 95% MC interval [-3.369,-0.331]) and -4.00 pp at r=20
([-5.638,-2.362]). Thus unequal observed spans are not necessary for the
coverage deficit under the declared random-count and endpoint conventions.
The covariance mechanism quantitatively explains a finite-moment linear
diagnostic and approximates the observed central spread and coverage; it does
not prove a fitted population-MSE comparison, a dependent-gap CLT or a robust
PFML interval.

## 6. Evidence and reproduction

The original 250-file supplement remains sealed with SHA-256
`c595df1d521767a7191510c6e2decb9111a886b1ca6d5e6d7de7a3aa51cc2c30`.
Its generated narrative and Draft 2 are historical. This review, the revised
[mechanism note](week10_mechanism.md), the boundary proof and Draft 3 supply
the corrected interpretation. New saved-data diagnostics live outside the
original seal in [the review package](../artifacts/week10_analytical_review/).

The [calibration manifest](../artifacts/week10_analytical_review/calibration/manifest.json)
binds its source, inputs and four tables. In the pinned Week 10 environment,
`ou_irregular.week10.analytical_review.verify(package, output)` rederives them,
rechecks all 250 original file hashes, and reads all 42,000 raw records to
verify the 14,000 PFML sources. See the exact commands in
[reproduction instructions](week10_reproduction.md). Deterministic tests verify
the pole example, Hessian transformation, conditional mixture probabilities,
pairing, availability and saved-data guards. They do not add research
replications or replace academic review.
