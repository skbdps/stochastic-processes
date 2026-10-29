# Week 6: registered full-grid estimands and reporting

Prepared 30 October 2026 (India time), after completed Week 5 publication
`d49fa4beff0800fc804afacbfdc7d7836311fdfe` and before headline E2 generation.
This document specifies implementation and presentation of the existing E2
registration. It does not amend the scientific protocol or numerical gates.

## Design and source law

E2 is the full Cartesian product

\[
q\in\{0.1,0.5,1,2\},\qquad
c\in\{0,0.5,1,2,3\},\qquad
n\in\{250,1000,4000\},\qquad R=2000.
\]

Here \(q=\theta\bar\Delta\), \(c\) is nominal gap CV, and n counts observations.
There are K=n-1 transitions. Truth remains theta=1, mu=0, sigma=0.5, and the
registered root seed is 20261002. The 60 cells contain 120,000 source paths and
360,000 exact/PFML/Euler fit records. Every path is retained regardless of fit
or interval validity.

For c=0 all gaps equal the nominal mean \(\bar\Delta=q/\theta\). Otherwise the
raw gaps are iid Gamma with shape \(c^{-2}\) and scale \(\bar\Delta c^2\).
For nominal c>=1.5 apply \(\max(\Delta,10^{-6}\bar\Delta)\), without subsequent
renormalization. Thus the E2 floor applies to c=2 and 3; its threshold remains
1.5 even though E2 contains no c=1.5 cell. Realized moments, minimum gaps, spans
and clips are saved for every source.

OU paths use exact transitions with stationary initialization. Fitted
likelihoods condition on the initial observation. Simulation truth is used to
generate paths and evaluate errors; observed-data numerical normalization in
the unchanged inference service does not receive truth.

Scientific cell IDs retain the registered convention. The source address is
`SeedSequence(root_seed, spawn_key=(4, *eight_32bit_cell_id_words, replication)).spawn(2)`.
Gap and OU randomness are separate, while all three methods share a source.
Overlapping scientific cells keep their IDs across E1 and E2 but use independent
root seeds. E2 validation uses scope `validation` and root seed 21261002; it
cannot replace headline observations.

## How the dimensions change the experiment

At fixed n, changing q changes the nominal observation span in relaxation-time
units: \(\theta T_0=(n-1)q\). Across the registered grid that quantity ranges
from 24.9 to 7998. Varying n changes both the number of transitions and nominal
span when q is fixed. These are the stated finite-sample designs; they do not
hold realized duration or a fixed-horizon asymptotic regime constant.

For an irregular cell, its regular reference has the same n, truth and nominal
mean gap. There are 12 distinct regular controls, one per (n,q), and 48
irregular cells. Realized spans are not matched. Different cells are independent;
contrasts using the same control are mutually correlated. A regular cell compared
with itself is an identity, not two independent samples.

## Unchanged estimands and inference

All definitions and derivations in [week5_mathematics.md](week5_mathematics.md)
continue unchanged. In particular:

- Bias and RMSE condition on valid point fits. Their counts and Monte Carlo
  uncertainty remain visible; missing valid fits are not silently resampled.
- Conditional coverage divides by available intervals. Operational coverage
  divides by all R planned paths, with unavailable intervals encoded as
  noncoverage. Conditional coverage is undefined when no interval is available;
  operational coverage is then zero with a nondegenerate Wilson uncertainty
  bound.
- Independent regular contrasts combine independent variances. Within-source
  estimator comparisons recompute both point metrics on the same common-valid
  subset and retain its count. Paired operational-coverage differences use all
  R paths and the covariance of the two coverage indicators.
- Mu has zero truth. Its signed bias, RMSE and standard errors use original
  state units; coverage proportions and SD/SE ratios are dimensionless.
  Relative-to-truth metrics and the relative-bias practical criterion are
  undefined. A mu RMSE ratio is well-defined when both RMSEs are available and
  the regular-reference RMSE is positive.
- Wald intervals use the full conditional-NLL Hessian in observed-data
  normalized coordinates, the registered full/half-step comparison and all
  positive-definiteness, condition, standardized-score and whitened-discrepancy
  gates. The covariance Jacobian remains
  `diag(theta_x, state_scale, sigma_x)`. Natural-scale intervals are symmetric
  and are not clipped. PFML/Euler inverse Hessians are model-based working
  covariances, not sandwich estimates.

## Pointwise classifications

For theta and sigma, a supported breach occurs when the operational-coverage
Wilson upper bound is below 0.80, or the entire approximate 95% relative-bias
Monte Carlo interval lies beyond +0.10 or -0.10. Within criteria requires the
coverage Wilson lower bound to be at least 0.80 and the entire relative-bias
interval inside [-0.10,0.10]. Other cases are unresolved. Missing bias
uncertainty prevents within-criteria classification, while a supported coverage
breach still counts. Mu is labelled `not_applicable_zero_truth`.

The E2 reporting layer copies these classifications from the independently
audited core table without rounding or rethresholding. Displayed percentages
may be rounded; labels always use the unrounded calculations. These labels
describe the investigated points. Connecting or filling grid squares does not
establish a continuous operating boundary, interpolation guarantee or
simultaneous coverage across the grid.

## Complete tables, figures and failures

The parameter-summary, independent-reference and paired-comparison tables each
have 540 rows. The spacing table has 60 rows, and numerical diagnostics have
180 cell/estimator rows. The outcome table has one row per observed status
composition; its length can exceed 180 when a cell/method has multiple statuses.
Every cell/method's status counts must sum to R.

Week 6 adds three derived companion CSVs under `figures/`, which is the frozen
workflow's derived-product namespace:

| File | Contract |
| --- | --- |
| classification_table.csv | All 540 parameter rows, grouped by n/parameter/estimator/q/CV |
| classification_counts.csv | All 108 n/parameter/estimator/class combinations, including zero counts; each group sums to 20 cells |
| failed_cells.csv | All 180 cell/estimator rows, including zero failures, with separate point-invalid and interval-unavailable counts and complete status composition |

Point-invalid paths are already among interval-unavailable paths, so those two
counts must not be added to claim distinct failed fits. Poor coverage or a
supported scientific breach is not a numerical failure. The complete raw failed
fit records remain in the core `failed_fits.csv.gz` export.

F5 shows operational coverage and F6 shows conditional coverage for every
sample size, parameter and estimator on the full 4-by-5 q/CV grid. Companion
availability maps retain the same cells. Undefined quantities are labelled NA,
not converted to zero. Figure captions distinguish denominators, the registered
pointwise criteria and model-based covariance. Full uncertainty is retained in
the source tables rather than discarded by the heatmap presentation.

The report and figure manifests bind every derived CSV/PNG/PDF to the saved
parent-table hashes and design. Immutable regeneration, full-key checks and
visual inspection precede completion. The numerical engine and independent
Week 5 audit remain unchanged, and explicitly listed Week 6 dependencies are
frozen before headline execution. Results are recorded separately after the full
run; configuration files alone do not constitute completion.
