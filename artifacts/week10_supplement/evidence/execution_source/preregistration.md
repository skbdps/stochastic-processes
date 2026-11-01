# OU irregular sampling preregistration

Version 1, 27 September 2026. This protocol is frozen after the disclosed Week 3
exploratory work and numerical repairs, and before the Week 4 checkpoint or any
E1/E2 headline runs. It adopts the [plan amendment](docs/plan_amendment_week3_validation.md).
The historical plan, exploratory artifacts and all later amendments remain
available. The Git commit containing this protocol and its configuration is
the freeze record; result manifests will record its SHA and file digests.

## Question and estimands

Measure finite-sample point-estimation error and interval coverage for exact,
realized-mean-gap PFML, and Euler OU estimation across the stated observation
designs. All likelihoods condition on the initial observation. Data use the
exact OU transition, with stationary initialization, theta=1, mu=0, sigma=0.5.
No method is required to win, and no operating boundary is required to appear.
The Week 3 smoke findings are prior exploratory information, not confirmatory
evidence under this registration.

## Frozen designs

| Run | Design | Paths per cell | Root seed |
| --- | --- | ---: | ---: |
| Week 4 checkpoint | Exact estimator; q=0.5, CV=0, n=1000 | 2000 | 20260927 |
| Week 4 bootstrap diagnostics | q=0.5; CV=0 and 2; n=1000; source replication 0; all three estimators | One fixed source path per cell; 500 bootstrap draws per estimator | Source 20260927; bootstrap 20260928 |
| E1 | q=0.5, n=1000; CV in {0,0.25,0.5,1,1.5,2,3}; all estimators | 2000 | 20261001 |
| E2 | q in {0.1,0.5,1,2}; CV in {0,0.5,1,2,3}; n in {250,1000,4000}; all estimators | 2000 | 20261002 |

Here n counts observations, K=n-1 counts transitions, q=theta times nominal
mean gap. CV=0 is regular spacing; CV>0 is Gamma with shape CV^-2 and scale
nominal_mean_gap times CV^2. Apply the literal floor max(gap,1e-6 times nominal
mean gap) when nominal CV>=1.5. Do not renormalize the floored gaps; report clip
fractions, realized moments and spans. This freezes epsilon for headline work;
the earlier smoke use was exploratory.

Scientific cell IDs are SHA256 hashes of n, q, nominal gap, CV, truth and floor
settings. Week 4 source streams use SeedSequence(root_seed,
spawn_key=(4, *eight_32bit_cell_id_words, replication)).spawn(2), separating
gap and path randomness. E1/E2 will retain this addressing rule with their own
root seeds. Each estimator receives the same source path within a cell. Different
CV cells have independent streams and are not paired. Bootstrap draw streams
are independent children of a deterministic cell/estimator seed; the seed is
recorded, and no draw is replaced after failure.

## Point fitting and inference procedure

Use the validated strict Week 3 point estimators through the adapter interface.
Keep their definitions and point-validity classification. Inference does not
turn a boundary or failed point estimate into an ordinary interior estimate.
All numerical scaling is derived from observed data, never simulation truth.

For y=(x-a)/b, u=t/m, use normalized working coordinates
eta=(log(theta_u),mu_u,log(sigma_u)). The observed information is the numerical
Hessian of the full conditional NLL, not an optimizer's approximate inverse
Hessian and not a profiled one-dimensional Hessian. Central differences use
h_j=0.001 max(1,abs(eta_j)). Compare the h and h/2 matrices; use the half-step
matrix only if finite and the relative Frobenius discrepancy is <=0.005.
Require a positive-definite symmetric matrix, normalized-coordinate condition
number <=1e10, and sqrt(gradient^T inverse(H) gradient)<=0.01. Failures are
explicit interval-invalid outcomes; do not regularize or replace failed Hessians
until a separately justified, disclosed procedure has been registered.

Also check stability in weak-curvature directions: if H_half=L L^T is its
Cholesky factorization, require the largest absolute eigenvalue of
L^-1 (H_full-H_half) L^-T to be <=0.01. Record this whitened discrepancy as well
as the Frobenius discrepancy. This additional gate was chosen during mathematical
review before the checkpoint; a small unweighted matrix error alone can mask
large covariance error along a weakly identified direction.

Map covariance to original natural parameters with
J=diag(theta_x,b,sigma_x). The primary intervals are symmetric natural-scale
Wald intervals estimate +/- z_(0.975) SE. Retain negative lower bounds for
positive parameters without clipping; record them. Alternative log-Wald,
profile, sandwich or bias-corrected intervals are additional methods, not silent
substitutes. Hessian-based PFML/Euler intervals are deliberately model-based
working intervals under possible misspecification, not robust sandwich intervals.

Bootstrap diagnostics fix the observed timestamps and observed initial state,
simulate exact OU paths at each fitted plug-in parameter vector, and refit the
same estimator. Report linear-interpolated 2.5%/97.5% percentile endpoints and
the sample SD of valid draws. Require at least 475 of 500 valid draws and at
least two; otherwise withhold endpoints and report the failure fraction.
The sample SD may still be reported as a descriptive conditional spread when at
least two draws are valid, even if the interval is withheld; it is not a certified
bootstrap SE in that case. With fewer than two valid draws it is missing.
If the source fit is invalid, report the check as skipped with zero draws
attempted. These two selected-path checks do not estimate bootstrap coverage.
In particular, PFML/Euler plug-in bootstrap is diagnostic and can reproduce
their misspecification bias; it is not an automatic bias correction.

## Denominators and Monte Carlo uncertainty

Retain one record for every planned source path/estimator, including simulation,
point-fit, Hessian and interval failures. Record planned, point-valid and
interval-valid counts separately. Point bias and RMSE condition on valid point
fits, with their availability and MCSEs. Use absolute metrics for mu=0;
relative metrics divide by absolute truth only when nonzero.

Report both interval coverage among valid intervals and operational coverage
over all planned paths. The operational numerator is the number of valid
intervals covering truth; an unavailable interval contributes zero to that
numerator, rather than disappearing from the denominator. These are distinct
estimands. Operational coverage is the primary practical criterion below;
conditional coverage and interval availability accompany it so numerical
exclusions are visible. Use plug-in binomial MCSE and 95% Wilson intervals for
both coverage proportions. MCSEs and Wilson intervals describe Monte Carlo
uncertainty, not uncertainty for an individual fitted path.

Matched CV=0 controls retain the same n, truth and nominal mean gap, not the
same realized or exactly expected span. Report same-estimator excess bias,
RMSE ratios and coverage differences with clear denominators. Within-cell
estimator contrasts use shared-path pairing and common-valid counts.

## Practical classification and checkpoint disposition

Preserve the planned thresholds for theta and sigma: operational coverage below
80% or absolute relative bias above 10%. Before headline data are examined,
define uncertainty-aware classification as follows:

- Supported breach: the 95% Wilson upper bound for operational coverage is
  below 0.80, or the approximate 95% MC interval for relative bias lies entirely
  above +0.10 or entirely below -0.10.
- Within the stated criteria: the Wilson lower bound is at least 0.80 and the
  whole bias MC interval is within [-0.10,+0.10].
- Otherwise unresolved. Missing bias uncertainty prevents a within-criteria
  claim, though a supported operational-coverage breach can still be reported.

These are pointwise classifications over the investigated grid, not simultaneous
familywise guarantees or proof of a sharply located boundary. Publish the
underlying estimates, MCSEs and counts. Mu is not subject to a relative-bias
threshold at zero truth.

The Week 4 exact regular checkpoint retains the historical 94-96% coverage
band as a diagnostic reference. Coverage outside it triggers investigation of
likelihood, optimizer, Hessian, transformations and Monte Carlo variability.
Unresolved implementation defects block headline runs. Once those checks pass,
persistent undercoverage is a property of the specified procedure and must be
reported; do not adjust the interval, seeds or repeats to force nominal coverage.
Coverage inside the band alone does not prove implementation correctness.

## No peeking and amendments

Do not run E1/E2 until Week 4 implementation validation and its disposition are
documented. This task runs only the registered checkpoint and bootstrap checks.
Run every planned headline cell once; do not increase replications, widen grids,
change epsilon, select seeds, drop failures or shift thresholds in response to a
preferred ranking or missing breach. A technical rerun must preserve inputs,
archive the prior output, and document the defect and changed code. Do not
silently reinterpret historical results under a later validity policy.

E3 selected-cell diagnostics, floor sensitivity and any 10,000-replication
follow-up remain labelled exploratory additions under the original plan; they
must not replace complete E1/E2 outputs. Record any compute-budget reduction or
scientific amendment before running affected headline cells, with rationale and
the prior version retained. Any material inference-procedure change after the
checkpoint is versioned, independently validated and disclosed before headlines.
Professor outreach remains deferred. Written theory notes do not certify the
owner's personal S4 learning session, and Tang-Chen full-text reading remains
outstanding.
