# Plan amendment: Week 3 validation and interpretation

**Version 1 — 27 September 2026.** This amendment records decisions made **after
the exploratory Week 3 smoke run**. It responds to W3-07 and the Week 4 handoff in
the *Week 3 Feedback and Validation Checklist*, which identifies baseline commit
`4513a4b`. It is not a retrospective preregistration and does not report new
experimental results or certify that a software correction passed.

The original [OU_master_plan_v2.pdf](../OU_master_plan_v2.pdf) remains unchanged
as the historical plan. For future implementation and interpretation, the
replacement clauses below take precedence over the specified passages of that
PDF. All other scientific design choices remain in force unless a separately
dated amendment changes them. The numerical fixes, tests, provenance and rerun
evidence are recorded separately in
[week3_validation_update.md](week3_validation_update.md). Historical results are
preserved rather than silently reinterpreted under a new validity policy.

## 1. Certification is an investigation checkpoint

**Replaces the interpretation in plan §1.5, the coverage-artifact response in §7,
and the certification wording in §8 and the Week 4 exit condition.**

Keep the planned exact-MLE equidistant checkpoint at
`theta * nominal_mean_gap = 0.5`, `n = 1000` observations, with 2,000 independent
replications. Investigate coverage discrepancies through checks of the
likelihood, optimization, fit validity, observed Hessian, coordinate
transformations, interval construction and Monte Carlo variability. Unresolved
implementation defects block headline experiments.

The former 94–96% band is a diagnostic reference around nominal 95% coverage,
not a target that the implementation must be tuned to attain. If the checks pass
and undercoverage persists, report it as a property of the stated finite-sample
interval procedure. Do not change seeds, delete runs, alter the method or repeat
the checkpoint until the observed coverage enters the band. A scientifically
motivated alternative interval procedure is a separately specified method with
its own documented evaluation, not a replacement result concealed as a repair.

For an actual coverage probability of 0.95 and 2,000 independent replications,
the binomial Monte Carlo standard error is

$$
\sqrt{0.95(1-0.95)/2000} \simeq 0.004873,
$$

or approximately **0.487 percentage points**. This describes simulation
variability conditional on that coverage probability. It does not prove that a
finite-sample Wald procedure has 95% coverage. Week 4 must state the interval
availability and validity denominators before interpreting this checkpoint.

The proposed below-80% exact-MLE coverage at a low-span corner remains a
**hypothesis**, not a required finding. A finite-span bias expansion alone does
not determine coverage: the sampling distribution and interval procedure also
matter. The checkpoint is complete when its implementation checks and scientific
disposition are documented; merely entering the nominal band is insufficient.

## 2. Name the regular-grid comparison precisely

**Replaces the decomposition and isolation claims in plan §§0, 0.5, 1.4 and the
paper-design language in §5.**

Each irregular cell has an equidistant control with the same true parameters,
observation count and **nominal** mean gap. The reported excess bias, RMSE ratio
and later coverage gap compare the same estimator and parameter across those
two designs. Call these *matched nominal-design comparisons* or regular-grid
comparisons.

They do not match every realized observation span. High-CV flooring can also
change expected span, even when its effect on the mean is small. The comparison
therefore measures the combined change in the sampling design, including its
random-span and floor effects; it does not isolate a pure effect of CV or
attribute a discrepancy to span alone. Regular and irregular cells have
independent simulation streams; equal replication indices do not make them
paired observations.

These comparisons are not Aït-Sahalia–Mykland's formal FIML/IOML/PFML information
decomposition. This study does not implement IOML. Retain PFML attribution and
the source-specific distinctions in
[week3_mathematics.md](week3_mathematics.md#2-primary-literature-reading-notes-and-limits).
An expected-span or realized-span matched experiment is optional future work
requiring a separately specified sampling design. It is not a required repair
to the existing smoke run.

## 3. Allow robustness, null and mixed findings

**Replaces the failure-presuming research framing in §0, the instruction in §5
to extend the grid until failure appears, and the outcome-dependent part of
§8's definition of done.**

Ask where the planned estimator and interval procedures meet or fail the stated
practical criteria within the investigated domain. Report all planned cells,
including robust exact-MLE behavior, cases favorable to PFML or Euler, and
mixed results. Do not extend the grid merely to produce a failure, an estimator
ranking, or a visible operating boundary. Retain the already planned E1/E2 grid;
any later extension needs an independent scientific rationale and a dated
designation as additional or exploratory work.

Keep the existing practical criteria for theta and sigma:

- coverage of the specified nominal 95% interval below **80%**; or
- absolute relative bias above **10%**.

These thresholds are unchanged. In Week 4, **before headline data are examined**,
specify how Monte Carlo uncertainty near either threshold affects
classification, including an unresolved category if appropriate. Specify the
fit/interval validity policy and denominators at the same time. Neither point
estimates from 20 smoke paths nor a numerical crossing alone establish a
well-localized operating boundary. Do not move a threshold to make one appear.

An acceptable completed conclusion is “no breach established within the
investigated domain,” accompanied by the tested domain and its Monte Carlo
uncertainty. It is not a claim of universal robustness. The final threshold
table and F7 should state supported classifications and unresolved regions; if
no boundary is established, they should say so. A required “quantitative
threshold statement” must not force the invention of a crossing. Investigating
and transparently reporting persistent finite-sample undercoverage can satisfy
the certification checkpoint once implementation defects have been resolved.

## 4. Align notation, model claims and evidence

**Clarifies plan §§1.2, 1.4, 1.5 and their later uses.**

- `n` means observations; `K = n - 1` means transitions. Nominal span is
  `K * nominal_mean_gap`, so at `n = 250` and `q = theta * nominal_mean_gap = 0.1`,
  nominal dimensionless span is **24.9**, not 25. Realized span is calculated
  from the actual timestamps.
- Replace “point estimates clean by design” with “correctly specified
  conditional transition likelihood; finite-sample bias remains possible.”
  Stationary initialization does not add the stationary initial density to the
  fitted likelihood. All three current estimators condition on the observed
  initial value.
- Use absolute bias and RMSE for `mu = 0`. Relative bias/RMSE and their relative
  MCSEs are undefined at zero truth and must remain missing. Any later
  scale-normalized mu metric needs an explicitly named denominator.
- Optimizer termination, statistical point-fit validity and interval validity
  are distinct. Retain all planned rows and numerical candidates. Summaries
  conditional on valid fits need valid and planned counts plus outcome counts;
  a boundary candidate must not disappear from the record or be replaced by
  truth. A negative Euler coefficient `1 - theta * delta` alone is not evidence
  of an invalid Euler fit.
- The raw-Gamma PFML population limit and a benchmark for floored Gamma gaps are
  different calculations. Finite-sample smoke means are not required to equal
  either asymptotic target, and truth proximity is not a numerical acceptance
  criterion.
- Tang–Chen 2009's **full-text working-depth reading remains outstanding**.
  The existing access limitation and source notes remain in force. Do not
  claim verified correction coefficients or attribute a coverage guarantee to
  that paper before checking its full assumptions and parameter conventions.

## 5. Preserve the exploratory smoke experiment

No wording correction in this amendment changes the data-generating experiment.
Where a verified fitting or validity defect requires a rerun, retain:

| Design component | Preserved choice |
| --- | --- |
| Truth | theta = 1, mu = 0, sigma = 0.5 |
| Estimators | Exact, realized-mean-gap PFML, Euler |
| Nominal grid | q in {0.1, 0.5, 1}; CV in {0, 1, 2} |
| Size | 250 observations, 249 transitions, 20 replications per cell |
| Random streams | Root seed 20260910, scientific cell IDs and gap/path stream mapping |
| Gap law | Existing regular/Gamma generation; literal floor of 1e-6 times nominal mean gap for CV >= 1.5 |
| Path and objective | Exact transition simulation, stationary initial draw, likelihood conditional on the initial observation |
| Planned records | 9 cells, 180 path/spacing records, 540 fit rows, 81 parameter-summary rows |

Preserve the historical `artifacts/week3_smoke/` outputs. Store corrected-code
results separately and reconcile them on the planned keys. Verify actual
timestamp and observation arrays, not only seeds, when asserting the same
inputs. Capture code/dependency provenance and disclose a deliberate replay
under changed code. Helper-only or wording-only changes do not themselves
require new simulation; relevant acceptance tests still apply.

The smoke floor value remains exploratory. Its preservation in a regression
run does not retroactively freeze it for E1/E2. Week 4 must explicitly freeze
epsilon with the other headline design choices.

## 6. Week 4 handoff: normalization and covariance

Data-derived normalization may improve optimization without changing the
statistical objective. Record the normalization and preserve its coordinate
maps for the later Hessian and delta-method work. Do not use the simulation
truth for numerical scaling.

For a fixed observed dataset, write normalized observations and time as

$$
y=(x-a)/b,\qquad u=t/h,\qquad b>0,\ h>0.
$$

If an optimizer fits the OU parameters in these coordinates, map back with

$$
\theta_x=\theta_u/h,\qquad
\mu_x=a+b\mu_u,\qquad
\sigma_x=b\sigma_u/\sqrt h,
\qquad L_x=L_u+K\log b.
$$

The NLL correction is the state-density Jacobian over the K conditional
transitions; conditioning on timestamps adds no time-density term. At a fixed
dataset the numerical normalization constants are held fixed when computing
parameter derivatives. They do not create extra estimated model parameters.

If the estimated covariance is in normalized working coordinates
`eta_u = (log(theta_u), mu_u, log(sigma_u))`, the map to original natural
parameters has Jacobian

$$
J=\operatorname{diag}(\widehat\theta_x,\ b,\ \widehat\sigma_x),
\qquad
\widehat{\operatorname{Cov}}(\widehat\theta_x,\widehat\mu_x,
\widehat\sigma_x)=J\widehat{\operatorname{Cov}}(\widehat\eta_u)J^T.
$$

If covariance is first transformed to the **original** working coordinates,
the working-to-natural Jacobian is instead
`diag(theta_x, 1, sigma_x)` as in the original plan. The two routes must agree;
using the latter matrix directly on a normalized-coordinate covariance omits
the mu scaling. Record the actual coordinates used by the implementation and
test state/time transformations of the covariance and resulting intervals.
These formulas are a handoff requirement, **not a claim that covariance or
interval construction is implemented or validated**.

Week 4 retains observed-Hessian conditioning checks, explicit interval validity
and coverage denominators, Wald/delta/bootstrap work, the certification
investigation, and a committed preregistration before E1/E2. Nonidentified or
unresolved boundary candidates must not acquire ordinary finite Wald intervals
merely because an optimizer terminated successfully.

Before Weeks 5–6, extend the smoke-only configuration/plot limits to support the
already planned E2 grid, including its four q values; do not reduce that grid to
fit the current interface. Keep detailed anomaly and epsilon-sensitivity work
in Week 7/E3, with selected reruns labelled as diagnostics and the complete
headline record preserved. Clustering remains behind the Week 8 gate.
Professor outreach remains deferred by the owner's decision.

## Adoption and evidence boundaries

This document adopts the W3-07 interpretation corrections and records their
post-exploratory timing. It does not close W3-01 through W3-06 by assertion,
change the original results, establish nominal coverage, complete the deferred
reading, or replace Week 4 preregistration. Consult the separate
[validation update](week3_validation_update.md) for the actual implementation
and acceptance evidence.
