# Week 3 feedback assessment and validation update

27 September 2026. Baseline: `4513a4b1cafa88560bffef5e267b53fae2ed3fc1`.

The feedback is valid. We reproduced its reported numerical defects and fixed
them without changing the scientific smoke design. The corrected run uses the
same 180 timestamp/observation pairs and retains all 540 planned fits. Every
corrected smoke fit passes the new interior-validity rule and agrees with an
independent reference within the tolerances specified before patch validation.
The original [smoke outputs](../artifacts/week3_smoke/) remain byte-for-byte
unchanged; [updated results](../artifacts/week3_validation_rerun/) and
[validation evidence](../artifacts/week3_validation_evidence/) are separate.

This update repairs numerical fitting and reporting. It does not establish a
coverage boundary or overturn the original smoke findings. Its effect on the
ordinary smoke estimates is small; the larger defects arise in the reported
unit transformations and out-of-grid stress fixtures.

## Feedback disposition

| Item | Decision and implemented change | Acceptance evidence |
| --- | --- | --- |
| W3-01 | Valid correctness defect. Fit in data-derived state/time units; restore all parameters and the conditional NLL Jacobian. | Fitted state scales 1e-6, 1, 1e6 with zero/substantial offsets; time factors 1e-3, 1, 1e3; regular and irregular paths including the reported path. |
| W3-02 | Valid classification defect. Separate solver termination, interior validity, boundary/suspected boundary, degeneracy and numerical failure; preserve candidates. | Negative-AR-slope stress fixture excluded for exact/PFML but remains valid for Euler; constant/lost-resolution/noiseless paths and solver exceptions tested. |
| W3-03 | Valid helper defect. Guard the finite-data and admissible AR1-to-OU domain. | Invalid slope example now raises ValueError; interior reference and bad-input tests pass. |
| W3-04 | Valid helper defect. Use log1p(z)/z with its continuous limit at zero for raw Gamma PFML targets. | CV=0, 1e-8, 1e-4, 0.5, 1, 2 and input/overflow cases tested. |
| W3-05 | Valid missing acceptance evidence. Add independent exact profile-score, PFML OLS and Euler WLS references. | All 540 corrected fits and 81 transformed-fixture comparisons pass fixed parameter/NLL gates. |
| W3-06 | Valid integration and reporting improvements. One validity rule across exports/counts/plots; strict replay provenance; input hashes; MCSE bars; shared-path estimator comparisons. | Every outcome injected into hand-checkable integration tests; all-invalid cells, mismatches, raw integrity and paired denominators tested. |
| W3-07 | Valid plan corrections, adopted after exploration. Certification is an investigation checkpoint; nominal twins and null/mixed findings are stated precisely. | [Versioned plan amendment](plan_amendment_week3_validation.md); original PDF preserved. |

All seven items are closed for their stated Week 3 scope. The unchanged Week 4
and later work is listed below. No item was accepted merely because a revised
estimate moved closer to the simulation truth.

## What was reproduced

On the original q=0.5, CV=2, replication 0 path, multiplying observations by
1e-6 changed the old exact theta estimate from 1.260668750 to 1.500927755 and
the old Euler estimate from 0.421122437 to 0.320483405, despite success flags.
Their mapped-back NLLs worsened by about 0.449175 and 0.348854 respectively.
PFML also showed a smaller NLL deterioration. See
[reproduced_defects.json](../artifacts/week3_validation_evidence/reproduced_defects.json).

The old 250-observation, gap-5, seed-1 stress fixture had unrestricted AR slope
-0.108111049. Exact/PFML returned a finite theta near 4.159644 with success,
although the regular OU likelihood supremum approached phi=0. That fixture is
outside the original smoke grid. It now receives `boundary`, with a finite
diagnostic candidate retained and `valid_for_point_summary=False`.

The supplied increasing-path AR1 helper example returned negative theta, and
`pfml_limit(1, .5, .5, cv=1e-8)` returned zeros. Both defects were reproduced.
The repaired helper rejects the inadmissible AR slope; the small-CV target is
approximately (1, .5).

## Solver and mathematical contract

The strict Week 3 interface now solves the same conditional objectives using
normalized regression/profile methods. It uses no simulation truth for scaling:
with a the observed state mean, b its positive RMS centered scale, and m the
realized mean gap, set y=(x-a)/b and u=t/m. Recover

```text
theta_x = theta_y / m
mu_x    = a + b * mu_y
sigma_x = b * sigma_y / sqrt(m)
NLL_x   = NLL_y + (n - 1) * log(b)
```

PFML and regular-grid exact fitting use conditional AR1 OLS with an admissible
OU map. Euler uses weighted regression of differences; its transition
coefficient may be negative without violating theta>0. Irregular exact fitting
profiles drift and diffusion variance at each theta, expands the dimensionless
search range, refines sampled minima, and compares endpoint and limiting
objectives. `profile_search_ranges`, expansion counts, objective evidence and
solver messages are saved. Failed polishing retains the best finite candidate
and marks it invalid. A finite profile search is not a proof of a global optimum
for every possible irregular path; unresolved tails remain `boundary_suspected`.

The old YAML start-count settings are accepted and recorded as requested
controls for compatibility. These new solvers do not perform the old 3+5
L-BFGS-B schedule: actual `attempted_starts` is zero and `retried` is false.
`solver` identifies the algorithm. The switch changes the numerical solution
method, not the estimator definition. Existing legacy `ou.fit`, its moment
initialization, gap-clamping objectives and simulator remain unchanged for
Week 2 reproduction. The strict Week 3 path never calls those legacy objectives;
tests enforce this separation. Only the two public analytic helpers in `ou.py`
were repaired. Historical notebook outputs were not regenerated.

The [mathematical notes](week3_mathematics.md) and code docstrings explain the
profile, normalization, boundary, regression and population-limit formulas.
The [SciPy 1.17 L-BFGS-B documentation](https://docs.scipy.org/doc/scipy-1.17.0/reference/optimize.minimize-lbfgsb.html)
confirms the absolute forward-difference step used when no Jacobian is supplied;
this supports investigating scale sensitivity, without establishing it as the
only cause. [NumPy log1p](https://numpy.org/doc/stable/reference/generated/numpy.log1p.html)
documents the numerically stable logarithm used by the repaired helper. The
raw-Gamma target is not a floored-law benchmark or a finite-sample correction.

## Validity and reporting

For new records, `success` aliases `optimizer_success`. Point metrics instead
require the explicitly classified `valid_interior` outcome, successful numerical
termination, finite parameters/NLL and positive theta/sigma. The corresponding
`valid_for_point_summary` flag must agree; contradictory records raise an error.
The policy is named `classified_interior_v2`.

Raw rows retain every candidate and reason. `failures.csv` contains all
point-summary exclusions, including boundaries; `boundary_fits.csv` isolates
confirmed/suspected boundaries. Manifest totals and summaries report planned,
valid, excluded, numerical-success and individual outcome counts. Historical
rows without classifications retain their original `legacy_success_finite_v1`
policy; they are not retroactively certified as interior estimates.

Bias and RMSE remain conditional on the stated validity policy. Approximate
plus/minus 1.96 MCSE bars show simulation uncertainty for these Monte Carlo
summaries, not a confidence interval for one path's OU parameter. With only
20 paths per cell they warrant cautious interpretation. Relative mu metrics
remain undefined at zero truth. All Week 3 coverage fields remain uncomputed.

The RMSE-ratio figure explicitly compares each estimator with its own CV=0
nominal-design control. New `estimator_comparisons.csv` instead compares methods
on their shared paths within a cell, reporting bias/MSE differences, paired
MCSEs, planned pairs and common-valid denominators. It does not pair independent
CV cells. All three estimator pairs and all three parameters are reported.

## Validation and effect on the original smoke results

The full suite passes **300 tests**. The unmodified baseline suite also passes
its original 141 tests in the present environment. Python 3.12.14 and the five
numerical/plot/YAML dependency versions match the original run. The recorded
Linux kernel differs; no cross-platform byte-identity promise is made.

Acceptance gates were fixed before patch validation: parameter rtol=1e-4 and
atol=1e-6 in observed-data normalized coordinates; absolute NLL discrepancy at
most 1e-7 + 1e-9 times the absolute reference NLL. They were not loosened to
preserve a success count. The exact reference uses roots of an analytic profile
score in stationary-variance coordinates, independently of production scalar
minimization. PFML/Euler references use separate least-squares implementations.

| Recorded check | Result |
| --- | --- |
| Original configuration and archived output hashes | Unchanged |
| Timestamp and observation array identities | All 180 pairs match |
| Planned/recorded fit keys | All 540 match |
| Corrected valid interior fits | 540; no excluded outcomes in this smoke grid |
| Parameter-summary / paired-comparison rows | 81 / 81 |
| Largest corrected NLL difference from independent reference | 1.60e-12 or less |
| Largest corrected parameter error as a fraction of its allowed tolerance | 0.00582 or less |
| Largest absolute before/after changes in theta / mu / sigma | 5.13e-6 / 2.03e-6 / 1.46e-7 |
| Fit-level transformed-reference checks | 81; maximum NLL discrepancy 3.84e-11 or less |

The full [before/after table](../artifacts/week3_validation_evidence/before_after.csv)
records every key, all parameter changes, NLLs, solver/outcome diagnostics and
reference gates. Two old rows exceeded the new normalized-parameter engineering
gate: q=0.1, CV=2, replication 5 for exact and Euler. Their gate fractions were
2.10 and 1.56, while their old objective gaps were only 3.78e-11 and 1.20e-10.
These are small parameter-solution discrepancies in a nearly flat objective,
not large likelihood errors or evidence that the original smoke conclusions
were invalid. Their corrected fits pass both gates. All other old smoke rows
already passed these comparison gates.

The canonical run is recorded in `run_metadata.json`. During figure QA, an
overlapping ratio-axis label was shortened and the same configuration was run
again. Raw fit records were byte-identical across those two executions; the
pre-QA manifest is retained in the evidence directory. No alternate seed, grid,
truth, validity threshold or preferred estimator ranking was selected. All four
final figures were visually checked.

## Reproduce and audit

From the repository root with the pinned requirements installed:

```bash
python -m pytest -q
python -m ou_irregular.runner run --config configs/week3_smoke.yaml --output artifacts/my_validation_run
python -m ou_irregular.runner replot --run-dir artifacts/my_validation_run
python scripts/validate_week3_update.py --baseline-run artifacts/week3_smoke --updated-run artifacts/my_validation_run --evidence artifacts/week3_validation_evidence
```

The audit command preserves baseline files, verifies the frozen 180 input hashes,
reconciles every fit against independent references and rewrites only the named
derived comparison evidence. Its canonical output is
[reconciliation.json](../artifacts/week3_validation_evidence/reconciliation.json).
Source hashes and dependency-version fingerprints accompany the run. Dependency
fingerprints identify versions, not every installed binary file.

Strict replay rejects changed scientific package source, Python or dependency
versions before fitting. Use a full cell ID from `cells.json`:

```bash
python -m ou_irregular.runner replay --run-dir artifacts/my_validation_run --cell-id FULL_CELL_ID --rep 0
```

A deliberate patched-code replay of historical data must be explicit and leave
a separate provenance record:

```bash
python -m ou_irregular.runner replay --run-dir artifacts/week3_smoke --cell-id FULL_CELL_ID --rep 0 --allow-provenance-mismatch --output artifacts/my_patched_replay
python -m ou_irregular.runner replot --run-dir artifacts/week3_smoke --output artifacts/my_historical_replot
```

The old run has no input-array hashes; patched replay labels that limitation
rather than claiming historical array identity. For this validation update,
arrays were reconstructed from the frozen baseline implementation, hashed before
patch fitting, and verified against the updated run. The baseline commit, source
and environment record is in `baseline_provenance.json`. Replotting an old run
preserves its old validity semantics and requires separate output when provenance
differs. Neither command overwrites the archived smoke results by default.

## Week 4 handoff and remaining limits

Proceed with observed-Hessian and coordinate-transformation checks, interval
validity, Wald/spot-check procedures, checkpoint investigation and preregistration.
The amendment supplies the normalization covariance map; intervals themselves
are not implemented here. Preserve the planned headline grid and thresholds,
allow robust or mixed findings, and define uncertainty handling near thresholds
before inspecting headline results. Smoke-only grid limits must be expanded for
E2 later, not used to shrink its planned design.

Tang–Chen full-text reading remains outstanding. No unverified correction
coefficient is used. Detailed anomaly/floor-sensitivity work remains E3; new
models, span-matched designs and clustering are not required Week 3 repairs.
Professor outreach remains deferred.
