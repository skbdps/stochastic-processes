# Sources for Week 10 Draft 1

This map supports the methods manuscript at the completed research baseline
`4bddb527f3c8bb88bc1ebec28c8818e03dc251aa`. Repository links below point to
unchanged source files in this release; the recorded hashes in
[validation.json](../artifacts/week10_draft1/validation.json) bind the exact
bytes checked for Draft 1. The dated preregistration and exploratory protocols
govern their respective experiments. Older aspirational master-plan wording
does not replace those adopted definitions.

| Draft passage | Authoritative source and scope |
| --- | --- |
| Introduction and attribution | [Week 9 literature comparison](../docs/week9_literature_comparison.md), plus the primary papers in the draft references. The initial citation check is recorded in the [review record](../artifacts/week10_draft1/review.md). This is not an exhaustive pre-submission novelty search. |
| 2.1 Model and simulation | [Week 3 mathematical definitions](../docs/week3_mathematics.md), [exact OU simulator](../ou_irregular/ou.py), [preregistration](../preregistration.md). |
| 2.2 Estimator objectives and computation | [Strict estimators](../ou_irregular/estimators.py), [Week 4 estimator adapter](../ou_irregular/inference/adapters.py). Inspect the implemented AR(1), Euler WLS and scalar profile methods, rather than the historical multistart optimizer proposal. |
| 2.3 Normalization and intervals | [Preregistration](../preregistration.md), [Week 4 mathematics](../docs/week4_mathematics.md), [coordinates](../ou_irregular/inference/coordinates.py), [Hessian checks](../ou_irregular/inference/hessian.py), [Wald intervals](../ou_irregular/inference/wald.py). |
| 2.4 Expected-objective targets and calibration | [Week 9 benchmark note](../docs/week9_benchmarks.md), [prospective E4 population targets](../docs/week8_population_targets.md). The used-gap law includes its floor; covariance theory is not inferred from targets alone. |
| 3.1 Core grid, floor and seeds | [Preregistration](../preregistration.md), [E1 mathematics](../docs/week5_mathematics.md), [E2 mathematics](../docs/week6_mathematics.md). |
| 3.2 Exploratory E3 selection and design | [E3 protocol](../docs/week7_exploratory_protocol.md), [Week 7 completion and accounting](../docs/week7_progress.md). The follow-up and paired floor studies have different independent units. |
| 3.3 MMPP and full-marginal controls | [E4 protocol](../docs/week8_exploratory_protocol.md), [E4 mathematics](../docs/week8_mathematics.md), [source contract](../docs/week8_source.md). The protocol also links the primary MAP/phase-type references behind the specialized construction. |
| 3.4 Metrics and Monte Carlo uncertainty | [E1 mathematics](../docs/week5_mathematics.md), [summary calculations](../ou_irregular/experiments/summaries.py), [contrasts](../ou_irregular/experiments/comparisons.py). Different-cell comparisons are independent; within-path methods and floor variants use their declared pairing. |
| 3.5 Practical classifications | [Preregistration](../preregistration.md), [adopted interpretation amendment](../docs/plan_amendment_week3_validation.md), [Week 9 protocol](../docs/week9_protocol.md). |
| 3.6 Audit, bootstrap and reproduction scope | [Week 4 progress](../docs/week4_progress.md), [Week 8 completion](../docs/week8_completion.md), [Week 9 completion](../docs/week9_completion.md), [Week 9 reproduction](../docs/week9_reproduction.md). These reports preserve original failures and the limits of their actual checks. |

## Boundary between methods and later results

Draft 1's counts, grids, seeds and thresholds describe fixed study designs or
retained accounting. Performance estimates for the later Results section must
come from the [generated claim ledger](../artifacts/week9_publication/claims.json)
and bound tables. No new number should be typed into a result statement merely
because it appears plausible or agrees with an anticipated outcome.

The [Week 9 package state](../artifacts/week9_publication/package_state.json)
and [manifest](../artifacts/week9_publication/manifest.json) retain their
original identities. The draft is a new document and is not inserted into the
sealed package or retrospectively described as preregistered.
