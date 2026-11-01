# Week 10 explanatory supplement

The completed analysis addresses both requested extensions: controlling the
observed span and explaining the PFML mean precision mechanism. It amends the
old figure-freeze schedule explicitly and preserves the completed original E4
and Week 9 packages.

Start with the [critical analytical review](week10_analytical_review.md),
the revised [mathematical derivation](week10_mechanism.md), and
[Draft 3](../manuscript/week10_draft3.md). The original
[generated results and interpretation](../artifacts/week10_supplement/interpretation/results.md)
and Draft 2 remain historical, unchanged sealed evidence; their moment and
coverage qualifications are corrected in the new review. The
[protocol](week10_supplement_protocol.md) was recorded before new supplemental
effect tables or fitted results. The study remains exploratory because its
questions were selected after observing original E4.

## Scientific outcome

At a common observed span T=499.5, empirical PFML mean RMSE is higher under clustered
sampling at r=5 and r=20, and its working Wald intervals have lower coverage
than the independent parent-marginal controls. The full pointwise uncertainty
is in the generated report and [matched contrast table](../artifacts/week10_supplement/analysis/matched_iid_comparisons.csv).
Span variation is therefore not necessary for the effect under this design.
The experiment changes N and applies a shared endpoint convention; it does not
identify a universal direct effect of dependence with every other feature held
fixed.

The new covariance calculation starts from the actual PFML regression mean,
including its endpoint correction. Matching one-gap marginals fixes the
adjacent-observation correlation and population working variance, but leaves
multi-gap OU correlations different. The resulting net covariance contribution
predicts the fixed-coefficient diagnostic's spread and approximates the saved
fitted central spread. Exact finite-design results for a fixed-
coefficient linear diagnostic are separated from approximations to the fitted
estimator and its coverage. The critical review identifies a boundary pole in
the ideal fitted ratio, so empirical RMSE and its historical moment-based
intervals are not justified ideal population-risk confidence statements.
Coverage indicators remain bounded and the coverage finding is retained.
No new robust interval or general dependent-gap
CLT is claimed.

Saved-record span standardization also retains the effect. Its reported
uncertainty conditions on support and weights, and its unequal tail exclusions
preclude a causal percentage-of-effect decomposition.

## Completed evidence

| Item | Actual scope |
| --- | --- |
| Fixed design and source checks | Seven cells; separate development/source-only/full seeds; 1,536 source-only schedules and 22 checks |
| Development fitted run | 28 paths / 84 fits; all valid; every source reconstructed and 42 endpoint fits replayed |
| Full fixed-horizon study | 14,000 paths / 42,000 fits; all point estimates and intervals available; zero pipeline errors |
| Full source/fit audit | Every schedule and OU path reconstructed; 14,000 independent PFML mean identities checked; 42 endpoint fits replayed exactly |
| Original-span analysis | All 14,000 original sources / 42,000 fits retained; 69 original files bound; 216 checks recover original comparisons |
| Mechanism validation | Independent matrix integration, finite-count covariance, dense conditional covariance and fitted-versus-fixed coefficient diagnostics |
| Independent result review | 1,092 numerical checks of raw arithmetic, contrasts, predictions and source counts; no discrepancy |
| Publication outputs | Complete tables, generated numerical claims and interpretation, three PNG/PDF figure pairs, Draft 2 |
| Completion gate | Final tests and source preservation, actual figure inspection, bound final validation, guarded seal and post-seal readback |

The authoritative actual test totals and checks are in
[final validation](../artifacts/week10_supplement/evidence/final_validation.json).
The [completion seal](../artifacts/week10_supplement/completion_seal.json) binds
the scientific payload. The
[visual review](../artifacts/week10_supplement/evidence/visual_review.json)
identifies the inspected original PNGs and PDF rasterizations. Source and
analysis freezes record actual file hashes and the uncommitted source status;
they do not claim that the current Git head already contained the new code.

## Preserved failures and limits

Initial test collection in the old renderer environment failed because it
lacked PyYAML. One source-policy syntax error was fixed before fitted
development observations. A new isolated pinned environment supplied the
required dependency; network access initially failed in the restricted
environment and the authorized dependency download then succeeded.

Independent guard tests exposed insufficient completion checks and an overly
broad filename exclusion. These were repaired before sealing, with the
original failed checks and dispositions retained. No estimator, scientific
threshold, seed or research replication count was changed. The source's
possible insufficient-observation outcome remains an explicit source failure
under the existing Sample contract; none occurred in the full study.

The mathematical results are specialized covariance calculations under the
specified stationary noise-free Gaussian OU and exogenous schedules. The
sampling origin is event-Palm, the horizon includes a deterministic endpoint,
and observation counts remain random. Inference remains pointwise with
model-based PFML/Euler standard errors. Hardware-independent bit replay,
general novelty, academic acceptance and personal theory mastery are not
certified by these artifacts.

## Manuscript handoff

Draft 3 revises the mechanism and focused results using the additive review;
the original Draft 2 and 250-file seal remain unchanged. Week 11 still
assembles the complete E1–E4 Results, practical guidance,
related work and discussion and obtains academic feedback. Week 12 handles
revisions, final abstract/layout, literature refresh and submission preparation.
The supplement strengthens the paper's explanation without completing those
later authoring stages. Follow [the reproduction instructions](week10_reproduction.md)
for saved-data checks before considering any new simulation.
