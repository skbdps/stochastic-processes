# Week 5: E1 implementation, validation and results record

Prepared 29 October 2026 UTC. This record separates implementation evidence from
completed headline computations. The governing protocol remains the published
preregistration and dated Week 3 amendment; the implementation plan makes no
scientific amendment.

**Completion:** the full registered E1 is complete, independently audited and
visually reviewed. It contains 14,000 source paths, 42,000 point/interval records
and all 63 parameter-summary rows. All 42,000 point fits and intervals are valid.
No E1 scientific rerun, resampling, threshold change or estimator substitution
was made. Week 6 implementation begins only after publication of this record.

## Status at the pre-E1 freeze

The published handoff `cf9ecb6043f029b679606d737b825740ad029c54` was verified
against the live `weeks5-8-sequential` branch before work. Main was
`4e35ba23d79ec6dc16f10035df27729817e23678`. The handoff was an unfinished
Week 5 checkpoint. E1/E2 headline data had not run. No Week 6–8 implementation
was started while completing Week 5.

The new modules were implemented and tested sequentially: comparison formulas;
explicit provenance; immutable per-cell storage; execution/failure retention;
table export; independent artifact audit; figures; and the week-specific CLI and
completion workflow. Independent reviews covered both formulas and failure
handling. The detailed contracts and derivations are in
[week5_architecture.md](week5_architecture.md) and
[week5_mathematics.md](week5_mathematics.md).

## Review dispositions

- The initial design/summary modules were independently reviewed before reuse.
  A generator-consumption defect was repaired, complete planned-key validation
  was strengthened, and nonfinite/Boolean candidate data and overflow fail closed.
- Independent regular-reference and paired method formulas were checked against
  hand calculations, including common-valid subsets and zero-truth mu handling.
- Resume checks bind design, procedure, source files and the numerical runtime.
  Crash-orphan recovery preserves the original partition boundary, even if a
  later requested partition size differs.
- The auditor recalculates every table independently, checks every saved fine
  Hessian/covariance and fixed numerical gates, and replays replications 0 and
  R-1 in every cell. Full paths and coarse Hessians are not stored for every
  replication; their exact recomputation is confined to the fixed replay subset.
- Final integration review found that a completion seal must bind scientific
  root metadata as well as raw/table files. Regression tests now reject changed
  procedure, pipeline-error counts, overrides, cell accounting, code-commit
  metadata and unexpected nested files. Legitimate execution/recovery history
  remains recorded without preventing byte-identical reuse of an existing audit.
- The mu figure label was clarified to “Bias (original units)”; it reports signed
  bias, not the magnitude of bias. Dense CV labels were rotated for readability.

## Development validation and historical preservation

The first complete seven-cell validation is retained at
`artifacts/week5_validation/`. It used R=8 and the separate validation seed
21261001: 56 paths and 168 fit records. It finished its numerical audit but was
not sealed as completed after final integration review identified completion-state
and figure-label changes. Its evidence is preserved unchanged. The second
validation, retained at `artifacts/week5_validation_layout_review/`, verified
the final guards but exposed a caption-margin problem in visual inspection.
`artifacts/week5_validation_final/` uses identical validation addresses with
the corrected figure margins. All 14 raw partitions are byte-identical across
these three development runs, as recorded in
`artifacts/week5_evidence/development_replay_comparison.json`. Changed development
source versions are preserved under that evidence directory. The final validation
passed the independent audit and all five figures were visually inspected before
its completion seal. These are development checks, not E1 headline
results or coverage certification. The subsequent headline run will use seed
20261001 and exactly R=2000 without changing the registered numerical procedure.

The original 540-test checkpoint suite passed before implementation. The complete
pre-E1 suite passed 624 tests; its output is saved at
`artifacts/week5_evidence/tests_pre_e1.txt`. The figure-margin-only change was
subsequently verified through focused deterministic-plot tests and visual review.
Fault-injection tests exercise source, estimator, inference, persistence and
post-audit figure failures without resampling a scientific failure.

The historical Week 4 audit was run in a detached checkout at its published
Week 4 revision `0a465d6b42195f4136a81c0b0cbd4d607d833893`, rather than
weakening its exact recursive source-map rule. It passed all 2,000 source replays,
3,000 saved bootstrap-draw checks and summary arithmetic. The separate report is
`artifacts/week5_evidence/historical_week4_audit.json`. The preservation manifest
records 72 baseline scientific source/artifact files checked against main, with
no changes. The original master-plan PDF was read as a historical attachment and
was not restored to the repository.

The runtime uses Python 3.12.14 and the pinned NumPy, SciPy, pandas, matplotlib
and PyYAML versions. Notebook-only packages are unavailable in the current
environment, and the configured shell proxy prevented installing them. Their
absence is recorded; notebook execution is not claimed.

## Published freeze and reproducibility detail

The exact code tree `e9751ad40a408870f00980e210653799800885c3` was frozen
locally as `da519f589281e343b1748f5e4f7a68a51ceb1d7c` and published by the
GitHub connector as `559bf903ecba1bf7ccbb72b78276a70e1c908b32`. Both commits
have the same parent and identical files; their author/committer metadata,
timestamps and messages differ. The guarded branch update returned success
before the E1 run invocation.

An initial attempt to import the connector commit locally failed because its
original timezone was +0530. Execution was inadvertently launched before that
local-ref synchronization finished, so its provenance faithfully records the
local freeze `da519...`. This bookkeeping issue did not alter any scientific
source or input. The root manifest has not been relabelled. All 35 recorded
dependency digests match both commits. The exact local commit bytes and mapping
are preserved in `artifacts/week5_evidence/`.

A focused independent fresh-clone check verified that
`scripts/restore_recorded_commit.py` restores the otherwise absent original
commit object, resolves all 35 dependency hashes and leaves refs and copied run
files unchanged. This preserves the existing auditor's exact commit-file check.
Run the following before auditing the published E1 in a fresh clone:

```bash
python scripts/restore_recorded_commit.py --mapping artifacts/week5_evidence/freeze_publication_mapping.json --run-dir artifacts/week5_e1
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m ou_irregular.week5.runner audit --run-dir artifacts/week5_e1 --output artifacts/my_e1_audit.json
```

## Headline E1 results

The design is q=0.5, n=1000 observations (999 transitions), truth (1,0,0.5),
CV={0,0.25,0.5,1,1.5,2,3}, R=2000 per cell, seed 20261001. Exact stationary OU
sources and conditional likelihoods follow the registration. The one raw
execution used seven workers and 250-path immutable partitions; it took 57.55
seconds, excluding postprocessing and audit. Its 56 partitions retain every
planned source/estimator record. The final status is `completed`.

The independent audit passed all planned-key checks, raw-to-table arithmetic,
fine-Hessian/covariance identities, recorded numerical gates and the fixed replay
of 14 source paths / 42 fits. It verified 63 parameter rows, 63 independent
regular-reference rows, 63 paired-comparison rows, seven spacing rows and 21
rows each for outcomes and numerical diagnostics. There are no pipeline errors,
invalid point fits or unavailable intervals. The failed-fit export is retained
with its schema and zero data rows. Operational and conditional coverage are
equal here because all intervals are available, not because their definitions
were combined. Each per-cell/method availability is 2000/2000; its 95% Wilson
interval is approximately [99.8083%,100%].

An independent mathematical/results review also recomputed 71 selected
raw-record arithmetic and classification identities, all matching. The largest
recorded condition number is 261.31, relative Hessian discrepancy 3.88e-7,
whitened discrepancy 3.02e-6 and standardized score 5.42e-6. The observed major
coverage failures therefore occur with accepted numerical inference, rather
than through unavailable intervals. This does not certify unrun E2 cells.

The full machine-readable evidence is under
[artifacts/week5_e1](../artifacts/week5_e1/). The following table reports coverage
and truth-relative bias in percent. Each coverage denominator is 2000. W means
within criteria, B supported breach, U unresolved; these pointwise labels use
the registered Monte Carlo uncertainty rules, not point estimates alone.

| CV | Estimator | Theta coverage | Sigma coverage | Theta relative bias | Sigma relative bias | Theta / sigma class |
| ---: | --- | ---: | ---: | ---: | ---: | --- |
| 0 | exact | 95.40% | 95.05% | +1.03% | +0.12% | W / W |
| 0 | pfml | 95.40% | 95.05% | +1.03% | +0.12% | W / W |
| 0 | euler | 2.30% | 0.00% | -20.78% | -20.57% | B / B |
| 0.25 | exact | 95.20% | 95.15% | +1.23% | +0.20% | W / W |
| 0.25 | pfml | 95.25% | 94.15% | -0.28% | -0.56% | W / W |
| 0.25 | euler | 1.70% | 0.00% | -21.58% | -20.20% | B / B |
| 0.5 | exact | 95.80% | 95.30% | +1.10% | +0.16% | W / W |
| 0.5 | pfml | 89.75% | 79.45% | -4.64% | -2.72% | W / U |
| 0.5 | euler | 0.25% | 0.00% | -24.30% | -19.27% | B / B |
| 1 | exact | 95.85% | 95.10% | +1.07% | +0.10% | W / W |
| 1 | pfml | 30.45% | 5.75% | -18.07% | -9.89% | B / B |
| 1 | euler | 0.00% | 0.00% | -32.86% | -16.38% | B / B |
| 1.5 | exact | 96.35% | 95.45% | +1.51% | +0.10% | W / W |
| 1.5 | pfml | 1.10% | 0.05% | -32.06% | -18.16% | B / B |
| 1.5 | euler | 0.00% | 0.00% | -42.43% | -13.39% | B / B |
| 2 | exact | 95.70% | 95.25% | +1.46% | +0.08% | W / W |
| 2 | pfml | 0.00% | 0.00% | -43.98% | -25.75% | B / B |
| 2 | euler | 0.00% | 0.10% | -51.22% | -10.87% | B / B |
| 3 | exact | 95.20% | 94.35% | +2.39% | +0.00% | W / W |
| 3 | pfml | 0.00% | 0.00% | -61.01% | -38.46% | B / B |
| 3 | euler | 0.00% | 6.90% | -64.47% | -7.50% | B / B |

### Interpretation and uncertainty

Exact theta/sigma estimates meet the registered practical criteria at all seven
CV values. Exact coverage ranges are 95.20–96.35% for theta, 94.15–95.60% for mu
and 94.35–95.45% for sigma. These are observed finite-sample outcomes, not proof
of exact 95% calibration. Exact theta bias is positive: approximately 1.03% at
regular spacing and 2.39% at CV=3. At CV=3 its relative-bias Monte Carlo interval
is [1.89%,2.90%], still inside the practical criterion.

PFML theta is within criteria through CV=0.5 and has supported breaches at
CV>=1. PFML sigma is within criteria at CV=0 and 0.25, unresolved at CV=0.5,
and has supported breaches at CV>=1. The unresolved sigma cell covers
1589/2000=79.45%; its Wilson interval [77.62%,81.16%] crosses 80%. It must not
be labelled a supported coverage breach. At CV=1, sigma coverage is only
115/2000=5.75%, so the supported breach follows from coverage even though the
bias interval straddles the -10% threshold.

Euler theta and sigma have supported breaches at every E1 CV, including the
regular reference. Deterioration is not monotone for every metric: Euler sigma's
relative bias improves from -20.57% at CV=0 to -7.50% at CV=3, but coverage at
CV=3 is only 138/2000=6.90%. A small bias alone does not establish satisfactory
interval performance. PFML/Euler intervals continue to use their registered
model-based inverse Hessians; no robust covariance or interval substitution
was introduced.

Zero observed coverage is retained with a Wilson upper bound of approximately
0.1917% at R=2000. Its zero plug-in coverage MCSE does not prove that the
population coverage probability is zero.

At CV=3, exact RMSE relative to its independent regular reference is 1.415 for
theta (95% log-delta Monte Carlo interval 1.349–1.485), 1.365 for mu
(1.307–1.427), and 0.863 for sigma (0.825–0.903). Thus this irregularity sweep
increases theta/mu error while decreasing exact sigma RMSE at this endpoint.
References match nominal n/truth/mean gap, not realized span. Ratios sharing a
regular reference are correlated and do not form simultaneous inference.
Euler sigma's ratio of 0.378 compares with its own poor regular baseline and
does not make its absolute performance acceptable.

Exact does not uniformly minimize point RMSE. At CV=0.25, PFML theta RMSE is
0.08464 versus exact 0.08652. On the same 2000 common-valid paths, the paired
exact/PFML RMSE ratio is 1.0222 (Monte Carlo interval 1.0106–1.0340).
At CV>=0.5, exact theta RMSE is lower than PFML's in these data. Such comparisons
use covariance from shared paths, not an independent-sample standard error.

Mu is retained in original state units because its truth is zero. At CV=3,
exact/PFML/Euler RMSEs are 0.03095/0.03702/0.03913 and operational coverages
94.15%/93.70%/99.50%. The paired exact-minus-PFML mu RMSE difference is
-0.00607 (Monte Carlo interval -0.00688 to -0.00526); their coverage difference
is 0.45 percentage points (interval -0.60 to 1.50). All relative-to-truth mu
metrics remain undefined. The increasing Euler mu coverage is not evidence that
its misspecified covariance is robust or universally calibrated.

### Spacing, figures and remaining scope

The registered floor was applied without renormalization. Mean clipped-gap
fractions at nominal CV 1.5, 2 and 3 are 0.1718%, 2.4597% and 17.8223%.
At CV=3 the across-path mean realized gap is 0.501919 and mean realized span
501.416925; the regular span is 499.5. The saved source summaries include
realized CV, minimum gaps, spans, clips and their Monte Carlo uncertainty.
These results describe the registered floored designs. Floor sensitivity and
possible numerical anomalies will be examined only under the later dated E3
diagnostic protocol, after the complete E2 run.

F2 presents both coverage denominators for all three parameters; its companion
shows availability. F3 gives theta/sigma relative bias and a companion for signed
mu bias and RMSE. F4 gives all three independent regular-reference RMSE ratios.
All five PNGs were visually inspected individually, and matching PDFs plus a
hashed figure manifest are retained. Figures are reproduced from saved tables.

E1 alone does not establish a universal smooth operating boundary, simultaneous
grid guarantees, behavior at other n/q values, or novelty. The next stage is
the full registered 60-cell E2 grid, implemented only after Week 5 publication.

Personal reading completion and theory mastery are not certified by this code
or its tests. The project makes no novelty or publication-readiness claim.
