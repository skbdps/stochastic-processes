# Week 5: E1 implementation, validation and results record

Prepared 29 October 2026 UTC. This record separates implementation evidence from
completed headline computations. The governing protocol remains the published
preregistration and dated Week 3 amendment; the implementation plan makes no
scientific amendment.

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

## Headline E1 results

Pending at this pre-execution code checkpoint. No headline numerical conclusion
is inferred from the small validation sample. This section will be completed only
after the registered 14,000 paths / 42,000 fits, independent audit, figure review
and completion seal have actually finished.

Personal reading completion and theory mastery are not certified by this code
or its tests. The project makes no novelty or publication-readiness claim.
