# Research implementation handoff

29 October 2026. Repository: https://github.com/skbdps/stochastic-processes

## Exact status

Main at `4e35ba23d79ec6dc16f10035df27729817e23678` is complete through the
Week 4 computational work. Branch `weeks5-8-sequential` contains the execution
plan and an unfinished Week 5 implementation checkpoint. No E1 or E2 headline
experiments have run. Weeks 6–8 implementation has not begun. Personal theory
mastery and outstanding reading are not certified by software completion.

## Work added in this checkpoint

- `docs/weeks5_8_execution_plan.md`: sequential strategy and weekly exit gates.
- `ou_irregular/experiments/design.py`: immutable cells/run designs, validated
  domains, duplicate rejection and the existing registered cell hash convention.
- `ou_irregular/week5/config.py` and `configs/week5_e1.yaml`: exact E1 design
  checking. Separate validation streams cannot be labelled headline runs.
- `ou_irregular/experiments/summaries.py`: complete Cartesian key verification,
  composition with Week 4 coverage/classification, RMSE Monte Carlo SE,
  availability Wilson intervals and MCSE, matched-subset reported-SE diagnostics,
  once-per-source spacing summaries and outcome counts.
- `tests/experiments/`: 20 new tests for design identity, frozen settings,
  hand-calculated accounting and missing/duplicate groups.

Validation at this checkpoint: **540 tests passed** with `python -m pytest -q`
(520 existing and 20 new). All preexisting scientific source and archived
artifacts are unchanged from main. The new modules have unit tests but still
require independent review, end-to-end integration and artifact validation.

## Next implementation step

Review the two new modules first. Then add and test a comparison module:
independent regular-reference contrasts, and paired within-source estimator
contrasts with common-valid point subsets. A CV=0 self-reference is not an
independent comparison. Operational paired coverage uses all planned paths.
Next implement durable execution, provenance/resume checks, E1 plotting, replay
and artifact auditing, mathematical comments and a progress report. Freeze
tested code before the full E1 run. Publish a completed Week 5 before Week 6.

## Governing scientific specifications

Read `preregistration.md`, `docs/plan_amendment_week3_validation.md`,
`docs/week3_validation_update.md`, and the Week 4 architecture, mathematics
and progress documents. The original master-plan PDF remains in history at
`088db0a4a05698af34b03d73aa0677b8e1033146`; it was deliberately removed from
the current checkout. Do not restore it into the repository.

E1 is seven CV cells at q=0.5, n=1000, R=2000, seed 20261001. E2 is the full
60-cell q/CV/n grid, R=2000, seed 20261002. Preserve exact numerical inference,
the literal high-CV floor, seed addressing, all failures and both coverage
denominators. No favorable-outcome tuning or silent compute reductions.

Reuse numerical services through composition. Do not relax the Week 3 smoke
schema or Week 4 validation runner. The historical Week 4 auditor compares its
entire recursive package source map, so additive modules invalidate that exact
historical environment comparison. Run it at its historical revision; new
audits should verify explicitly recorded dependency files. This is not license
to edit old manifests or weaken numerical checks.

## Remaining weeks

Week 6 implements full E2, coverage heatmaps and the failure log. Week 7 specifies
exploratory selection rules before anomaly/floor-sensitivity runs, preserves
headline outputs and produces the one-page memo. Week 8 records a core-completion
gate; MMPP is optional and needs a new precise protocol, including event-stationary
initialization and a full-marginal independent-gap control for dependence claims.
Prepare documents but do not send professor emails or require professor approval.

The registered numerical checkpoint already passed; its two selected bootstrap
source paths are diagnostic only, not a bootstrap coverage study or automatic
bias correction. Proposal refinements and newly derived population benchmarks
are not implemented or independently certified by this checkpoint.
