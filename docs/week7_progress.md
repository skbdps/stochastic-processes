# Week 7 complete: exploratory diagnostics and fixed sensitivity studies

Completed 30 October 2026. Scientific freeze:
`6c97d939b4a9902e4127751ff263b5e03ebe7188`, published before the full
computations. This report accompanies the completed Week 7 publication on
`weeks5-8-sequential`; it does not itself constitute the later Week 8 gate.

**All three fixed parts of the dated exploratory protocol have run.** The
selected-path archive is complete, the two new studies have passed their
numerical audits and actual figure reviews, and both new study archives are
sealed `completed`. No headline estimate, interval, seed, replication count,
threshold or scientific dependency was changed. The [one-page memo](week7_findings_memo.pdf)
and [detailed selected-path disposition](week7_numerical_disposition.md) give
complementary summaries.

## Completed evidence and accounting

| Component | Sources or physical realizations | Fit records | Valid point records | Valid interval records | Pipeline errors |
|---|---:|---:|---:|---:|---:|
| Fixed selected E1/E2 diagnostics | 53 selected sources | 159 original replays | Original decisions retained | Original decisions retained | 0 failed selected-source evaluations |
| Independent follow-up | 20,000 | 60,000 | 59,579 | 59,556 | 0 |
| Paired floor sensitivity | 24,000 | 72,000 | 72,000 | 72,000 | 0 |

The floor realizations represent **8,000 independent (CV, replication)
coupling groups**, with three epsilon variants per group. The follow-up has
20,000 independent source paths. Estimators within a path are paired; fit
records are not independent experimental units. Sources:
[selected audit](../artifacts/week7_selected/audit.json),
[follow-up audit](../artifacts/week7_followup/audit.json),
[floor audit](../artifacts/week7_floor/audit.json), and
[completion verification](../artifacts/week7_evidence/full_run_completion_verification.json).

The original 1,494 tracked files were checked against their recovered
published-byte baseline before this completion report's README update; all
matched. That includes all 86 tracked E1 archive files and 335 tracked E2
archive files. Both headline experiments, their six unresolved parameter
classifications in total, and their retained failures remain unchanged. The
[baseline map](../artifacts/week7_evidence/resume_baseline_files.json) and
completion verification preserve the check. New exploratory results are
neither pooled into nor substituted for E1/E2.

## Selected paths: boundaries and weak curvature

All 53 fixed sources reproduce their original arrays and all 159 original
estimator records. The 46 PFML boundary cases have nonpositive independently
calculated OLS coefficients. For 39 of those sources the observed gaps are
exactly equal, supporting the same conditional AR(1) boundary argument for
the exact estimator. Those 46 PFML plus 39 exact records account for the
original 85 point-invalid outcomes. A PFML boundary under unequal gaps does
not prove a boundary of the exact OU likelihood.

The 14 originally valid exact points received the fixed profile diagnostics.
All stored full fitted objectives agree with the best saved finite candidate
within the prescribed resolution. All five starts and their outcomes remain
recorded, including three converged local searches that stop above the best
candidate. A finite search does not prove a global optimum.

The rejected E2 exact interval at n=250, q=2, CV=0.5, replication 465 remains
`unstable_hessian`. Its first whitened discrepancy is 0.0956515 in the
registered calculation and 0.0956538 in the independent longdouble ladder,
both above the unchanged 0.01 gate. Step halvings give approximately fourfold
decay. This supports finite-difference truncation amplified by weak curvature;
extra precision does not remove the original discrepancy. At h/2, the full
matrix's first diagonal is about 124.508, but its algebraic Schur complement
is only 0.00748436. The diagonal cannot substitute for profile information.

The regular control at replication 1424 remains accepted: its first
discrepancy is 0.00680701. Both cases are weakly separated from the iid
limiting likelihood, and both original decisions follow the same numerical
rule. No finer-step interval replaces either decision. The actual longdouble
format has 64 significand bits in 16 storage bytes, not arbitrary precision.
The [numerical disposition](week7_numerical_disposition.md) gives exact values,
boundary qualifications and the scope of independent checking.

## Independent follow-up: availability affects operational coverage

The fixed follow-up uses n=250, q=2, CV=0 and 0.5, R=10,000 per cell and seed
20261007. All 444 unavailable-interval records are retained, on 231 distinct
source paths. There are 421 invalid points and 23 additional interval
rejections among otherwise valid points:

| CV | Estimator | Invalid points | Valid points with unstable Hessian | Available intervals |
|---|---|---:|---:|---:|
| 0 | Exact | 202 boundary | 10 | 9,788 |
| 0 | PFML | 202 boundary | 10 | 9,788 |
| 0 | Euler | 0 | 0 | 10,000 |
| 0.5 | Exact | 1 boundary suspected | 0 | 9,999 |
| 0.5 | PFML | 16 boundary | 3 | 9,981 |
| 0.5 | Euler | 0 | 0 | 10,000 |

The irregular exact record is replication 4369. Its saved reason is that an
expanded profile endpoint or limiting objective is unresolved and no finite
optimum is certified. Its finite candidate and optimizer-success flag are
retained, but neither makes it an eligible point estimate. This remains an
unresolved statistical/numerical outcome under the original rule; it is not
promoted to a proven boundary or silently refit. See the
[complete failure export](../artifacts/week7_followup/failed_fits.csv.gz) and
[outcome table](../artifacts/week7_followup/outcome_counts.csv).

The following coverage percentages use **all 10,000 planned paths** for each
method and cell; relative bias uses valid point fits only.

| CV | Estimator | Theta coverage | Sigma coverage | Theta relative bias | Sigma relative bias |
|---|---|---:|---:|---:|---:|
| 0 | Exact / PFML | 92.63% | 92.71% | +7.55% | +2.37% |
| 0 | Euler | 0.00% | 0.00% | -56.52% | -50.70% |
| 0.5 | Exact | 94.77% | 94.97% | +5.47% | +1.49% |
| 0.5 | PFML | 72.29% | 73.96% | -13.90% | -8.26% |
| 0.5 | Euler | 0.00% | 0.00% | -59.50% | -44.73% |

The regular exact/PFML theta coverage among the 9,788 available intervals is
94.636%, compared with 92.63% operational coverage. This difference is a real
consequence of keeping unavailable intervals in the planned denominator.
The operational theta Wilson interval is 92.101%-93.126% for regular exact,
94.316%-95.189% for irregular exact, and 71.404%-73.159% for irregular PFML.
Euler's zero observed theta/sigma coverage has a 95% Wilson upper limit of
0.0384%, rather than establishing a population probability of zero.

Both exact theta/sigma cells meet the adopted practical criteria in this
separate study. Irregular PFML theta/sigma and Euler theta/sigma show supported
breaches. These labels do not overwrite any headline label. Mu operational
coverage ranges from 92.90% to 96.25% across these six method/cell combinations;
its relative bias and relative RMSE remain undefined because its truth is zero.
All parameter-specific counts, absolute metrics, MCSEs, Wilson intervals,
paired-method and independent regular-design comparisons remain in the
[follow-up tables](../artifacts/week7_followup/parameter_summary.csv).

## Paired floor sensitivity: small effects at the tested settings

The complete fixed design is n=1000, q=0.5, CV in {0,1.5,2,3}, epsilon in
{1e-8,1e-6,1e-4}, R=2000 per variant and seed 20261009. Each CV/replication
shares its underlying Gamma draws, stationary initial-state draw and
standardized OU innovations across epsilon. These are coupled exact OU paths,
not a single continuous path observed on different grids. No gap
renormalization or realized-span matching is used. CV=0 gives identical
arrays and numerical results across all three variants and is one shared
negative control.

All 72,000 fits and intervals are valid. All 72 nonreference epsilon contrasts
therefore use the full common-valid subset of 2,000 pairs. The largest absolute
estimated bias difference across methods/parameters is 0.0000341530 (PFML
theta, CV=3, epsilon=1e-4 versus 1e-6). The largest RMSE-ratio departure from
one is exact theta at CV=3, epsilon=1e-4: **1.00023261**, with pointwise 95%
paired Monte Carlo interval **[0.99974475, 1.00072070]**. This is an estimated
RMSE change of about 0.0233%, with uncertainty spanning no change.

Numerical availability does not imply accurate coverage. At epsilon=1e-6,
exact theta coverage ranges from 94.70% to 95.30% across the four CVs. PFML
theta coverage is 1.25%, 0% and 0% at CV=1.5, 2 and 3, respectively, and PFML
sigma coverage is 0% at all three irregular settings. Small floor effects
therefore coexist with large working-model coverage failures.

Some very small effects are distinguishable from zero under the pointwise
Monte Carlo calculations. For example, PFML sigma at CV=3 and epsilon=1e-4
has a bias difference of 0.0000123503 (MCSE 0.0000025648) and RMSE ratio
0.99993672 [0.99991083, 0.99996261]. These findings are retained; the study is
not reported as exact invariance or as a formal equivalence test. No
simultaneous inference across contrasts is claimed.

The largest observed absolute operational-coverage difference is 0.05
percentage points. Several contrasts have no observed paired coverage
variation and therefore collapsed empirical normal-Monte-Carlo intervals;
these do not prove population equality. At CV=3, the mean clipped-gap fraction
changes from 10.698% to 17.832% to 29.714% across the three floor settings.
Many extremely short gaps can be clipped while the changes in realized mean
gap and the reported estimator metrics remain small. This result concerns
the fixed grid and coupling only. Sources:
[all paired contrasts](../artifacts/week7_floor/epsilon_comparisons.csv),
[per-variant summaries](../artifacts/week7_floor/parameter_summary.csv), and
[realized-spacing summaries](../artifacts/week7_floor/spacing_summary.csv).

## Validation, audit scope and visual review

The runner's focused final gate passed 82 tests with three documented
floor-only skips; this is not a new certificate for the entire historical
test suite. The discovered sealing defect, initial failing tests/source,
repair and same-input integration evidence remain in the
[runner validation record](week7_runner_validation.md). The fixed scientific
code was then published before the full studies.

The selected audit replays every selected source and all 159 original fits.
Its explicit `independent_recalculation_of_all_diagnostic_math=false` remains
unchanged. An additional read-only review passes 4,396 targeted saved-data
checks, including independent binary64 algebra on 352 saved Hessian matrices;
it does not reconstruct Hessian stencils from observations.

The follow-up audit checks every raw record and table and replays four fixed
endpoint sources / 12 fits. The floor audit additionally reconstructs all
24,000 source realizations, verifies all 48,000 observation/time hashes, and
replays 24 fixed endpoint sources / 72 fits. The separate
[stdlib-only review](../artifacts/week7_evidence/independent_week7_results_review_result.json)
independently reproduces all 126 parameter-summary rows and 72 paired epsilon
rows from saved raw data, with unchanged archive hash maps before/after. It
generates no observations or fits and explicitly does not independently
certify every other table field or Hessian calculation.

All eight new PNGs were actually opened and visually inspected: two selected
diagnostic figures, two follow-up figures and four floor figures. Their exact
figure-manifest bindings and review notes are retained. Both study seals were
read back successfully after review. Full invocation logs, including empty
stderr files, are in [full invocations](../artifacts/week7_evidence/full_invocations/).

No new unexpected root-manifest transition or pipeline exception was observed
in these full runs. The historical unattributed state regressions and lost
forensic evidence remain documented in the [runtime record](week7_runtime_evidence.md)
and [recovery record](week7_recovery_20261030.md). Current nonrecurrence and
verified hashes establish current-run integrity within the checks performed;
they do not identify or prove a causal repair of the earlier incidents.

## Reproduction and readback

Use the pinned scientific environment and thread settings, from the repository
root. Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0, pandas 2.2.3 and matplotlib
3.10.8 were recorded; exact provenance is embedded in each archive. Use fresh
directories for distinct numerical runs. Never write over published archives.

```bash
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MPLCONFIGDIR=/tmp/ou_recovery_mpl
python -m ou_irregular.week7.runner verify --kind followup --run-dir artifacts/week7_followup --require-completed
python -m ou_irregular.week7.runner verify --kind floor --run-dir artifacts/week7_floor --require-completed
python artifacts/week7_evidence/independent_week7_results_review.py --followup artifacts/week7_followup --floor artifacts/week7_floor --output /tmp/my_week7_scalar_review.json
```

The exact full-execution commands were:

```bash
python -m ou_irregular.week7.runner run --kind followup --output artifacts/week7_followup --workers 2 --partition-size 500
python -m ou_irregular.week7.runner run --kind floor --output artifacts/week7_floor --workers 12 --partition-size 500
```

For a new run change both output paths. Selected diagnostics used
`path_diagnostics.run_selected` with the committed
`artifacts/week7_evidence/selection.json` and sealed E1/E2 roots, followed by
`plots.write_selected_figures` into a separate directory. The standalone
[runner contract](week7_runner.md) describes resumption, read-only audits and
the mandatory actual visual review before sealing. The preserved saved-matrix
review script records its original absolute input and output paths. For
another checkout, use a separate copy with those two path constants changed
to the intended archive and a fresh external output; retain the published
script and its original review record unchanged.

## Limits and next gate

This work establishes no universal global-optimization guarantee, new estimator,
sandwich-covariance correction, personal theory mastery, Tang-Chen full-text
reading, or publication-readiness claim. The new follow-up boundary-suspected
record remains unresolved under its original eligibility rule. All
exploratory interpretation remains distinct from registration.

Week 7 is complete for publication. After this checkpoint is published, Week 8
must commit its gate verdict separating current scientific completeness,
known implementation defects, current archive integrity and the remaining
historical operational uncertainty. No Week 8 implementation or optional
MMPP results are included in this Week 7 completion checkpoint.
