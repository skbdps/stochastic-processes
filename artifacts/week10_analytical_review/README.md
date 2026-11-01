# Week 10 analytical correction evidence

This additive package supports the [critical analytical review](../../docs/week10_analytical_review.md)
and [Draft 3](../../manuscript/week10_draft3.md). It contains no new research
replications. The original 250-file supplement, Draft 2, estimators and
recorded intervals remain unchanged.

| Evidence | Scope |
| --- | --- |
| [before_review.json](before_review.json), `before_review/` | Pre-correction working-document hashes and exact copies; changes remain distinguishable from the original draft. |
| [calibration/manifest.json](calibration/manifest.json) | Current four saved-data tables, source/input hashes, exact conditional-mixture calculation, SD comparisons and conservative paired coverage bounds. |
| [calibration_verify.json](calibration_verify.json) | Recalculated tables; all 250 original files verified; all 42,000 raw fit records read and 14,000 PFML centers, SEs and coverage indicators checked. |
| `calibration_initial/`, `initial_source/` | Preserved first-pass diagnostics and source. The first version's normal paired intervals degenerated at zero discordance; current output adds conservative finite-sample bounds. Original approximations remain explicitly labeled. |
| [numerical_validation/critical_numerical_audit.py](numerical_validation/critical_numerical_audit.py), [result](numerical_validation/critical_numerical_audit.json) | Independent 70-digit finite-count checks on the adopted grid and a separately identified off-grid cancellation example. |
| `math_validation/`, [Week 10 test output](week10_tests.stdout), [JUnit](week10_tests.xml) | Deterministic pole/Hessian checks, diagnostic probability/pairing tests and the complete Week 10 suite. |
| [original_verify.stdout](original_verify.stdout) | Original CLI readback of the unchanged seal, five original summary tables, six saved-span exports and generated claim bindings. It verifies historical bytes without endorsing the superseded interpretation. |
| [final_validation.json](final_validation.json) | Final test, preservation, document-link and numerical-statement checks, bound to current source and document hashes. |

The current paired confidence limits subtract simultaneous 97.5% two-sided
Clopper–Pearson bounds for the two discordance probabilities, giving at least
95% coverage for each paired difference. Across-cell comparisons remain
pointwise. SD-prediction intervals are approximate schedule-MC intervals;
the same schedules contribute to the observed experiment and the conditional
benchmark. Fitted RMSE summaries do not establish finite moments of the ideal
ratio estimator.

Reproduce the current calibration with the commands in
[Week 10 reproduction](../../docs/week10_reproduction.md). The numerical audit
can also be rerun from the repository root using the pinned Week 10 environment
and `PYTHONPATH=.`. Preserve outputs separately when repeating a check.
No current review claim transfers to changed source bytes or a changed draft.
