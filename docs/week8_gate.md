# Week 8 core-completion gate

Dated 30 October 2026, after publication and readback of completed Week 7 at
`d8e921df5a395cd17b62be51c6b8b5fabd7f9221`. No Week 8 scientific implementation
or validation/full-study observations preceded this verdict.

## Verdict: proceed with the guarded optional E4 extension

Within the completed checks, no known unresolved defect in the core
likelihood, estimator, inference procedure or current archive integrity
requires Week 8 buffer work. Proceed with the separately specified
[exploratory MMPP protocol](week8_exploratory_protocol.md), preserving every
earlier source file, archive and scientific decision. This is a bounded
implementation decision, not a universal correctness or publication-readiness
certificate.

| Question | Observed evidence | Disposition |
|---|---|---|
| Is the registered core complete? | Published E1: 14,000 paths / 42,000 fits; E2: 120,000 paths / 360,000 fits. All prior tracked files matched the recovery baseline before the completion README update. | Complete; preserve existing results and all six unresolved E1/E2 classifications. |
| Did the fixed Week 7 studies actually finish? | 53 selected sources / 159 original-fit replays; 20,000 independent follow-up paths / 60,000 fits; 24,000 floor realizations / 72,000 fits in 8,000 coupling groups. Both new studies audited, visually reviewed and sealed. | Complete and published before this gate. |
| Did the selected diagnostics establish a core numerical defect? | All 85 original invalid points have the stated boundary evidence. All 14 valid exact profiles match their stored fits within the fixed diagnostic resolution. The rejected Hessian discrepancy persists in longdouble and decreases approximately quadratically with step. | No newly established estimator/inference defect. Retain the original rejection and finite-search limitations. |
| Are current failure records complete? | Follow-up retains 421 invalid points and 23 additional interval rejections, including one unresolved `boundary_suspected` point. Floor study has all points/intervals valid. Both have zero pipeline exceptions. | These are retained statistical/numerical outcomes, not concealed implementation errors. |
| Is the repaired runner defect resolved within tested scope? | Initial failing seal tests/source/outputs retained; repaired fault-injection tests and same-input integrated runs passed. Full study completion readbacks and independent saved-record reviews passed. | No known unresolved seal defect in the tested contract; do not imply immunity to every external writer. |
| Do historical runtime incidents remain uncertain? | Earlier unexpected manifest-state changes were observed. The writer/cause remains unattributed; some original forensic files were lost. Current full runs had no such observed incident and all declared hashes verify. | Retain this operational limitation. Current safeguards and nonrecurrence do not establish a causal repair of the past. |

The scientific completeness and current byte/identity checks support guarded
extension despite the separately disclosed historical uncertainty. The
remaining boundary-suspected follow-up record does not certify an incorrect
accepted estimate: the procedure conservatively withholds eligibility. Small
floor effects, poor working-model coverage and unresolved Monte Carlo labels
are possible findings, not grounds for changing the experiment or calling the
software defective.

## Conditions of the extension

Publish this gate and the dated protocol before implementation. Implement,
review and validate one module at a time. Publish the final reviewed dependency
freeze before full E4 observations. Use the isolated checkout, guarded root
transitions, durable immutable partitions, explicit source policy, full
accounting, source-specific audits and actual figure inspection. Any new
pipeline exception, unexplained state mutation, failed audit or independently
demonstrated numerical defect must be preserved and investigated before an
E4 completion claim. No failed source may be silently replaced or resampled.

The optional extension is justified by the original clustered-sampling stage
and the feasible full-marginal independent control, not by a desired result.
Its fixed full design contains seven cells and 14,000 new independent source
paths. Neither E4 nor this gate amends E1/E2 or retroactively registers E3.

Evidence: [Week 7 completion](week7_progress.md),
[selected numerical disposition](week7_numerical_disposition.md),
[completion verification](../artifacts/week7_evidence/full_run_completion_verification.json),
[prepublication readback](../artifacts/week7_evidence/prepublication_readback.json),
[independent scalar review](../artifacts/week7_evidence/independent_week7_results_review_result.json),
[runner repair validation](week7_runner_validation.md), and
[historical runtime evidence](week7_runtime_evidence.md).
