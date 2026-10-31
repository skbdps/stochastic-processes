# Week 8 progress and Weeks 5–8 exit review

Prepared 31 October 2026. Current release disposition:
**completed scientific study; release evidence and guarded main publication recorded below**. This is the current operational record; older
module checkpoint documents retain their dated scope and frozen scientific
bytes. The [completion report](week8_completion.md) records the original run
identity and later operations; the [artifact index](../artifacts/README.md)
distinguishes headline results, development validation and recovery evidence.

## Prior weeks remain complete

The [sequential execution plan](weeks5_8_execution_plan.md) required publishing
one completed week before implementing the next. Its old Week 5 handoff
paragraph is historical and does not describe the current checkpoint.

| Stage | Governing exit requirement | Existing evidence and disposition |
| --- | --- | --- |
| Week 5 E1 | Execution-plan Week 5 steps 1–8: immutable registered design, complete uncertainty/comparisons, deterministic durable execution, all 14,000 paths / 42,000 fits, complete audit, actual figures and publication before E2. | [Week 5 completion](week5_progress.md), [sealed archive](../artifacts/week5_e1/) and [evidence](../artifacts/week5_evidence/). Complete; preserve its unresolved PFML-sigma classification and identical-tree freeze mapping. |
| Week 6 E2 | Execution-plan Week 6: every registered 60-cell grid point, 120,000 paths / 360,000 fits, parameter/sample-size figures, complete failure/uncertainty record, audit and publication before E3. | [Week 6 completion](week6_progress.md), [sealed archive](../artifacts/week6_e2/) and [independent reviews](../artifacts/week6_evidence/). Complete; retain 86 unavailable intervals and five unresolved classifications. |
| Week 7 E3 | Execution-plan Week 7 and [dated protocol](week7_exploratory_protocol.md): fixed post-headline selection, preserved headline outputs, numerical/floor investigation, findings memo, dispositions and reproduction evidence. | [Week 7 completion](week7_progress.md), [numerical disposition](week7_numerical_disposition.md), [runtime evidence](week7_runtime_evidence.md), [findings memo](week7_findings_memo.pdf). 53 sources / 159 original-fit replays, 20,000 independent follow-up paths / 60,000 fits, 24,000 floor realizations / 72,000 fits in 8,000 coupling groups. Both new studies sealed; 444 unavailable follow-up intervals retained. |
| Week 8 permission to extend | Execution-plan Week 8 and [gate](week8_gate.md): publish core-completion verdict after E3, document unresolved operational limitations, specify optional clustered-sampling design before affected observations. | Gate dated 30 October follows Week 7 publication `d8e921df5a395cd17b62be51c6b8b5fabd7f9221`. Dated [E4 protocol](week8_exploratory_protocol.md) and [mathematics](week8_mathematics.md) published before implementation/data. Proceed verdict has bounded scope. |

The prior completion records do not erase historical failures. The Week 7
boundary-suspected follow-up point remains unresolved under its eligibility
rule; earlier unexpected manifest-state transitions remain unattributed.
Current preservation checks are required in addition to those original
completion records: **[preservation and clean-checkout release verification](../artifacts/week8_evidence/completion/release_candidate.json)**.

## E4 exit criteria matrix

| Governing requirement | Actual evidence or remaining verification |
| --- | --- |
| Protocol §§1–3: exactly seven independent cells, event-stationary Palm initialization, full-marginal IID reset, exogenous exact OU, fixed truth/q/n/R/seed/namespace, r=20-only floor without renormalization/padding. | Original [resolved design](../artifacts/week8_e4/resolved_design.json), [metadata](../artifacts/week8_e4/run_metadata.json), complete raw partitions and source audit. Fixed full run has 14,000 independent sources / 42,000 fits; no scientific override. |
| Protocol §4: analytic/limiting identities and fixed source-only validation before fitted full observations; fixed complete R=4 pipeline and published implementation freeze. | [Source-only validation](../artifacts/week8_source_validation/validation.json): 1,536 sequences, 108 passed checks. [Development gate](../artifacts/week8_evidence/validation_gate.json): sealed seven-cell R=4 validation, 28 paths / 84 fits / 42 endpoint replays. Scientific freeze `0ff57919f5d8345f52d8689732ff94448ca3ca1d` was the original full-run checkout. |
| Protocol §5 and preregistration denominators: retain planned failures, valid-point bias/RMSE, operational versus conditional coverage, undefined relative mu, paired within-source versus independent between-cell contrasts. | All nine saved tables and complete-header zero-row failure export in [full archive](../artifacts/week8_e4/). All 42,000 points/intervals available, zero pipeline errors. [Interpretation](week8_results_interpretation.md) retains denominator and pointwise limitations. |
| Protocol §6 / [audit contract](week8_audit.md): complete keys/schemas, source/dependency bindings, independent arithmetic and all-source reconstruction, endpoint original-fit replay. | Original saved [full audit](../artifacts/week8_e4/audit.json) passed: 56 partitions, 65 dependencies, 14,000 sources, 28,000 array hashes, 42,000 record trace hashes, 14 endpoint sources / 42 fits. New recovered-evidence verification/re-audit: **[passed independent full re-audit, Actions 41214670935](../artifacts/week8_evidence/completion/independent_reaudit.json)**. Full fitting/coarse-Hessian reconstruction is endpoint-only. |
| Protocol §6 / [plot contract](week8_plots.md): render complete figures from saved tables and inspect every current original PNG. | All 18 original PNGs actually opened individually by two Codex AI agents, nine each; axes, legends, labels, uncertainty and saved-table consistency checked. Recorded 2026-10-31 at 15:19:11.466414 UTC and 15:19:55.951874 UTC. Manifest SHA-256 `5d05f882a19e74c4d903b19aff55a575cb531468d33a3cfbb993ab0eac8e5734`; no blocking issue. Review explicitly does not claim human inspection or use R=4/synthetic attestations. Combined record: **[all 18 actual AI image reviews bound to the guarded seal](../artifacts/week8_e4/visual_review.json)**. |
| [Runner contract](week8_runner.md), stages/failure handling: seal only verified evidence after actual review; verify completed state; retain original evidence and renderer state. | Original artifact remains preserved at its pending-review state. Guarded seal of verified working copy and completed readback: **[completed archive; original pending-review metadata retained separately](../artifacts/week8_e4/run_metadata.json)**. |
| Protocol §6 / runner replot contract: regenerate all original figure pairs externally from saved audited records without changing the sealed archive. | 36-file full-study replot and complete before/after archive hash comparison: **[all 36 external PNG/PDF files byte-identical; source archive unchanged](../artifacts/week8_evidence/completion/replot_verification.json)**. Sibling renderer state is retained; validation reproduction alone cannot satisfy this full-study requirement. |
| Protocol §6 / execution-plan evidence: publish completion report, mathematical interpretation, exact reproduction commands and limitations. | New [completion](week8_completion.md), [interpretation](week8_results_interpretation.md), [reproduction](week8_reproduction.md), [publication diagnosis](week8_publication_recovery.md), updated README and artifact index. Scientific contracts remain byte-identical. The completed records bind the original experiment and subsequent verification operations. |
| Execution-plan preservation and complete validation: old source/output bytes retained, no weaker assertions, failures preserved, final applicable regression/historical replay checks with complete accounting. | Pre-full development gate: 2,061 accounted / 2,055 passed / six design-inapplicable skips / zero failures. New clean-checkout split-environment release regression: **[composite gate: 2,061 identities, 2,055 passes, six design-inapplicable skips](../artifacts/week8_evidence/completion/regression_assembly.json)**. Current prior-source/archive preservation: **[preservation and clean-checkout release verification](../artifacts/week8_evidence/completion/release_candidate.json)**. |
| User-authorized release scope: guarded scientific publication and merge/fast-forward to main, obey protections, inspect newer work, read back exact main commit and key artifacts. | Publication-only repair preserves the failed workflow/log and uses recovered original outputs. Expected-head publication, remote main SHA and key-artifact readback: **[guarded main publication and remote readback](https://github.com/skbdps/stochastic-processes/blob/codex/week8-release-readback-20261031/artifacts/week8_evidence/completion/remote_main_readback.json)**. |

## Results and scope disposition

The saved full study supports 16 within-criteria and 26 breach theta/sigma
labels, with 21 mu labels inapplicable and no unresolved E4 labels. All exact
theta/sigma cases are within the fixed criteria. Matched PFML mu at r=5/20
shows increased RMSE and reduced coverage with pointwise MC intervals
excluding their nulls; matched bias contrasts do not establish mean shifts.
Full-marginal comparisons include realized-span consequences. PFML/Euler
covariance remains model-based without sandwich correction. The
[interpretation](week8_results_interpretation.md) traces every displayed
number to saved tables and retains mixed/null comparisons.

The exit review is supported by the bound experimental, visual, numerical and validation records above. Publication is read back independently at the recorded remote commit. No professor
message, external approval gate, paper submission, personal theory certificate,
novelty claim or Weeks 9–12 manuscript completion is included.
