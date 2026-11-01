# Week 8 E4 completion and publication record

Prepared 31 October 2026. Release disposition: **completed scientific study; release evidence and guarded main publication recorded below**.
This record follows the dated [gate](week8_gate.md),
[exploratory protocol](week8_exploratory_protocol.md) and completed
[development validation](week8_runner_validation.md). It supplements those
frozen documents without changing their bytes or scientific requirements.

## Original execution and recovery

The full E4 study was executed once in Actions run
[41200456847](https://github.com/skbdps/stochastic-processes/actions/runs/41200456847).
Its workflow commit was `4deaa94d5735bcfff6c6cb47f8cff88df2fb2eb3`; the scientific
checkout actually executed was the published freeze
`0ff57919f5d8345f52d8689732ff94448ca3ca1d`. The scientific CLI returned zero
with empty stderr. The saved state before recovery was
`audited_pending_visual_review`.

The subsequent Git publication request failed with HTTP 422:
`tree.sha 40af4433b49d2cef073384b77a2d2213963c5cb6 is not a valid blob`.
The workflow nevertheless uploaded the original 124-file artifact, ID
`13042844395`, named `week8-full-e4-41200456847`. Its original ZIP is
35,558,419 bytes with SHA-256
`e9a780206925c4085ff4bb47cf988c7ccec66ce40db038f3fdef82703fe52513`.
Recovery used that bundled original; it did not rerun the full experiment.
The original ZIP, per-member hashes, logs, invocation, dependency capture and
saved numerical audit remain distinguishable from later review and sealing.

Publication repair and precise failure disposition:
the [native-Git publication recovery](week8_publication_recovery.md), verified in [Actions 41218787133](https://github.com/skbdps/stochastic-processes/actions/runs/41218787133), which committed the recovered checkpoint as `f1ad8cc5a7cd2831c6888f2840961178eeeceba7`. A successful later publication does
not by itself establish the exact cause of the original API failure.

## Fixed scientific accounting

The design contains event-stationary MMPP cells at r=1,5,20, independently
generated full-Palm-marginal IID controls at the same ratios, and one independent
regular control. Each cell contains 2,000 paths with 1,000 observations and
999 transitions. Truth is (theta,mu,sigma)=(1,0,0.5), q=0.5, seed=20261011 and
source namespace=8. Only r=20 receives the specified posthoc floor; gaps are
neither renormalized nor padded by extra latent-chain evolution.

| Saved evidence | Observed accounting |
| --- | ---: |
| Complete cells / immutable partitions | 7 / 56 |
| Independent source paths | 14,000 |
| Exact, PFML and Euler fit records | 42,000 |
| Valid point estimates / available intervals | 42,000 / 42,000 |
| Retained unavailable-outcome records / pipeline errors | 0 / 0 |
| Parameter / paired-estimator / matched-IID / regular-comparison rows | 63 / 63 / 27 / 54 |
| Source / spacing / numerical / outcome table rows | 7 / 7 / 21 / 21 |
| Original figure pairs | 18 PNG/PDF pairs |

The zero-row failure export retains its complete required header. These are
saved-result counts, not a claim that unavailable estimates cannot occur in
other settings. Bias and RMSE retain their valid-point conditioning;
operational coverage uses all planned paths, and conditional coverage uses
available intervals. These two coverage values coincide here because all
intervals are available. Relative mu metrics remain undefined.

The original full audit reconstructed every source, verified 28,000
observation/time-array hashes and 42,000 schedule-trace hashes, checked all
nine tables independently, and replayed all three original fits at replications
0 and 1999 in each cell: 14 sources / 42 fits. Its saved audit SHA-256 is
`9fa6c9974ede84412b42c209c0d004741fbf9afa3da921955f986688766e482f`.
It verified 65 declared source/dependency files. All-record coordinate,
covariance, validity and endpoint checks do not imply that every original
fit or coarse Hessian was recomputed: complete fitting replay is endpoint-only.
The exact OU primitive, calibrated cell definitions and analytical gap-moment
routines are shared dependencies, and replay uses the original inference
services. See the [audit contract](week8_audit.md) for independence boundaries.

Saved-evidence verification and new numerical readback/re-audit:
**[passed independent full re-audit, Actions 41214670935](../artifacts/week8_evidence/completion/independent_reaudit.json)**.

## Actual full-study figure inspection and sealing

Two review components cover individual inspection of all 18 original full-study
PNGs, checking axes, labels, legends, complete series, captions, uncertainty,
undefined quantities and agreement with saved tables. The first nine-image
review was recorded at 2026-10-31T15:19:11.466414+00:00; the other nine-image
review was recorded at 2026-10-31T15:19:55.951874+00:00. Both found all
panels legible and no blocking issue. Their records bind individual PNG hashes
and the same figure-manifest SHA-256:
`5d05f882a19e74c4d903b19aff55a575cb531468d33a3cfbb993ab0eac8e5734`.
The development R=4 review and synthetic unit-test attestations were not used
as full-study inspection evidence.

The review retained limitations visible in the figures: near-zero coverage
bars are small on the 0–100% overview scale; nearly identical raw/transformed
series overlap because the floor effect is tiny; zero empirical-normal
contrast intervals follow the saved calculation and do not prove certainty.
Saved tables supply full precision. Raster legibility review does not replace
numerical validation.

Combined review and guarded completion seal:
**[all 18 image reviews bound to the guarded seal](../artifacts/week8_e4/visual_review.json)**. Completed-state verification:
**[completed archive; original pending-review metadata retained separately](../artifacts/week8_e4/run_metadata.json)**. External figure reproduction:
**[all 36 external PNG/PDF files byte-identical; source archive unchanged](../artifacts/week8_evidence/completion/replot_verification.json)**. The reproduction check compared all 36 PNG/PDF files and the archive's complete before/after byte map. The sibling renderer state and original execution logs remain retained.

## Results and bounded interpretation

All 14 exact-estimator theta/sigma combinations meet the adopted pointwise
practical criteria. Regular PFML adds two within-criteria combinations; the
remaining 26 theta/sigma combinations are supported breaches. All 21 mu
classifications are inapplicable at zero truth. E4 has no unresolved practical
classification. Exact coverage over all parameters and designs ranges from
94.10% to 95.55%. The criteria concern the fixed 80% coverage / 10% relative
bias thresholds, rather than a claim of exact nominal calibration.

The clearest matched-marginal difference is PFML mu: MMPP/IID RMSE ratios are
1.083874 at r=5 and 1.157863 at r=20; operational coverage is lower by 2.70 and
3.85 percentage points, respectively, with pointwise MC intervals excluding
the null. No matched bias interval excludes zero. These contrasts hold the
entire one-gap marginal fixed and include dependence-induced changes in the
realized-span distribution. The ratio curve alone changes both marginal law
and dependence. Full values, uncertainties, derivations and caveats appear in
the [results interpretation](week8_results_interpretation.md) and the saved
[parameter](../artifacts/week8_e4/parameter_summary.csv) and
[matched-IID](../artifacts/week8_e4/matched_iid_comparisons.csv) tables.

## Development gate, preservation and final release checks

The pre-full-run gate accounts for 2,061 collected cases: 2,055 passed, six
design-inapplicable paired-floor checks on non-paired fixtures skipped, and
zero failures. It also retains 1,536 source-only sequences / 108 passed checks
and the completed R=4 validation: 28 paths / 84 fits / 42 endpoint-fit replays.
All 36 validation replot files reproduced exactly. These are development
checks and remain separate from the full-study evidence above.

Final clean-checkout regression/replay checks, complete case accounting and
retained failures: **[composite gate: 2,061 identities, 2,055 passes, six design-inapplicable skips](../artifacts/week8_evidence/completion/regression_assembly.json)**. Prior tracked source
and experimental-evidence preservation: **[preservation and clean-checkout release verification](../artifacts/week8_evidence/completion/release_candidate.json)**.
Guarded scientific-branch publication and remote-main readback:
**[guarded main publication and remote readback](https://github.com/skbdps/stochastic-processes/blob/ddf6046b76e942b526244fe9bb67d895d6f7b987/artifacts/week8_evidence/completion/remote_main_readback.json)**. Exact commands and runtime requirements
are in [reproduction](week8_reproduction.md).

Weeks 5–7 retain their previous completion and failure dispositions. E2's
86 unavailable intervals and the Week 7 follow-up's 444 unavailable intervals
remain scientific outcomes; the floor study has 24,000 realizations in
8,000 independent coupling groups. Historical unattributed manifest changes
and lost forensic evidence remain disclosed. CPU-dispatch sensitivity was
demonstrated, but the earlier combined-process replay failure's exact cause
was not established. No seed, threshold, floor, inference gate or archived
outcome is changed to obtain a preferred result. The work does not certify
personal reading mastery, novelty, paper acceptance or Weeks 9–12 completion.


## Final validation across recorded runtimes

The fresh native invocation [41215391840](https://github.com/skbdps/stochastic-processes/actions/runs/41215391840) accounted for every one of 2,061 cases and remains recorded as failed: 2,054 passed, six existing design-inapplicable skips, and one strict historical Week 3 input-hash failure on the AVX2-only AMD EPYC 7763 host. Its complete 271-file output is preserved unchanged. Four complete groups passed: core 653; current Week 7 769 plus six skips; Week 8 492; historical optional-package fixtures 137. These groups total 2,057 cases, 2,051 passes and six skips.

The entire unchanged four-case Week 3 group then passed in [41219755639](https://github.com/skbdps/stochastic-processes/actions/runs/41219755639), including exact observation/time hashes and all fitted values. This used Intel SDE 10.13.1 with the documented `-no-follow-child -icl` driver profile, on the same committed scientific source. Effective NumPy dispatch, the actual Git subprocess smoke, complete runtime and pytest shared process ID were retained. Emulator bytes were checked against Intel's published SHA-256. No scientific code, test, seed, hash comparison or tolerance changed. Emulation is an additional explicitly recorded execution environment, not a claim of native AVX-512 hardware.

The separate assembly record reconciles the four complete native groups and this entire four-case historical group against the original 2,061-node collection exactly once: 2,055 passed, six design-inapplicable skips, zero failures in the selected complete groups. This is a composite gate across recorded runtimes, not one successful invocation. The original failed native summary is never rewritten. Earlier unavailable-host probes and unsuccessful emulator launchers are operational failures retained by artifact identity, digest, source commit and diagnostic disposition. Their tests did not run. The successful result establishes this particular replay environment; it does not uniquely establish the cause of the AVX2-host mismatch or promise bitwise portability everywhere.

CI orchestration is not a research result or an additional release requirement. No new active CI workflow is installed on main. Reproduction uses the preserved commands, guarded verifiers, test identities and actual numerical evidence. Historical workflow branches remain available for provenance. The release verifier rechecks the original 124 files, the complete E4 seal, both regression attempts, all group identities and runtime origins, every pre-existing Git object except the deliberately updated README, and documentation links from a clean checkout.
