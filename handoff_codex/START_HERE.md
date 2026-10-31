# OU research: Codex handoff

Prepared 31 October 2026 from live GitHub reads, not from the earlier interrupted chat summaries.
Repository: https://github.com/skbdps/stochastic-processes

## The most important correction

**Do not restart at Week 5, rebuild the Week 8 runner, or launch a new full E4 simulation by default.** The live repository advanced beyond the last conversational status report. Weeks 5–7 are published; the complete Week 8 implementation and integrated R=4 validation are published; a full R=2000-per-cell E4 execution subsequently succeeded, including its numerical audit and figure generation. Its final GitHub publication-tree operation failed. Recover and verify that saved full-study output first.

This handoff does not declare the project complete. The full E4 archive remains `audited_pending_visual_review`; its new figures need actual inspection, a truthful review record, guarded sealing, result interpretation, durable publication, README updates and final merge to main.

## How to use the bundle

Read this file, `CODEX_PROMPT.md`, and `DOCUMENT_INDEX.md`. The ZIP contains the exact original master-plan PDF as a historical reference, a pinned source/document snapshot, recorded validation evidence, CI workflow snapshots and logs, and the original full E4 Actions artifact ZIP. `MANIFEST.json` and `SHA256SUMS.txt` bind the bundle's files. Run `python3 verify_bundle.py` before using it. The packager performs integrity/readback checks, not a new numerical audit or visual review.

The lightweight `repository_snapshot/` is not a Git checkout and omits earlier large experimental archives. Use `python3 setup_repository.py --destination ../stochastic-processes-codex` to make a fresh full-history clone at the pinned scientific freeze; it refuses an existing destination, never overwrites local work, and changes no remote refs. A full clone supplies earlier experiments, historical commit objects and archived CI evidence. Do not run numerical provenance checks from the lightweight snapshot.

The complete E4 output is in `checkpoints/week8_full_e4_41200456847.zip`. Extract it into a new evidence directory outside the repository, preserve that original ZIP, and work on a verified copy. Its roots are `week8_e4/`, `week8_e4-figure-renderer/`, and `evidence/`. Keep the renderer state beside the archive when copying. Do not overwrite any existing evidence directory.

## Exact current checkpoints

| Role | Branch | Commit at handoff |
|---|---|---|
| Scientific implementation freeze and completed validation | `weeks5-8-sequential` | `0ff57919f5d8345f52d8689732ff94448ca3ca1d` |
| Final publication destination, still old baseline | `main` | `4e35ba23d79ec6dc16f10035df27729817e23678` |
| Successful runner CI and retained previous failures | `validation/week8-runner-20261031` | `afcec2b8c82e6f4e1ba7fd58fbde56582a57e1ee` |
| Full-E4 orchestration and failed-publication diagnostic | `validation/week8-completion-20261031` | `2a6f5dbf798171f16a23bbd2a47f7ecda378e98b` |
| Historical replay investigation | `validation/week8-replay-diagnostic-20261031` | `1ac7b78fb0217dcf4efd138c43871fa1904b21e9` |

Always fetch and compare current remote heads before writes. These pins are a reproducible handoff, not permission to reset a newer branch. The handoff branch itself is packaging only; do not merge it wholesale into the scientific branch or main.

## Progress and evidence

### Weeks 1–4

The simulator, estimator and modular inference foundations, mathematical notes, registered regular validation and historical archives are preserved. The original Week 3 smoke is historical; corrected evidence is in `artifacts/week3_validation_rerun/`. Week 4's regular exact checkpoint has 2,000 paths and its registered coverage evidence. Six selected-path bootstrap checks are not a bootstrap coverage study. Personal theory mastery is not certified by these files.

### Week 5: E1 complete

Seven cells; 14,000 independent source paths; 42,000 estimator-fit records; complete summaries, audits and figures. See `docs/week5_progress.md` and `artifacts/week5_e1/`. Preserve all old bytes and recorded failures/classifications.

### Week 6: E2 complete

Sixty cells; 120,000 paths; 360,000 fits. The published record retains 359,915 valid point estimates, 359,914 valid intervals, and 86 unavailable intervals, with zero pipeline exceptions. See `docs/week6_progress.md`. Poor coverage/bias and unresolved Monte Carlo classifications are findings, not reasons to alter seeds or gates.

### Week 7: E3 complete

Published completion at `d8e921df5a395cd17b62be51c6b8b5fabd7f9221` precedes the Week 8 gate. The fixed selected-path investigation includes 53 source paths/159 original-fit replays. The follow-up contains 20,000 independent paths/60,000 fits and retains 444 unavailable-interval records. The floor study has 24,000 realizations/72,000 fits but only 8,000 independent coupling groups. Both new studies are audited, visually reviewed and sealed. Read `docs/week7_progress.md`, `docs/week7_numerical_disposition.md`, `docs/week7_runtime_evidence.md` and the findings memo. Historical unexpected manifest-state changes were not causally attributed; do not erase that limitation.

### Week 8: implementation and complete development validation published

The adopted E4 protocol was specified after the core-completion gate and before affected data. All design, source, execution, reporting, independent-audit, plotting and runner modules are published. The source-only validation retains 1,536 sequences and 108 passed checks. The complete R=4 validation contains 28 independent paths/84 fits with 42 endpoint-fit replays, zero pipeline errors, and a completed seal.

`artifacts/week8_evidence/validation_gate.json` records 2,061 cases accounted for, 2,055 passed, six design-inapplicable skips and zero failures. The six skips concern paired-floor checks instantiated for non-paired follow-up fixtures; none bypasses strict historical replay. The successful runner CI is run `41198942670`; the final validation gate is run `41200166912`. Read `docs/week8_runner_validation.md` for exact scope.

The 18 validation PNGs were byte-identical to earlier actually inspected images. The validation review explicitly transfers that prior judgment, with its original reviewer/time and hashes; it does not claim a new inspection. All 36 external PNG/PDF replot files reproduced byte-for-byte without altering the sealed validation archive. This validation review must NOT be reused for the new full-study figures.

### Week 8: full E4 execution exists, publication not completed

Actions run: https://github.com/skbdps/stochastic-processes/actions/runs/41200456847
Workflow execution commit: `4deaa94d5735bcfff6c6cb47f8cff88df2fb2eb3`
Scientific checkout actually executed: `0ff57919f5d8345f52d8689732ff94448ca3ca1d`
Job: `124535057563`
Artifact ID: `13042844395`
Artifact name: `week8-full-e4-41200456847`
Original ZIP size: 35,558,419 bytes
Original ZIP SHA-256: `e9a780206925c4085ff4bb47cf988c7ccec66ce40db038f3fdef82703fe52513`
Recorded GitHub expiry: 30 December 2026 (the bundle preserves a copy independently).

The original CLI exit code is 0 and stderr is empty. The saved numerical audit reports seven complete cells, 56 partitions, 14,000 source paths, 42,000 fit records, 42,000 valid points and intervals, zero retained failure records and zero pipeline errors. It reports all 14,000 sources reconstructed, 28,000 observation/time array hashes, 42,000 schedule-trace hashes, and 42 endpoint fits replayed on 14 endpoint sources. The run generated 18 PNG/PDF pairs. The lifecycle is `audited_pending_visual_review`, not completed.

The workflow then passed its saved-evidence verification and wrote `evidence/full_summary.json`, but failed while creating the Git publication tree:

```
gh: tree.sha 40af4433b49d2cef073384b77a2d2213963c5cb6 is not a valid blob (HTTP 422)
```

The precise cause of that Git object failure has NOT been established here. Do not assume it was a scientific failure, a permission failure, an empty-lock-file issue, or a repaired defect. The Actions artifact upload succeeded afterward and retained 124 files. The publication step did not advance the scientific branch. Recover `evidence/full_summary.json`, `week8_e4/audit.json`, all raw partitions, all tables, figures and renderer state from the artifact; do not attempt to manufacture missing results from the console excerpt.

## CI history: what is resolved and what is not

Earlier CI failures included predecessor-artifact timing, historical/current environment mixing, spawn-unsafe test startup, and a strict Week 3 byte-replay mismatch. The later release gate supersedes the old 'one test still failing' status. The exact historical assertions passed in isolation; no assertion or scientific gate was removed. A separate diagnostic showed CPU-dispatch-dependent floating-point results: disabling AVX-512 changed some observation hashes on that host. This does not prove the cause of the earlier combined-process failure, whose CPU state was not captured.

Retain CPU dispatch/flags, NumPy/BLAS runtime, thread settings, packages and source hashes. A matching Python/package version is not by itself evidence of exact replay across macOS/ARM and Linux/x86 or across different CPU dispatch. Do not modify old manifests or weaken hash comparisons to make a new host pass. The unresolved current blocker is the publication-tree failure plus genuinely unfinished visual review, interpretation and release work.

## Immediate continuation sequence

1. Verify the handoff manifest; fetch live refs; inspect newer work without overwriting it. Read governing documents before execution or edits.
2. Recover the original E4 ZIP and verify its recorded digest. Extract to a new directory, hash all member bytes, preserve an untouched original, and read the complete metadata/audit/full summary. Inspect the job log and exact workflow at the pinned completion branch.
3. Recreate the recorded scientific checkout/runtime. Run the saved-evidence verifier first, then the appropriate independent numerical/readback checks in a compatible environment. Diagnose any discrepancy before further mutation. Preserve every failure.
4. Repair the publication-only workflow or use authenticated Git to commit the exact recovered artifacts, with path-by-path/hash verification and expected-head guards. Do not rerun 14,000 paths just to work around a publication API problem. Do not automatically merge orchestration/payload branches.
5. Inspect all 18 NEW full-study PNGs individually for legibility, complete series/labels, correct uncertainty/axes/captions and consistency with saved tables. Use the original images, not validation figures or only file signatures. Record actual reviewer/time and the precise current figure-manifest digest. When actual image inspection is unavailable, stop the completion claim and request that narrowly scoped review; never invent it.
6. Seal only after successful review and scientific readback. Re-verify the seal and reproduce figures to a separate directory without modifying the original archive. Retain renderer provenance and original execution logs.
7. Write the Week 8 completion report and mathematical interpretation from saved tables, including matched-marginal MMPP/IID contrasts, Monte Carlo uncertainty, model-based PFML/Euler limitations, raw/transformed spacing diagnostics and honest pointwise scope. Update README, progress, artifact indexes and concrete reproduction commands. The current README is stale about plotting/runner progress and must be corrected.
8. Check every Week 5–8 exit requirement, preserve earlier source/archive bytes, run final clean-checkout validation and applicable review checks. Publish reviewed evidence to the research branch and merge/fast-forward to main only after completion. Read back the remote main commit and key artifacts. Report exact SHAs and checks, not merely that code was written.

## Scope and authorization carried from the user

The user authorized completing the remaining Weeks 5–8 work, debugging CI (or removing/replacing broken orchestration if needed), updating README/docs and publishing the completed work to main. That is not authorization to delete failed scientific evidence, weaken checks, tune results, force-push over other work or bypass branch protections. Ask only for a genuinely missing credential/permission or required actual visual inspection that your environment cannot perform. No professor messages, external approval gate, arXiv submission or Weeks 9–12 manuscript work is included in this continuation.

## Scientific constraints that remain binding

Preserve exact OU simulation, stationary initialization, conditional likelihoods, registered exact/PFML/Euler definitions, full observed-data-normalized 3x3 Hessian and natural-scale unclipped Wald intervals. Never use simulation truth for numerical scaling. Keep all failures and planned keys. Operational coverage uses all planned paths, unavailable intervals count as noncoverage; conditional coverage uses valid intervals. Bias/RMSE condition on valid points. Relative mu metrics remain undefined. Within-cell estimator comparisons are paired; different E4 cells are independent even at equal replication labels. The full-marginal MMPP/IID comparison isolates dependence within this fixed design including realized-span consequences; an intensity-ratio curve or regular comparison alone does not.

E4 is seven cells: MMPP r=1,5,20; independent full-Palm-marginal IID r=1,5,20; one regular control. Truth=(1,0,0.5), q=0.5, n=1000, R=2000, full seed=20261011, namespace=8. Only r=20 has the adopted floor, with no renormalization or latent-clock padding. The R=4 and source-only validations use separate fixed seeds. The dated protocol, preregistration and amendments have priority over conflicting aspirational wording in the historical master plan. Do not silently enlarge scope, reduce replications or re-register existing observations.

## Evidence sources

- Scientific freeze: https://github.com/skbdps/stochastic-processes/commit/0ff57919f5d8345f52d8689732ff94448ca3ca1d
- Validation gate: https://github.com/skbdps/stochastic-processes/blob/0ff57919f5d8345f52d8689732ff94448ca3ca1d/artifacts/week8_evidence/validation_gate.json
- Integrated validation report: https://github.com/skbdps/stochastic-processes/blob/0ff57919f5d8345f52d8689732ff94448ca3ca1d/docs/week8_runner_validation.md
- Full-study diagnostic: https://github.com/skbdps/stochastic-processes/blob/2a6f5dbf798171f16a23bbd2a47f7ecda378e98b/e4-full-diagnostic-41200456847.json
- Full-study workflow: https://github.com/skbdps/stochastic-processes/blob/2a6f5dbf798171f16a23bbd2a47f7ecda378e98b/.github/workflows/week8-full-e4.yml
- Full-study artifact: https://github.com/skbdps/stochastic-processes/actions/runs/41200456847/artifacts/13042844395

The builder records live branch metadata, workflow metadata, original logs, file hashes and the artifact readback in `evidence/`. New handoff prose describes observed stored evidence; it is not a fresh independent scientific certification.
