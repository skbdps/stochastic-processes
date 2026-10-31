# Prompt for Codex

Continue and finish my OU irregular-sampling research project through Week 8. Read `START_HERE.md`, `DOCUMENT_INDEX.md`, the pinned scientific documents, the actual validation evidence and the retained full-E4 artifact before editing. Accuracy, mathematical correctness, reproducibility and preserving the protocol matter more than speed. Carry out the work; do not stop at a plan or call the task complete because code exists.

## Scope and authorization

You are authorized to inspect, debug, test, document and publish the remaining work in `skbdps/stochastic-processes`, update README/progress/reproduction documentation, and merge the fully verified result to main. You may fix or replace broken CI orchestration, but not remove substantive tests or audit requirements to hide failures. Preserve original failure evidence. Do not force-push, reset newer work, delete research artifacts, modify old scientific outputs or bypass repository protections. Ask only for genuinely missing execution/network/Git credentials or a narrowly scoped review that cannot be performed in your available environment. Do not send professor emails or submit a paper. Weeks 9–12 are not part of this task.

## Starting checkpoint — verify live refs

Repository: https://github.com/skbdps/stochastic-processes
Scientific branch: `weeks5-8-sequential`
Published full implementation/validation freeze: `0ff57919f5d8345f52d8689732ff94448ca3ca1d`
Main at handoff: `4e35ba23d79ec6dc16f10035df27729817e23678`
Runner CI branch: `validation/week8-runner-20261031`, at `afcec2b8c82e6f4e1ba7fd58fbde56582a57e1ee`
Completion/orchestration branch: `validation/week8-completion-20261031`, at `2a6f5dbf798171f16a23bbd2a47f7ecda378e98b`
Replay-diagnostic branch: `validation/week8-replay-diagnostic-20261031`, at `1ac7b78fb0217dcf4efd138c43871fa1904b21e9`

Fetch before making claims or writes. Check for newer work and preserve it. Start from the scientific branch, not main and not an old partially uploaded plotting branch. The `handoff/codex-20261031` branch is a packaging aid, not a scientific release branch; do not merge it wholesale.

## Critical saved progress

Weeks 5, 6 and 7 are complete and published. Week 8 source-only validation, all modules, the complete runner and sealed fixed R=4 validation are also published. The final development gate accounts for 2,061 cases: 2,055 pass, six existing design-inapplicable skips, zero failures. The old strict Week3 replay failure was subsequently investigated and passed unchanged in an isolated process. Read the recorded CPU-dispatch limitations; do not claim its earlier exact cause was proven.

MOST IMPORTANT: full E4 has already run. Actions run `41200456847` executed all 14,000 paths and 42,000 fits at the published freeze. Its scientific CLI returned 0, generated tables and figures, and reported a passed full numerical audit with zero pipeline errors. Its saved state is `audited_pending_visual_review`.

The workflow failed AFTER that during Git tree creation:
`tree.sha 40af4433b49d2cef073384b77a2d2213963c5cb6 is not a valid blob (HTTP 422)`.
Do not assume a cause without inspecting the precise Git object/path and publication request. Do not waste or change the experiment by rerunning it to avoid this publishing failure.

The original artifact ZIP is bundled as `checkpoints/week8_full_e4_41200456847.zip`:
- GitHub artifact ID `13042844395`, name `week8-full-e4-41200456847`.
- Size 35,558,419 bytes.
- SHA-256 `e9a780206925c4085ff4bb47cf988c7ccec66ce40db038f3fdef82703fe52513`.
- It contains `week8_e4/`, `week8_e4-figure-renderer/` and `evidence/`, including the full summary, original logs, raw records, all manifests and figure files.

## Execute in this order

1. Verify the bundle manifest. Make or inspect a full-history repository checkout; the lightweight source snapshot is reference material, not a working Git clone. Re-read the authority hierarchy and completion criteria. Confirm exact branches and artifacts against live GitHub.
2. Recover the full E4 artifact without altering it. Verify the ZIP digest, extract safely to a new directory, preserve an untouched original, and inspect every required file and recorded hash. Recreate the source/dependency/thread/CPU environment needed for the saved audit. Use `runner.verify_run` before any expensive computation. Never relabel regenerated data as the original run.
3. Diagnose and repair the publication-only failure, or publish exact recovered evidence through authenticated Git with guarded heads. Do not edit the frozen scientific modules merely to repair a Git upload. Preserve the failed workflow/log and document what was actually fixed. A success in file upload is not a completion certification.
4. Check full planned source/fit accounting, failure export, parameter/contrast/source tables, source reconstruction and endpoint replay requirements. Read the saved independent audit and perform appropriate readback/re-audit without modifying it. Diagnose numerical differences on another host instead of relaxing hash/tolerance gates.
5. Visually inspect all 18 new full-study PNGs. Check axes, legends, labels, captions, error bars, complete cell/method display, undefined values and agreement with the saved numerical tables. Record actual inspection and bind every PNG plus the figure-manifest digest. Do not reuse the R=4 review, invent an inspection, or use the synthetic unit-test review as evidence.
6. Seal with the existing guarded runner only after actual review. Verify the completed archive, regenerate all figures externally from saved records, and verify the original archive remains unchanged. Preserve the sibling renderer-state directory and command logs.
7. Write the Week 8 completion report, result interpretation, limitations and exact reproduction commands. Every numerical claim must trace to a saved table. Distinguish MMPP versus independent full-marginal controls from r-curves or regular comparisons. Explicitly retain working-model covariance and pointwise-inference limitations. Update README, progress, mathematical notes where needed and artifact indexes so they describe actual state rather than old pending-module text.
8. Verify that all Week 5–8 exit conditions are met and prior source/evidence is preserved. Run the applicable clean-checkout regression/replay checks using the documented environment split, with complete case accounting and failure logs. Publish completed evidence and docs, then merge/fast-forward the reviewed scientific branch into main according to repository rules. Read back main and verify the final SHA and key artifacts before reporting completion.

## Governing scientific rules

The dated amendment and preregistration supersede conflicting historical-plan wording. E3/E4 have their own dated exploratory protocols. No outcome-driven seed, threshold, floor, sample-count or grid changes. No repeated runs to achieve 94–96% coverage. No discarded failed paths. Do not overwrite historical archives. The original master-plan PDF is included solely as a historical reference and must not be restored to the repository.

Keep exact OU data generation, stationary initialization, conditional fitted likelihood, registered estimators and full-Hessian natural-scale unclipped 95% Wald inference. Preserve availability/count denominators and uncertainty formulas. Relative mu metrics are undefined. Within-path estimator contrasts are paired; different cells are independent; equal replication labels do not make pairs. Do not conflate 24,000 Week7 floor realizations with 24,000 independent groups.

E4 is exactly seven cells, theta=1, mu=0, sigma=0.5, q=0.5, n=1000, R=2000, seed=20261011, namespace 8. MMPP/IID r=1,5,20 plus one regular control. Maintain event-stationary phase initialization, full-marginal IID controls and the specified r=20 floor transform. Do not infer dependence-only effects from an r-curve alone or hold realized span fixed after the fact.

## Required final response

Report the final remote main SHA; the scientific freeze and full-run provenance; actual tests/audits and their scope; fixed experiment counts; any failures and dispositions; actual figure review and reproduction evidence; documentation updates; and remaining limitations. Do not certify personal reading mastery, paper acceptance, novelty or Weeks 9–12 completion. Keep durable progress checkpoints throughout so an interrupted session does not lose work. An unavailable tool or credential is a specific limitation to report, not grounds to declare incomplete work done.
