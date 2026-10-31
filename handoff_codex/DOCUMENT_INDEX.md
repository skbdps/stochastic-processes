# Document and evidence index

All `repository_snapshot/` files are byte copies from scientific freeze `0ff57919f5d8345f52d8689732ff94448ca3ca1d`. This snapshot contains code and governing documentation for reading, not the full historical data or a `.git` directory. The manifest records every copied file's source. Use a full Git clone for execution and old experiment audits.

## Read first

1. `START_HERE.md`: current verified checkpoint, corrected progress, remaining work and authorization.
2. `CODEX_PROMPT.md`: continuation instructions ready for Codex.
3. `original_docs/OU_master_plan_v2.pdf`: original historical master plan. Recovered from historical Git blob `2dd700ced7fdb44526e42f477e8f4bbd93b3af41`, present at historical commit `088db0a4a05698af34b03d73aa0677b8e1033146`. Kept only inside the downloadable handoff, not restored to the repository. Its historical schedule, aspirational novelty and professor-related text do not override subsequent instructions.
4. `repository_snapshot/preregistration.md` and `repository_snapshot/docs/plan_amendment_week3_validation.md`: governing registered protocol and amendment. Read these rather than executing conflicting older plan wording.
5. `repository_snapshot/docs/weeks5_8_execution_plan.md` and `handoff_week5_start.md`: historical sequencing and original handoff context, not the current resumption point.

## Completed core and mathematical foundations

Read the progress reports and mathematical notes for Weeks 3–7 under `repository_snapshot/docs/`. In particular:
- `week3_mathematics.md`, `week3_validation_update.md`.
- `week4_architecture.md`, `week4_mathematics.md`, `week4_progress.md`.
- `week5_architecture.md`, `week5_mathematics.md`, `week5_progress.md`.
- `week6_mathematics.md`, `week6_progress.md`.
- `week7_exploratory_protocol.md`, `week7_progress.md`, `week7_numerical_disposition.md`, `week7_runtime_evidence.md`, and the findings memo.

All files currently under the pinned `docs/` directory are included, so supplementary module contracts and historical incident records are available without guessing names. Mathematical derivations, prior-art notes and original progress evidence are supporting material, not certificates of the user's personal reading or novelty.

## Week 8 authority and implementation contracts

Read in this order under `repository_snapshot/docs/`:
- `week8_gate.md` and `week8_exploratory_protocol.md`.
- `week8_mathematics.md` and `week8_population_targets.md`.
- `week8_implementation.md`, `week8_source.md`, `week8_execution.md`, `week8_reporting.md`, `week8_audit.md`, `week8_audit_progress.md`, `week8_plots.md`, `week8_runner.md`.
- `week8_runner_validation.md`, followed by the actual `repository_snapshot/artifacts/week8_evidence/validation_gate.json`.

The builder checks these paths exist. Module-level and older progress files can contain historical 'next step' language; use the later validation gate and current full-study artifact for current status. The snapshot README also contains stale status prose that Codex must update before the final release.

## Implementation files supplied

`repository_snapshot/ou_irregular/`, `tests/`, `scripts/`, `configs/`, the requirements files and pytest configuration are included. Read the existing Week 8 runner rather than writing a replacement automatically. This is a CPU numerical project; no GPU setup is needed for the existing pipeline. The production and historical replay environments are recorded in evidence; do not replace them with arbitrary latest dependencies.

## Critical full-study checkpoint

`checkpoints/week8_full_e4_41200456847.zip` is the original, digest-verified GitHub Actions artifact, not a regenerated experiment.

After safe extraction into a new directory:
- `week8_e4/resolved_design.json` and `run_metadata.json`: exact scientific identity, original commit, runtime, source hashes and lifecycle.
- `week8_e4/cells/`: all raw immutable partitions and cell manifests.
- `week8_e4/audit.json`: saved passed source-specific numerical audit.
- `week8_e4/parameter_summary.csv`, `matched_iid_comparisons.csv`, `regular_comparisons.csv`, `paired_comparisons.csv`, source/spacing/numerical summaries and complete failure export.
- `week8_e4/figures/`: all 18 new full-study PNG/PDF pairs and manifest, awaiting actual visual review.
- `week8_e4-figure-renderer/renderer_state.json`: preserve beside the archive.
- `evidence/invocation.json`, `dependencies.txt`, `command.txt`, `stdout.txt`, `stderr.txt`, `exit.json`, `full_summary.json`, `tracked_source_preservation.json`: original execution and publication preparation evidence.

`evidence/E4_READBACK_SUMMARY.json` is a new handoff-time hash/count readback only; it does not replace the original numerical audit. `evidence/E4_ARTIFACT_METADATA.json` records the original GitHub artifact ID, SHA-256, size and expiry. The original full-study artifact retains the published-freeze provenance unchanged.

## CI failure and recovery documentation

`ci_snapshots/completion/`, `ci_snapshots/runner/`, and `ci_snapshots/replay_diagnostic/` contain exact pinned workflow definitions and diagnostic readbacks from the corresponding staging branches. `evidence/github/` includes workflow/job metadata and logs for the full E4 execution, integrated validation and final release gate.

The current full-study failure is the Git publication-tree HTTP 422 error, not an unexecuted experiment. Earlier failures and the strict historical replay investigation are retained. `docs/week8_runner_validation.md` distinguishes observed hardware-dispatch sensitivity from an unproven explanation of the earlier failure.

## Files deliberately not duplicated here

The bundle is not a full 400+ MB Git history export. Earlier large E1/E2/E3 raw archives and all historical branches remain in the repository and arrive with the fresh full-history clone. Their progress reports and key manifests are included for orientation. Never infer a missing earlier archive was deleted merely because it is not duplicated in this handoff ZIP.

The 29 October Word handoff and continuation prompt described a much earlier partial Week 5 checkpoint. This handoff replaces their operational status while retaining the governing plan, amendments and repository handoff documents. Do not resume at that obsolete point.

## Completion boundary

The remaining task is to recover/verify the existing E4 result, fix publication, inspect the new full-study figures, seal honestly, interpret the results, update documentation, and publish the fully validated Weeks 5–8 work to main. No professor correspondence, paper submission or Weeks 9–12 manuscript/figure-freeze expansion is authorized by this handoff.
