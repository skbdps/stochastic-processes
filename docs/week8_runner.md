# Week 8 E4 runner and completion contract

The standalone `ou_irregular.week8.runner` composes the previously reviewed E4
source, durable executor, reporting, independent auditor and renderer. It changes
none of their mathematical definitions, scientific designs, seeds, fit eligibility,
Monte Carlo formulas or inference thresholds. Earlier week runners are unchanged.

## Fixed plans and prerequisites

`run --validation` means all seven cells, R=4, seed 21261011: 28 independent
source paths and 84 fit records. `run` without that flag means all seven cells,
R=2000, seed 20261011: 14,000 independent paths and 42,000 records. No command
accepts a seed, replication-count, grid, estimator or numerical-threshold override.
Worker count and partition size affect scheduling/storage only, not source identities.
The root stores the complete resolved design, explicit dependency closure, original
code commit, environment, execution attempts and completed-cell accounting.

Before full E4 data, the dated core-completion gate and exploratory protocol must
be published, the source-only validation must pass, and the complete fixed R=4
pipeline must be checked, visually inspected and sealed. Publish the reviewed full
implementation freeze before executing R=2000. The runner enforces the fixed plans,
committed full-run dependencies, source identities, immutable evidence and stage
semantics. It does not replace the documented scientific review by a software flag.

## Commands

From a clean repository checkout, use Python 3.12.14 with the pinned
`requirements-lock.txt`. Use the same environment and thread variables for execution,
resumption, verification, numerical audit and rendering:

```bash
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
python -m ou_irregular.week8.runner run --validation --workers 2 --partition-size 2 --output /tmp/ou_e4_validation
python -m ou_irregular.week8.runner verify --run-dir /tmp/ou_e4_validation
```

Successful `run` finishes at **audited_pending_visual_review**, not completed.
Inspect every one of the 18 PNGs in `figures/`. Create an external review JSON
record naming the reviewer, all current PNG filenames exactly once, the SHA-256
of `figures/figure_manifest.json`, and `all_panels_legible: true`, with factual
inspection notes. That Boolean attests actual inspection; no test or template may
be passed off as a real visual review. Then:

```bash
python -m ou_irregular.week8.runner seal --run-dir /tmp/ou_e4_validation --review /tmp/e4_validation_visual_review.json
python -m ou_irregular.week8.runner verify --run-dir /tmp/ou_e4_validation --require-completed
# Only after the reviewed full implementation freeze has been published:
python -m ou_irregular.week8.runner run --workers 4 --partition-size 250 --output /tmp/ou_e4_full
```

The full run requires its own inspection and review record. R=4 figures are
prominently labelled validation only and cannot provide E4 scientific findings.

```bash
# Retain original inputs and completed partitions after an interruption.
python -m ou_irregular.week8.runner run --workers 4 --partition-size 250 --output /tmp/ou_e4_full --resume
# Continue derived stages without entering raw execution.
python -m ou_irregular.week8.runner postprocess --run-dir /tmp/ou_e4_full
# Repeat the independent numerical audit without modifying the archive.
python -m ou_irregular.week8.runner audit --run-dir /tmp/ou_e4_full > /tmp/e4_repeat_audit.json
# Regenerate all 18 PNG/PDF pairs from saved audited tables, outside the archive.
python -m ou_irregular.week8.runner replot --run-dir /tmp/ou_e4_full --output /tmp/e4_external_replot
```

A repeated audit reconstructs every source and replays the original fits at fixed
replications 0 and R-1 in every cell. `audit --no-fit-replay` still reconstructs every
source but omits those fits and cannot certify completion. `verify` and closed-stage
`run --resume` generate neither sources nor fits. Replotting likewise never samples,
refits or changes the source archive. Retain the sibling `*-figure-renderer` state
as execution evidence. On another machine, exact rendering additionally requires
the recorded renderer versions/settings; do not rewrite an old state to hide drift.

## Stages and failure handling

The sequence is raw execution -> immutable tables -> full source-specific audit ->
immutable figures -> audited_pending_visual_review -> actual visual review ->
completed. Reporting and rendering own their locks; a root-byte guard spans stage
boundaries. Audit publication and lifecycle changes hold exclusive run locks.
Partial tables/figures are adopted only when they exactly equal the corresponding
original-evidence candidates. A saved full audit is verified rather than silently
replaced. Rendering interruption cannot trigger new scientific observations/fits.

An incomplete raw run can resume using the same identities and deterministic
partitions. Source/pipeline exceptions block completion and are retained for
investigation. Boundary or unavailable-interval outcomes remain in scientific
accounting under the existing rules; they are not replaced by successful draws.
A closed run is read-only. An incompatible requested validation/full mode is rejected.
Unknown artifacts, malformed schemas, changed source/runtime dependencies, symlinks,
conflicting immutable outputs and unexplained root changes are not repaired silently.

The final seal binds every artifact except the lock file and root manifest, which
contains the seal. The review binds the exact figure manifest and complete PNG set.
After writing the review, the runner rechecks numerical/figure bindings and compares
artifact maps before and immediately at the completed transition. A later external
write invalidates the seal; the software does not claim immunity to arbitrary writers.

## Verification scope

Unit tests may use explicitly labelled synthetic review records on disposable copies
to exercise sealing and corruption rejection. Such records are not actual visual
reviews and are never used to mark the retained native validation or full E4 complete.
Native CLI validation retains its original stdout/stderr, input plan, all raw records,
tables, independent audit and figures. Test counts and validation results are reported
only after observed execution. A successful implementation checkpoint is not a
completed E4 study. Completion also requires the fixed full experiment, genuine
figure review, results interpretation, limitations and final publication.
