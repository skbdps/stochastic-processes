# Week 7 standalone runner and completion contract

Prepared 30 October 2026, after the published plotting checkpoint and before
the full Week 7 studies. This is an implementation contract. The dated
[exploratory protocol](week7_exploratory_protocol.md), unchanged registration
and adopted amendment retain their scientific authority.

## Fixed studies and dependencies

`ou_irregular.week7.runner` composes the published guarded executor, saved-record
reporter, independent auditors and deterministic renderer. It changes no source
generator, estimator, inference rule, denominator, uncertainty formula, seed,
selection or replication count. Full `followup` means the two specified cells
with R=10,000 and seed 20261007. Full `floor` means all twelve specified variants
with R=2,000 and seed 20261009. Validation uses the complete same cells, an
explicit positive integer R and the previously specified seed+1,000,000.

One explicit dependency union includes the entire plotting dependency set,
`ou_irregular/week7/runner.py` and this note. The same union is supplied to raw
execution, reporting, the floor auditor and rendering. All runner entry points
require exact membership, not only matching hashes for whichever files happen
to be present. Full exploratory execution requires committed bytes; publication
of that reviewed freeze precedes affected scientific observations. Later
unrelated Git commits do not invalidate the unchanged recorded dependency set.

The registered-source follow-up uses the original numerical audit. The floor
study uses the independently reviewed paired-source audit, including all source
array hashes and the fixed endpoint fit replay. The runner additionally checks
the exact ordered failure-export schema, including zero-row exports. Duplicate,
missing, surplus and reordered CSV headers are refused without changing the
historical auditor or the source record format.

## Lifecycle and resumption

The ordinary sequence is raw execution, complete tables, independent full
audit, immutable audit publication, saved-evidence figures, and a guarded
transition to `audited_pending_visual_review`. Successful numerical execution
and rendering alone do not complete a run. Actual inspection of every PNG
must precede a separate `seal` call with a reviewer and the exact current
figure-manifest hash.

The new runner never calls the historical unguarded workflow writers. Reporting
and rendering own their execution locks; the runner takes separate locks for
audit publication and lifecycle transitions. Expected root-manifest bytes are
checked across the intervening stages and before/after root writes. Unexpected
root drift preserves the observed bytes in an incident directory and blocks
automatic mutation. An advisory lock cannot exclude an external writer that
ignores it. These safeguards do not establish the cause of the historical
execution-state incidents.

An ordinary postprocessing exception leaves the true `raw_complete_pending_audit`
state and every partial artifact unchanged; the exception is raised to the
caller and must be retained in the invocation log. It does not reset execution
history or fabricate a completed state. Existing `implementation_review_required`
evidence requires a documented disposition and is not automatically repaired
or retried by this runner. Statistical boundary outcomes and unavailable
intervals remain ordinary retained results; a pipeline exception prevents
completion.

If `audit.json` already exists, postprocessing verifies its full replay scope,
endpoint identities, exact input hashes, planned counts, source policy and
numerical tolerances before reusing it. It does not regenerate raw observations
or replay fits. This is necessary for partial plotting: the numerical floor
auditor intentionally rejects an unfinished figure directory, whereas the
saved-audit renderer can verify and adopt each byte-identical partial figure.
Changed declared files or an invalid saved audit are rejected, never replaced.

`run --resume` on an audited or completed archive performs read-only verification
and leaves timestamps, execution attempts, resume history and bytes unchanged.
Repeated sealing of a completed archive likewise verifies and returns the
existing state, requiring the same visual-review record. Raw-complete resumption
uses only postprocessing; incomplete raw work uses the guarded cell executor.

## Verification and sealing

`verify` rechecks saved source/dependency identities, design, complete accounting,
audit bindings, exact table and figure sets, hashes and the completed seal where
present. It generates no sources or fits. `audit` is a separately requested
independent numerical recheck and prints a report without changing archive
files. Its `--no-fit-replay` option omits endpoint fits; for the floor study it
still reconstructs every paired source. Such a report cannot certify full
completion. Saved audit publication always uses the complete replay mode.

Sealing requires a passed full unchanged audit, its semantic metadata and
scientific-artifact bindings, complete renderer verification, and a real visual
review naming every PNG exactly once. The root seal binds all archive files
except the root manifest itself and the execution-lock file. A review file
left by an interrupted seal can be adopted only with identical bytes. The
numerical and figure bindings are reverified after review publication, enclosed
by equal artifact hash maps; the map is checked again before the root write and
the completed seal is read back before success is returned. This closes the
review-write gap found by independent fault injection. A later external write
invalidates the exact seal and is refused on readback, never folded into a new
silently accepted artifact map. The
renderer sibling `<run-name>-figure-renderer/renderer_state.json` is preserved
alongside the run for partial-render recovery; completed figure manifests also
embed its identity and can be verified without that sibling.

## Commands

Use the recorded Python/scientific environment. Set `OPENBLAS_NUM_THREADS=1`
and `OMP_NUM_THREADS=1` before Python, including during resumption and readback.
Mutable outputs must be outside the synchronized workspace and disjoint from
all previous experiment roots. For example, run the following from the isolated
repository checkout, with fresh output directories:

```bash
python -m ou_irregular.week7.runner run --kind followup --output /tmp/ou_w7_followup_validation --validation-replications 2 --workers 2 --partition-size 1
python -m ou_irregular.week7.runner run --kind floor --output /tmp/ou_w7_floor_validation --validation-replications 2 --workers 2 --partition-size 1
```

After validation/review and publication of the exact full code freeze:

```bash
python -m ou_irregular.week7.runner run --kind followup --output /tmp/ou_w7_followup --workers 2 --partition-size 500
python -m ou_irregular.week7.runner run --kind floor --output /tmp/ou_w7_floor --workers 12 --partition-size 500
```

Append `--resume` with the identical design/environment only when continuing
an existing invocation. For an existing validation include its same validation
replication count. Readback and sealing examples are:

```bash
python -m ou_irregular.week7.runner verify --kind floor --run-dir /tmp/ou_w7_floor
python -m ou_irregular.week7.runner seal --kind floor --run-dir /tmp/ou_w7_floor --review /tmp/actual_floor_visual_review.json
python -m ou_irregular.week7.runner verify --kind floor --run-dir /tmp/ou_w7_floor --require-completed
```

The review must be authored after inspection, not inferred from successful
rendering. It contains `reviewer`, `reviewed_pngs`, `all_panels_legible: true`,
and `figure_manifest_sha256`; additional explanatory review evidence is allowed.

## Selected-path diagnostics remain a separate reviewed entry point

The already published `path_diagnostics.run_selected` accepts the fixed
repository selection and sealed E1/E2 roots. It has its own lifecycle,
immutable source payloads, frozen dependencies and complete-selection checks.
It is called directly; the new study runner does not wrap it, monkeypatch its
dependency list or claim that its raw provenance includes the new runner.
Use `plots.write_selected_figures` afterward with a separate figure directory.
Preserve the diagnostic evidence even when its recorded status requires review.

The full scientific workflow still consists of all fixed selected diagnostics,
both complete exploratory studies, numerical and visual review, a one-page
findings memo and detailed dispositions, then publication of completed Week 7.
None of those results is certified merely by this module or its tests existing.
