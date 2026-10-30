# Guarded Week 7 raw execution

This module is a technical reconstruction refinement after the interrupted
unpublished Week 7 implementation. One guarded orchestrator handles both fixed
E3 follow-up and paired-floor studies through the unchanged durable cell engine.
It changes no scientific design, seed namespace, estimator, inference gate or
partition format. It does not supply reporting, audit, figures or a completed
week. The governing scientific definitions remain in
`week7_exploratory_protocol.md` and `week7_mathematics.md`.

## Interface and design authority

```python
require_design(kind, design)
source_policy(kind, design)
preflight(kind, design, output_dir, *, extra_files=(), resume=False)
run_raw(kind, design, output_dir, *, extra_files=(), workers=1,
        partition_size=500, resume=False, progress=True)
```

`kind` is `followup` or `floor`. The only full designs are the unchanged
`selection.followup_design()` and `selection.floor_design()`. A development
design must equal the corresponding `.validation(R)` object in every serialized
field, including name, seed, scope, all cells, their order and replication count.
Validation retains the full physical grid and uses the already specified
seed+1,000,000 namespace with a small fixed R. There is no arbitrary design,
source, inference-service or unfrozen-full-run override in this interface.

Full exploratory execution always requires committed dependency bytes. Before
any directory creation or source call, provenance capture verifies the explicit
dependency set and numerical environment. Resume also verifies every recorded
file against its actual bytes in the recorded Git freeze, without requiring
the recorded commit to equal the current HEAD when unrelated later commits
have been added. The file manifest includes the composed cell engine, storage,
provenance, source adapter, selection protocol and reused root-state guard,
including its imported dependencies. Extra entry-point dependencies must be
supplied explicitly and must match on resume.

`preflight` is read-only preflight of a proposed **raw execution mutation**:
it creates no directory, lock, source or manifest and does not alter lifecycle
history. It is not a generic archive auditor; closed-state and synchronized-
output policies intentionally still apply. Its `Preflight` result contains `output`,
`metadata`, `expected_metadata_bytes`, current captured `provenance`, verified
`cells`, `cached_completion` and `require_frozen`. On a resume, the recorded run
provenance remains in `metadata`; the fresh capture checks the current exact
dependency membership and runtime.

## Scientific metadata and source composition

Follow-up runs retain project `ou-irregular-experiments`, no scientific
overrides, and the exact historical `PROCEDURE`. They remain compatible with
the registered-source artifact auditor after their separate reporting stage.
The cell worker continues to construct `OuSampleSource(design.seed)` and the
unchanged registered Wald services.

Floor runs use project `ou-week7-paired-floor` and the explicit override
`paired_floor_source`. Every inference setting remains equal to `PROCEDURE`;
only `source_address` changes to
`SeedSequence(root,spawn_key=(4,*8_registered_floor_anchor_cell_words,rep)).spawn(2)`.
The worker receives `PairedFloorSource(design.seed)` and otherwise uses the
same cell execution and inference code. These outputs require the separate
paired-source auditor; they must not pass as ordinary registered-source runs.

The new `execution_kind`, `execution_contract`, `require_frozen` and
`source_policy` metadata are not excluded by the existing
`semantic_metadata_digest`. Consequently the audit's semantic binding includes
them. The contract is `week7_guarded_cell_execution_v1`. Source policy records
the policy name, root seed, complete actual-to-anchor ID map, anchor floor,
stationary initialization, absence of floor renormalization, and planned counts:

| Kind | Actual source realizations | Independent coupling groups | Fit records |
|---|---:|---:|---:|
| Follow-up | 2R | 2R | 6R |
| Paired floor | 12R | 4R | 36R |

The floor policy name is the source adapter's fixed `COUPLING_POLICY`; the
follow-up name is `registered_ou_source_v1`. A follow-up cell maps to itself;
a floor cell maps to its same-CV epsilon=1e-6 anchor. Planned counts are not
reported as completed counts. Actual completion is separately read back from
every durable cell.

## Pure preflight and durable evidence

Before any mutation, preflight checks the resolved design's exact canonical
JSON, root schema, procedure, source policy, freeze policy, explicit dependency
membership, environment and lifecycle state. Duplicate JSON object keys,
nonstandard NaN/Infinity JSON, Boolean integer counts and noncanonical manifest
representations are rejected. Root completed-cell claims must match complete
durable evidence. A root map may lag already completed cells after a crash;
that is recoverable only in an incomplete execution state.

Known root files include raw metadata/design, the established derived tables,
epsilon comparisons and their separate `epsilon_table_manifest.json`, the
ordinary table manifest, audit and review files, and the
`cells`/`figures` directories. Unknown artifacts, unsafe directory objects,
symlinks, unknown cells and unknown files inside cells block automatic
execution. Declared figure/table manifests, when present, bind their existing
file sets and hashes. A planned empty cell directory with no evidence can
start normally; a cell containing evidence without its manifest cannot.
If a figure manifest has not yet been written, ordinary PNG/PDF files can be
retained as an unfinished derived set; other file types, directories and
symlinks are rejected. This permits postprocessing recovery without claiming
figure completion. The later plotting/audit stage must validate exact planned
filenames and unchanged bytes before adopting those files. A present figure
manifest still requires its complete recorded file set and hashes.

Each existing cell is opened with `CellStore(..., resume=True)` so its
constructor verifies rather than creates it. A strict outer manifest check
adds exact field and integer-type checks. The historical storage component
then verifies the contiguous committed replication prefix, partition digests,
row/key accounting and any single next orphan partition. Partition size can
change on resume, but an orphan's original range cannot change. An orphan can
be adopted only when the unchanged worker regenerates byte-identical content;
different evidence is never overwritten.

Preflight scans **every durable partition**, including the orphan and cells
missing from the root's completed map, for the explicit string `pipeline_error`
field. Missing or malformed error evidence is not interpreted as zero errors.
Any saved nonempty pipeline error blocks automatic resumption before a source
is called. The raw evidence remains intact for implementation review. Ordinary
boundary, point-fit or Hessian/interval statistical failures with an empty
pipeline-error field remain planned outcomes and do not trigger resampling.

## Lifecycle rules

| Saved state | Raw entry-point behavior after successful verification |
|---|---|
| `running_incomplete`, `incomplete` | Resume the durable prefixes under the exclusive root lock. |
| `raw_complete_pending_audit` | Return verified cached metadata, without creating a lock or changing history. |
| `implementation_review_required` | Return cached metadata only if every cell is complete and all pipeline-error counts are zero; this permits postprocessing recovery while preserving its original error fields. |
| `audited_pending_visual_review`, `completed` | Reject entry into raw execution. |

Unknown states and evidence without a root manifest are rejected. Raw cached
completion requires the full completed-cell map and zero pipeline errors; it
does not reset the status, append an attempt or call the cell engine. A
completed cell missing from an incomplete root map is instead recovered through
the unchanged worker's cached complete-cell return, without source generation.

An executing attempt appends its worker/partition settings and timestamps to
the preserved attempt/resume histories. These orchestration settings never
enter scientific seed addresses. Completion can produce only
`raw_complete_pending_audit` or `implementation_review_required` when newly
saved pipeline errors require review. Exceptions retain an honest `incomplete`
state and their diagnostic fields. No raw result is labelled `completed`.

## Exclusive orchestration and root-state drift

The orchestrator acquires the existing nonblocking `run_lock` and rechecks
preflight after locking. A spawn-based process pool assigns disjoint cells to
workers; only the parent writes root metadata. Each worker uses the unchanged
`_execute_cell`, `CellStore` and registered services. Already committed source
outcomes are never resampled. Pending worker tasks are cancelled on an observed
orchestration exception where cancellation is possible; already running cell
workers can finish their durable work before process-pool shutdown returns.

The existing `_ManifestGuard` compares expected root bytes before and after
root replacements, before submitting/executing cells, before completion
verification and in the error handler. Unexpected state is preserved as
expected/observed incident evidence in an adjacent directory. The handler does
not overwrite the changed root with a stale cached copy. A bound incident
blocks subsequent automated mutation even if someone restores an earlier root
file. Advisory locks and byte checks cannot prove the absence of every
noncooperating external writer or explain the historical runtime incidents.

Output is rejected inside another saved experiment or through a symbolic-link
path. If `/usr/local/scripts/sync_share.sh` is present, mutable outputs under
`/workspace/scratch` are rejected before creation. Isolated `/tmp` outputs are
used instead. No platform synchronization service is modified. These are
operational containment measures, not a claim that the earlier incidents have
been causally diagnosed or repaired.

Every final raw-completion decision rereads all cell manifests/partitions and
rechecks dependencies before the guarded root commit. Further scientific
completion still requires the separately implemented reporting, independent
audit, reproducible figures, actual visual review and publication steps.
