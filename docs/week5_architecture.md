# Week 5 execution and reproducibility architecture

Prepared before E1 execution on 30 October 2026 (India time). Week 5 adds a
separate experiment layer around the unchanged Week 3 estimators and Week 4
inference. No historical runner is extended to accept a later phase.

## Responsibilities

| Module | Responsibility |
| --- | --- |
| experiments/design.py | Immutable scientific cells and run identities; registered cell hashing; separate development streams |
| week5/config.py | Reject any change to the frozen E1 configuration |
| experiments/summaries.py | Complete planned-key accounting, availability, coverage, classification and matched-subset SE diagnostics |
| experiments/comparisons.py | Independent regular references and paired estimator comparisons |
| experiments/provenance.py | Explicit dependency hashes, configuration identity, runtime versions and committed-code freeze guards |
| experiments/storage.py | Flushed immutable partitions, atomic manifests, execution lock and verified crash recovery |
| experiments/execution.py | Per-cell scheduling and composition of the unchanged numerical services |
| experiments/analysis.py | Reproduce all tables from saved raw partitions, processing one cell at a time |
| experiments/audit.py | Independently check complete keys, arithmetic, numerical records and predetermined replay |
| week5/plots.py | F2--F4 and companion availability/mu figures from saved tables |
| week5/runner.py | Week-specific command-line lifecycle and completion gate |

All paths are relative to ou_irregular/. The mathematical definitions are in
[week5_mathematics.md](week5_mathematics.md). Components listed here are the
intended final Week 5 architecture; implementation/test status and actual
headline results are recorded separately in the progress report.

## Source addressing and failure retention

Scientific cell IDs retain the published Week 4 SHA256 convention. For every
replication, OuSampleSource uses the registered root seed and spawn key
(4, *eight_32bit_cell_id_words, replication), then spawns separate gap and OU
path streams. Worker count, execution order and partition boundaries do not
enter that address. A validation design uses a different root seed and a
validation scope; it cannot be labelled a headline design.

Every source yields exact, PFML and Euler records, including failed sources,
failed point fits, invalid intervals and unexpected dependency exceptions.
Inference exceptions preserve the available point candidate. The trace adapter
retains the derivative object returned by the unchanged finite-difference
provider. It neither recomputes the objective nor changes the returned object.

## Persistence and resume

Raw rows are canonical JSON lines in deterministic gzip files, partitioned by
cell and contiguous replication range. JSON null represents unavailable numeric
information. Complete point and interval diagnostics, candidate parameters,
NLL, covariance, fine Hessian/gradient, normalization, source hashes and gap
diagnostics remain available for subsequent E3 investigation.

Each partition is written to a temporary file, flushed with fsync, renamed,
and followed by a flushed atomic manifest update. Existing different partition
bytes are never overwritten. An operating-system advisory lock prevents two
orchestrators from writing the same run concurrently. Workers write different
cell directories and use spawn to avoid inheriting that lock or BLAS state.

Resumption requires the same design, registered procedure, explicit dependency
hashes and recorded numerical environment. Every committed partition's hash,
range and planned keys are checked. A saved statistical failure is an outcome,
not a reason to draw another path. A previously committed partition is never
refitted during ordinary resume.

A crash between partition rename and manifest update can leave one orphan.
Recovery discovers its original range even if the requested future partition
size has changed. Deterministic replay must reproduce its exact bytes before
it is adopted. A mismatch, overlap, unexplained raw partition, multiple orphans or orphan
inside a completed cell causes failure. Recovery events are recorded. An
interruption inside an uncommitted in-memory block can require its deterministic
re-execution; that repeats the same inputs and does not resample failures.

The run manifest distinguishes running_incomplete, incomplete,
raw_complete_pending_audit, audited_pending_visual_review,
implementation_review_required and completed.
Raw completion never certifies the scientific result by itself. Unexpected
source/fit/inference dependency exceptions block completed status and remain
visible even if other scheduled paths finish.

## Explicit provenance and historical preservation

A run records exact source/configuration/dependency file hashes, the published
preregistration commit and digest, the resolved design, its code commit, numerical
package versions and Python/platform information. Headline execution requires
every listed dependency to have been committed before execution. Worker count,
partition size, elapsed time and resume attempts are separate execution records.

The dependency manifest names actual files rather than recursively requiring
the entire future package tree to remain unchanged. Added unrelated later
packages therefore do not invalidate a past run; changes to its dependencies
still do. Explicit source verification is repeated at completion and during
audit/replot. Package-version fingerprints do not identify every binary on a
machine, so cross-platform bitwise reproducibility is not promised.

The historical Week 4 auditor retains its original strict recursive-map rule.
It is run in a detached checkout of the published Week 4 revision, with a
separate report. Its archived manifest and numerical gates are not rewritten.
All original Week 3/Week 4 scientific source and archived evidence remain
byte-for-byte preserved in the continuation branch.

## Artifact and completion checks

The analysis layer checks the exact set of cell directories and every planned
cell/replication/estimator key. It exports parameter summaries, independent
regular-reference contrasts, paired contrasts, spacing summaries, numerical
diagnostics, status counts and the complete failed-fit log. Raw JSON partitions
remain authoritative; CSV tables and figures are derived outputs with hashes.
Table calculations process one cell at a time, while the current exporter
retains the accumulated failed-fit log until writing its compressed CSV. Memory
can therefore increase with widespread failures; it is not strictly bounded
by one cell in that case.

The independent auditor recalculates arithmetic from raw records and checks
validity/status agreement, source identity, numerical gate records, covariance
transforms and endpoint construction. It replays fixed replication IDs in
every cell. Figure review checks labels, uncertainty conventions, denominators,
missing outcomes and readability. The completion seal binds the numerical
artifacts and scientific root metadata to the independent audit, then binds
the reviewed figure manifest to the same design. Only explicit lifecycle
fields may change between audit and seal. The progress report then records the actual
findings and limitations. GitHub publication advances the named continuation
branch only after checking its expected head; it never force-pushes over newer
work. Each completed week is published before the next week's implementation.

## Runtime qualification

The available Python is 3.12.14. The five numerical/plotting/YAML dependencies
match the pinned versions, and pytest is 9.1.1. The runtime initially lacks the
notebook-only nbformat, IPython and ipykernel requirements; package installation
is unavailable through the current shell proxy. Those missing packages do not
enter this experiment's numerical execution or tests. Their absence is recorded
in the environment manifest. No notebook execution is claimed. Historical
numerical replay uses the same pinned numerical package versions.
