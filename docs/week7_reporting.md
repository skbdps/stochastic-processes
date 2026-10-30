# Week 7 reporting from saved raw records

The reporting module derives tables from a complete, verified raw run. It does
not simulate observations, fit models, retry failures, alter inference gates,
write root lifecycle metadata, audit a completed study, create figures or seal
an archive. It composes the historical standard table writer unchanged and
adds the already reviewed paired epsilon table for the floor study.

## API and prerequisites

```python
preflight_reporting(kind, design, output_dir, *, extra_files=())
raw_input_binding(output, design, metadata)
iter_epsilon_tables(output, design, identity)
write_tables(kind, design, output_dir, *, extra_files=())
```

`kind` and `design` follow the guarded Week 7 execution contract. The exact
full follow-up/floor design or its exact separate-namespace validation design
is required. Reporting accepts only `raw_complete_pending_audit` or a complete,
zero-pipeline-error `implementation_review_required` state. The latter permits
recovery from a postprocessing failure without erasing its original error
fields. Incomplete, audited and completed archives cannot enter this writing
interface. There is no automatic recreation of a declared missing table in a
closed archive.

`FILES` is an explicit union of the execution dependencies, historical analysis
and comparison dependencies, paired-floor comparison module, this reporting
module and this technical note. A run must have captured this complete set,
plus any supplied later entry-point files, before raw execution. The reporter
requires exact dependency membership and current/recorded file bytes and the
same numerical environment. Full exploratory runs require the recorded Git
freeze. An earlier validation run lacking reporting dependencies is preserved;
its provenance is not rewritten to make it appear eligible for a new reporter.

The read-only preflight composes execution preflight, requiring complete
Cartesian raw keys, complete cell manifests, zero saved pipeline errors,
canonical design, source policy and freeze identity. It then snapshots raw
input identity and validates any existing table declarations. It creates no
output files or locks. The writing function subsequently acquires `run_lock`
and repeats preflight before calculations or writes.

## Standard tables and failure retention

`experiments.analysis.write_tables` remains unchanged. It produces these exact
seven artifacts and their existing `table_manifest.json` format:

- `parameter_summary.csv`
- `regular_comparisons.csv`
- `paired_comparisons.csv`
- `spacing_summary.csv`
- `outcome_counts.csv`
- `numerical_diagnostics.csv`
- `failed_fits.csv.gz`

Every original invalid point, unavailable interval and numerical candidate
remains in its raw record and, when appropriate, in the failed-fit table. The
original planned-denominator, conditional/operational coverage, Monte Carlo
uncertainty and practical-classification rules are unchanged.

When a standard manifest exists, its membership must be exactly these seven
files, with exact `rows`/`sha256` descriptors. Every declared file must exist
and match its digest, CSV structure and row count. Fixed-shape tables also
retain their full design counts. A missing or changed declared file blocks
reporting before any new output is written. A partially written *set* of
unmanifested standard files may be adopted only when the unchanged writer
reproduces their exact bytes; it never overwrites different prior evidence.

## Bounded epsilon calculations

For floor runs, the reporter additionally computes `epsilon_comparisons.csv`
with the exact ordered 71-field `floor_comparisons.COLUMNS` contract and all
72 planned comparison rows. `iter_epsilon_tables` reads one complete CV triple
at a time, calls the reviewed comparison function on that complete subdesign,
releases the raw combined DataFrame, and yields only its 18 summary rows. Raw
records from the prior CV are not retained in the suspended generator frame.
The standard writer retains its own historical memory behavior; this bounded
claim concerns the additional epsilon-comparison pass.

All four epsilon groups are calculated and validated **before any table
output**. Consequently malformed coupling metadata or incomplete common-key
accounting cannot first create standard outputs and only later fail the paired
preflight. Only the resulting 72 summary rows are concatenated, in deterministic
CV/epsilon/estimator/parameter order. Their complete key set is verified before
CSV serialization. No fit, statistical failure or cell is resampled or omitted.
Follow-up runs reject epsilon CSV/manifest artifacts rather than silently
ignoring an unexpected paired-floor result.

The unchanged `csv_bytes` serializer fixes 17-digit numeric formatting and
newline representation. The epsilon table uses immutable atomic publication;
an existing unmanifested table is adopted only if its bytes match exactly.
Previously declared epsilon artifacts must already exist unchanged.

## Raw-input binding and epsilon manifest

`raw_input_binding` records:

- The saved provenance run-identity hash and complete design hash.
- The resolved-design file's byte hash.
- The exact source policy, including the actual-to-anchor map and planned counts.
- Every cell-manifest byte hash and its complete ordered partition map, including
  filenames, ranges, row counts and hashes verified against actual partition bytes.

It excludes root lifecycle bytes such as raw/audited/completed status and later
review timestamps. Those fields legitimately change after reporting; they
must not make the table's scientific input identity unstable. Complete cell
manifests and raw partitions do not change during later reporting/audit/sealing.

`epsilon_table_manifest.json` has the exact schema:

```text
schema_version = 1
project = ou-week7-paired-floor-reporting
reporting_contract = week7_saved_raw_tables_v1
run_identity
design_sha256
raw_inputs
raw_input_sha256
table = {filename, rows, columns, sha256}
difference_direction = epsilon_minus_registered_floor
ratio_direction = epsilon_over_registered_floor
```

`raw_input_sha256` hashes the canonical serialized raw-input descriptor. The
table descriptor fixes `epsilon_comparisons.csv`, all 72 rows, the ordered
71-column contract and its exact CSV hash. The manifest contains no newly
generated timestamp, so resuming unchanged evidence reproduces identical bytes.
This manifest binds table provenance; it is not a substitute for the later
independent numerical/source audit.

## Root and raw drift, partial outputs and completion boundaries

The reporter never writes `run_metadata.json`. A composed `_ManifestGuard`
checks its exact expected bytes before and after calculation/write stages and
in the exception handler. An unexpected root change preserves the expected
and observed bytes as an adjacent incident and blocks automatic continuation.
There is no stale `finally` write that could erase foreign state.

An advisory root lock cannot rule out a noncooperating raw-file writer. The
reporter therefore recomputes and compares the entire raw-input binding around
table production, including actual partition hashes, and rechecks dependencies
and layout. A changed partition plus a consistently updated cell manifest is
still detected because the initial binding has changed. The epsilon manifest
is not published after such a change. Existing raw files and any immutable
partial derived outputs remain available for review; the reporter does not
roll them back, delete them or rewrite the root to conceal the failed attempt.

`write_tables` returns `{"standard": ..., "epsilon": ...}`; epsilon is `None`
for a follow-up run. A successful return means the requested derived tables
and their manifests were verified. Root lifecycle status remains byte-for-byte
unchanged. Independent auditing, figures, actual visual review, scientific
disposition and publication are separate later steps.
