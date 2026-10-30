# Week 8 guarded fitted execution

The execution module composes the reviewed E4 source with the unchanged
Week 3/4 estimators, Hessian/Wald services and durable partition store. It
accepts only the complete fixed seven-cell full R=2000 or fitted-validation
R=4 design. It does not alter the source law, root seeds, observation count,
point eligibility or interval gates.

## Public boundary

`preflight(design, output, extra_files=(), resume=False)` inspects a proposed
raw execution without creating directories, locks, sources or attempts.
`run_raw` adds worker count, partition size and progress controls; those are
execution settings and never enter scientific source addresses.
`readback(design, output, extra_files=())` performs read-only identity,
provenance, layout and cell inspection, including later audited/completed
states. Readback does not authorize mutation or independently certify the
numerical summaries, audit or visual seal.

Full exploratory execution requires committed dependencies. Provenance binds
the exact source-file set, numerical environment, complete resolved design and
recorded Git freeze. Later stages supply their additional dependencies before
execution; the full-run freeze therefore covers the complete pipeline.
E4 declares its own source policy and namespace 8 while retaining every
existing inference setting. Its registered source is not a caller-supplied
scientific override.

## Preservation and resumption

Mutable output must use an isolated ordinary path. Symlinks, scientific-archive
descendants, preserved incident blocks and synchronized workspace locations
are rejected by the mutation boundary. Read-only verification of closed
archives remains available.

The root manifest is guarded against byte changes at every lifecycle write.
An unexpected outside write preserves incident evidence and blocks automatic
continuation instead of overwriting the observed root. Only the parent process
writes the root; each worker owns a disjoint scientific cell. Durable
partitions retain every planned estimator record, including statistical
failures and pipeline errors.

Resume verifies committed prefixes and any single identifiable orphan
partition through the existing CellStore contract. It never resamples a saved
pipeline error. Raw-complete evidence is returned without a new attempt, lock
or source allocation. Audited or completed evidence cannot be reset through
raw execution; later read-only stages must verify it in place.

Raw completion is `raw_complete_pending_audit` only when every planned cell is
durable and no pipeline error is present. Statistical point/interval failures
remain in the planned records and do not become pipeline errors merely because
their estimates or intervals are unavailable. Pipeline errors require review
and cannot receive scientific completion.

The known derived-artifact layout is declared in advance for the later table,
audit, figure and review modules. This raw module does not create those outputs
or mark the study completed. Its focused tests use literal records and mocked
workers; the native R=4 validation and full E4 study are separate later gates.
