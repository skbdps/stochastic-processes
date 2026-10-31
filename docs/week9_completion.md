# Week 9 completion: publication evidence and figure freeze

Completed 31 October 2026. The renderer is frozen at
`dad63eced437736cc32f7a481b687e226cb4bdaf`, based on the completed research release
`b6edb85d1e7a21ed19fb54b5c71213ed04703309`.

**Status: the saved-data publication package is implemented, verified,
actually visually reviewed, externally reproduced and sealed.** The work uses
the completed original studies. It generated zero new scientific source paths
and zero model fits. Manuscript drafting, human academic review and submission
preparation remain later stages.

## Deliverables and exit checks

| Requirement | Actual completed evidence |
| --- | --- |
| Pin and verify the completed sources | [Catalog](../artifacts/week9_publication/source_catalog.json) binds 831 files to original baseline Git blobs, SHA-256, size and source membership; saved seals/audits and complete keys are checked. |
| Complete classifications and bounded guidance | [T1](../artifacts/week9_publication/tables/T1.md): 360 E2 and 42 independent E1 theta/sigma rows, an 18-row E2 group summary, and all six unresolved cases. Labels use unrounded original bounds. |
| Traceable numerical claims | [24 generated records](../artifacts/week9_publication/claims.md), with [exact bindings](../artifacts/week9_publication/claims.json) and stable CSV; all human-readable views are rederived by verification. |
| Complete publication figures and supporting evidence | [Figure manifest](../artifacts/week9_publication/figure_manifest.json): 11 new PNG/PDF pairs and 40 byte-identical original pairs, with original review and figure manifests copied alongside them. [Captions](../artifacts/week9_publication/captions.md) identify scope and uncertainty. |
| F1 distributions from saved fits | All 42,000 E1 fit keys are read; [126,000 parameter-specific ECDF rows](../artifacts/week9_publication/tables/f1_ecdf.csv.gz) retain source partition digests, replication identity and eligibility. Main selection is post-result presentation; supporting pages show all seven CV settings. |
| Mathematical and literature interpretation | [Benchmark note](week9_benchmarks.md) adopts population targets plus empirical calibration; [dated primary-source comparison](week9_literature_comparison.md) records versions and access limits, including relevant 2026 coverage work. |
| Meaningful tests and independent saved-record references | [92 tests passed, zero failures/errors/skips](../artifacts/week9_evidence/week9_test_summary.json). [Actual log](../artifacts/week9_evidence/final_week9_tests/stdout.txt) and [JUnit](../artifacts/week9_evidence/final_week9_tests/junit.xml) are retained. |
| Clean-checkout reproduction | A detached clean checkout of the renderer commit reproduced all 141 immutable files exactly. [Comparison](../artifacts/week9_evidence/external_reproduction.json) checks PNG/PDF/table/claim/caption/note bytes and records an empty clean-checkout status. |
| Actual visual review | Three AI reviewers individually inspected all 11 new PNGs and 11 Poppler-rendered PDF pages. [Actual review](../artifacts/week9_evidence/actual_visual_review.json), original component records and [byte-transfer proof](../artifacts/week9_evidence/visual_inspection/byte_identical_transfer.json) bind the final figures. This is not human review. |
| Guarded freeze and completed readback | [Finalization result](../artifacts/week9_evidence/finalization_result.json), [completed lifecycle](../artifacts/week9_publication/package_state.json), [final readback](../artifacts/week9_evidence/final_readback.json), and an identical-review idempotence check. Renderer files must match their recorded committed bytes. |
| Preserve earlier science and failed evidence | [Baseline preservation](../artifacts/week9_evidence/source_preservation.json) and [failure dispositions](../artifacts/week9_evidence/README.md). The only intentional prior-file changes are the project README and artifact index. |

Immutable package-manifest SHA-256:
`370f30e1628821d90cfa42568a60d51dc2c80eb8da71c7036f93c4979d03efbe`.
The package manifest excludes only its own digest carrier, the lifecycle and
the separately hashed actual visual review. The lifecycle records the renderer
commit, review hash and completion time. The verifier detects changed or
unexpected files, missing PNG/PDF roles, altered accounting, changed captions,
and rehashed publication statements that disagree with pinned source data.

## What the frozen results say

E1 retains 19 within-criteria, 22 supported-breach and one unresolved theta/sigma
classification. E2 retains 160 within, 195 supported breaches and five unresolved.
These are parameter/design classifications, not an overall method score or
nominal-95%-coverage certificate. E1 and E2 keep independent study identities at
overlapping nominal settings.

Exact theta's five E2 breaches occur at n=250,q=0.1 and are driven by positive
relative bias, not the historical prediction of coverage below 80%. Coverage
there is 93.6–95.3%, while relative bias is 17.1–22.6%. Irregularity does not
uniformly worsen all metrics, and exact is not uniformly lowest in observed
finite-sample RMSE. Nonmonotone classifications require discrete F7 maps rather
than a maximum-safe-CV curve. These statements and counterexamples trace to
their saved rows in the generated ledger.

The clearest E4 matched-marginal difference remains PFML mu: MMPP/IID RMSE
ratios 1.083874 at r=5 and 1.157863 at r=20, with coverage lower by 2.70 and
3.85 percentage points. Their pointwise Monte Carlo intervals exclude the null;
the comparisons do not establish a mean-bias shift. The full set of methods and
parameters remains available. Matching one-gap marginals does not fix realized
span or establish equal variance. At r=20, MMPP-derived and IID gaps share the
declared floor transformation.

Empirical SD versus RMS reported SE uses the same interval-valid subset and is
descriptive. PFML/Euler intervals retain the registered model-based covariance.
No robust covariance method, general dependent-gap CLT or departure from a
fully matched first-order coverage theory is established. The proposal's
population-target/empirical-calibration fallback is explicitly adopted.

## Actual validation scope and preserved failures

The new renderer runs in a separate CPython 3.12.10 macOS ARM64 environment.
Its complete installed closure and thread settings are recorded in the
[build record](../artifacts/week9_publication/build_record.json) and
[renderer lock](../requirements-week9-lock.txt). This environment is not relabelled
as the original Linux CPython 3.12.14 scientific environment.

The 92 passing Week 9 cases cover all 402 headline theta/sigma classifications,
all unresolved identities, denominator accounting, original ECDF samples,
source-bound claims, actual plotting artists and asymmetric uncertainty,
source/output corruption, exact inventories, portable note links, and rejected
incomplete reviews. An independent reference calculation also recomputes the
PFML mu r=5/r=20 contrasts and SD/RMS-SE quantities directly from saved E4 fit
records. The tests forbid scientific source generation and fitting during the
publication build. They do not create a positive synthetic publication review.

Development failures are retained: the first publication build changed the
source catalog's JSON serialization and was rejected by its exact digest;
an initial CLI test harness unnecessarily imported a historical YAML dependency;
and a later test found insertion-order-dependent claims CSV columns. These were
reporting/test-harness defects. Their original observed records, available logs,
exact source snapshots and clearly labelled replay evidence remain in
[week9_evidence](../artifacts/week9_evidence/). The final 92-case run follows the
repairs and additional inventory/accounting checks. No scientific inputs were
replaced or repeated to repair these issues.

An additional historical test probe on the pre-existing Python 3.14 Mac
environment remains **failed**: 78 passed, seven failed, 13 fixture errors and
three existing design-inapplicable skips. Eighteen affected cases reject the
archived Linux numerical environment; two assume extended longdouble precision,
while this host supplies 64-bit longdouble. The [disposition](../artifacts/week9_evidence/legacy_saved_tests/disposition.json),
[full log](../artifacts/week9_evidence/legacy_saved_tests/stdout.txt) and
[JUnit](../artifacts/week9_evidence/legacy_saved_tests/junit.xml) retain those
outcomes. Old source and test bytes are unchanged. This probe is not called a
successful historical replay, and its tests/thresholds were not weakened.
The completed Linux release evidence remains documented in
[Week 8](week8_completion.md). The narrowly scoped Week 9 data-only freeze is
supported by its own passing tests, source checks and exact external rebuild.

## Reproduction and Week 10 handoff

Follow [the exact commands](week9_reproduction.md) from a full-history checkout.
The external rebuild checks all immutable derivative bytes in the recorded
renderer environment; byte identity is not promised across different font,
platform or numerical-library runtimes. The original audit/replay limitations,
historical operational incidents, all unavailable records and unresolved
statistical findings remain intact.

The manuscript handoff consists of T1, the F1–F8 figure families and complete
supporting figures, generated captions/claims, mathematical definitions,
benchmark scope, literature comparison and reproducibility evidence. Human
academic assessment, final journal-specific layout, further full-text reading
where access was limited, manuscript drafting and pre-submission literature
refresh remain later tasks. No ongoing CI workflow is introduced.
