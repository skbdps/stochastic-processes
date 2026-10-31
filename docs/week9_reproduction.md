# Week 9 saved-evidence reproduction

The Week 9 command creates a derivative publication package from the already
completed E1, E2, E3 and E4 archives. It verifies their pinned baseline bytes,
saved seals, planned keys and reporting contracts. It does not reconstruct OU
paths, rerun a numerical audit or fit an estimator. The original audits and
their endpoint-only fitting replay limits remain those in the source reports.

Use a full-history checkout containing baseline
`b6edb85d1e7a21ed19fb54b5c71213ed04703309` and the Week 9 code. The source catalog
is `configs/week9_publication.json`; its own digest is bound in the input module.
Each archived file is checked against its original Git blob, SHA-256 and size.
The original scientific dependency files and required review records are also
bound. The verification is a new data-only contract, without changing the old
runtime-dependent verifiers or metadata.

## Renderer environment

The actual Week 9 renderer uses CPython **3.12.10**, macOS ARM64. This is
separate from the original Linux CPython 3.12.14 scientific execution/replay
environment. It is suitable for the saved-data operations claimed here; it is
not a cross-platform numerical-replay certificate.

`requirements-week9.txt` lists direct dependencies;
`requirements-week9-lock.txt` records the installed package closure. The build
also records the complete installed distribution versions, Python, OS,
architecture and relevant thread settings. Recreate this environment for
byte-identical figure reproduction. Other hosts can verify pinned data, but
byte-identical PNG/PDF rendering is not promised across different runtimes,
font libraries or platforms.

From the repository root, create an isolated environment outside the tracked
source tree:

```bash
python3.12 -m venv /tmp/ou-week9-render
/tmp/ou-week9-render/bin/python -m pip install -r requirements-week9-lock.txt
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
unset MKL_NUM_THREADS VECLIB_MAXIMUM_THREADS
```

Check the interpreter version and capture it with the validation evidence.
The original experiment requirements and provenance are retained unchanged.

## Build and read-only verification

Choose destinations that do not exist. The build rejects a conflicting output,
input overlap, symlinks or a changed source catalog. All historical archives
are hashed before and after generation.

```bash
/tmp/ou-week9-render/bin/python -m ou_irregular.week9.cli build \
  --sources configs/week9_publication.json --output /tmp/ou-week9-first
/tmp/ou-week9-render/bin/python -m ou_irregular.week9.cli verify \
  --run-dir /tmp/ou-week9-first
```

The first build remains `pending_visual_review`. The package includes complete
T1 tables, the six unresolved E1/E2 classifications, full plotting data, claims,
captions, new figures and identical copies of the original figures. The F1
companion table uses deterministic gzip and includes all 126,000
parameter-specific values from 42,000 retained E1 fit records. It is not a new
fit archive. The source partition digest and fit key are retained for every
sample.

The read-only verifier checks exact file membership and hashes, source
preservation, renderer source identity, table arithmetic, machine-readable and
readable claims, caption definitions and transferred source-review identities.
It leaves the package and source archives unchanged.

## Independent external rendering

```bash
/tmp/ou-week9-render/bin/python -m ou_irregular.week9.cli build \
  --sources configs/week9_publication.json --output /tmp/ou-week9-second
/tmp/ou-week9-render/bin/python -m ou_irregular.week9.cli compare \
  --left /tmp/ou-week9-first --right /tmp/ou-week9-second
```

The comparison first verifies both packages, then compares every immutable
output, including PNG, PDF, tables, numerical claims, captions and notes.
Lifecycle timestamps and the separately bound visual review are excluded from
the deterministic content manifest. Preserve the comparison result outside the
sealed package as validation evidence. Failed builds and failed checks are
retained with their error records; use a new destination after a code change.

## Actual visual review and finalization

Inspect every new PNG individually and render every new PDF page for layout
inspection. Check legends, all intended series, axes, uncertainty bars, cell
annotations, captions and source-table agreement at full size and intended
publication viewing size. The review must say that it is AI review when
performed by AI; it is not human academic review.

The review JSON identifies `reviewer`, an actual timezone-aware `reviewed_utc`,
`review_type: actual_ai_raster_and_pdf_review`, and the exact current
`figure_manifest_sha256`. `reviewed_png_sha256` maps every new PNG's package
path to its digest. `reviewed_pdf_pages` maps every new PDF to its page number
(one), `source_pdf_sha256`, `rendered_page_sha256`, `render_command` and renderer
version. It also includes `all_panels_legible: true`, `blocking_issues: []`, and
actual inspection notes. The completion evidence retains the rendered pages
and their hashes. No template or unit-test review supplies publication evidence.

Original figures copied without byte changes retain their explicitly identified
original review records. A new layout always requires its own actual review.

After actual review, finalize only the new derivative package:

The renderer code, input catalog, requirement files and copied source notes must
first match a committed Git revision. Finalization verifies those bytes against
that commit and records its SHA in the lifecycle record. The moving checkout
HEAD is excluded from deterministic figure content, so later commits that add
the sealed artifacts do not invalidate otherwise identical reproduction.

```bash
/tmp/ou-week9-render/bin/python -m ou_irregular.week9.cli finalize \
  --run-dir /tmp/ou-week9-first --review-file /path/to/actual-review.json
/tmp/ou-week9-render/bin/python -m ou_irregular.week9.cli verify \
  --run-dir /tmp/ou-week9-first --require-completed
```

Repeating finalization with the identical review is read-only and idempotent.
A completed review cannot be replaced. Changed source, renderer or package
bytes block verification rather than silently updating the recorded hashes.

## Validation scope and later writing

Run the Week 9 tests and the affected saved-table/plotting regression tests in
the recorded environment, preserving collection, outcomes, JUnit and logs.
Tests use saved records and explicit forbidden-call guards; no new full study
is necessary for this presentation layer. Existing numerical regression and
historical CPU-sensitive replay evidence remain documented in the
[Week 8 completion report](week8_completion.md).

The Week 9 package freezes discrete, pointwise performance statements, not a
continuous operating boundary, robust interval method or general dependent-gap
limit theorem. The [benchmark note](week9_benchmarks.md),
[figure specification](week9_figure_plan.md),
[protocol](week9_protocol.md) and
[literature comparison](week9_literature_comparison.md) state the exact scope.
Final manuscript drafting, academic review and the pre-submission literature
refresh remain later stages. Ongoing GitHub CI is not required; the commands,
saved checks and evidence provide the reproducibility record.
