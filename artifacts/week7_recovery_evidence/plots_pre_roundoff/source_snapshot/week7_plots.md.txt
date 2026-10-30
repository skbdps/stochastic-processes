# Week 7 figures from saved evidence

This module renders the dated exploratory protocol's follow-up, floor
sensitivity and two selected-case views. It changes no source addressing,
fit, Hessian gate, interval, classification or archived headline result.
It does not generate observations, refit any path, recompute a profile or
Hessian, or invoke a numerical audit. A figure is a presentation of the
saved evidence and its limitations, not an additional scientific experiment.

The renderer is `ou_irregular/week7/plots.py`. The numerical input and its
audit must already exist. The executable Week 7 runner is a separate module;
these APIs never advance a run's lifecycle or certify visual inspection.

## APIs and eight planned figures

`write_figures(kind, output_dir, design, *, extra_files=())` writes into
`output_dir/figures`, where `kind` is `followup` or `floor`. The exact full
design, or its exact separate-seed `.validation(R)` design, is required.
`verify_figures` has the same arguments and performs a read-only complete
readback.

`write_selected_figures(selected_dir, output_dir, *, extra_files=())` writes
into an explicit external output directory. It reads closed selected-case
diagnostics; both its output and its sibling renderer state must be disjoint
from that input archive. `verify_selected_figures` has the same arguments and
is read-only.

Each following stem produces a PNG and a PDF, with no additional image names:

| Stem | Saved evidence shown |
| --- | --- |
| `W7_followup_coverage` | Operational and conditional coverage for theta, mu and sigma, both CV cells and all three estimators |
| `W7_followup_availability` | Joint point/interval availability and two disjoint failure counts |
| `W7_floor_paired_bias` | Both nonreference floors versus 1e-6, all four CVs, three parameters and estimators |
| `W7_floor_paired_coverage` | Paired operational coverage differences over all planned sources |
| `W7_floor_rmse_ratio` | Common-valid paired RMSE ratios and their saved asymmetric log-delta limits |
| `W7_floor_spacing` | Per-source mean gap, clipping fraction, minimum gap and span summaries at all three floors |
| `W7_selected_E2_cv0p5_rep000465` | Saved profile and Hessian diagnostics for the fixed irregular E2 case |
| `W7_selected_E2_cv0_rep001424` | Saved profile and Hessian diagnostics for the fixed regular E2 case |

The pure builders `build_followup_figures(tables, design)`,
`build_floor_figures(tables, epsilon, design)` and
`build_selected_figure(payload, *, validation)` return Matplotlib figures.
Table dictionaries are keyed by the canonical CSV filenames. These builders
make the plotted values inspectable without a durable output or a source
generator. `profile_series(payload)` returns the exact display arrays and
a deterministically derived zoom interval; `figure_bytes(figure, stem)` is the deterministic
PNG/PDF serialization boundary.

Pure builder fixtures may be literal numerical examples. A selected payload
with `literal_validation_fixture=True` requires `validation=True` and
visibly states that it is a **literal validation fixture, not an original
E2 replay**. The full archive adapter separately requires complete closed
payloads, original-record bindings and their saved readback audit. A successful
literal plot test therefore does not certify an original selected path.

## Statistical meaning of the panels

Follow-up coverage includes all 18 parameter-summary rows. The upper row
uses all R planned sources, including unavailable intervals as noncoverage.
The lower row uses each method's interval-valid denominator. Vertical bars
are the saved pointwise 95% Wilson Monte Carlo intervals; the dotted line
is nominal 95% coverage. Missing conditional coverage is a visible NA gap.
PFML and Euler intervals retain their model-based working covariance
interpretation.

Point and interval validity are joint across the three parameters. The
availability plot takes one checked parameter copy of those joint counts.
Its two failure categories are R minus point-valid, and point-valid minus
interval-valid. They are disjoint, sum to all unavailable intervals, and
are never counted three times. Zero bars are explicitly labelled zero;
count axes use integer ticks. Availability bars use saved Wilson limits.

All 72 paired-floor rows and all 71 fields remain in the input contract.
CV axes identify their equally spaced design cells explicitly; visual spacing
does not imply equal numerical distance between the nominal CV values.
The displayed difference is the metric under epsilon minus the metric under
the registered 1e-6 floor; the RMSE ratio is RMSE(epsilon)/RMSE(1e-6). Positive coverage difference
means higher coverage under the displayed epsilon; it does not imply that
epsilon itself is numerically larger. The reference floor has no self
contrast.

The bias and ratio panels use the common-valid point subset for each
estimator, CV and floor pair, and print those S counts in estimator order.
Both underlying metrics were recomputed on the same subset by the reviewed
comparison module. Theta and sigma bias differences are truth-relative
percentage points; mu uses absolute state units. The practical absolute
relative-bias threshold is not a threshold for a bias contrast.

Operational coverage differences use every one of the R planned paired
sources, not S. Their saved uncertainty is the paired normal approximation
and may extend beyond plus or minus 100 percentage points; the plot does
not clip it. A collapsed normal-MC interval reflects zero observed
variation, which by itself does not prove population equality.

RMSE ratios use a logarithmic ordinate and the saved asymmetric log-delta
endpoints. The renderer does not replace them with a symmetric interval.
An exactly zero defined numerator ratio cannot be drawn on a log axis:
it is explicitly labelled 0 at the panel baseline and omitted from the
positive log curve. Undefined ratios are labelled NA instead. An available
estimate with unavailable MC uncertainty is drawn without a bar. Compact
tick formatting preserves legibility when all ratios lie very close to one.

The source-spacing panels use R source summaries per physical variant, not
3R estimator rows. Each plotted field must report R available source values.
The minimum-gap panel shows the mean of each path's minimum, not a global
minimum. Its signed-log scale preserves small-sample normal-MC limits that
cross zero. The table retains realized CV although it is not a fifth panel.
There are n observations and K=n-1 transitions: n=1000 means K=999. Floors
are not followed by renormalization or realized-span matching. All three
variants within a CV/replication share the base Gamma, initial-state and
transition-innovation randomness. Their curves are paired; the CV=0
control has exactly identical physical sources.

## Selected profile and Hessian views

The two selected identities are fixed by the protocol, with n=250 and q=2.
Each figure has the complete saved profile grid, a local view around the
stored theta and best finite candidate, exact-Hessian eigenvalues, condition
numbers, standardized scores and three adjacent-step discrepancies.

For every profile grid value the renderer subtracts the saved
`best_finite.nll_extended` from `grid_nll_extended` in `numpy.longdouble`
**before** converting the completed difference to float64 for Matplotlib.
Subtracting the rounded float64 absolute NLL values can erase the small
differences that motivate these diagnostics. Saved extended differences
also supply the stored-profile, stored-full and limiting-likelihood markers.
Negative differences are retained. No new likelihood values are evaluated or
interpolated; plotted line segments connect saved grid values. No optimizer
or likelihood is called. The zoom includes the
stored and best locations with 0.15 units of log-theta padding, clipped
to the saved grid.

The figure displays all four recorded relative steps, 0.001, 0.0005,
0.00025 and 0.000125, and checks the identities of all three adjacent
pairs. The registered procedure uses its original h/h2 comparison to
determine interval validity. Finer diagnostic steps do not replace that
decision. Signed eigenvalues and zero diagnostics remain visible; missing
inverse-based quantities on non-PD steps remain NA. The original accepted
or unavailable interval status is printed.

Float64 and independent longdouble exact-Hessian ladders use the same
recorded float64-normalized input representation. The saved longdouble
profile separately normalizes original observations and gaps using the
recorded observed-data scales. Their representation differences are
printed rather than silently treated as mathematical disagreement.
Precision labels report the saved significant binary digits, fraction bits
and storage bytes. Sixteen storage bytes do not imply IEEE quadruple
precision or arbitrary precision.

The saved resolution, tail-target status and finite-search disposition are
printed. Finite search does not prove a global optimum. Neither a better
diagnostic objective value nor a different diagnostic Hessian supplies a
replacement fit, interval or covariance. Closed `diagnostic_review_required`
evidence is allowed if its retained failures and saved audit agree. An
unavailable fixed case produces an explicitly unavailable panel; an
evaluated fixed case may still be rendered when another selected source
failed. The figure manifest preserves the archive's review-required scope.

## Saved input identity and lifecycle

Run rendering requires complete zero-pipeline-error raw accounting and a
saved **passed full numerical audit** with the expected source contract,
fixed endpoint replay counts, run identity, design, semantic metadata and
actual raw/table byte map. Floor inputs must disclose the paired-source
override and all-source array verification. The renderer does not invoke
that audit again: a legitimate interrupted figure stage can contain
partial PNG/PDF files that a full numerical archive audit rejects until
the figure manifest is complete.

The standard table manifest must declare exactly its six CSV summaries and
the complete failed-fit export. Declared files cannot be missing or changed.
Floor rendering separately validates the 71-column/72-row epsilon table
and its explicit raw-input manifest. CSVs use pandas' round-trip float
parser to preserve the canonical scientific floating-point keys.
Pure builders additionally reject incomplete key sets, nonintegral or
Boolean counts, out-of-range denominators and inconsistent own/common/
coverage contingency accounting.
Before renderer-state creation, the archive adapter checks the complete keys
of all six summary tables, canonical cell coordinates, per-source spacing
availability, and complete outcome/numerical count accounting. Fixed saved
audit endpoint descriptors are checked against saved raw endpoint triplets,
their source hashes, estimator membership, selection and tolerance declaration;
this is a byte/identity readback, not a repeated numerical replay.

The normal writable state is `raw_complete_pending_audit` with a saved full
audit. A root audit seal has not necessarily been added at this stage.
`audited_pending_visual_review` and a complete, zero-pipeline-error
`implementation_review_required` postprocessing retry are also supported.
If the root declares an audit hash, it must match; the audited/completed
states require it. `completed` is verification-only: missing or changed
figures cannot be recreated inside a sealed archive.

Persistent run bindings use the semantic metadata digest, run identity,
design and saved audit/artifact hashes. Legitimate later audited/completed
transitions therefore do not invalidate a figure. Exact root bytes are
instead guarded during each rendering call. For selected diagnostics the
input is already closed, so its exact root bytes and every source payload
are also part of the persistent input binding.

## Separate immutable renderer provenance

Earlier saved raw archives may predate the plotting module. Their dependency
records are verified as recorded and never upgraded or rebound. The renderer
captures a separate explicit provenance record covering `plots.FILES`,
the historical core dependency set and exactly the caller's `extra_files`.
Its identity includes its actual input descriptor and a complete render
environment. A full exploratory rendering requires all those bytes to be
committed before any output is created. Verification checks the actual
recorded Git commit's dependency bytes; it does not demand equality to
the current HEAD, which may advance through unrelated artifact commits.

The full effective rc configuration begins with the installed Matplotlib
defaults and the module's specified style. Ambient notebook changes to
TeX, transparency, font family or DPI cannot affect output. Renderer
provenance also binds Matplotlib, FreeType, Pillow, FontTools, zlib and the
bundled TrueType font bytes. PNG metadata is static; PDF creation and
modification dates are suppressed. No capture timestamp enters figure
bytes.

Before writing any images, the writer creates one immutable sibling:

```text
<run-name>-figure-renderer/renderer_state.json
```

For external selected figures the same naming rule uses the external output
directory name. The sibling directory contains exactly that one descriptor.
It is outside the strict scientific archive layout, guarded against symbolic
links, surplus files and nesting in another saved experiment. It records the
first renderer provenance, descriptor digest and canonical expected names.
Resume reuses that exact record, including its original capture time.
Preexisting partial figures with no original renderer state are rejected.
The sibling is required resumability evidence and should be preserved and
published alongside its run/output.

The final `figure_manifest.json` embeds the complete renderer state, its
digest, input bindings and the exact PNG/PDF hash map. Complete verification
is self-contained and deliberately does not require or inspect an optional
external sibling. The historical `input_table_sha256` compatibility field
contains only the seven standard table artifacts. Epsilon CSV/manifest
hashes are kept in a separate explicit field and verified by this module.
Selected manifests instead bind the closed diagnostic descriptor, serialized
selection, all payload hashes, both selected source identities and audit.

## Publication, interruption and readback

Read-only preflight rejects unknown filenames, symbolic links, invalid input
identities, changed declared products and nonfrozen full-scope dependencies
before creating output or lock files. An exclusive run/output lock guards
mutation. Input raw/table/audit bytes and exact root bytes are rechecked
before and after rendering stages. A foreign root change is preserved by
the reviewed manifest guard as adjacent expected/observed incident evidence
and blocks automatic mutation; no stale root write or finally-update exists.

Images are committed immutably as complete files. A later retry may adopt
an existing partial file only if re-rendering produces exactly the same
bytes. A declared complete manifest must have every expected file and no
surplus files; declared missing products are not regenerated. State and
input drift stop publication while preserving partial evidence. The figure
manifest is the final rendering product, after all expected PNG/PDF pairs
are present and rechecked.

The module writes `visual_review_performed=False`. Root lifecycle bytes,
audits and raw records remain untouched. An actual human/agent visual review
of every exported PNG, the separately guarded runner transition and the
completion seal remain distinct required steps.
