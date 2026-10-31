# Week 8 saved-evidence figure contract

This module creates 18 separate single-axis PNG/PDF pairs from saved E4 tables.
It does not generate schedules or observations, reconstruct sources, fit models,
calculate new summary statistics from raw estimates, or replace prior results.
The builder's simple rescaling of stored Monte Carlo errors is documented below.

## Complete figure plan

Three coverage figures show theta, mu and sigma for all seven independent cells
and all three estimators. Operational coverage uses all planned paths, including
unavailable intervals as noncoverage; bars are pointwise Wilson Monte Carlo
intervals. Availability shows both point and interval availability, with the same
planned denominators and Wilson limits. Optimizer success is not substituted for
validity. Full integer counts and conditional coverage remain in saved tables.

Nine matched-control figures show bias differences, operational-coverage
differences and RMSE ratios separately for each parameter. Differences are MMPP
minus independently generated matched IID; ratios are MMPP over IID. Each cell's
own valid point subset supplies bias/RMSE; no matching of replication labels is
invented across independent cells. Bias uses parameter units, including mu
(which has zero truth and no relative metric). Ratios use a logarithmic axis and
the saved log-delta uncertainty. Difference intervals are the saved normal Monte
Carlo intervals. All intervals are pointwise, not simultaneous. NA marks an
undefined point; a dagger marks unavailable uncertainty without deleting a
finite point. A zero ratio is explicitly marked at the bottom of the log axis,
not misrepresented as missing or dropped. The matched law includes the entire one-gap marginal and floor
transform, but does not fix realized observation span.

Five source figures show mean gaps, adjacent-gap cross-products, average
within-path CV, realized-span standard deviation and observed clipping. Raw and
posthoc-transformed moments remain distinct. The analytical transformed-law
references on mean/product figures are saved targets, not new source draws.
Within-path CV is explicitly not equated to population marginal CV. Moment
whiskers are 1.96 times the saved across-path MCSE, with paths as independent
units. The span SD is the stored span MCSE times sqrt(R), using ddof=1; no interval
for SD is invented. Clipping is descriptive, with integer counts; dependent gaps
are not treated as independent Bernoulli trials. The clipping display is
nonnegative and has an upper limit of at least one gap per million, so an
all-zero sample does not acquire an artificial negative plotting range. Only r=20 is floored. The plots
retain the two independent exponential r=1 cells and the distinct regular
reference; they are not mislabeled as three regular controls.

## Evidence and renderer identity

Pure readback verifies the full saved audit against the complete raw/table
artifact map, scientific metadata, counts, all source-replay declarations and
fixed endpoint identities. Raw records are read and validated, not regenerated.
The saved numerical audit is not replaced by this binding check, and the plotter
does not claim to rerun its numerical reconstruction. Completed archives must
retain the root, artifact and visual seals. Partial figures are allowed only in
open precompletion stages and only with their original renderer state.

Rendering has separate provenance with an explicit dependency set, effective
Matplotlib settings, font-file hashes, compression/font-library versions and
fixed PDF metadata. No font files are distributed. The immutable renderer state
is a sibling `<run>-figure-renderer/renderer_state.json`; it is also embedded in
the figure manifest, so closed readback does not depend on the external sibling.
Full exploratory rendering requires committed renderer dependencies. Validation
may record dirty development files honestly. Raw provenance is never rewritten
to add a later renderer. Default colors are retained, with no custom palette.

A root-manifest guard spans rendering. Every candidate pair is generated and all
existing partial bytes compared before any missing figure is published. Writes
are immutable; a changed table, audit, root, renderer or partial figure fails
rather than overwriting history. The figure manifest is last and binds all36
image/PDF files. Repeated complete calls are read-only. An actual visual review
of every PNG is still required by the later runner; no rendering success claims
human visual inspection or full Week8 completion.

## External regeneration

`write_figures(source_run, design, destination=fresh_output)` uses the same
bindings, renderer identity, immutable publication and readback checks while
writing only the separate destination and its sibling renderer state. The source
archive is not copied, refitted, unlocked or rewritten. Destination ancestors
and descendants of the source are rejected, as are unrelated files at the
destination. This supports later `replot` commands on sealed runs. An external
replot never replaces or inherits a claim of actual visual review.
