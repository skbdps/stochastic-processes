# Week 8 saved-record reporting

The reporting module consumes complete saved E4 partitions. It generates no
sampling schedules, OU observations, point fits or intervals. The existing
summary and paired-estimator calculations are composed with explicit E4
reference maps and a source-diagnostic summary.

## Complete planned outputs

| Artifact | Planned rows | Purpose |
|---|---:|---|
| `parameter_summary.csv` | 63 | Seven cells, three estimators, three parameters |
| `paired_comparisons.csv` | 63 | Three estimator pairs per cell and parameter |
| `regular_comparisons.csv` | 54 | Six stochastic cells against the independent regular reference |
| `matched_iid_comparisons.csv` | 27 | Three MMPP cells against their same-ratio IID controls |
| `spacing_summary.csv` | 7 | Existing realized-spacing summary |
| `source_summary.csv` | 7 | Raw, transformed and stored-time source diagnostics |
| `numerical_diagnostics.csv` | 21 | Three estimators per cell |
| `outcome_counts.csv` | Variable | Every observed eligibility/status category |
| `failed_fits.csv.gz` | Variable, including zero | All unavailable point/interval or pipeline outcomes, with the complete raw header |

Mechanism and ratio remain explicit in the reported scientific metadata.
The regular reference has its declared ratio-zero sentinel and its own cell
identity. References are selected by the fixed scientific maps, never by
replacing a cell's CV or matching realized spans.

## Eligibility, denominators and comparisons

Bias and RMSE condition on valid point estimates. Operational coverage uses
all planned paths, with unavailable intervals counted as noncoverage.
Conditional coverage uses available intervals. Optimizer success, point
validity, interval validity and pipeline success remain separate. Counts,
Monte Carlo uncertainty, Wilson limits and existing practical classifications
are retained. Relative mu metrics are undefined because its truth is zero.

Estimator comparisons within a cell are paired because methods share the
same source arrays. Their point contrasts use the intersection of eligible
fits, and their coverage contrasts retain paired discordance. Existing
paired formulas remain unchanged.

Different cells are independent. Regular-reference and matched-IID contrasts
therefore use each cell's own eligible point subset, with component counts
reported separately. Their bias and coverage difference standard errors use
the square root of the sum of component variances. Positive RMSE ratios use
the existing log-delta uncertainty calculation. Undefined component metrics
remain undefined; zero derivatives are not invented at singular cases.
Differences are the stochastic cell minus its reference, and ratios are that
cell divided by its reference. Equal replication labels never establish a
paired subset across cells. Comparisons sharing the regular sample are not
mutually independent.

The matched MMPP/IID comparison holds the entire one-gap marginal fixed,
including the floor transform. It permits dependence-induced differences in
realized span. The [prospective target calculation](week8_population_targets.md)
helps interpret the results; every simulation accuracy/coverage summary still
uses the true OU parameters.

## Source diagnostics

Source diagnostics must agree across all three method records at each
replication. The source summary uses one copy per source path. It retains
raw and posthoc-transformed gap moments, adjacent products, CDF measurements,
clipping and phase/switch diagnostics as declared by the module's metric
schema, together with actual stored-time spacing diagnostics and analytical
expectations. Across-path Monte Carlo errors use independent source paths as
the units. Neither repeated method records nor dependent gaps are counted as
additional independent paths.

Raw and transformed vectors describe the generated gaps before cumulative
roundoff. Realized spacing describes differences of the stored times actually
used for OU observations. The floor is not a renormalization or latent-chain
time evolution.

## Evidence binding and publication

Pure preflight verifies the fixed design, recorded numerical environment,
exact dependency set and durable cell identities. Table writes require
complete raw evidence with no pipeline error, use a run lock and the root
manifest guard, and preserve the raw lifecycle state for the later auditor.
Partial publication can be adopted only when independently regenerated table
bytes are identical to the saved bytes.

`table_manifest.json` records exactly the nine table files, their row counts
and SHA-256 hashes. `table_binding.json` separately binds the resolved design,
run identity, source policy and every complete cell manifest/partition. The
binding excludes root lifecycle bytes that legitimately change at later
audit and completion stages. Raw bindings, provenance and root bytes are
rechecked around publication; an altered input cannot silently become a new
table history.

Read-only table verification checks the saved binding, declared files, exact
schemas and complete planned keys, including the complete sorted failure
header when no failures occur. It does not claim independent numerical
recalculation. The next source-specific audit supplies that separate check.
The fitted R=4 experiment and full E4 study follow the complete pipeline freeze.
