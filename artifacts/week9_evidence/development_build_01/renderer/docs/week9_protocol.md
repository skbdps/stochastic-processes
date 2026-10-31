# Week 9: publication evidence and figure-freeze protocol

Dated 31 October 2026, after completion of E1, E2, the exploratory E3 studies,
and E4. Week 9 is a saved-evidence reporting stage. Its figure selection,
grouping, captions and cross-study synthesis are post-result presentation
decisions, not a new preregistration or independent confirmatory experiment.

## Authority and scope

The [preregistration](../preregistration.md), its adopted
[Week 3 amendment](plan_amendment_week3_validation.md), and the dated
[E3](week7_exploratory_protocol.md) and [E4](week8_exploratory_protocol.md)
protocols govern the scientific definitions and reporting rules. The original
master plan supplies the Week 9 figure-freeze and threshold-table objectives.
Later outcome-neutral interpretation replaces its failure-presuming language.
The research proposal also permits population-target analysis and empirical
calibration when a fully matched covariance derivation is beyond scope.

The input release is `b6edb85d1e7a21ed19fb54b5c71213ed04703309`.
Week 9 adds a verified source catalog, numerical claims, publication tables,
figure assemblies, captions, a benchmark note and a scoped literature review.
It preserves the completed experiment archives and the scientific methods that
produced them. The derivative renderer has its own source/environment identity.
Its completion report must distinguish implemented products, actual checks and
limitations; this protocol alone does not establish their completion.

## Source roles and independent units

| Source | Role in the publication | Accounting and independent unit |
| --- | --- | --- |
| E1 | Registered variability curves | Seven cells, 14,000 independent source paths, 42,000 fits |
| E2 | Registered joint-grid maps and primary practical-classification table | 60 cells, 120,000 independent source paths, 360,000 fits |
| E3 selected diagnostics | Explanation of selected E1/E2 outcomes | 53 selected existing source paths, 159 original-fit replays; not a prevalence sample or new headline data |
| E3 independent follow-up | Exploratory check of the fixed selected setting | 20,000 independent source paths, 60,000 fits |
| E3 floor sensitivity | Exploratory paired sensitivity analysis | 24,000 realizations, 72,000 fits in 8,000 independent coupling groups |
| E4 | Exploratory matched-marginal dependent-gap comparison | Seven independent cells, 14,000 paths, 42,000 fits |

Archive identity is part of every cross-study key. E1 and E2 can share a
scientific cell ID while having independent root seeds. Their outcomes are not
pooled or selected between to obtain a preferred classification. Estimators
within a source are paired; distinct E1/E2/E4 cells are independent. E3 floor
variants use the declared coupling. Contrasts sharing a regular reference are
not mutually independent.

## Fixed reporting rules

Point bias/RMSE use the point-valid subset and retain its count. Operational
coverage uses every planned path, treating unavailable intervals as noncoverage.
Conditional coverage uses interval-valid paths and is undefined when none are
available. Point-invalid outcomes are already included among unavailable
intervals; the two counts must not be added as disjoint failures.

Keep all original numerical candidates, statuses, thresholds, gates and negative
natural-scale Wald lower bounds. Report zero observed coverage with its Wilson
bound; zero plug-in MCSE is not certainty about the underlying probability.
Empirical SD and reported SE comparisons use the same interval-valid subset.
The RMS reported SE is distinct from the mean reported SE.

For theta and sigma, a supported practical breach is established by operational
coverage's Wilson upper limit below 0.80, or a relative-bias MC interval wholly
beyond +0.10 or -0.10. Within criteria requires the Wilson lower limit at least
0.80 and the full bias interval inside [-0.10,0.10]. Other cases remain
unresolved. The existing unrounded classification is authoritative; an
independent check of its logic must agree before export. Mu has zero truth, so
relative-to-truth metrics and this relative-bias classification remain
inapplicable. Within criteria does not mean precisely 95% coverage.

The main E2 T1 inventory contains 360 theta/sigma rows; E1 has a separate 42-row
inventory. Retain all six unresolved E1/E2 parameter classifications. F7 uses
the complete discrete grid, with separate methods, parameters and sample sizes.
Nonmonotone classifications prohibit a single maximum-safe-CV interpretation.
No interpolation supplies evidence at untested coordinates.

All research uncertainty statements are pointwise. Main-text selection does not
convert exploratory E3/E4 comparisons into confirmatory results or provide
simultaneous error control. Publish complete supporting tables and figure
inventories alongside the selected displays.

## Claims and theoretical scope

Each quantitative publication claim must bind its exact source row keys,
columns, file digests, denominators, units, uncertainty method and displayed
rounding rule. Distinguish supported, descriptive, unresolved, theoretical under
stated assumptions, and unsupported statements. Generate numerical captions
from these bound values rather than editing numbers independently.

The [benchmark note](week9_benchmarks.md) connects population targets to empirical
bias and interval calibration. A full nuisance-matched covariance benchmark,
dependent-gap central limit theorem, and robust PFML/Euler interval procedure
are not established by the present experiments. The manuscript may describe
measured finite-sample performance; it must not claim a demonstrated departure
from first-order coverage theory on the strength of the current comparisons.

The [literature comparison](week9_literature_comparison.md) records the checked
sources, version dates and access limits. It supports a bounded contribution
statement, not a certificate that no overlapping work exists. The final
pre-submission search remains a later task.

## Build, review and freeze

The build reads the pinned saved evidence and writes a fresh derivative output
directory. Before and after generation, verify the consumed inputs against the
source catalog. Reject a destination within or overlapping the input archives,
unsafe path traversal, unexpected symlinks and conflicting existing output.
Verify exact planned keys in addition to aggregate counts.

The publication package binds its tables, plotting data, captions, claim values,
figures, source catalog, renderer sources and environment. Its read-only verifier
must detect missing, unexpected or changed evidence. Rebuild to a separate
directory in the recorded renderer environment and compare deterministic output
bytes. Saved-data regeneration is distinct from numerical source reconstruction
and fit replay.

Newly rendered or composited layouts require actual image review at manuscript
viewing size. Record the actual reviewer, time, inspected files, relevant hashes,
and legibility/consistency findings. A prior review transfers only to explicitly
identified identical bytes; it is not a new inspection of changed layouts.
Finalize the derivative package only after the numerical checks and actual
review pass. A source mismatch, failed accounting check or independently
demonstrated defect blocks that completion claim until investigated.

The [figure specification](week9_figure_plan.md) gives the publication mapping.
Week 9 completion supplies a figure/table order, accepted claims, uncertainty and
limitations to manuscript preparation. It does not certify personal theory
mastery, academic review, journal acceptance or completion of the manuscript.
