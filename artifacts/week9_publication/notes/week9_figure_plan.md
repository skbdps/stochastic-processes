# Week 9 publication figures, tables and captions

Dated 31 October 2026. This post-result presentation specification follows the
[Week 9 protocol](week9_protocol.md). Publication figure IDs are a mapping onto
verified saved evidence; they do not rename or alter historical archive files.
Original figures and complete source tables remain available behind the final
publication selection.

## Figure families

| ID | Data and assembly | Required interpretation |
| --- | --- | --- |
| F1 | ECDFs of saved E1 estimates, with rows for regular, CV=1 and CV=3 and columns for theta/mu/sigma; supporting ECDFs expose all seven CV settings | Post-result editorial selection of design anchors; all three methods; truth marker; valid and planned counts; no density-smoothing inference |
| F2 | Existing E1 operational/conditional coverage curves and companion availability | Wilson Monte Carlo intervals; 95% nominal reference; separate 80% practical reference for applicable parameters; model-based working PFML/Euler inference |
| F3 | Existing E1 relative theta/sigma bias and separate absolute mu display | MC uncertainty in original definitions; no truth-relative mu values or classification |
| F4 | Existing E1 RMSE ratios against the regular design | Same estimator/parameter; independent-cell uncertainty; nominal n/mean-gap match, not a realized-span match |
| F5 | Complete E2 operational-coverage maps for each n | Every registered q/CV point; unavailable intervals counted as noncoverage; common color scale |
| F6 | E2 conditional coverage and availability companions for each n | Distinguish the available-interval denominator from all planned paths; keep low availability visible |
| F7 | Discrete theta/sigma practical classifications: three methods by two parameters for each of three n values | All 18 complete 4-by-5 grids; within/breach/unresolved; no smoothed boundary or single safe-CV threshold |
| F8 | E4 matched MMPP/IID mu RMSE ratios and coverage differences, with every method and r=1,5,20 | Independent-cell pointwise MC uncertainty; exploratory status; full one-gap marginal/floor match; realized-span consequences included |
| Supporting E3 | Original selected-case diagnostics, independent follow-up coverage/availability and paired floor sensitivity | Selected diagnostics are not prevalence estimates; retain boundary-suspected and Hessian limitations; 8,000 floor coupling groups |
| Supporting E4 | Complete original 18 figure pairs and source/spacing diagnostics; calibration comparison from saved summaries | Preserve all parameter/method contrasts, null/mixed results, tiny floor effects and mean-within-path versus marginal CV distinctions |

F1 includes every point-valid retained estimate within its displayed designs.
The full planned count and invalid count accompany each series even when no
estimate can enter its ECDF. Mu uses natural units; theta/sigma axes identify
whether values or truth-relative errors are shown. A marker at truth is not a
confidence interval. Repeated fits or regenerated paths cannot supply an F1
sample. Every plotted sample binds its study/cell/replication/method and original
partition digest.

F7 is the practical-classification display, not a nominal-calibration guarantee.
E2 PFML theta at n=250,q=0.1 has B,B,B,W,U across CV=0,0.5,1,2,3. Euler sigma
at that n/q has B,B,B,U,W. These patterns require discrete cells and explicit
unresolved encoding. Main and supporting layouts together show every panel.

F8's emphasis on mu is an editorial choice after the E4 results. The full
theta/sigma contrasts remain accessible. The 27 matched bias intervals include
zero; that does not prove equivalence. The isolated r=20 exact-sigma coverage
contrast is not a general beneficial-clustering result. A corrected-inference
F7b is outside this package because the existing studies do not evaluate such
a method.

## T1 and supporting data

T1 has a full 360-row E2 theta/sigma inventory, a compact 18-row
n/method/parameter count table, and all five unresolved E2 rows. The independent
42-row E1 inventory and its one unresolved classification remain separate.
Full parameter data retain mu with its inapplicable relative-bias criterion.

Each detailed table row identifies study, scientific cell, n, q, nominal CV,
nominal dimensionless span, method, parameter, source digest and classification.
Include bias with MC limits, operational/conditional coverage with counts and
Wilson limits, point/interval availability and the breach driver
(bias/coverage/both). Derive the driver from the unrounded registered bounds.
The class itself must agree with the saved audited class.

Expose the numerical data used by each final figure in companion machine-readable
tables. Keep paired estimator and floor comparisons distinct from independent
regular-reference and MMPP/IID comparisons. E1 and E2 retain their independent
experiment identities at overlapping coordinates.

## Caption and visual contract

Every caption states the study's registered/exploratory role, truth and design,
metric, eligible/planned denominator, uncertainty method and nominal level,
independent unit, and contrast direction. State that inference is pointwise.
Coverage differences use percentage points; RMSE ratios are dimensionless.
Label missing values explicitly. Preserve asymmetric log-delta ratio limits.
Zero observations of coverage do not establish zero population coverage.

Use stable method colors plus distinct markers, legible axes and shared scales
where comparison requires them. Evaluate at the intended manuscript width as
well as full resolution. A compact main display can refer to full supporting
pages without hiding methods, uncertain classifications or unfavorable results.
PDF font embedding and PNG/PDF text placement are part of actual visual review.

## Preservation and review

Reusing an original figure means copying exact bytes and recording their source
and prior review. Restyling or compositing creates a new derivative whose
current bytes require actual review. The publication manifest binds the
original data and new rendering provenance separately.

Review every new layout for correct series, uncertainty, scales, labels,
denominators, complete intended display and caption agreement with source data.
Record actual reviewer/time and file hashes, and identify AI review accurately.
Freeze only after independent table checks, deterministic external reproduction
and the actual visual review pass. This figure specification is not itself a
visual-review record.
