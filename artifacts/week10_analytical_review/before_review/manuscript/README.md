# Manuscript drafts

The current manuscript is [Week 10 Draft 2](week10_draft2.md), dated
1 November 2026. It adds a PFML mean covariance mechanism, span controls and
focused supplementary results to the original methods baseline.
It is an editable Markdown manuscript with GitHub-rendered mathematics.

The original [Draft 1](week10_draft1.md) describes the completed studies at
release `4bddb527f3c8bb88bc1ebec28c8818e03dc251aa`. Draft 2 additionally uses
the separately audited [Week 10 supplement](../docs/week10_supplement_completion.md).
It introduces a new experimental design and analytical benchmark while
preserving the original estimators, intervals, archives and Week 9 package.

Read the [source map](week10_source_map.md) when checking a method statement,
and the [draft status and validation record](../docs/week10_draft1_status.md)
for the historical baseline review scope; the supplement has its own linked
validation and review records. Draft 2 has not received human academic review and is
not a submission-ready manuscript. No final PDF or journal layout is claimed.

## Author review of Week 10

Work through the draft in this order. Record corrections and questions in a
subsequent revision rather than changing frozen study records.

| Passage | Question for the author |
| --- | --- |
| Introduction | Is this the research question and contribution you intend to defend? Does the financial motivation fit the scope without implying a real-data study? |
| Model and observation conventions | Can you derive the exact OU transition and explain stationary simulation versus a likelihood conditional on the initial observation? |
| Exact, PFML and Euler | Can you explain which information or moments each procedure replaces, and why PFML mean estimation can respond to sampling dependence despite its regression formula? |
| Numerical coordinates and intervals | Can you explain the full Hessian, covariance transformation, natural-scale intervals and separate point/interval validity rules? |
| Population targets | Can you distinguish an expected-objective minimizer from finite-sample mean, variance and coverage? Can you explain why full-marginal matching does not imply equal precision? |
| Mean covariance mechanism | Can you derive the PFML endpoint correction, distinguish fixed-coefficient variance from fitted-estimator coverage, and explain why matching lag one leaves higher-lag covariance different? |
| Span supplement | Can you explain common-support exclusions, conditional MC uncertainty, event-Palm initialization, random N and the deterministic terminal observation? |
| Experimental methods | Can you explain the independent units, controls, floors, coverage denominators and uncertainty-aware classifications without treating all fit records as independent paths? |

Author review includes checking derivations, revising prose and verifying
calculations. The author remains responsible for understanding the methods
and accepting the claims before independent academic assessment.

## Subsequent manuscript work

Week 11 integrates the core Results with the new focused supplement, practical
guidance, related work and limitations, using the existing [claims ledger](../artifacts/week9_publication/claims.md),
[T1](../artifacts/week9_publication/tables/T1.md),
[figure plan](../docs/week9_figure_plan.md) and
[literature comparison](../docs/week9_literature_comparison.md).
Every reported result must retain its source, denominator and uncertainty.
An introductory claim about the study's aims is not evidence for an empirical
conclusion.

Week 12 addresses academic feedback, refreshes the literature search, writes
the final abstract and checks the assembled manuscript, figures and
reproducibility package. Actual academic review, outreach and submission are
not completed by publishing this initial draft in the repository.
