# Week 6: E2 full-grid progress and results

Prepared 30 October 2026 (India time). Week 5 was completed and published as
`d49fa4beff0800fc804afacbfdc7d7836311fdfe` before Week 6 implementation began.

## Code-freeze checkpoint

The strict configuration, complete-grid reporting, heatmaps and independent Week 6
entry point are implemented and reviewed. The complete suite passes **657 tests**
(`artifacts/week6_evidence/tests_pre_e2.txt`). No E2 headline path has been generated
at this checkpoint. All frozen Week 5 engine and inference files are unchanged.

E2 retains all 60 cells, 2000 paths per cell, seed 20261002, truth (1,0,0.5),
all sample sizes and the original floor rule. Reporting preserves 540 parameter
rows, 108 classification-count rows and all 180 cell/estimator failure-log rows,
including zero failures. It rejects missing complete groups and duplicate-plus-
missing substitutions, and distinguishes noncoverage from numerical failure.

The separate validation stream uses seed 21261002 and R=4 across all 60 cells:
240 paths, 720 fits and no pipeline errors. The original and fresh validation
archives have byte-identical raw partitions and all 23 derived figure-directory
files. All nine PNGs were visually inspected across three reviewers. These small
validation results establish implementation checks, not scientific conclusions.
The explicit 37-file dependency manifest includes the new Week 6 package without
altering prior sources or archives.

### Preserved development execution-state incident

The first validation attempt completed all raw cells, then failed the audit's
root-state gate. A traced resume observed `raw_complete_pending_audit` before table
export and `running_incomplete` afterward. The audit correctly refused completion.
The original exception logs and metadata snapshots are retained under
`artifacts/week6_evidence/`; no failed or inconvenient source paths were replaced.

Review found no root-metadata write in table export, no source/bytecode mismatch,
no unexpected file link, and no identified competing writer. Subsequent traced
table export, full workflow, fresh full validation and unchanged CLI resume passed.
The CLI resume generated zero new raw sources and replayed exactly the 120 fixed
endpoint paths required by the audit; all 60 raw partition hashes stayed unchanged.
The successful checks do **not** establish the cause of the earlier state change.
It remains an unresolved execution-state incident, not a verified scientific or
inference defect and not a claimed repair.

Disposition before headline E2: keep the original logs and snapshots; use the
unchanged fail-closed state, provenance, arithmetic, replay and visual-seal checks;
accept no incomplete or altered archive as complete. If a headline state failure
recurs, preserve it and investigate without replacing paths or weakening the gate.
This operational uncertainty must be revisited in Weeks 7–8. It does not authorize
changing the scientific protocol or treating passing reruns as evidence of a fix.

See [week6_mathematics.md](week6_mathematics.md) for the estimands and complete-grid
reporting contract. The code freeze is published before the full run, with its
exact published Git commit recorded in the run provenance.

## Headline E2 results

Pending. The full 120,000-path / 360,000-fit run, audit, figures and completion
seal must finish before results or Week 6 completion are claimed. Week 7
diagnostic selection and reruns have not begun.
