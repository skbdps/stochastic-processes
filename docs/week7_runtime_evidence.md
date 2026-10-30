# Week 7 execution-state evidence and limits

Reviewed 30 October 2026 (India time and UTC), during reconstruction of the
selected-path diagnostics module. This document distinguishes observations
preserved in the published repository from historical observations whose local
files were lost. It does not identify a causal repair.

## Published Week 6 evidence

The completed Week 6 publication
`7ed521b4567da5c973c409e47e1b1d9906fa3b47` retains the original development
incident records in [week6_evidence](../artifacts/week6_evidence/).
The [recovery trace](../artifacts/week6_evidence/validation_recovery_trace.log)
records a raw-complete state before table export, a running-incomplete state
afterward, and the auditor's refusal to proceed. The trace establishes that an
unexpected state change was observed; it does not identify the writer.

The files named `validation_before_tables.json` and
`validation_after_tables.json` are byte-identical **successful later-check
snapshots**, both with SHA-256
`70972b0009539aa2a5752c913fa0d3d0de5db501a22801c29ee44948ceaf5ae0`.
They are not before/after byte evidence of the failing transition. The analogous
workflow snapshots have that same digest. The fresh-validation before/after
snapshots also agree, with SHA-256
`40998b59a15924b7e875dd576384a3df1918d8a04c1a3839ac32439855bc6d31`.
These comparisons qualify the scope of the preserved evidence; the archived
files and Week 6 report are left unchanged.

The [unchanged CLI resume trace](../artifacts/week6_evidence/unchanged_cli_resume_trace.log)
records zero newly generated raw sources, 120 fixed audit source replays and
60 unchanged partition hashes. The
[workflow trace](../artifacts/week6_evidence/workflow_write_trace.log) and
[fresh pipeline trace](../artifacts/week6_evidence/fresh_pipeline_write_trace.log)
record later successful state transitions. Success in those checks demonstrates
that verified resumption worked on those attempts. It does not explain the
earlier unexpected change.

## Lost unpublished Week 7 evidence

The interrupted session had also observed a root-state regression in a small
Week 7 development preview, then a successful isolated temporary-directory
preview. The original local traces, snapshots and validation artifacts were
lost when the runtime restarted. Their existence and reported outcomes remain
historical conversation evidence only; the reconstruction does not fabricate
equivalent files or count those checks as fresh verification. See the
[recovery record](week7_recovery_20261030.md).

A workspace synchronization script was observed in the runtime. That is a
possible operational explanation, not an attributed writer or established
cause. No platform service is modified or disabled for this project.

## Current safeguards and interpretation

Mutable recovery work is performed in an isolated temporary checkout. New
selected-path execution checks expected manifest bytes before updates and
preserves unexpected state as separate incident evidence instead of writing
cached state over it. Immutable source records, exact resume identities,
exclusive execution locks, fixed source addresses, dependency verification,
readback audits and explicit completion states remain necessary independently
of the execution location. An advisory lock constrains cooperating writers;
it cannot exclude an external process that ignores that lock.

The new safeguards are validated on synthetic fixtures before scientific
execution. Their observed test outcomes are recorded with the actual current
code. A passing test or later run means that the stated check passed on that
attempt; it does not establish that an unexplained historical runtime problem
has been repaired. The Week 8 core-completion decision must retain this
operational uncertainty and the loss of unpublished forensic evidence.
