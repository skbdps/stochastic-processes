# Sequential execution plan: Weeks 5–8

Prepared 29 October 2026, before E1/E2 execution. Baseline: `4e35ba2`.
The original plan in Git history, the dated Week 3 amendment, and
`preregistration.md` govern this work. This implementation plan does not amend
the registered headline designs, random streams, inference or thresholds.

## Architecture and order

Preserve every existing Week 3 and Week 4 module and archived result. Add a
small `experiments` package for experiment contracts, execution, persistence,
summaries and comparisons. Separate week-specific entry points compose these
services. Inference receives observations only, never simulation truth.
Implement and check one module at a time, then complete one week before writing
the following week's implementation. Independent read-only reviews may overlap.

### Week 5: E1 irregularity sweep

1. Validate an immutable scientific design and exact registered configuration.
2. Add complete-record summaries and Monte Carlo uncertainty; verify formulas
   against hand calculations, including unavailable fits/intervals.
3. Add independent regular-reference and paired within-path comparisons.
4. Add deterministic execution with per-cell durable outputs, explicit
   provenance, failure retention and safe resume checks.
5. Add E1 figures, replay/artifact checks and mathematical documentation.
6. Freeze tested code, then execute all 14,000 paths / 42,000 fits once.
7. Audit complete planned keys, raw-to-summary arithmetic, source preservation
   and sample replay. Inspect figures; document actual results and limitations.
8. Publish Week 5 before beginning Week 6.

### Week 6: E2 coarseness, irregularity and sample size

Reuse the reviewed engine with the full registered 60-cell design: 120,000
paths and 360,000 fits. Add parameter/sample-size-specific heatmaps and a
failure log. Preserve all cells, including unresolved and null findings. Repeat
the complete-record, provenance, arithmetic and replay checks; document and
publish before Week 7.

### Week 7: E3 exploratory diagnostics

After examining complete E1/E2 results, write a dated selection rule and
diagnostic protocol before selected reruns. Investigate boundary, optimization,
curvature and minimum-gap-floor issues. Preserve headline outputs; separately
label any sensitivity runs and technical reruns. Produce the one-page findings
memo, documented numerical disposition and executable reproduction commands.
Professor contact remains deferred. Resolve or document anomalies before Week 8.

### Week 8: completion gate and optional clustered sampling

Commit a gate verdict. If core implementation defects remain unresolved, use
Week 8 as buffer and report the limitation. Otherwise specify an optional MMPP
protocol before its results: time versus event stationary initialization,
intensities/sojourns, gap-law calibration, independent OU source paths, controls,
replication count, seeds, floor rule and contrasts. A control matching only
mean/CV cannot isolate serial dependence; use an independent-gap control with
the same full event-stationary marginal law for that claim. Validate simulation
against analytical moments and limiting cases before experimental execution.

## Evidence and interpretation

Keep conditional point metrics and interval coverage distinct from operational
coverage over every planned path. Report Monte Carlo uncertainty and availability.
Different CV cells are independent; estimators within a source path are paired.
Regular references match nominal design, not each realized observation span.
Use absolute metrics for zero-truth mu. No numerical repair may change the
registered inference silently; defects require preserved inputs and an explicit
versioned disposition. The work is complete only when claimed results have
actually run and been audited, not when configuration files merely exist.

Historical Week 4's strict audit compares its entire recursive source map; run
that historical audit at its recorded revision. New audits verify their explicit
recorded dependency files, allowing unrelated later additions without falsely
changing old scientific provenance.

## Handoff checkpoint

Implementation paused at the owner's request on 29 October 2026 to continue in
a new chat. The E1 design and initial summary components are implemented and
unit tested. No E1/E2 headline runs have occurred. See `handoff_week5_start.md`
for the exact checkpoint and remaining work.
