# Week 7 findings memo

**30 October 2026 | Exploratory E3 | Scientific freeze 6c97d939**

## Completion

The fixed 53-path diagnostic selection reproduced all 159 original estimator
records. The independent follow-up completed 20,000 paths / 60,000 fits; the
floor sensitivity completed 24,000 physical realizations / 72,000 fits,
representing 8,000 independent coupling groups. Both new studies passed full
audits, actual figure review and completion seals, with zero pipeline errors.
E1/E2 files and decisions remain unchanged. [1]

## Main findings

**Boundaries and curvature.** All 85 original point-invalid records are
consistent with the saved conditional-regression boundary evidence. All 14
evaluated exact profiles match their stored fits within the fixed diagnostic
resolution. The rejected irregular exact interval reproduces its whitened
Hessian discrepancy in longdouble: 0.095654 versus the unchanged 0.01 gate.
Approximately fourfold decay on step halving supports finite-difference
truncation amplified by weak curvature. No original interval is replaced;
finite searches do not prove global optimality. [2]

**Independent follow-up, n=250 and q=2.** The study retains 421 invalid points
and 23 further unavailable intervals: 444 records on 231 source paths. One
irregular exact candidate remains `boundary_suspected`; no finite optimum
is certified. Operational theta coverage is 92.63% for regular exact/PFML,
94.77% for irregular exact, 72.29% for irregular PFML and 0% observed for
Euler. Regular exact/PFML conditional coverage is 94.64%; excluding unavailable
intervals would hide that denominator difference. Irregular exact theta's
Wilson interval is 94.32%-95.19%; irregular PFML's is 71.40%-73.16%. Euler's
zero count has a 0.0384% Wilson upper limit. New results do not overwrite
headline classifications. [3]

**Floor sensitivity, n=1000 and q=0.5.** All 72,000 points and intervals are
valid. Across the fixed CV/epsilon contrasts, the largest estimated RMSE-ratio
change is exact theta at CV=3 and epsilon=1e-4: 1.00023261, with paired 95%
Monte Carlo interval [0.99974475, 1.00072070]. The largest absolute observed
coverage difference is 0.05 percentage points. Small nonzero bias/RMSE effects
are retained. These results support modest sensitivity on this grid, not
universal invariance or formal equivalence. CV=0 is one shared exact control;
other variants share innovations, not one continuous trajectory. [4]

## Assurance and limits

Independent saved-data review reproduced all 126 parameter-summary rows and
72 epsilon contrasts; a separate review passed 4,396 checks on selected
diagnostics, including algebra on 352 saved Hessians. All eight PNGs were
visually inspected. Current archive hashes verify, but the earlier
unattributed manifest-state incidents remain a historical operational
limitation. Nonrecurrence does not establish their cause. Week 8 requires a
separate committed gate after Week 7 publication. No new core estimator or
inference defect was established by these checks. [1,2]

### Evidence

1. [Completion report and verification](week7_progress.md).
2. [Selected numerical disposition and exact review scope](week7_numerical_disposition.md).
3. [Follow-up summaries](../artifacts/week7_followup/parameter_summary.csv) and [retained outcomes](../artifacts/week7_followup/outcome_counts.csv).
4. [Paired floor contrasts](../artifacts/week7_floor/epsilon_comparisons.csv) and [source audit](../artifacts/week7_floor/audit.json).
