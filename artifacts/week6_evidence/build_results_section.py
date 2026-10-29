"""Generate the numerical tables used in the Week 6 narrative from audited CSVs."""
from pathlib import Path
import json
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
run = ROOT / 'artifacts/week6_e2'
audit = json.loads((run / 'audit.json').read_text())
assert audit['status'] == 'passed' and audit['full_completion_audit']
s = pd.read_csv(run / 'parameter_summary.csv')

def markdown(headers, rows):
    return '\n'.join(['| ' + ' | '.join(headers) + ' |', '| ' + ' | '.join(['---'] * len(headers)) + ' |'] +
                     ['| ' + ' | '.join(map(str, row)) + ' |' for row in rows])

rows = []
for n in [250, 1000, 4000]:
    for method in ['exact', 'pfml', 'euler']:
        for parameter in ['theta', 'sigma']:
            group = s[(s.n == n) & (s.estimator == method) & (s.parameter == parameter)]
            counts = group.practical_classification.value_counts()
            rows.append([n, method, parameter, *[int(counts.get(k, 0)) for k in
                         ['within_criteria', 'supported_breach', 'unresolved']]])
classes = markdown(['n', 'Method', 'Parameter', 'Within', 'Supported breach', 'Unresolved'], rows)
group = s[(s.n == 250) & (s.theta_mean_gap == .1) & (s.estimator == 'exact') & (s.parameter == 'theta')]
short_span = markdown(['CV', 'Relative bias', '95% MC interval for bias', 'Operational coverage', '95% Wilson interval'],
    [[f'{r.cv:g}', f'{100*r.relative_bias:.3f}%', f'[{100*r.relative_bias_mc_lower:.3f}, {100*r.relative_bias_mc_upper:.3f}]%',
      f'{100*r.coverage_operational:.2f}%', f'[{100*r.operational_wilson_lower:.3f}, {100*r.operational_wilson_upper:.3f}]%']
     for r in group.itertuples()])
group = s[s.practical_classification == 'unresolved']
unresolved = markdown(['n', 'q', 'CV', 'Method', 'Parameter', 'Relative bias', '95% bias MC interval', 'Coverage', '95% coverage Wilson interval'],
    [[r.n, f'{r.theta_mean_gap:g}', f'{r.cv:g}', r.estimator, r.parameter, f'{100*r.relative_bias:.3f}%',
      f'[{100*r.relative_bias_mc_lower:.3f}, {100*r.relative_bias_mc_upper:.3f}]%', f'{100*r.coverage_operational:.2f}%',
      f'[{100*r.operational_wilson_lower:.3f}, {100*r.operational_wilson_upper:.3f}]%'] for r in group.itertuples()])
section = '''## Headline E2 results

The complete registered E2 run used frozen published commit
`960436615a710d21de3439d50257f50af72a45e4`. It began at
2026-10-29 21:49:11 UTC, finished raw execution at 21:56:24 UTC, and finished
the full audit and figure generation at 22:04:37 UTC. Raw execution took
433.71 seconds with 12 workers and 500-path durable partitions. These times
exclude publication and independent review. There were no technical reruns,
resampled source paths, changed numerical gates or changed scientific settings.

All **60 cells, 120,000 paths and 360,000 planned fit records** are retained in
240 partitions. The full audit passed: 359,915 point-valid records, 359,914
interval-valid records, 86 retained records with unavailable intervals, and
zero pipeline exceptions. It independently checked every summary and contrast
table, the full key set and saved Hessian/covariance identities, and replayed
replications 0 and 1999 in every cell (120 source paths, 360 fitted results).
The replay covers a fixed subset; it is not a claim that all source arrays or
all coarse Hessians were independently regenerated. See the scope in
`artifacts/week6_e2/audit.json`.

The core archive contains 540 parameter summaries, 540 independent regular
contrasts, 540 paired estimator contrasts, 60 spacing rows, 180 numerical
diagnostic rows and 184 observed status-composition rows. The Week 6 reports
retain all 540 classification rows, all 108 classification-count combinations,
and all 180 cell/method failure rows including zero-failure cells. All nine
PNG figures and their PDF counterparts show the full grid without truncation.

### Pointwise practical classifications

Each n/method/parameter group below contains 20 cells. Labels combine the
registered operational-coverage Wilson bounds and relative-bias Monte Carlo
interval, without rounding before classification. They are pointwise grid
statements, not simultaneous guarantees or interpolated boundaries. Mu's zero
truth makes the relative-bias criterion inapplicable.

''' + classes + '''

Exact theta is within criteria in 55 of 60 cells; its five supported breaches
are all at n=250,q=0.1. Exact sigma is within criteria in every cell. PFML has
40 theta and 41 sigma supported breaches; Euler has 51 theta and 58 sigma
supported breaches. Five parameter-specific E2 cells remain unresolved.
These are separate parameter statements, not a method-wide certification.

### Short-span exact-theta bias, without the predicted coverage collapse

At n=250,q=0.1 the nominal dimensionless span is 24.9. Exact theta has a
supported positive relative-bias breach at every CV, including regular spacing:

''' + short_span + '''

The earlier master plan anticipated an exact-coverage collapse in this corner.
That outcome did **not** occur in E2: exact-theta operational coverage is
93.6–95.3% here, while relative bias is 17.1–22.6%. The supported breaches are
driven by bias, not coverage below 80%. The regular reference already has
substantial bias, so the entire effect cannot be assigned to irregularity.
Across all E2 cells, exact theta coverage ranges from 91.6% to 96.4%, and exact
sigma coverage from 91.75% to 96.35%. These observations do not prove a general
coverage theorem or exact nominal calibration at every design.

### Numerical availability is separate from noncoverage

All unavailable intervals occur at n=250,q=2:

| CV | Method | Point-invalid | Point-valid without interval | Total unavailable intervals |
| --- | --- | --- | --- | --- |
| 0 | exact | 39 | 0 | 39 |
| 0 | PFML | 39 | 0 | 39 |
| 0.5 | exact | 0 | 1 | 1 |
| 0.5 | PFML | 7 | 0 | 7 |

The point-invalid records are labelled `boundary`; the one additional exact
record is `unstable_hessian`. All other cell/method combinations have 2000 valid
point fits and intervals, including every Euler fit. No unavailable record is
discarded from operational coverage. For example, exact theta at n=250,q=2,CV=0
has 1832 covering intervals: operational coverage 1832/2000=91.60%, versus
conditional coverage 1832/1961=93.42%. Conditional point metrics likewise describe
the 1961 valid fits, not the rejected boundary candidates.
Each of exact and PFML has 97 negative theta lower limits among its 1961 valid
intervals (194 estimator records in total); these limits remain unclipped, as
registered. The 86 unavailable fit records occur on 46 distinct
source paths; estimator failures on a shared path are not independent events.

The largest recorded condition number is about 23,216.74, below the 1e10 gate.
The one rejected exact Hessian has whitened full/half-step discrepancy about
0.095652, exceeding 0.01; a small ordinary Frobenius discrepancy does not override
that separate gate. These are preserved diagnostic classifications. Week 7 will
investigate whether a software defect or an actual finite-sample/numerical issue
explains them before any stronger causal claim is made.

### Irregularity does not uniformly worsen every metric

The regular-reference comparisons match nominal n, truth and mean gap, and use
independent source paths. They do not match realized span. At n=4000,q=2,CV=3,
exact theta RMSE is 0.0472633 versus 0.0579830 for the regular reference: ratio
0.815123, approximate 95% Monte Carlo interval [0.778016,0.854000]. Exact sigma's
ratio is 0.445275 [0.426063,0.465352], whereas mu's is 1.681608
[1.608710,1.757810]. A blanket monotone deterioration statement is unsupported.

The paired estimator table also retains common-valid counts. At n=250,q=2,CV=0.5,
exact versus PFML theta point comparisons use 1993 common-valid paths; the
paired RMSE ratio is 0.979086 [0.917394,1.044928]. Their operational-coverage
difference uses all 2000 planned paths and is +22.15 percentage points, with a
paired MCSE of 0.9315 percentage points.
The two comparisons answer different questions and use different denominators.
Their uncertainty has not been replaced with independent-sample formulas.

Nor is exact estimation uniformly lowest in finite-sample RMSE. At n=250,q=0.1,
CV=3, theta RMSE is 0.474796 for exact, 0.353144 for PFML and 0.315932 for Euler.
Exact's positive bias is 22.65%, PFML's is -10.11% and Euler's is -17.54%.
Different mixtures of bias and dispersion can reverse the RMSE ranking; this
does not remove each method's separately reported practical classification.

### Absolute mu results and unresolved classifications

Mu relative-to-truth metrics remain undefined. For the regular n=250,q=0.1 cell,
all three methods have signed bias -0.004096 in original state units (MCSE
0.002240), RMSE 0.100242 (MCSE 0.001586), and coverage 89.50%.
Exact mu at CV=3 in that same n/q design
has bias -0.001487 (MCSE 0.002412), RMSE 0.107849 (MCSE 0.001799), and coverage
89.15%. Small mean bias alone therefore does not certify interval calibration.
Across E2, exact mu coverage ranges from 89.15% to 95.65%; Euler mu reaches 100%
in one cell. High coverage can accompany broad working intervals and is not a
claim of superior estimation accuracy.

The five unresolved parameter-specific classifications are retained explicitly:

''' + unresolved + '''

These labels are not forced into either category. In particular, the E2 PFML
sigma n=1000,q=0.5,CV=0.5 result is separate from the independent E1 result at the
same nominal design. The two headline experiments are not pooled or selected
between to obtain a preferred classification.

### Reproducibility, review and remaining work

Use the commands in the README with the recorded numerical environment. E2's
recorded freeze is the actual published commit; unlike the earlier E1 local-
metadata mapping, E2 needs no commit-object restoration helper. The 657-test
pre-run suite and two complete development validations precede E2. All 72
historical scientific source/artifact files remain identical to the main handoff
baseline, and the frozen Week 5 engine/inference sources remain unchanged.

Independent review additionally streamed all 360,000 saved records to reconcile
source sharing, status composition, all failure exports and provenance. A second
reviewer recomputed raw metrics for six preselected cells (12,000 paths, 36,000
fit records), checked all 540 classifications and 108 class totals, and passed
43,944 recorded checks. See `artifacts/week6_evidence/independent_n1000_review.json`
and `independent_n4000_review.json`, with their standalone scripts. All nine
headline PNGs were visually inspected across three reviewers before the archive
was sealed `completed`.

The full headline run passed on its first attempt without the development
manifest-state incident recurring. That success does not determine the incident's
cause. Its evidence remains open for the Week 7/8 operational review. Headline
outputs will remain immutable during diagnostic follow-up. No personal theory
mastery, completion of readings, novelty or publication readiness is certified
by these computational results. No professor email or approval was sought.
'''
path = ROOT / 'docs/week6_progress.md'
old = path.read_text()
path.write_text(old[:old.index('## Headline E2 results')] + section)
