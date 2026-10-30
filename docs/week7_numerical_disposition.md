# Independent disposition of the completed selected-path diagnostics

Reviewed 30 October 2026 from the closed `artifacts/week7_selected` archive,
recording published code freeze `6c97d939b4a9902e4127751ff263b5e03ebe7188`.
This review generated no observations, fits, objective evaluations, profile
searches or finite-difference stencils. It reads the already saved results.

Evidence: the [selected archive audit](../artifacts/week7_selected/audit.json),
[independent saved-data review](../artifacts/week7_evidence/selected_saved_evidence_independent_review.json),
[review script](../artifacts/week7_evidence/review_selected_saved_evidence.py),
and [actual figure review](../artifacts/week7_evidence/selected_visual_review.json).
The separate full-study results are in the [Week 7 completion report](week7_progress.md).

## Outcome and verification scope

All 53 fixed selected sources are retained and evaluated. The archived replay
audit reports 159 verified original estimator records, no failed sources, and
no replacement of headline results. Its explicit scope remains source/fit
replay followed by immutable evidence readback; it does not independently
recalculate all diagnostic mathematics.

The additional script `review_selected_saved_evidence.py` performs **4,396
independent saved-data checks**, all passing. It verifies file bindings,
boundary-decision logic, agreement of the independently calculated OLS slopes
with archived PFML unrestricted coefficients, saved profile comparisons and
selection rules, and algebra on **352 saved Hessian matrices**. The whitened
discrepancy is independently reconstructed as a symmetric generalized
eigenvalue problem, rather than by repeating the producer's Cholesky formula.
All archive file bytes remain unchanged. The JSON review records the exact
input hashes, checks, tolerances, counts, limitations and both detailed cases.

Specifically, matrix verification uses binary64 NumPy/SciPy linear algebra on
the saved numerical entries: 296 float64 matrices and the binary64 displays of
56 independent longdouble matrices. The extended decimal entries are checked
to round identically to their saved displays. Spectral algebra is not repeated
at longdouble or arbitrary precision, and the Hessians themselves are not
independently reconstructed from objective evaluations.

The matrix-comparison tolerance is only a reviewer arithmetic tolerance; no
registered fit, interval or classification threshold changes. Observed arrays
are not stored in these payloads, so this review does not independently
recompute OLS or the likelihood directly from observations.

## Boundary outcomes

The 46 selected sources with boundary outcomes all have a nonpositive saved
OLS slope and positive regression residual variance. The independent OLS
diagnostic agrees with the archived unrestricted PFML coefficient. Consequently
all **46 PFML boundary records** are consistent with the conditional AR(1)
likelihood supremum at phi=0: theta tends to infinity, rather than attaining a
finite interior estimate.

For **39 regular sources**, the observed gaps are exactly equal. The same
quadratic-SSE argument therefore applies to their **39 exact-estimator
boundary records**. Together these account for all **85 point-invalid records**.
Their finite stored boundary candidates remain ineligible for ordinary point
or interval summaries. The invalid original points were not promoted into
Hessian or profile targets.

For the remaining seven CV=0.5 sources, a nonpositive OLS slope proves a
boundary only for the PFML working likelihood. Their unequal observed gaps
make the exact OU likelihood a different function. All seven exact fits remain
point-valid; one has the retained unavailable interval discussed below.

## Fixed profile searches

There are 14 valid exact fits and therefore 14 evaluated exact profiles;
the other 39 profiles are explicitly skipped because their original exact
point was invalid. Every evaluated profile reached both prescribed tail
targets and retained all five fixed local starts. All 70 local searches
reported convergence, but convergence is not a guarantee of the best objective.

All 14 original full fitted objectives agree with the best saved finite search
candidate within the fixed diagnostic resolution. The largest reported
stored-full-minus-best difference is about **4.77e-12**, far smaller than the
corresponding resolution. No optimizer failure or erroneous original optimum
is established by these comparisons. A finite search still does not prove a
global optimum.

Three local searches finish above the best objective by more than the fixed
resolution. Two stop on essentially iid upper-tail plateaus: the +2 start for
the regular control and the +2 start for the rejected irregular case. The -2
start for the irregular case finds a separate poorer finite local minimum.
The complete fixed grid and other starts retain the better solution; none of
these local outcomes is silently omitted.

## The rejected exact interval: E2 n=250, q=2, CV=0.5, replication 465

The original fitted natural parameters are theta=16.9118447112,
mu=0.0371230570 and sigma=1.8346580317. The point remains valid; the original
interval remains `unstable_hessian`.

| Diagnostic | Registered float64 | Independent longdouble |
|---|---:|---:|
| Whitened discrepancy, h=.001 versus .0005 | 0.0956515203 | 0.0956538193 |
| Whitened discrepancy, h=.0005 versus .00025 | 0.0233311472 | 0.0233549034 |
| Whitened discrepancy, h=.00025 versus .000125 | 0.0057427564 | 0.0058047901 |
| Condition number at h=.0005 | 103,967.42 | 103,966.96 |
| First-coordinate Hessian diagonal at h=.0005 | 124.5076766 | 124.5076766 |
| Algebraic Schur complement at h=.0005 | 0.0074843582 | 0.0074843909 |

The full matrix is positive definite, but it contains a very weak curvature
direction. At the original comparison the Frobenius discrepancy is only
2.19e-6, while the whitened discrepancy is about 0.09565, above the registered
0.01 gate. This is precisely why the directional gate matters.

**Additional precision does not remove this discrepancy.** The independent
longdouble result reproduces it, and successive step halvings reduce it by
approximately fourfold. This pattern supports **finite-difference truncation
amplified by weak curvature** as the main explanation. At finer steps float64
roundoff also becomes visible, but float64 cancellation alone does not explain
the original rejection. The registered rejection is consistent with its
specified numerical procedure and must not be reclassified using a finer step.

The stored full objective is only 2.97e-15 above the best finite diagnostic
candidate, against a resolution of 4.534e-7. The iid limit is 0.0001483221
higher in NLL: separated at the diagnostic numerical resolution, yet nearly
flat in practical likelihood terms. This supports weak information about the
theta/sigma direction. The lower-start finite local minimum is approximately
0.04658313 above the best objective; the stored estimator reached the better
solution. The negative PFML OLS slope, -0.0010189155, is not evidence of an exact
likelihood boundary under these irregular gaps.

## Accepted regular control: E2 n=250, q=2, CV=0, replication 1424

The original theta is 3.1476157946, with positive OLS phi=0.00184508194.
The original interval remains valid. Its first whitened discrepancy is
**0.00680701495**, below the registered 0.01 threshold. Independent longdouble
gives **0.00680726171**, then 0.00169891928 and 0.000424543715 at the finer
adjacent pairs. The same approximately fourfold decay occurs.

At h=.0005 the condition number is about 23,216.7 and the Schur complement is
0.0335153, compared with a theta-coordinate diagonal of 124.5229. The point is
interior but weakly separated from the iid closure: the latter is only
0.0004238120 higher in NLL. The stored full objective matches the finite-search
best within about 5.55e-17. Its accepted status and the irregular case's rejected
status are both consistent with the unchanged numerical gate.

## Interpretation limits and disposition

The saved longdouble precision is **64 significand bits in 16 storage bytes**
(18 reported decimal digits), not IEEE 128-bit arithmetic or arbitrary
precision. The independent Hessian ladder lifts the same float64-normalized
data as registered inference. The profile uses separately recorded longdouble
normalization of the original observations. Representation differences are
tiny and explicitly retained; they are not a substitute covariance procedure.

These selected diagnostics establish no new defect in the registered OU
likelihood, estimator mapping or inference gate. They support the boundary
interpretation, correct retention of the unstable interval, and weak-curvature
finite-difference sensitivity. The algebraic Schur complements are at the
stored working point and are not unconditionally claimed to equal profile
curvature away from nuisance stationarity.

These outcomes are exploratory selected-path evidence. They do not estimate
population failure frequencies or replace the independent fixed follow-up,
the complete floor-sensitivity experiment, or the preserved E1/E2 results.
The historical unattributed manifest-state incidents remain a separate
operational limitation for the Week 8 gate.
