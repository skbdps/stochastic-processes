# Week 7 E3 exploratory protocol

Dated 30 October 2026 (India time), 29 October 2026 UTC. Prepared after completed
Week 6 publication `7ed521b4567da5c973c409e47e1b1d9906fa3b47`, before any affected
diagnostic replay or new E3/floor-sensitivity observations. The implementation
and this protocol will be committed before those computations begin.

This is a new, outcome-informed exploratory protocol, not prior registration of
E3 methods and not an amendment to E1/E2. The original preregistration and all
headline archives remain unchanged. Computational results do not certify the
owner's theory mastery or reading completion. No email or outside approval gate
is part of this work.

## 1. Selection fixed before reruns

Verify both completed headline seals, their audit bindings and all artifact
hashes. Retain every source path with any point-invalid, interval-unavailable or
pipeline-error record. In the published inputs this selects 46 E2 source paths
and 86 affected estimator records; E1 has none. Additionally, for each headline
experiment and estimator, select the interval-valid record with the largest
finite whitened Hessian discrepancy and with the largest standardized score.
Break ties by lexicographically smallest scientific cell ID then smallest
replication number, and take the union of sources and selection reasons.
Inspect all three estimators for each selected source. Selection uses saved
records only and is fully serialized before any selected replay.

Retain all six unresolved parameter classifications across E1/E2 in the selection
manifest. They are descriptive follow-up issues, not a reason to change a
threshold, choose between E1/E2, pool their records or repeat until a label changes.

## 2. Boundary, optimizer and Hessian diagnostics

Recreate selected headline sources using their original root seed and scientific
address. Verify exact observation/time hashes and the original three-method
records. Independently calculate the intercept-and-slope least-squares PFML
coefficient. For regular exact/PFML data, a nonpositive fitted AR coefficient
places the constrained positive-correlation likelihood optimum at the phi=0
boundary, corresponding to theta tending to infinity; it does not supply a
finite interior OU estimate.

For selected valid exact fits, independently evaluate the conditional profile
likelihood in observed-data normalized coordinates, using extended precision
where available, a fixed log-theta grid and five deterministic starts. Compare
with the stored fit and the iid limiting likelihood. Record objective values,
local solutions and numerical resolution; finite searches do not prove a
universal global optimum. Also calculate the theta-to-zero Brownian-with-drift
limit with alpha=theta*mu allowed to remain finite. The 1D profile curvature is
a Schur complement, not the theta diagonal of the full information matrix.

For normalized observed gaps d, fix the log-theta range
L=max(-40,min(-30,log(1e-10)-log(max(d)))) and
U=min(40,max(12,log(50)-log(min(d)))). Use a grid with spacing at most 0.05 and
polish every finite sampled local-minimum bracket. Record whether the target
tail inequalities theta_L*max(d)<=1e-10 and theta_U*min(d)>=50 are reached;
cap-limited tails remain unresolved. Evaluate the stored fit as a candidate.
Five genuine local-search initializations use its normalized log-theta plus
{-2,-1,0,1,2}, clipped to [L,U], and are passed to an optimizer that accepts x0.
The reviewed implementation records optimizer settings, convergence and every
candidate. Extended-precision evaluation records the platform's actual mantissa;
it is not claimed to be arbitrary precision. Objective differences within
1e-7+1e-9*abs(best NLL) are numerically unresolved for disposition, with unrounded
differences retained. This diagnostic resolution rule does not replace any
headline fitting or inference gate.

For selected valid points, inspect the full Hessian at relative steps 0.001,
0.0005, 0.00025 and 0.000125, reporting eigenvalues, condition, directional
stability and score where defined. These are diagnostic step comparisons only.
No new Hessian, interval or covariance replaces the registered headline result.
An additional exact-likelihood ladder uses the same original fitted eta and the
same float64 normalized sample/gaps as registered inference, but evaluates an
independent full conditional NLL and its central differences in longdouble.
Only completed Hessian/gradient arrays are converted to float64 for matrix
diagnostics. All four steps remain recorded; no replacement intervals are
constructed. The profile scan instead normalizes original observations/gaps
in longdouble with the recorded observed scales. Record these representation
differences and both full-NLL values so precision changes are not conflated.
In particular, inspect the rejected E2 exact fit at (n=250,q=2,CV=0.5,rep=465)
and the accepted high-discrepancy regular control at rep=1424. Diagnose any
verified software defect separately from a weak-information sample or working-
model failure. Preserve unresolved numerical uncertainty when evidence is weak.

## 3. Fixed independent 10,000-replication follow-ups

Both and only E2 cells with numerical unavailability are selected:
(n=250,q=2,CV=0) and (n=250,q=2,CV=0.5). Run exactly R=10,000 new paths per cell,
root seed 20261007, using the unchanged registered OU/Gamma source addressing
`SeedSequence(seed, spawn_key=(4,*eight_cell_id_words,rep)).spawn(2)` and the
unchanged three estimators and Wald gates. Truth remains (1,0,0.5). The new seed
is independent of E1/E2. This is justified by the boundary rates and the rare
Hessian failure, not a desire to force coverage into a preferred band.

Retain every planned record, including failures and candidates, with durable
partitions and identity-verified resumption. Report all original summary,
uncertainty, paired/independent-contrast and failure fields. These 20,000 independent new source paths and 60,000 fit records are separately
labelled exploratory follow-ups; headline counts and labels are never replaced
or pooled. Complete the fixed R regardless of observed interim results.

## 4. Paired epsilon sensitivity

Use q=0.5,n=1000,truth=(1,0,0.5), CV in {0,1.5,2,3}, and floor fractions
epsilon in {1e-8,1e-6,1e-4}. R=2000 per variant, root seed 20261009. The twelve
variants contain 24,000 retained path realizations and 72,000 fitted records,
but only 8,000 independent (CV,rep) coupling groups.
This is a 100-fold bracket on each side of the registered floor, fixed before
sensitivity data are generated. Apply the floor for nominal CV>=1.5 with no
renormalization. Retain clip counts, moments, minimum gaps and realized spans.

For each (CV,rep), anchor the source's scientific address to that CV's cell ID
with the registered epsilon=1e-6. Compose `OuSampleSource` using this anchor ID
but the variant's actual floor argument. Thus all three epsilon variants share
the underlying Gamma draws and OU normal innovations, including stationary X0;
different CVs remain independent. Each stored variant retains its actual cell
ID, anchor ID, floor fraction and full spawn keys. This deliberate common-random-
number construction is newly specified here and is not a claim about E1/E2
pairing. CV=0 is an exact negative control: all three floor variants must produce
the same times, values and numerical results. It is one shared control sample,
not three independent control samples.

The shared standardized innovations couple different exact OU paths. They are
not the same continuous OU trajectory observed at different times. Changing the
floor may change gaps, realized span and states. Numerical normalization is
recomputed from each variant's own observed data, without simulation truth.

Compare each nonreference epsilon against 1e-6 within each CV/method/parameter.
Recompute point metrics on the common-valid paired subset, retain its count,
and use paired differences and delta-method RMSE-ratio uncertainty. Operational
coverage differences use all planned pairs and unavailable intervals count as
noncoverage. Zero-truth mu relative metrics remain undefined. Full per-variant
summaries and availability are also retained. Separate regular-reference
contrasts, if produced, match nominal design; multiple contrasts sharing the
same paired regular control are not mutually independent.

A separate audit explicitly validates the new anchor addressing and pairing;
these source overrides must not pass or be represented as the registered-source
audit. It reuses reviewed pure inference/arithmetic checks through composition,
checks complete keys and source sharing, verifies all source hashes by replay,
and replays full fits at replications 0 and R-1 in each variant. Figures regenerate
from saved tables. Report small, null, robust or adverse findings as observed.

## 5. Validation, execution and reporting boundaries

Development runs use the same definitions with seed+1,000,000 and small fixed R.
They are distinct from exploratory evidence and cannot supply its reported
results. Freeze reviewed dependencies before full runs. Preserve E1/E2 bytes,
the original development manifest-state incident, diagnostic inputs and any
failed attempt. Do not weaken acceptance gates or change a failed manifest by
hand. Record whether the incident recurs; nonrecurrence alone does not identify
its cause. Week 8 will decide whether the remaining operational uncertainty
requires buffer work rather than optional E4.

Complete and publish a one-page findings memo and detailed numerical disposition
after audits and visual review. No professor message will be sent. A subsequent
protocol change must be dated before its newly affected observations, with all
previous outputs retained.
