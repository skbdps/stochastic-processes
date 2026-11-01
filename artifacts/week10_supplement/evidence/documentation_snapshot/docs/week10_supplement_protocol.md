# Week 10 supplement: observation span and the PFML mean mechanism

Dated 1 November 2026, after the original E4 results and Week 10 first draft.
The author requested this extension after observing E4. It is an exploratory,
post-result supplement, not a retroactive preregistration or a replacement of
E4. It explicitly amends the historical Week 9 experiment-freeze schedule.
Original experiments, failures, figures and seals remain unchanged. The
question is explanatory; neither persistence nor disappearance of the effect
is required for success.

## Questions and sequence

1. How do the original MMPP/IID differences change after descriptive
   standardization to a common realized-span distribution?
2. Do differences remain in a separately generated design with exactly equal
   observed spans and random observation counts?
3. Can OU covariance and the actual PFML mean formula quantitatively explain
   the precision loss and the working standard error's failure to follow it?

This protocol and its machine-readable settings are recorded before producing
new span-effect tables or fitted supplement results. Development uses separate
seeds. The complete implementation and its dependencies are hash-bound before
full execution. Later repairs and failed attempts must remain visible.

## A. Saved-record span analysis

Use every original E4 fitted record, retaining planned, point-valid and
interval-valid counts. Compare independently generated MMPP and IID cells at
r=1,5,20, for all three estimators and parameters. Mu is the motivating focus;
other outcomes are retained as context. Original data are already observed.

Use the fixed, left-closed/right-open span bins
`[0,400,450,475,500,525,550,600,infinity]`. Publish all bins, including empty
ones. For each ratio, common support includes bins containing at least 30
source paths in each arm. This selection uses spans and counts, never fitted
errors or coverage. Shared weights are each eligible bin's pooled source
count divided by the total pooled count on common support. Report exclusions
separately for each arm and retain the unadjusted comparison.

Standardize operational coverage, availability, bias and mean squared error
to those weights; take the square root only after averaging MSE. Point
metrics condition on valid points within each bin. MC uncertainty uses
independent paths and independent arms, conditional on realized bin counts,
support and weights. It does not account for estimating that standardization
distribution. Within-bin residual span variation remains. Conditioning on T
also changes the distribution of schedules; this is a descriptive diagnostic,
not a causal direct-effect identification or a fixed-horizon experiment.

## B. Fixed-horizon experiment

Retain truth `(theta,mu,sigma)=(1,0,0.5)`, nominal mean gap m=0.5, the same
calibrated arrival/switch rates, and the same parent event-gap floor rule as
[E4](week8_exploratory_protocol.md). Use **T=499.5**, the original nominal
999-transition span. Run seven independent cells: MMPP and full-parent-
marginal IID at r=1,5,20, plus regular spacing m. Use exactly **2,000 paths per
cell**, 14,000 paths and 42,000 fits. Methods share a path; cells are independent.

Time zero is an observation and an event; initialize the phase from the same
Palm distribution alpha as E4. Retain phase between MMPP gaps and reset from
alpha for every IID gap. Generate each complete parent gap by the competing
arrival/switch construction. Apply the original posthoc floor to r=20 gaps,
without latent-clock padding. Keep arrivals strictly before T. Observe X at
the deterministic endpoint **T**, so every actual observation span is exactly
499.5. A parent arrival exactly at T supplies that endpoint once. Record the
first excluded parent gap, the terminal interval and whether the endpoint was
an arrival. The terminal interval is a boundary observation interval, not a
new draw from the parent gap distribution; it is not floored or normalized.

N, the number of observations including both endpoints, is random. A field
named `n` retained for compatibility with historical cell identities denotes
the nominal reference 1000; every source must also record `actual_n`, and
reporting must label the distinction. The regular cell has N=1000.

The two stochastic mechanisms share the parent complete-gap marginal and
the same boundary rule. Stopping and endpoint insertion mean that the finite
window's empirical retained-gap distributions need not be identical. This
contrast tests performance at a common observed span; it does not hold N
fixed or identify a universal dependence-only effect. Initializing at an
event is also distinct from an arbitrary calendar-time observation window.

Full root seed: **20261101**. Development root seed: **21261101**, four paths
per cell. Source-only validation seed: **22261101**. New namespace: **10**.
SeedSequence addresses bind the complete scientific cell identity and
replication, with separate schedule and OU streams; worker count and partition
size do not change results. No new seed or additional replications are chosen
in response to observed effects.

Generate exact stationary OU transitions. Reuse unchanged exact, realized-
mean-gap PFML and actual-gap Euler fits, full 3x3 observed Hessians and
natural-scale unclipped 95% Wald intervals. No robust-SE replacement, bias
correction, truth-based numerical scaling or outcome-based exclusion is added.
Operational coverage includes all planned paths; unavailable intervals count
as noncoverage. Point errors condition on valid estimates. All inference is
pointwise; mu relative-to-truth error remains undefined.

Report all cells/parameters, independent MMPP-minus-IID coverage differences
and RMSE ratios with MC intervals, and N/terminal-gap/floor diagnostics.
Compare the fixed-horizon contrasts with original E4 descriptively; their
different seeds, random N and endpoint convention prohibit paired inference.

## C. Quantitative mechanism

Start with the interior PFML AR(1) regression identity, including its endpoint
correction. Use the stationary OU covariance
`Cov(X_s,X_t)=sigma^2/(2 theta) exp(-theta |t-s|)` to derive the conditional
variance of a linear mean estimator at a fixed working AR coefficient.
Compute this in linear time and check it against dense covariance arithmetic.

For the original fixed-N event-stationary design, integrate the joint
phase-type gap kernel, including the floor, to calculate lag correlations and
finite-N variance predictions. Compare MMPP with its matched IID kernel.
Evaluate the original nominal working-variance and Gaussian-coverage
approximations separately from exact finite-N linear-estimator variance.
Reconstruct each new fixed-horizon schedule for conditional variance predictions
and compare their averages with the empirical fitted-mean MSE/variance.

These calculations use the model's known generating law as an analytical
benchmark, not as an estimator input. Distinguish exact identities from the
fixed-rho approximation to a fitted-rho estimator and from a normal-coverage
approximation. Fitted nuisance parameters, endpoint correction, random N and
studentization require explicit limitations. No general dependent-gap CLT or
new robust confidence-interval procedure is claimed.

## Validation, evidence and completion

Before the full run: test sampling/address independence, endpoint and floor
rules, r=1 Poisson limit, regular limit, complete accounting, refusal to
overwrite existing evidence, deterministic resume, retained failure behavior,
and independent analytical covariance calculations. Audit the complete
development run, including all-source replay and selected full-fit replay.

After the full run: verify all planned keys and partition hashes, replay every
source and endpoint fits at replications 0 and R-1 in each cell, and rederive
reported summaries from raw records. Preserve source/runtime/thread/package
identity and all failed checks. Inspect every new figure as an actual image;
bind the inspected bytes and numerical tables before sealing. Recheck the
completed archive and original-source preservation. Update the Week 10
manuscript, source map, status and reproduction instructions with outcomes and
limitations, including null or contrary findings. Completion requires actual
results and validation, not merely implemented code.
