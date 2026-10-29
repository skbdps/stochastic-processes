# Week 5 mathematics and reporting contract

Prepared 30 October 2026 (India time), before the E1 headline run. This is an
implementation and derivation record, not an amendment to the frozen scientific
protocol. The [preregistration](../preregistration.md) and
[dated amendment](plan_amendment_week3_validation.md) govern the experiment.
The [Week 3](week3_mathematics.md) and [Week 4](week4_mathematics.md) notes define
the unchanged model, estimators, numerical normalization and interval procedure.
Actual results and publication checkpoints belong in the separate progress report.

## Design and inferential scope

E1 uses theta=1, mu=0, sigma=0.5, q=0.5, n=1000 observations and R=2000
independent source paths in each of seven nominal-CV cells:
0, 0.25, 0.5, 1, 1.5, 2 and 3. There are 999 transitions per path, 14,000 paths,
42,000 estimator fits and 63 parameter-summary rows. All three estimators see
the same path within a cell. Different cells have independent source streams.

The exact transition simulator starts from the stationary distribution; all
fitted conditional objectives omit the initial-state density. At nonzero CV,
Gamma gaps have shape CV^-2 and scale nominal_mean_gap*CV^2. At nominal CV>=1.5,
the literal floor is 1e-6*nominal_mean_gap. There is no subsequent rescaling.
The stored source diagnostics retain raw and realized moments, floor/clipping
information, minimum actual gap and realized span. Neither a nominal CV label
nor a matched regular reference asserts equal realized spans.

The original numerical services are composed without altering their source.
The full 3-by-3 conditional-NLL Hessian is evaluated in observed-data normalized
coordinates. The registered steps, full/half-step comparison, whitened
comparison, positive definiteness, condition-number and standardized-score
gates remain unchanged. A trace adapter records the provider's returned fine
Hessian and gradient without changing the calculations or results. Natural
covariance uses J=diag(theta_x,state_scale,sigma_x); the symmetric 95% Wald
endpoints are retained even when a positive parameter's lower endpoint is
negative. PFML/Euler curvature covariances remain model-based working estimates.

## Complete accounting

The planned Cartesian product is cell x replication x estimator. A complete
experiment must contain exactly that set, not merely complete groups among
those that happen to be present. Simulation, point-fit and interval failures
retain their rows, numerical candidates, status, reasons and diagnostics.
Unexpected dependency exceptions block a completed implementation claim.

For one parameter, let S be the point-valid count, V the interval-valid count,
C the number of valid intervals covering truth, and R the planned count.
Point bias/RMSE condition on S valid fits. Conditional coverage is C/V, missing
when V=0; operational coverage is C/R, with unavailable intervals contributing
zero. Point and interval availability are S/R and V/R. Every proportion has a
plug-in binomial Monte Carlo SE and a 95% Wilson interval. All-zero or all-one
plug-in SEs are zero, but their Wilson intervals still express uncertainty.

For errors e on the point-valid subset, bias is mean(e), RMSE is sqrt(mean(e^2)),
and the sample SD convention uses denominator S-1. Therefore

    MCSE(bias) = sample_sd(e) / sqrt(S)
    MCSE(RMSE) = sample_sd(e^2) / (2 sqrt(S) RMSE).

These uncertainties are unavailable when S<2; the RMSE delta approximation
also requires positive RMSE. Relative bias/RMSE and their SEs divide by the
absolute truth only when truth is nonzero. Mu truth is zero here, so those
relative quantities remain missing; its absolute bias/RMSE are fully reported.

Empirical SD versus reported model SE is a descriptive comparison. Both use
the same interval-valid subset: sample SD of its point estimates, mean of its
reported SEs, and square root of the mean squared reported SE. Using the
point-valid subset for the first number and the interval-valid subset for the
second would introduce a hidden selection difference.

## Independent regular-reference comparisons

For every cell, find its CV=0 control with exactly the same n, truth, nominal
mean gap and floor settings. Missing references make the complete-grid analysis
fail explicitly. This guards against silently producing a partial result.

Let A denote an irregular cell and B its regular control. Each uses its own
point-valid subset, with sizes S_A and S_B. Write b_A,b_B for bias, r_A,r_B for
RMSE, and u_A,u_B for their RMSE MCSEs. Independence gives

    excess_bias = b_A - b_B
    MCSE(excess_bias)^2 = MCSE(b_A)^2 + MCSE(b_B)^2.

For positive RMSEs, Q=r_A/r_B and first-order propagation gives

    MCSE(log Q)^2 = (u_A/r_A)^2 + (u_B/r_B)^2
    MCSE(Q) = Q * MCSE(log Q).

The ratio's approximate 95% Monte Carlo interval is
Q*exp(+/- z_0.975*MCSE(log Q)), where z_0.975=1.959963984540054. This positive
log-ratio reporting convention concerns simulation-summary uncertainty; it is
not a substitute for the registered natural-scale parameter Wald intervals.
A zero numerator and positive denominator yield ratio zero but missing delta
uncertainty. A zero or missing denominator yields a missing ratio. Mu's RMSE
ratio is well-defined when its denominator is positive: it divides two RMSEs,
not by the zero true mean.

Operational coverage differences subtract C_A/R_A and C_B/R_B. Their independent
variance estimates add the component plug-in binomial variances. The optional
conditional-coverage difference instead uses V_A and V_B, comparing the two
selected interval-valid populations; its denominators are retained explicitly.
Difference MC intervals use the normal approximation and are not constructed
by subtracting Wilson endpoints. Each component's Wilson interval stays in the
main parameter summary.

A regular row compared with itself is an identity: defined differences are zero
with zero MC uncertainty, and a defined positive-RMSE ratio is one with zero MC
uncertainty. Undefined inputs remain undefined. Treating that row as two
independent Monte Carlo samples would fabricate uncertainty.

Each irregular-versus-regular contrast individually uses independent cell
samples. Multiple contrasts reuse the same regular baseline and are therefore
correlated with one another. Their pointwise intervals are not simultaneous
grid guarantees. These nominal-design comparisons include random-span and
floor effects; they do not isolate a pure CV effect or implement the formal
Ait-Sahalia--Mykland information decomposition.

## Paired estimator comparisons

Within a cell, take the intersection of the point-valid paths for estimators
A and B, of size S. Recompute both methods' metrics on precisely that subset.
Subtracting their separately summarized survivor means would estimate a
different quantity. Report R, both individual point-valid counts, S, both
common-subset biases/RMSEs, and the direction A-minus-B or A-over-B.

For paired errors a_i,b_i, define U_i=a_i^2, V_i=b_i^2, M_A=mean(U),
M_B=mean(V), R_A=sqrt(M_A), R_B=sqrt(M_B). Then

    bias_difference = mean(a-b)
    MCSE(bias_difference) = sample_sd(a-b)/sqrt(S)
    mse_difference = mean(U-V)
    MCSE(mse_difference) = sample_sd(U-V)/sqrt(S).

The RMSE-difference delta influence is U_i/(2 R_A)-V_i/(2 R_B). Its sample SD
divided by sqrt(S) estimates the MCSE. For the log-RMSE ratio, use influence
(U_i/M_A-V_i/M_B)/2, then multiply its MCSE by Q=R_A/R_B for the ratio MCSE.
Centering each influence changes no sample variance. This formulation retains
the squared-error covariance and avoids cancellation from subtracting large
nearly equal variance terms. The saved squared-error covariance supports
independent arithmetic checks. Delta derivatives require positive RMSEs;
estimating their sample Monte Carlo variance requires at least two common paths.

For operational coverage, define I_A and I_B as valid-and-covering indicators
on every planned path. Unavailable intervals are zero, even if finite diagnostic
candidate endpoints remain. Set D=I_A-I_B. The paired difference is mean(D)
and its MCSE is sample_sd(D)/sqrt(R). Four counts (both cover, A only, B only,
neither) expose its pairing. With discordant counts n10 and n01, the equivalent
estimated variance is

    [(n10+n01)/R - ((n10-n01)/R)^2] / (R-1).

This sample-variance convention differs from a plug-in population variance
by R/(R-1). It is documented and matches the other paired mean contrasts.
For R=1 the difference exists but this MCSE does not. Shared identical coverage
indicators produce zero paired uncertainty even when each separate coverage
estimate has nonzero binomial uncertainty.

## Practical criteria and interpretation

The existing ThresholdPolicy is reused. For theta and sigma, a supported
breach requires the operational-coverage Wilson upper bound below 0.80, or
the entire relative-bias MC interval beyond +0.10 or -0.10. Within-criteria
requires the coverage lower bound at least 0.80 and the entire relative-bias
interval within [-0.10,+0.10]. Otherwise the classification is unresolved.
Missing bias uncertainty cannot support within-criteria, but cannot erase a
supported coverage breach. Mu has no zero-truth relative-bias classification.

These are pointwise statements on the registered grid, conditional on the
specified fit/interval procedures and explicit availability. A correct exact
transition likelihood does not imply finite-sample unbiasedness or 95% coverage.
No finding is required to show collapse, a preferred ranking or a smooth
boundary. Written mathematical notes do not certify personal theory mastery.

## Validation evidence and audit limits

Hand-calculated fixtures include independent errors [0,2] versus [0,1], giving
excess bias 1/2 with MCSE sqrt(5)/2 and RMSE ratio 2 with MCSE sqrt(2). The
paired fixture [1,2,3] versus [0,1,1] gives RMSE ratio sqrt(7) with delta MCSE
sqrt(3)/4. Different valid subsets, disjoint subsets, zero RMSE, unavailable
intervals, exact coverage endpoints, malformed flags and shuffled row order
are tested explicitly. Independent review verifies these formulas separately
from the implementation.

The artifact audit checks every planned saved outcome and its derived
arithmetic. The stored fine Hessian/gradient/covariances support direct numerical
checks; the coarse Hessian is not stored. Predetermined source replays and
refits at replications 0 and R-1 of every cell verify actual array hashes and
the complete procedure, including its full/half-step calculations. This is a
bounded reproducibility audit, not a proof that every optimizer found a global
maximum or every floating-point derivative equals its analytic derivative.
