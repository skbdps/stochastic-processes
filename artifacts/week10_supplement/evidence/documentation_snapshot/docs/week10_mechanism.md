# Week 10: a quantitative mechanism for PFML mean uncertainty

Dated 1 November 2026. This is a post-result analytical supplement under the
[Week 10 protocol](week10_supplement_protocol.md). It preserves the original
estimators, intervals and experiment records. The calculation distinguishes
an exact estimator identity, exact variances of fixed-coefficient linear
diagnostics, and approximations to the fitted estimator and its coverage.

## 1. Begin with the implemented mean estimator

Write the observations as X_0,...,X_K, where K=N-1. At an interior PFML
solution, its conditional Gaussian AR(1) regression has fitted slope
rho_hat and intercept a_hat. The intercept normal equation gives

$$
\widehat\mu_P
=\frac{\widehat a}{1-\widehat\rho}
=\frac{\sum_{i=1}^K X_i-\widehat\rho\sum_{i=0}^{K-1}X_i}
       {K(1-\widehat\rho)}
=\overline X_{0:K-1}
 +\frac{X_K-X_0}{K(1-\widehat\rho)}. \tag{1}
$$

This is an algebraic identity for the estimator used in the study. The
realized mean gap maps rho_hat into theta_hat but does not enter (1) directly.
The mean is an event-index average plus an endpoint correction. Interior
observations receive equal weights even when many occur in a short burst;
the first observation has a negative weight in the equivalent linear
expression when rho is held fixed. It is therefore imprecise to identify
PFML mu exactly with the ordinary average of all N observations.

Equation (1) explains the relevant weighting issue. A cluster of closely
spaced observations contains many highly correlated measurements of nearly
the same OU state. Counting those observations equally does not give each
unit of calendar time equal influence. This is a precision mechanism under
exogenous sampling, not evidence that bursts select systematically high or
low values of X. The stationary conditional mean of each observation remains
mu_0. The original E4 results likewise did not establish a mean-bias shift.

Hold a coefficient rho in (0,1) fixed and call the resulting diagnostic
mu_tilde(rho). Its weights are

$$
w_0=-\frac{\rho}{K(1-\rho)},\qquad
w_i=\frac1K\ (1\le i<K),\qquad
w_K=\frac1{K(1-\rho)},\qquad \sum_i w_i=1. \tag{2}
$$

At the parent-law target phi=E[exp(-theta_0 Y)], the exact discrepancy is

$$
\widehat\mu_P-\widetilde\mu(\phi)
=\frac{(X_K-X_0)(\widehat\rho-\phi)}
 {K(1-\widehat\rho)(1-\phi)}. \tag{3}
$$

For the fixed-count design, if rho_hat converges to phi<1 and the stationary
endpoints remain bounded in probability, this difference is o_p(K^-1).
If rho_hat is also root-K consistent, it is O_p(K^-3/2). These are conditional
implications of stated regularity, not a new proof of consistency or a
dependent-gap central limit theorem. Boundary solutions and finite-sample
near-unit slopes remain governed by the existing validity rules.

## 2. Exact variance conditional on an exogenous schedule

For stationary OU data, with V=sigma_0^2/(2 theta_0),

$$
\operatorname{Cov}(X_{t_i},X_{t_j}\mid t_0,\ldots,t_K)
=V e^{-\theta_0|t_i-t_j|}.
$$

Consequently, any weights fixed given the schedule have conditional variance

$$
v_w(t)=V\left\{\sum_iw_i^2+
2\sum_{i<j}w_iw_j e^{-\theta_0(t_j-t_i)}\right\}. \tag{4}
$$

This expression is exact at finite N. It integrates over the stationary X_0
draw, as the simulation summaries do; the conditional likelihood's exclusion
of the stationary initial density does not mean that repeated-sample variance
conditions on a particular numerical value of X_0.

The implementation evaluates (4) through the OU innovations representation
without allocating a dense covariance matrix. Put
a_i=exp[-theta_0(t_i-t_(i-1))]. Set b_K=w_K and work backwards with
b_i=w_i+a_(i+1)b_(i+1). Then

$$
v_w(t)=V\left\{b_0^2+\sum_{i=1}^K(1-a_i^2)b_i^2\right\}. \tag{5}
$$

The squared innovations form is stable even with signed endpoint weights.
It requires O(N) work. Dense covariance products independently check it.
For the ordinary sample average set every weight to 1/N; for the PFML
diagnostic use (2) at the fixed parent-law coefficient phi.

Both (4) and (5) apply to the newly generated fixed-horizon schedules,
including their deterministic endpoint and random N. Conditional unbiasedness
implies that averaging these conditional variances over the schedule law
gives the unconditional MSE of the fixed-coefficient diagnostic. The source
replication is the independent Monte Carlo unit for this average. Plugging
in rho_hat from the same X values would invalidate this interpretation because
the weights would then be random given the schedule.

## 3. Exact event-index correlation under the original E4 design

Use the [existing phase-type notation](week8_mathematics.md): A=Lambda-Q,
Palm row vector alpha, and joint gap/ending-phase density
exp(-A t)Lambda. A raw gap D is transformed to Y=max(D,f), with the adopted
floor f (zero when inactive). Define the matrix Laplace kernel

$$
B_s=\int_0^\infty e^{-s\max(t,f)}e^{-At}\Lambda\,dt
=e^{-sf}\left\{A^{-1}(I-e^{-Af})\Lambda
       +e^{-Af}(A+sI)^{-1}\Lambda\right\}. \tag{6}
$$

At f=0 this reduces to (A+sI)^-1 Lambda. At s=0 it is the event-phase
transition matrix P. The atom induced by flooring is included in the first
term; the hidden phase is not advanced through artificial floor padding.

For the stationary MMPP-derived sequence, integrating its ordered joint
kernel gives

$$
c_h=\frac{\operatorname{Cov}(X_i,X_{i+h})}{V}
=E\left[e^{-\theta_0(Y_{i+1}+\cdots+Y_{i+h})}\right]
=\alpha B_{\theta_0}^{\,h}\mathbf1. \tag{7}
$$

The IID control instead has c_h=phi^h, where phi=alpha B_(theta_0) 1.
Thus matching the complete one-gap marginal matches c_1, but does not match
the covariance at lags two and beyond. At r=1 the Poisson limit makes both
sequences coincide. At r=5 and r=20 in this calibration, the MMPP correlations
at lags 2 through 20 exceed those of the matched IID sequence, as explicitly
checked by the tests. All 999 lags enter the finite-count variance calculation
below, yielding the larger variance in the prediction table. No universal
correlation ordering for all dependent gap processes is asserted.

For fixed N, a deterministic weight vector therefore has exact unconditional
variance V times the quadratic form in the Toeplitz matrix with entries
c_|i-j|. In particular, (2) yields

$$
\frac{\operatorname{Var}\{\widetilde\mu(\phi)\}}V
=\frac{K-1}{K^2}+\frac{1+\phi^2}{K^2(1-\phi)^2}
 +\frac{2}{K^2}\sum_{h=1}^{K-1}(K-h)c_h
 -\frac{2\phi c_K}{K^2(1-\phi)^2}. \tag{8}
$$

This is the exact finite-count variance of a fixed-coefficient diagnostic.
It automatically includes the original design's random spans and their
dependence on its gaps. It is not a span-adjusted result.

The corresponding ordinary-average variance is

$$
\operatorname{Var}(\overline X)
=\frac V{N^2}\left\{N+2\sum_{h=1}^{N-1}(N-h)c_h\right\}. \tag{9}
$$

Since B_(theta_0) is substochastic with spectral radius below one in this
calibration, its correlation sum is available without truncating a series:

$$
G_M=1+2\alpha B_{\theta_0}(I-B_{\theta_0})^{-1}\mathbf1,
\qquad G_I=\frac{1+\phi}{1-\phi}. \tag{10}
$$

The leading variance coefficient is VG/K. This is an exact limit of the
linear-estimator variance times K. It is not, by itself, a normal-limit or
coverage theorem for the fitted estimator.

## 4. Why the original working standard error misses the increase

The PFML population residual variance is V(1-phi^2). In its stationary
Gaussian AR(1) working model the mean coordinate's population information
gives

$$
v_{\mathrm{work}}=\frac{V(1-\phi^2)}{K(1-\phi)^2}
=\frac VK\frac{1+\phi}{1-\phi}. \tag{11}
$$

The population mean-versus-slope information cross term is zero after
centering, so the jointly estimated slope does not alter this leading
working coefficient. The observed full Hessian can still have nonzero cross
terms in a finite sample. Equation (11) is a benchmark, not a replacement for
the actual full-Hessian interval used in the study.

The matched IID sequence has c_h=phi^h. Substitution in (8) gives exactly
v_work for the fixed-coefficient mean, even though random-gap transitions
do not make the entire unconditional observation process a Gaussian AR(1).
The MMPP has the same phi and thus the same working variance, while its
higher-lag covariances increase (8). This explains why the original reported
SEs can remain similar while the empirical spread increases.

A deliberately approximate coverage prediction is

$$
C_{\mathrm{normal}}
=2\Phi\!\left(1.9599639845
\sqrt{\frac{v_{\mathrm{work}}}
 {\operatorname{Var}\{\widetilde\mu(\phi)\}}}\right)-1. \tag{12}
$$

It substitutes a normal law and a fixed population interval width. A random
schedule creates a mixture of conditional Gaussian laws even for the linear
diagnostic; fitted slopes and studentization introduce further differences.
Consequently, (12) is not an exact coverage prediction or a new calibrated CI.

## 5. Deterministic magnitude predictions

The following values use the original N=1000, truth (1,0,0.5), nominal mean
gap .5 and the used floor law. They are generated by
`population_predictions` in [the theory module](../ou_irregular/week10/theory.py).
No coefficient is fitted to observed errors or coverage.

| r | IID fixed-rho SD | MMPP fixed-rho SD | SD ratio, MMPP/IID | MMPP normal working-coverage approximation |
|---|---:|---:|---:|---:|
| 1 | 0.02501251 | 0.02501251 | 1.000000 | 95.0000% |
| 5 | 0.02635641 | 0.02840156 | 1.077596 | 93.1063% |
| 20 | 0.03076606 | 0.03576256 | 1.162403 | 90.8229% |

These predictions concern SD, while the motivating empirical headline uses
RMSE. A numerical comparison should retain that distinction and report mean
bias and empirical SD alongside RMSE. The supplement's generated results
tables bind the comparison to the original data and report the fixed-horizon
diagnostics separately.

For the fixed-horizon experiment, (7)--(10) cannot be reused by merely replacing
N with its mean. Stopping, random N and the terminal boundary alter the joint
law. Instead, use each saved schedule in (5) with its actual N, then average
the conditional variances. Phi remains a specified parent-law reference
coefficient there; it is not claimed to be the exact finite-window expected-
objective minimizer. This conditional calculation is the mechanism check
that accompanies the span-controlled experiment.

## 6. Relationship to existing theory and limits

The Markovian arrival kernel and distinction between event and time
stationarity are standard. [Asanjarani and Nazarathy (2022), *Stationary
Markovian Arrival Processes: Results and Open Problems*](https://doi.org/10.1017/S1446181122000013),
Section 2, supplies the event-stationary phase-type framework. Equations
(1)--(12) above are explicit project derivations from that framework and the
OU covariance; no priority claim is made for these identities.

[Ait-Sahalia and Mykland (2004), *Estimators of Diffusions with Randomly
Spaced Discrete Observations: A General Theory*](https://www.princeton.edu/~yacine/gmmrandom.pdf),
Assumption 2 on page 2190, uses IID sampling intervals. Its observation-window
setup permits a random count, and equations (14) use covariance sums for
estimating equations. Those results motivate keeping covariance beyond the
Hessian visible, but do not supply a ready-made theorem for our dependent,
jointly fitted, finite-sample design. Sources were checked on 1 November 2026.

This supplement establishes finite-design covariance identities and tests a
quantitative approximation. It does not establish a general dependent-gap
CLT, a robust PFML confidence interval, or performance outside the specified
OU model, sampling laws and finite designs. The fitted-mean diagnostics and
span-controlled experiment must be assessed using their actual saved results,
including any failures or contrary findings.
