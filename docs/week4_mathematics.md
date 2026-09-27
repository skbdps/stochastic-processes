# Week 4 mathematics: uncertainty and the validation checkpoint

This document explains the inference addition to the OU project. Its governing
design is the original Week 4 row together with the
[dated Week 3 plan amendment](plan_amendment_week3_validation.md), which supersedes
outcome-dependent certification language. Numerical settings and experiment
selection are fixed in [preregistration.md](../preregistration.md); execution
evidence belongs in the Week 4 progress report, not in these derivations.

The code implements a statistical procedure. Producing these notes does not
establish that the owner completed the plan's personal S4 theory session or its
exercises. Tang–Chen's full-text working-depth reading also remains outstanding;
no unverified coefficient or coverage guarantee from that paper is used.

## 1. Mathematical and software scope

| Mathematical task | Implementation responsibility |
| --- | --- |
| Distinguish observed information, expected Fisher information and model misspecification | Evaluate the selected conditional objective; label its inverse-Hessian covariance accurately |
| Differentiate in stable coordinates and transform the full covariance | Separate normalization and coordinate maps from Hessian evaluation |
| Construct a specified 95% interval with explicit availability criteria | Preserve point-fit validity, interval validity and their separate reasons |
| Check the finite-sample Wald procedure at the planned regular cell | Retain all 2,000 replication outcomes and quantify Monte Carlo uncertainty |
| Compare with a specified conditional parametric bootstrap | Fix timestamps and initial value; keep all 500 scheduled bootstrap draws |
| Interpret coverage and bias without hiding numerical failures | Report planned and valid denominators, uncertainty and unresolved threshold classifications |
| Freeze future experimental choices before E1/E2 | Record grids, seeds, replications, floor, interval policy and no-peeking rules |

The Week 3 simulation, spacing, objectives, point estimators and saved experiments
retain their definitions. Week 4 adds uncertainty on top of those APIs.

## 2. Likelihood, score and information

Write the natural parameters as $\beta=(\theta,\mu,\sigma)$ and the working
parameters as $\eta=(\log\theta,\mu,\log\sigma)$. With $n$ observations there are
$K=n-1$ transitions. For each estimator, the total conditional negative
log-likelihood or Gaussian contrast is

$$
L(\eta)=\frac12\sum_{i=1}^{K}
\left[\log(2\pi)+\log v_i(\eta)
+\frac{\{x_i-m_i(\eta)\}^2}{v_i(\eta)}\right].
$$

For exact OU, $m_i=\mu+(x_{i-1}-\mu)e^{-\theta\Delta_i}$ and
$v_i=\sigma^2(1-e^{-2\theta\Delta_i})/(2\theta)$. PFML replaces every gap in
these expressions by the observed mean gap. Euler uses
$m_i=x_{i-1}+\theta(\mu-x_{i-1})\Delta_i$ and $v_i=\sigma^2\Delta_i$.
Every fit conditions on observed $x_0$ and timestamps. The stationary initial
density is absent even when simulation initially draws $X_0$ from stationarity.

The log-likelihood score is $s=-\nabla L$; observed information is
$H=\nabla^2 L$. The model-based covariance approximation is

$$
C_\eta=H(\widehat\eta)^{-1}.
$$

The objective is a **sum**, so its Hessian already accumulates the information
from all $K$ transitions. No further division by $K$ belongs in this covariance.
Using the Hessian of an averaged objective would require a corresponding factor.

Expected Fisher information is instead $I(\eta)=E_\eta[H(\eta)]$; under suitable
regularity conditions it also equals $E_\eta[ss^T]$. These are expectations under
the specified model, rather than identities asserting that a realized Hessian
equals a realized score outer product. Stanford's
[Fisher-information lecture](https://web.stanford.edu/class/archive/stats/stats200/stats200.1172/Lecture15.pdf)
states the regular likelihood identities and large-sample normal approximation.
Its displayed iid theorem is background, not a proof of OU finite-sample
normality or an iid assumption about the observed path.

For orientation, direct differentiation of one Gaussian transition, with its
lagged observation held fixed, gives the following **repository derivation**:

$$
I_{i,ab}(\eta\mid x_{i-1},t)
=\frac{\partial_a m_i\,\partial_b m_i}{v_i}
+\frac{\partial_a v_i\,\partial_b v_i}{2v_i^2}.
$$

To see this, set $r=x_i-m_i$. The score component is
$r\,\partial_a m_i/v_i+(r^2-v_i)\partial_a v_i/(2v_i^2)$; take products using
$E[r]=E[r^3]=0$, $E[r^2]=v_i$ and $E[r^4]=3v_i^2$. Expected information for the
entire conditional path additionally averages over the random lagged states.
Week 4 implements the **observed** Hessian, not a separate expected-information
estimator; the equation above clarifies the distinction.

### Misspecified contrasts require a different interpretation

Exact OU uses the data-generating transition model. PFML at irregular times and
Euler at nonzero gaps generally do not. Their inverse Hessians describe local
curvature of those contrasts and are labelled **model-based**, not robust
covariance estimates. White's
[original misspecification paper](https://faculty.utrgv.edu/diego.escobari/teaching/Econ8370/Papers/White(1982)Econometrica-Corey.pdf)
establishes why the usual information equality need not survive model
misspecification. The generic asymptotic form is
$A^{-1}BA^{-1}$, with $A$ the limiting curvature and $B$ the score covariance.
For dependent data, $B$ may require a long-run covariance rather than the iid
outer-product formula. A sandwich estimator is not implemented here.

Even a correct sandwich around a pseudo-true parameter would not automatically
remove its discrepancy from the true OU parameter. We measure coverage of the
specified intervals against the simulation truth; we do not assume PFML or
Euler intervals are nominally calibrated. All finite-sample coverage, including
exact-MLE coverage, is an empirical question.

## 3. Normalization and the covariance map

For an observed dataset, choose $a$ as its observed state mean, $b>0$ as the
centered root-mean-square state deviation (denominator $n$), and $h>0$ as its
realized arithmetic mean gap. Write $y=(x-a)/b$, $u=(t-t_0)/h$; subtracting a
fixed time origin is immaterial because only gaps enter. These constants are
held fixed throughout differentiation.
They are numerical scaling choices, not additional OU parameters and not
functions of a candidate parameter. Simulation truth is never used to choose
them.

The parameter map and conditional objective identity are

$$
\theta_x=\theta_u/h,\qquad
\mu_x=a+b\mu_u,\qquad
\sigma_x=b\sigma_u/\sqrt h,\qquad
L_x=L_u+K\log b.
$$

The constant objective shift does not affect derivatives. There is no
time-density Jacobian: timestamps are conditioned on, not fitted as random
observations of a parametric spacing model.

Let $C_{\eta_u}$ be covariance in normalized working coordinates. Applying the
first-order delta method to the natural-parameter map gives

$$
J=\operatorname{diag}(\widehat\theta_x,b,\widehat\sigma_x),
\qquad C_{\beta_x}=J C_{\eta_u}J^T.
$$

The full matrix is transformed, including off-diagonal terms. The
[delta-method lecture](https://web.stanford.edu/class/archive/stats/stats200/stats200.1172/Lecture17.pdf)
explains the first-order Taylor argument; the multivariate map here is derived
for this repository's coordinates. Equivalently, first transform to original
working coordinates with $A=\operatorname{diag}(1,b,1)$, then apply
$D=\operatorname{diag}(\widehat\theta_x,1,\widehat\sigma_x)$. Since $DA=J$, both
routes agree. Applying $D$ directly to $C_{\eta_u}$ would omit the mean's scaling.

These equations also specify a useful invariance test. Under state change
$x'=c+d x$ and time change $t'=k t$, where $d,k>0$,

$$
\beta'=(\theta/k,\ c+d\mu,\ d\sigma/\sqrt k),\qquad
C_{\beta'}=M C_\beta M^T,\quad
M=\operatorname{diag}(1/k,d,d/\sqrt k).
$$

Interval endpoints and standard errors must follow the same affine maps.
They must not depend on which physical measurement units were supplied.

## 4. Numerical Hessian and interval validity

Week 4 uses central differences at the fitted normalized working coordinates.
Writing $f=L_u$ and step $d_j=10^{-3}\max(1,|\widehat\eta_j|)$,

$$
H_{jj}\approx\frac{f(\eta+d_j e_j)-2f(\eta)+f(\eta-d_j e_j)}{d_j^2},
$$

$$
H_{jk}\approx\frac{f(\eta+d_j e_j+d_k e_k)
-f(\eta+d_j e_j-d_k e_k)
-f(\eta-d_j e_j+d_k e_k)
+f(\eta-d_j e_j-d_k e_k)}{4d_jd_k},\quad j\ne k.
$$

The implementation repeats the calculation at half the steps, records
agreement and selects the finer Hessian and gradient. Its acceptance settings
are project decisions fixed before the
checkpoint, not theoretical constants taken from a paper. They require a valid
interior point fit, finite derivatives, positive Hessian eigenvalues, a
normalized-coordinate condition number at most $10^{10}$ and relative
Frobenius Hessian disagreement at most $0.005$. The standardized score diagnostic
must satisfy

$$
\sqrt{g^T H^{-1}g}\le 0.01,\qquad g=\nabla L_u(\widehat\eta_u).
$$

This is a dimensionless local measure of remaining displacement from a stationary
point. Optimizer termination alone is insufficient. A Frobenius check can hide
error in a weak-curvature direction, so a further pre-frozen gate compares the
matrices after whitening by the fine Hessian:

$$
H_{\rm fine}=LL^T,\qquad
E=L^{-1}(H_{\rm coarse}-H_{\rm fine})L^{-T},\qquad
\max_j|\lambda_j(E)|\le0.01.
$$

For every nonzero direction $v$, this bounds the magnitude of the coarse/fine
quadratic-form discrepancy relative to $v^T H_{\rm fine}v$ by 1%. It directly
addresses directions that can dominate the inverse. Agreement of two numerical
steps still does not prove agreement with the exact derivative; independent
derivative and coordinate-invariance tests remain necessary. None of these
numerical gates proves global optimality, asymptotic normality or statistical
identification, and no failed Hessian is silently regularized.

After valid covariance transformation, each natural parameter receives the
specified symmetric Wald interval

$$
[\widehat\beta_j-z_{0.975}\sqrt{C_{\beta,jj}},\;
 \widehat\beta_j+z_{0.975}\sqrt{C_{\beta,jj}}],
\qquad z_{0.975}\simeq1.959964.
$$

This is a natural-scale interval using delta-method standard errors; it is not
an exponentiated log-scale interval. A negative lower endpoint for positive
$\theta$ or $\sigma$ is retained. Clipping it or rejecting it would silently
change the procedure. Nonfinite covariance, nonpositive variances, boundary or
nonidentified fits, unstable derivatives and failed numerical gates instead
produce an unavailable interval with its reason recorded. The candidate point
estimate remains in the raw record.

## 5. Conditional parametric bootstrap

The
[Stanford bootstrap lecture](https://web.stanford.edu/class/archive/stats/stats200/stats200.1172/Lecture19.pdf)
describes parametric resimulation from a fitted model. Our OU adaptation is
explicitly conditional: for a selected observed path and selected estimator,
take its fitted $(\widehat\theta,\widehat\mu,\widehat\sigma)$, hold the original
timestamps and observed $x_0$ fixed, simulate a fresh exact OU path, and refit the
same estimator. Repeat the scheduled $B=500$ draws using independent recorded
child seeds. New spacings and stationary initial values are not drawn.

This preserves dependence through the OU recursion; independently resampling
individual observations would not implement this procedure. Numerical
normalization may be recomputed by each bootstrap refit, then mapped back to
the original units. That does not change the fixed timestamps or initial
condition.

The interval is the empirical 2.5% and 97.5% quantiles of valid bootstrap point
estimates, using NumPy's explicitly selected
[`method="linear"`](https://numpy.org/doc/stable/reference/generated/numpy.quantile.html).
For sorted values $v_0,\ldots,v_{m-1}$, write $(m-1)p=j+w$; the quantile is
$(1-w)v_j+w v_{j+1}$, with the endpoint convention when $j=m-1$.
At least $\max(2,\lceil0.95B\rceil)=475$ valid draws are required; otherwise the
interval is unavailable. Validity requires an interior point fit, successful
optimization, finite objective and parameters, and positive theta and sigma.
All attempted draws, including failures, remain recorded.
There is no redraw-until-success rule. The 95% availability cutoff is a
prespecified operational choice, not a theorem that discarding up to 5% failed
fits leaves inference unbiased. If the original point fit is invalid, report
the bootstrap as skipped with 500 planned and zero attempted draws. With at
least two valid draws, their sample standard deviation (denominator $m-1$) is
retained as a descriptive conditional spread even when endpoints are withheld;
with fewer than two, it is missing.

The fixed spot checks select replication 0 of two cells: $q=\theta\bar\Delta=0.5$,
$n=1000$, CV 0 and 2. Each checks exact, PFML and Euler. A pair of individual
datasets can expose disagreement in intervals but cannot estimate bootstrap
coverage. With 500 draws, a nominal 2.5% tail contains only about 12.5 draws, so
endpoint Monte Carlo uncertainty is material.

For PFML and Euler, simulating exact OU at an estimator-specific plug-in and
refitting that estimator may reproduce its pseudo-parameter or discretization
bias. These percentile intervals are diagnostic comparisons, not automatic
bias corrections, sandwich replacements or calibrated intervals for the true
OU parameters. Basic, studentized, BCa and bias-corrected bootstrap intervals
would be separately specified methods; none is silently substituted.

## 6. Coverage denominators and simulation uncertainty

For one estimator, interval method and parameter, let $R$ be the planned number
of independent path replications, $V$ the number with available valid intervals,
and $C$ the number of those intervals containing the true parameter. Report all
three counts, together with

$$
\widehat a=V/R,\qquad
\widehat p_{\rm valid}=C/V,\qquad
\widehat p_{\rm procedure}=C/R.
$$

Conditional coverage is missing when $V=0$. The planned-denominator value is
the rate of successfully producing a covering interval, with unavailable
intervals counted as procedure failures. It is not mislabeled as coverage of
the available intervals. A conditional value can look good even when a method
fails frequently; the companion availability rate makes that visible.

For a Bernoulli proportion $\widehat p$ over $N$ independent trials, its plug-in
Monte Carlo standard error is

$$
\operatorname{MCSE}(\widehat p)
=\sqrt{\widehat p(1-\widehat p)/N}.
$$

Use $N=R$ for procedure success and $N=V$ for the conditional proportion.
Independent simulation replications are the trials, not individual time-series
observations. A zero plug-in MCSE when all observed trials agree does not imply
the probability is known exactly.

Wilson intervals provide a bounded uncertainty summary even at 0 or 1. With
$z=z_{0.975}$, the center and half-width are

$$
c=\frac{\widehat p+z^2/(2N)}{1+z^2/N},\qquad
w=\frac{z\sqrt{\widehat p(1-\widehat p)/N+z^2/(4N^2)}}{1+z^2/N}.
$$

The interval is $[c-w,c+w]$. The
[NIST statistical handbook](https://www.itl.nist.gov/div898/handbook/prc/section2/prc241.htm)
gives this score-inversion formula and distinguishes it from the Wald interval
for a proportion. This Monte Carlo interval describes uncertainty in a measured
coverage probability; it is a different interval from the OU parameter's Wald
interval being evaluated.

Bias uses valid point estimates and its own denominator $S$, which can exceed
$V$. For nonzero truth $\beta_{0j}$,

$$
\widehat r_j=\frac{\overline{\widehat\beta_j}-\beta_{0j}}{|\beta_{0j}|},\qquad
\operatorname{MCSE}(\widehat r_j)
=\frac{\operatorname{sd}(\widehat\beta_j-\beta_{0j})}
{|\beta_{0j}|\sqrt S}.
$$

An approximate Monte Carlo interval is $\widehat r_j\pm z\operatorname{MCSE}$.
At $\mu_0=0$, relative bias is undefined and remains missing; mean inference is
still reported in absolute units. Invalid point fits must not be replaced by
truth or silently removed from planned counts.

## 7. Checkpoint interpretation and future threshold claims

The planned exact-MLE regular-grid checkpoint is $q=0.5$, $n=1000$, $R=2000$,
with truth $(1,0,0.5)$. If actual coverage were 0.95, its Monte Carlo SE would be
$\sqrt{0.95\cdot0.05/2000}\simeq0.004873$, or 0.487 percentage points. The
historical 94–96% band is a diagnostic reference, not a tuning target or a
required result. It is neither sufficient nor necessary evidence that every
implementation detail is correct. Each parameter's availability, coverage,
uncertainty and relevant numerical diagnostics must be inspected.

The classification policy is frozen before headline E1/E2. A point estimate crossing
coverage 0.80 or absolute relative bias 0.10 is not by itself a statistically
resolved boundary. For each of theta and sigma:

- **Supported breach:** the Wilson upper limit for planned-denominator
  operational coverage is below 0.80, or the relative-bias Monte Carlo interval
  lies wholly above +0.10 or wholly below -0.10.
- **Within the stated criteria:** the Wilson lower limit is at least 0.80 and
  the entire relative-bias Monte Carlo interval lies within [-0.10,+0.10].
- **Unresolved:** neither preceding condition holds. Missing bias uncertainty
  prevents a within-criteria claim, but cannot erase a supported coverage breach.

Retain fit and interval failure rates alongside every classification. These
rules evaluate the planned-denominator procedure-success rate, rather than
allowing a small successful subset to certify the whole procedure.
Pointwise Monte Carlo intervals are not simultaneous guarantees across all
grid cells and parameters.

If implementation checks pass and undercoverage persists, report that result
for the stated finite-sample procedure. Do not change seeds, omit failures,
adjust intervals or repeat the checkpoint until its coverage enters a desired
band. A later alternative method requires a dated separate specification.
Likewise, a well-tested grid may show no established breach; no operating
boundary needs to be invented. The two nominal-design spacing controls still
do not isolate a pure CV effect or constitute Aït-Sahalia–Mykland's formal
information decomposition.
