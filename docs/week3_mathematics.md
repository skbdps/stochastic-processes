# Week 3 mathematics and implementation contract

This document turns the Week 3 **harness v1** row of `OU_master_plan_v2.pdf` into mathematical definitions and implementation requirements. The accompanying progress record reports execution results. The derivations below were worked through for this repository's jointly estimated three-parameter OU model; they are not transcriptions of a paper's coefficient tables. This is development and smoke validation, not a certified coverage study or frozen preregistration.

Updated 27 September 2026 for the [numerical validation update](week3_validation_update.md).
The [plan amendment](plan_amendment_week3_validation.md) records post-exploratory
interpretation decisions and the Week 4 covariance-coordinate handoff.

## 1. What Week 3 requires

| Mathematical task | Programmatic task | Location |
|---|---|---|
| Specify equidistant and Gamma gaps, nominal CV, flooring, and the resulting moment changes | Generate positive gaps and retain clipping and span diagnostics | `ou_irregular/spacing.py` |
| Define conditional exact, PFML, and Euler Gaussian objectives consistently | Expose the three estimators with explicit convergence diagnostics | `ou_irregular/estimators.py`; shared likelihoods in `ou_irregular/ou.py` |
| Define bias, RMSE, their relative versions, and Monte Carlo uncertainty | Aggregate raw estimates without hiding failed fits | `ou_irregular/metrics.py` |
| Separate nominal grid coordinates, actual gap law, and realized observation span | Validate one YAML, seed each replication, simulate exact paths, fit each method, save raw output and configuration | `ou_irregular/runner.py` |
| Define precisely what a regular-grid twin controls | Compare each irregular cell with its regular counterpart and plot transparent labels | `ou_irregular/metrics.py`, `ou_irregular/plots.py` |
| Read the required literature at the depth needed to implement and interpret these choices | Record source conventions, applicability, and unresolved claims | This document |

The smoke design uses three values of $q=\theta\bar\Delta$, three CV values, and all three estimators. In the supplied exploratory configuration: $q\in\{0.1,0.5,1\}$, $\mathrm{CV}\in\{0,1,2\}$, $n=250$ **observations**, hence $K=n-1=249$ transitions, 20 independent replications per cell, and floor fraction $\varepsilon=10^{-6}$. Throughout this document, configuration field `n` means observation count; $K$ means transition count. The nominal observation span is $(n-1)\bar\Delta$, not $n\bar\Delta$. These settings test the harness; 20 replications cannot support a publication's operating-region or coverage claim.

## 2. Primary-literature reading notes and limits

### Aït-Sahalia and Mykland (2003)

The paper distinguishes FIML, which uses observed intervals; IOML, which integrates unobserved intervals out; and PFML, which treats intervals as a fixed mean. Its formal **cost of randomness** compares IOML with FIML; **cost of ignoring randomness** concerns PFML. These are not synonyms. Section 3's sampling assumptions distinguish progressively stronger exogeneity restrictions; the PFML analysis uses independent identically distributed intervals. Section 7 specializes to $dX_t=-\theta X_tdt+\sigma dW_t$, with the long-run mean fixed at zero. Results distinguish joint estimation of drift and diffusion variance from estimation of one parameter with the other known. Its asymptotic expansions concern long observation spans and small-gap expansions, not a numerical guarantee of finite-sample CI coverage. Relevant reading: Sections 2–4, 6, and 7. [Author-hosted published article: *The Effects of Random and Discrete Sampling When Estimating Continuous-Time Diffusions*, Econometrica 71, 483–549](https://www.princeton.edu/~yacine/random.pdf).

**Repository consequence:** retain the PFML attribution, but describe our twins as an empirical comparison inspired by that distinction. We do not implement IOML, and our twin contrasts do not reproduce the paper's formal decomposition.

### Aït-Sahalia and Mykland (2004)

The paper develops asymptotic analysis for estimating equations under random sampling, including exact likelihood and Euler approximations. Section 4.3 uses the zero-mean OU process; Table 1 reports estimation of $\theta$ with $\sigma^2$ known and of $\sigma^2$ with $\theta$ known. Section 4.4 shows why higher moments of the interval distribution matter. An estimator can have smaller variance than exact MLE while being inconsistent, so variance alone is insufficient. Read Sections 2 and 4, especially the Euler objective and the qualifications surrounding Table 1. The journal year is **2004**; the arXiv deposit is dated 2005. [Primary author reprint: *Estimators of Diffusions with Randomly Spaced Discrete Observations: A General Theory*, Annals of Statistics 32, 2186–2222](https://arxiv.org/pdf/math/0503679).

**Repository consequence:** estimate $(\theta,\mu,\sigma)$ jointly and do not transplant a known-nuisance asymptotic variance or bias coefficient into this different estimation problem. Setting the true $\mu$ to zero does not mean its estimate is fixed at zero.

### Tang and Chen (2009): verification status

The publisher-provided bibliographic record identifies the required paper as *Parameter Estimation and Bias Correction for Diffusion Processes*, Journal of Econometrics 149, 65–81. Its abstract describes bias/variance expansions for Vasicek and CIR, approximate likelihood for linear drift, and parametric bootstrap bias correction. [Publisher article](https://www.sciencedirect.com/science/article/pii/S030440760800208X); [publisher-supplied abstract and bibliographic record](https://ideas.repec.org/a/eee/econom/v149y2009i1p65-81.html).

**Access limitation:** the primary publisher full text returned HTTP 403 during this review. The full theorem assumptions and exact nuisance-parameter convention therefore have not been independently verified from Tang–Chen's text. Do not mark that paper's full working-depth reading complete or implement a coefficient attributed to it from memory. The executable harness does not depend on such a coefficient.

A related primary research article by Calderon discusses finite-sample OU bias and extends a Tang–Chen correction to measurement-noise settings; its equation (5) contains a bias contribution proportional to $1/(N\Delta t)$. Its model and estimation setup differ from ours. This supports investigating finite-span effects, but it does not substitute for reading Tang–Chen or prove any coverage threshold in our grid. [Calderon, *Correcting for Bias of Molecular Confinement Parameters Induced by Small-Time-Series Sample Sizes in Single-Molecule Trajectories Containing Measurement Noise*, equation (5)](https://arxiv.org/pdf/1304.4196).

**Correction to the plan's interpretation:** finite-span bias is a reason to measure coverage, not proof that exact-MLE coverage must fall below 80% at $(q,n)=(0.1,250)$. Even a bias formula cannot determine interval coverage without the sampling distribution and interval construction. Retain that corner as a hypothesis to test in Week 4 onward. No Week 3 output establishes a coverage breach.

## 3. Model, observations, and exact simulation

**Original derivation and repository convention.** Let

$$
dX_t=\theta(\mu-X_t)\,dt+\sigma\,dW_t,
\qquad \theta>0,\quad \sigma>0.
$$

There are $n$ observations $x_0,\ldots,x_{n-1}$ and $K=n-1$ transitions, with $\Delta_i=t_i-t_{i-1}>0$ for $i=1,\ldots,K$. All sums over transition index $i$ below run from 1 to $K$. Centering $Y_t=X_t-\mu$ and integrating the linear equation gives

$$
X_{t+\Delta}=\mu+(X_t-\mu)e^{-\theta\Delta}
 +\sigma\int_t^{t+\Delta}e^{-\theta(t+\Delta-u)}\,dW_u.
$$

The deterministic-kernel Wiener integral has mean zero and variance equal to the integral of its squared kernel. Thus

$$
\phi_i=e^{-\theta\Delta_i},\qquad
m_i=\mu+(x_{i-1}-\mu)\phi_i,\qquad
v_i=\frac{\sigma^2}{2\theta}(1-e^{-2\theta\Delta_i}).
$$

Given the gaps, independent $Z_i\sim N(0,1)$ produce an exact discrete observation path through $x_i=m_i+\sqrt{v_i}Z_i$. Start at stationarity,

$$
X_0\sim N\!\left(\mu,\frac{\sigma^2}{2\theta}\right),
$$

independently of future innovations and gaps. This removes a transient caused by a chosen initial value. Estimation nevertheless conditions on the observed $x_0$: the stationary density is **not** added to the fitted objective. An unconditional likelihood would be a different estimator.

Compute $1-e^{-2\theta\Delta}$ as `-expm1(-2 * theta * delta)`. For a tiny product, subtracting two nearly equal floating-point values would lose precision. The exact simulation recursion is used for every estimator comparison; Euler is only a fitting contrast.

## 4. Gap distribution and the floor

**Original derivation.** Use shape–scale Gamma notation, with raw gap $D\sim\Gamma(k,s)$ and density

$$
f(d)=\frac{d^{k-1}e^{-d/s}}{\Gamma(k)s^k},\qquad d>0.
$$

Integrating powers of $d$ and using $\Gamma(k+1)=k\Gamma(k)$ gives $E[D]=ks$, $\operatorname{Var}(D)=ks^2$, and $\mathrm{CV}=1/\sqrt{k}$. Hence for nominal mean gap $\bar\Delta$ and requested CV $c>0$,

$$
k=c^{-2},\qquad s=\bar\Delta c^2.
$$

$c=1$ gives exponential gaps. $c=0$ is a separate equidistant generator, $D_i=\bar\Delta$; passing $c=0$ to the Gamma formula is undefined.

For nominal $c\ge1.5$, the plan applies

$$
D'=\max(D,a),\qquad a=\varepsilon\bar\Delta.
$$

This is flooring, **not** rejection sampling from a truncated Gamma distribution. It creates a point mass at $a$ and leaves the upper tail unchanged. With regularized lower/upper incomplete Gamma functions $P(k,z)$ and $Q(k,z)=1-P(k,z)$,

$$
\Pr(D<a)=P(k,a/s),
$$

$$
E[D']=aP(k,a/s)+ks\,Q(k+1,a/s),
$$

$$
E[(D')^2]=a^2P(k,a/s)+k(k+1)s^2Q(k+2,a/s).
$$

These follow by splitting $E[\max(D,a)^r]$ at $a$ and integrating the Gamma upper tail. The implied floored population CV is

$$
c_{\mathrm{floor}}=
\frac{\sqrt{E[(D')^2]-E[D']^2}}{E[D']}.
$$

**Implementation boundary:** `spacing.py` computes the expected floored mean and expected clipping fraction, and reports empirical raw/floored CVs. The second-moment and population-CV equations above are reference derivations; dedicated routines evaluating them are not part of the Week 3 implementation. Empirical CV uses the population standard deviation (`ddof=0`) of the finite gap vector divided by its empirical mean.

Consequently, flooring increases the expected mean and changes the CV. A nominal label such as `CV=2` identifies the **input Gamma law**, not the exact CV after flooring or the realized sample CV. Small $\varepsilon$ need not imply a small clipping *fraction* when $k<1$: many raw draws can be extremely close to zero even when the total change in elapsed time is tiny.

Record the floor used, clipped count/fraction, raw and floored empirical mean/CV, and actual minimum gap. Do not secretly rescale a generated vector back to its nominal mean: that would introduce dependence among its entries and create another sampling law.

With fixed observation count $n$ and transition count $K=n-1$, distinguish:

$$
T_{\rm nominal}=K\bar\Delta=(n-1)\bar\Delta,\qquad
E[T_{\rm observed}]=K E[D'],\qquad
T_{\rm observed}=\sum_{i=1}^K D'_i.
$$

For regular gaps these coincide. For raw Gamma gaps the first two coincide but the last is random. After flooring even the first two differ. Report dimensionless span $\theta T$ when comparing cells. In particular, with `n=250` and $q=0.1$, $\theta T_{\rm nominal}=249\times0.1=24.9$ under this observation-count convention.

Sampling is exogenous: the gap RNG is independent of the OU innovations. A hitting-time or price-triggered observation rule is outside this harness's model.

## 5. The three conditional objectives

**Original derivation from Gaussian densities.** For any chosen mean/variance pair $(m_i,v_i)$, minimize the conditional negative log likelihood

$$
L(\theta,\mu,\sigma)=\frac12\sum_{i=1}^K
 \left[\log(2\pi)+\log v_i+\frac{(x_i-m_i)^2}{v_i}\right].
$$

| Estimator | Mean | Variance | Role of observed gaps |
|---|---|---|---|
| E-exact | $\mu+(x_{i-1}-\mu)e^{-\theta\Delta_i}$ | $\sigma^2(1-e^{-2\theta\Delta_i})/(2\theta)$ | Uses every actual gap |
| E-PFML | $\mu+(x_{i-1}-\mu)e^{-\theta\widehat\Delta}$ | $\sigma^2(1-e^{-2\theta\widehat\Delta})/(2\theta)$ | Replaces all gaps by $\widehat\Delta=K^{-1}\sum_{i=1}^{K}\Delta_i$ |
| E-euler | $x_{i-1}+\theta(\mu-x_{i-1})\Delta_i$ | $\sigma^2\Delta_i$ | Uses every actual gap in a first-order approximation |

PFML uses the **realized mean of the supplied, already floored gaps**. It is not allowed to use the simulation's hidden nominal mean. This is a feasible plug-in implementation of the fixed-mean likelihood. At CV zero, E-exact and E-PFML must return the same objective for every parameter vector.

Euler comes from $e^{-\theta\Delta}=1-\theta\Delta+O((\theta\Delta)^2)$ and $v=\sigma^2\Delta+O(\Delta^2)$. A small *mean* gap does not make every Gamma draw small; large gaps can invalidate the approximation locally. The Euler mean can have coefficient $1-\theta\Delta<0$ for large gaps. This is permitted in the contrast, but is not exact OU dynamics.

Use working parameters $(\log\theta,\mu,\log\sigma)$ and deterministic multiple starts. Retain success, objective, messages, and boundary information. A finite estimate is not automatically a converged optimum; a numerical bound is not evidence that the scientific parameter lies at a genuine model boundary.

### Profiling reference derivation

The 27 September update implements these profiling identities in the strict
Week 3 fitter. It profiles `alpha = theta * mu` rather than `mu` in the low-theta
tail to avoid an unnecessary unstable division. The original 10 September
implementation used joint L-BFGS-B optimization; its outputs remain archived.

For the exact objective with fixed $\theta$, let $b_i=1-\phi_i$, $y_i=x_i-\phi_i x_{i-1}$, and $w_i=(1-\phi_i^2)/(2\theta)$, so $y_i=b_i\mu+\eta_i$ and $\operatorname{Var}(\eta_i)=\sigma^2w_i$. Differentiating gives

$$
\widehat\mu(\theta)=
\frac{\sum_i b_iy_i/w_i}{\sum_i b_i^2/w_i},\qquad
\widehat\sigma^2(\theta)=\frac1K\sum_{i=1}^K
\frac{(y_i-b_i\widehat\mu(\theta))^2}{w_i}.
$$

Substituting these formulas into $L$ gives a one-dimensional profile objective.
The production search expands a dimensionless range, refines sampled local
minima, and compares endpoint and Brownian-drift/iid limiting objectives. This
finite numerical search does not prove a global maximum for every irregular
path. A separate test reference finds roots of the analytic envelope score in
stationary-variance coordinates; it does not call the production solver.

For PFML, conditional Gaussian regression with an intercept gives OLS slope
`phi` and residual variance `v`. An interior OU map requires `0 < phi < 1`:
`theta = -log(phi)/mean_gap`, `mu = intercept/(1-phi)`, and
`sigma^2 = 2*theta*v/(1-phi^2)`. The profiled SSE is a convex quadratic in phi,
so an unrestricted OLS slope outside this interval puts the constrained
optimum on its closure. A finite numerical approximation to that boundary is
retained but excluded from point summaries.

Euler uses weighted regression of `diff(x)/sqrt(gap)` on `sqrt(gap)` and
`-x_previous*sqrt(gap)`. The coefficients are `theta*mu` and `theta`; diffusion
variance is the residual sum of squares divided by K. Euler requires theta>0,
but does not require `1-theta*gap` to be positive. Degenerate residual variance
and rank-deficient data receive explicit classifications.

All strict fits now use observed-data coordinates `y=(x-a)/b`, `u=t/m`, where
a is the observed state mean, b its positive RMS centered scale, and m the
realized mean gap. Original units are recovered by `theta_x=theta_y/m`,
`mu_x=a+b*mu_y`, `sigma_x=b*sigma_y/sqrt(m)`, and `NLL_x=NLL_y+K*log(b)`.
No truth parameter enters these numerical scales. The affine state Jacobian
has one factor per conditional transition, and conditioning on timestamps
introduces no time-density Jacobian.

### Population PFML diagnostic for the joint model

This is an **original reference derivation for the model fitted by our implementation**, conditional on stationarity, iid exogenous gaps, finite moments, and an interior population optimum. It is a limiting mathematical sanity check, not a finite-sample correction. The raw, unfloored Gamma special case is implemented by `ou.pfml_limit` using `log1p(z)/z`; the floored-Laplace formulas below remain reference derivations. These are not generated simulation results.

Write $V=\sigma^2/(2\theta)$, $m=E[D']$, and $A=E[e^{-\theta D'}]$. The lag-one covariance of centered observations is $VA$. The best population Gaussian AR(1) regression therefore has coefficient $A$, mean $\mu$, and residual variance $V(1-A^2)$. Mapping those moments back into a fixed-gap OU transition yields

$$
\theta_*=-\frac{\log A}{m},\qquad
\mu_*=\mu,\qquad
\sigma_*^2=\sigma^2\frac{\theta_*}{\theta}.
$$

For an **unfloored** Gamma law its Laplace transform is $A=(1+\theta s)^{-k}$, so

$$
\frac{\theta_*}{\theta}
=\frac{\log(1+qc^2)}{qc^2},\qquad q=\theta\bar\Delta.
$$

For floored data substitute the actual $m$ above and

$$
A=e^{-\theta a}P(k,a/s)
 +(1+\theta s)^{-k}Q\!\left(k,(s^{-1}+\theta)a\right).
$$

The raw Gamma expression must not be labelled an exact benchmark for floored simulations. Population bias and finite-span estimation bias may act in different directions, so a 20-replication smoke average need not equal this limit.

## 6. Metrics, missing values, and Monte Carlo precision

**Original definitions and derivations.** For a scalar parameter $\alpha$ with true value $\alpha_0$, let $e_r=\widehat\alpha_r-\alpha_0$ over $R_s$ accepted finite fits:

$$
\widehat B=\frac1{R_s}\sum_re_r,\qquad
\widehat{\operatorname{RMSE}}=
\sqrt{\frac1{R_s}\sum_re_r^2}.
$$

When $\alpha_0\ne0$, the repository defines relative bias as $\widehat B/|\alpha_0|$ and relative RMSE as $\widehat{\operatorname{RMSE}}/|\alpha_0|$. This is normalization by the truth's **magnitude**: relative bias retains the sign of the signed error even when the true parameter is negative. At $\mu_0=0$, both relative quantities are undefined; retain missing values and use bias/RMSE in the original units. Do not divide by a tiny epsilon and call the result relative error. Optional scale normalization would require an explicitly named denominator, such as stationary standard deviation, and would be a separate metric.

Replications, not observations within one correlated path, are the independent Monte Carlo units. The estimated Monte Carlo standard error of bias is

$$
\operatorname{MCSE}(\widehat B)=s_e/\sqrt{R_s},
$$

where $s_e$ uses the sample-variance denominator $R_s-1$. For RMSE, set $u_r=e_r^2$; a delta-method approximation is

$$
\operatorname{MCSE}(\widehat{\operatorname{RMSE}})
\approx\frac{s_u}{2\sqrt{R_s}\,\widehat{\operatorname{RMSE}}}.
$$

It is undefined with insufficient successful replicates or at a degenerate zero RMSE; record that explicitly. Divide these MCSEs by the same absolute normalization factor for relative metrics.

Report both the requested/attempted count $R$ and accepted count $R_s$, failure count, and failure rate. Summaries over accepted fits estimate **performance conditional on success**. A method that fails on difficult paths can appear artificially accurate on its survivors; plotting failure rates alongside its errors makes this selection visible. Keep every failure in raw output and never replace it by the true parameter.

CI coverage is a Week 4 extension. Its future definition is the fraction of valid intervals containing the truth, accompanied by interval-failure counts and a stated denominator. The binomial MCSE formula $\sqrt{\widehat p(1-\widehat p)/R_{\rm valid}}$ does not generate intervals or certify them. Week 3 must not fill uncomputed coverage with zero, 95%, or a point-estimate success rate.

## 7. What the equidistant twin does and does not control

**Repository design definition.** For every irregular cell, its regular twin shares true parameters, observation count $n$ (therefore transition count $K=n-1$), and **nominal** mean gap $\bar\Delta$. Compute the same estimator and parameter on each side:

$$
\text{excess bias}=\widehat B_{\rm irregular}-\widehat B_{\rm regular},
\qquad
\text{RMSE ratio}=\frac{\widehat{\operatorname{RMSE}}_{\rm irregular}}
{\widehat{\operatorname{RMSE}}_{\rm regular}}.
$$

An absent twin or zero/undefined denominator produces a missing contrast, not an invented baseline. Since different cells use independent streams, these are contrasts of independent Monte Carlo summaries; a common replication index does not make the paths statistically paired. If uncertainty on the bias contrast is needed, independent-cell variance estimates add.

This controls $n$ and the nominal time scale. It does not equalize each path's observed span. Flooring also changes expected span, as derived above. Consequently the contrast measures the combined change in the sampling design, including random-span and floor effects; it cannot prove a pure isolated causal effect of CV. Explicit expected-mean or realized-span matched twins would be an additional design choice for preregistration, not a silent adjustment during this smoke run.

## 8. Reproducibility and the Week 4 boundary

NumPy documents `SeedSequence` spawning as a way to produce reproducible child streams with very high probability of independence. [NumPy, parallel random number generation](https://numpy.org/doc/stable/reference/random/parallel.html).

The repository's design derives a stable cell identifier from serialized scientific inputs and uses `SeedSequence(root_seed, spawn_key=(*cell_key_words, replication)).spawn(2)` for gap and path children. The estimator list does not determine data seeds, so estimator order does not change a path. Adding another grid cell should not reshuffle existing cells' random draws. Record the root seed, cell identity, and spawn keys; do not use Python's process-randomized `hash()` as a persistent scientific seed.

Save the resolved configuration, one raw record per estimator/replication, summaries, and figures. The raw record must allow a reader to distinguish a rejected fit, a missing mathematical quantity, and an actual zero. Figures identify the smoke run, nominal CV, estimator names, the versioned point-validity condition, and uncomputed coverage. Approximate MCSE bars quantify uncertainty in Monte Carlo bias/RMSE, not intervals for an individual path's estimate. Strict replay checks source/runtime provenance as well as hashes of actual input arrays when available.

The following remain **Week 4 or later**, irrespective of a successful Week 3 smoke run:

- Full observed-Hessian construction, conditioning checks, and interval validity policy.
- Wald intervals and delta-method covariance transformation for $(\theta,\mu,\sigma)$. The Jacobian from original log working coordinates is $\operatorname{diag}(\widehat\theta,1,\widehat\sigma)$; from the normalized working coordinates it is $\operatorname{diag}(\widehat\theta,b,\widehat\sigma)$. The amendment specifies the coordinate maps; interval implementation and validation remain outstanding.
- The certification cell at $q=0.5,n=1000$ with sufficient replications; smoke point estimates cannot certify coverage.
- Bootstrap spot checks and any finite-span bias correction.
- A committed frozen preregistration of grids, thresholds, seeds, replications, floor fraction, failure policy, and no-peeking rules.
- E1/E2 headline experiments, coverage operating boundaries, and a final comparison with asymptotic predictions.

Professor outreach remains deferred by the owner's decision. No mathematical or programming dependency requires sending that email.
