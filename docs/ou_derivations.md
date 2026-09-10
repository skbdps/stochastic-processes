# OU derivations and estimator definitions

This document supplies the mathematical reference for the Week 2 simulator and
likelihood notebook. **D1–D5 below are this repository's organization.** The
master plan refers to a separate `w1_w2_execution_doc.md`, which was not supplied;
we do not claim that its derivation numbering is reproduced here.

The notes are a completed written deliverable. They do **not** certify that the
owner has completed the plan's personal exercises or passed its timed,
closed-notes theory exit test. An optional self-check appears at the end.

## Setup and assumptions

Let $W_t$ be standard Brownian motion and consider

$$
dX_t=\theta(\mu-X_t)\,dt+\sigma\,dW_t,
\qquad \theta>0,\quad\sigma>0,\quad\mu\in\mathbb R.
$$

The SDE denotes the integral equation

$$
X_t=X_0+\theta\int_0^t(\mu-X_u)\,du+\sigma W_t.
$$

Assume that $X_0$ is fixed or independent of future Brownian increments.
The filtration contains $X_0$ and the Brownian past. Observation times are
strictly increasing. They are either deterministic or generated independently
of the entire OU path, including its initial value. For a random observation
grid, all transition and likelihood statements below condition on the realized
grid. This exogeneity assumption excludes state-dependent observation times;
a hitting-time sample, for example, need not have the transition law below
after conditioning on its observation times.

The only stochastic ingredients needed for D1–D5 are Brownian motion's defining
properties, the Gaussian limit/closure facts used in D1, and interchange of the
ordinary and Wiener integrals for the square-integrable deterministic kernel in
D2. No general Itô formula for nonlinear functions is needed.

## D1 — The Wiener integral for deterministic integrands

For a deterministic step function

$$
f(u)=\sum_{j=0}^{k-1}c_j\mathbf 1_{(u_j,u_{j+1}]}(u),
\qquad a=u_0<\cdots<u_k=b,
$$

define

$$
I(f)=\int_a^b f(u)\,dW_u
=\sum_{j=0}^{k-1}c_j(W_{u_{j+1}}-W_{u_j}).
$$

Brownian increments on disjoint intervals are independent centered Gaussians,
so

$$
\mathbb E[I(f)]=0,
\qquad
\mathbb E[I(f)^2]=\sum_jc_j^2(u_{j+1}-u_j)
=\int_a^b f(u)^2\,du.
$$

Linearity and the same argument on a common partition give

$$
\mathbb E[I(f)I(g)]=\int_a^b f(u)g(u)\,du,
\qquad
\mathbb E[(I(f)-I(g))^2]=\int_a^b(f-g)^2\,du.
$$

For any deterministic $f\in L^2([a,b])$, take step functions $f_k$
converging to $f$ in $L^2$. The final identity makes $I(f_k)$ Cauchy in
mean square. Define $I(f)$ as its mean-square limit. The result is independent
of the approximating sequence, and the identities above pass to the limit.
In particular, the **Itô isometry** is

$$
\boxed{\mathbb E\!\left[\left(\int_a^b f(u)\,dW_u\right)^2\right]
=\int_a^b f(u)^2\,du.}
$$

Each approximating integral is Gaussian; its variance converges to
$\int_a^b f^2$. The limit is therefore

$$
\boxed{\int_a^b f(u)\,dW_u\sim
N\!\left(0,\int_a^b f(u)^2\,du\right).}
$$

Integrals with deterministic integrands supported on disjoint time intervals
are independent. An integral over $(s,t]$ is independent of the filtration at
$s$. These assertions follow first for step sums and then for their joint
limits. The Gaussian conclusion here relies on **deterministic** integrands;
it is not a claim about arbitrary adapted random integrands.

## D2 — Integrating-factor solution of the OU equation

Set $Y_t=X_t-\mu$. Then

$$
dY_t=-\theta Y_t\,dt+\sigma\,dW_t.
$$

Multiplying by the deterministic integrating factor $e^{\theta t}$ suggests

$$
d(e^{\theta t}Y_t)=\sigma e^{\theta t}\,dW_t,
$$

and hence

$$
\boxed{X_t=\mu+(X_0-\mu)e^{-\theta t}
+\sigma\int_0^t e^{-\theta(t-u)}\,dW_u.}
$$

Here is a verification directly in the integral equation. Define

$$
Z_t=\int_0^t e^{-\theta(t-u)}\,dW_u.
$$

Interchanging the ordinary and Wiener integrals for this bounded deterministic
kernel on a finite time triangle yields

$$
\begin{aligned}
\theta\int_0^t Z_s\,ds
&=\int_0^t\left[\theta\int_u^t e^{-\theta(s-u)}\,ds\right]dW_u\\
&=\int_0^t(1-e^{-\theta(t-u)})\,dW_u
=W_t-Z_t.
\end{aligned}
$$

Also

$$
\theta\int_0^tY_0e^{-\theta s}\,ds=Y_0(1-e^{-\theta t}).
$$

Substitution therefore shows that
$Y_t=Y_0e^{-\theta t}+\sigma Z_t$ satisfies
$Y_t=Y_0-\theta\int_0^tY_s\,ds+\sigma W_t$.
For uniqueness, the difference $D_t$ of two solutions driven by the same
Brownian path and initial value obeys
$D_t=-\theta\int_0^tD_s\,ds$; the elementary integral inequality
(or the resulting ordinary differential equation) gives $D_t=0$.

## D3 — Exact conditional Gaussian transition

Split the solution at a fixed time $s$. For $d>0$,

$$
X_{s+d}=\mu+(X_s-\mu)e^{-\theta d}
+\sigma\int_s^{s+d}e^{-\theta(s+d-u)}\,dW_u.
$$

By D1, the final integral is Gaussian, centered, and independent of the past
at $s$. Its variance, including the factor $\sigma$, is

$$
\sigma^2\int_s^{s+d}e^{-2\theta(s+d-u)}\,du
=\frac{\sigma^2}{2\theta}(1-e^{-2\theta d}).
$$

Consequently,

$$
\boxed{X_{s+d}\mid X_s=x\sim
N\!\left(\mu+(x-\mu)e^{-\theta d},
\frac{\sigma^2}{2\theta}(1-e^{-2\theta d})\right).}
$$

The conditional law given the full past depends only on $X_s$, establishing
the Markov property needed for the likelihood factorization. For an exogenous
random grid, fix its realization first and apply this deterministic-time
argument.

For small positive $\theta d$, compute
$1-e^{-2\theta d}$ as `-expm1(-2 * theta * d)` to avoid subtracting
two nearly equal floating-point numbers. This changes only the numerical
evaluation, not the transition law.

## D4 — Stationarity, autocorrelation, and irregular AR(1)

Write $V=\sigma^2/(2\theta)$. If $X_0=x_0$ is fixed, D2 gives

$$
\mathbb E[X_t]=\mu+(x_0-\mu)e^{-\theta t},
\qquad
\operatorname{Var}(X_t)=V(1-e^{-2\theta t}).
$$

Thus the distribution tends to $N(\mu,V)$. More generally, the decaying
initial-value term in D2 tends to zero in probability for any almost surely
finite initial value, giving the same limiting distribution.

If instead $X_0\sim N(\mu,V)$, independently of future Brownian increments,
then D3 gives

$$
\mathbb E[X_t]=\mu,
\qquad
\operatorname{Var}(X_t)=e^{-2\theta t}V+V(1-e^{-2\theta t})=V.
$$

For $t\ge s$, write the D3 innovation as $\eta$, independent of $X_s$:

$$
\operatorname{Cov}(X_s,X_t)
=\operatorname{Cov}(X_s,e^{-\theta(t-s)}X_s+\eta)
=V e^{-\theta(t-s)}.
$$

The process is jointly Gaussian, has constant mean, and has covariance depending
only on time separation. It is therefore strictly stationary under this
initialization, with

$$
\boxed{X_t\sim N\!\left(\mu,\frac{\sigma^2}{2\theta}\right),
\qquad \rho(h)=e^{-\theta|h|}.}
$$

The autocorrelation formula assumes stationarity; for a fixed nonstationary
initial value, the covariance instead uses the time-dependent
$\operatorname{Var}(X_s)$.

Now observe $x_i=X_{t_i}$ at $t_0<\cdots<t_n$, with
$d_i=t_i-t_{i-1}>0$. There are $n$ transitions and $n+1$ observations.
Conditional on this grid,

$$
\boxed{x_i=\mu+(x_{i-1}-\mu)\phi_i+\varepsilon_i,\qquad
\phi_i=e^{-\theta d_i},\qquad
\varepsilon_i\sim N(0,v_i),\qquad
v_i=V(1-\phi_i^2).}
$$

The innovations are independent across disjoint intervals and independent of
the preceding states. This is a Gaussian AR(1) recursion with gap-dependent
coefficients and innovation variances. Exact simulation uses
$\varepsilon_i=\sqrt{v_i}z_i$ with independent $z_i\sim N(0,1)$.
It has no Euler discretization error at the observation times.

For constant gaps $d_i=d$, this becomes the usual AR(1) with
$\phi=e^{-\theta d}\in(0,1)$, intercept $a=\mu(1-\phi)$, and
innovation variance $v=V(1-\phi^2)$.

## D5 — Exact, PFML, and Euler conditional likelihoods

Let $\beta=(\theta,\mu,\sigma)$. By the Markov property and chain rule,

$$
L_c(\beta;x\mid x_0,t)
=\prod_{i=1}^n p_\beta(x_i\mid x_{i-1},d_i).
$$

This is the likelihood **conditional on the observed initial value** $x_0$
and on the observation grid. A parameter-free exogenous grid has no additional
factor relevant to maximizing over $\beta$.

### Exact conditional MLE

Define

$$
m_i(\beta)=\mu+(x_{i-1}-\mu)e^{-\theta d_i},
\qquad
v_i(\beta)=\frac{\sigma^2}{2\theta}(1-e^{-2\theta d_i}).
$$

The negative conditional log-likelihood is

$$
\boxed{Q_{\rm exact}(\beta)=\frac12\sum_{i=1}^n
\left[\log(2\pi v_i)+\frac{(x_i-m_i)^2}{v_i}\right].}
$$

Minimize this over $\theta>0,\sigma>0,\mu\in\mathbb R$. The code uses
working coordinates $(\log\theta,\mu,\log\sigma)$, enforcing positivity
by transformation. Reparameterizing the unknown parameters adds no Jacobian
term to this likelihood: the observed data have not been transformed.

If the **stationary initial density** is included, the objective instead becomes

$$
Q_{\rm joint}(\beta)=Q_{\rm exact}(\beta)
+\frac12\left[\log(2\pi V)+\frac{(x_0-\mu)^2}{V}\right].
$$

These are different finite-sample estimators. The notebook simulates a
stationary initial observation but fits the **conditional** likelihood; this
is a deliberate conditioning choice, not an omitted implementation term.
An exact transition likelihood does not imply unbiased finite-sample parameter
estimates or automatically valid Wald confidence intervals.

### Equidistant AR(1) analytic cross-check

At constant gap $d$, regress $y_i=x_i$ on $z_i=x_{i-1}$ with an
intercept. If $\sum_i(z_i-\bar z)^2>0$, the unconstrained Gaussian AR(1)
conditional MLE is

$$
\widehat\phi=\frac{\sum_i(z_i-\bar z)(y_i-\bar y)}
{\sum_i(z_i-\bar z)^2},\qquad
\widehat a=\bar y-\widehat\phi\bar z,\qquad
\widehat v=\frac1n\sum_i(y_i-\widehat a-\widehat\phi z_i)^2.
$$

The divisor is $n$, not the degrees-of-freedom correction used for an
unbiased regression variance estimate. When $0<\widehat\phi<1$ and
$\widehat v>0$, the mapping is

$$
\boxed{\widehat\theta=-\frac{\log\widehat\phi}{d},\quad
\widehat\mu=\frac{\widehat a}{1-\widehat\phi},\quad
\widehat\sigma^2=\frac{2\widehat\theta\widehat v}
{1-\widehat\phi^2}.}
$$

This is an interior OU solution only under the stated conditions. If the OLS
slope lies outside $(0,1)$, blindly applying this map does not give an
admissible OU estimate; a constrained fit or boundary analysis is required.

### PFML: replace every observed gap by its sample mean

Set

$$
\bar d=\frac1n\sum_{i=1}^n d_i
=\frac{t_n-t_0}{n}.
$$

Use the same objective as $Q_{\rm exact}$, but replace every $d_i$
inside $m_i$ and $v_i$ by $\bar d$. This defines this project's
PFML estimator (called `naive` in some function names). It keeps the sequence
of observations and their mean spacing while discarding variation among gaps.
It coincides with the exact conditional likelihood on an equidistant grid.

For irregular gaps it is a misspecified Gaussian contrast. It is **not** the
likelihood obtained by integrating the exact density over a gap distribution;
that integration generally gives a mixture of Gaussian transitions.

An independent population check explains the direction of its distortion.
Assume stationary sampling with i.i.d. positive gaps $D_i$, independent of
the OU path, finite mean $\delta=\mathbb E[D]$, and the long-run averaging
conditions needed for the sample regression moments to converge. Then

$$
q=\frac{\operatorname{Cov}(X_{t_i},X_{t_{i-1}})}{V}
=\mathbb E[e^{-\theta D}],
$$

so the limiting regression has slope $q$, intercept $\mu(1-q)$,
and residual variance $V(1-q^2)$. Mapping these to OU parameters gives

$$
\theta_*=-\frac{\log q}{\delta},\qquad
\mu_*=\mu,\qquad
\sigma_*^2=2\theta_*V=\sigma^2\frac{\theta_*}{\theta}.
$$

Since $d\mapsto e^{-\theta d}$ is strictly convex, Jensen's inequality
gives $q\ge e^{-\theta\delta}$, hence $\theta_*\le\theta$, with
strict inequality for a nondegenerate gap distribution. This is a population
limit, not a claim that each individual fitted value is too small.

For unfloored Gamma gaps with shape $k=c^{-2}$, scale
$s=\delta c^2$, and coefficient of variation $c>0$, direct integration
of the Gamma density gives

$$
q=(1+\theta s)^{-k},\qquad
\boxed{\theta_*=
\frac{\log(1+\theta\delta c^2)}{\delta c^2},\qquad
\sigma_*=\sigma\sqrt{\theta_*/\theta}.}
$$

The continuous $c\to0$ limit is $(\theta,\sigma)$.
For exponential gaps, $c=1$. These Gamma formulas do not apply unchanged
after flooring or renormalizing realized gaps: use the actual transformed gap
law and its mean. Likewise, a finite simulated sample need not have exactly
the configured population mean gap.

### Euler Gaussian contrast

For $\theta d_i$ small, expand the exact mean and variance:

$$
e^{-\theta d_i}=1-\theta d_i+O((\theta d_i)^2),
\qquad
v_i=\sigma^2d_i-\sigma^2\theta d_i^2+O(d_i^3)
$$

for fixed parameters as $d_i\to0$. Retaining the leading terms gives

$$
m_i^{\rm Euler}=x_{i-1}+\theta(\mu-x_{i-1})d_i,
\qquad
v_i^{\rm Euler}=\sigma^2d_i,
$$

and the contrast

$$
\boxed{Q_{\rm Euler}(\beta)=\frac12\sum_{i=1}^n
\left[\log(2\pi\sigma^2d_i)
+\frac{(x_i-x_{i-1}-\theta(\mu-x_{i-1})d_i)^2}{\sigma^2d_i}\right].}
$$

It uses the true observed gaps but approximates each transition. Large gaps
can violate the local small-step approximation even when the mean gap is
small. Data in this project are generated with D4's exact recursion, so Euler
estimation error is not confounded with Euler simulation error.

For stationary data at constant spacing $d$, another population regression
check gives

$$
\boxed{\theta_*^{\rm Euler}=\frac{1-e^{-\theta d}}{d},\qquad
\mu_*^{\rm Euler}=\mu,\qquad
(\sigma_*^{\rm Euler})^2
=\frac{\sigma^2}{2\theta d}(1-e^{-2\theta d}).}
$$

Indeed, matching the regression slope gives
$1-\theta_*^{\rm Euler}d=e^{-\theta d}$, while the residual variance is
the exact innovation variance and must equal $(\sigma_*^{\rm Euler})^2d$.
Both limits approach the true values as $d\to0$.

## Optional 30-minute blank-paper self-check

This exercise is supplied for later personal verification. **No attempt or
pass is recorded by the existence of this document.** The unavailable companion
document may specify additional or differently numbered tasks.

Without consulting the formulas above, use only Brownian motion's properties,
Gaussian limit/closure facts, and the deterministic-kernel integral interchange
to do the following:

1. **0–5 minutes:** Construct the deterministic Wiener integral from step sums;
   derive its mean, variance, and independence from the preceding filtration.
2. **5–11 minutes:** Derive the integrating-factor OU solution and verify it in
   the integral equation.
3. **11–16 minutes:** Derive the transition mean, variance, and Gaussian law.
4. **16–22 minutes:** Derive the stationary variance, autocorrelation, and the
   AR(1) recursion for unequal gaps. State the initialization and sampling
   assumptions.
5. **22–30 minutes:** Write the exact conditional negative log-likelihood,
   explain the stationary initial-density alternative, and construct PFML and
   Euler objectives. State why exact transitions alone do not guarantee
   unbiased estimates or calibrated confidence intervals.

Check the work against D1–D5 afterward and record any gaps honestly. The
population-limit extensions in D5 and numerical optimization details are useful
follow-up work rather than requirements of this suggested timed exercise.
