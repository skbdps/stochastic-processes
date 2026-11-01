# Parameter estimation for the Ornstein–Uhlenbeck process under irregular sampling

Suryansh Kumar

**Draft 2 · 1 November 2026 · Week 10 methods and explanatory supplement**

This revision preserves the original methods baseline and adds a quantitative
PFML mean mechanism, a saved-span sensitivity analysis, and an independent
fixed-horizon supplement. The original E1–E4 research release remains
`4bddb527f3c8bb88bc1ebec28c8818e03dc251aa`; supplemental evidence is separately
bound under `artifacts/week10_supplement`. Section 4 below covers only this
focused supplement. The integrated core Results, practical guidance, complete
related work, general discussion and final abstract remain Week 11–12 work.
This is a working draft for author review, not a submitted or academically
reviewed manuscript.

## 1 Introduction

The Ornstein–Uhlenbeck (OU) process is a basic continuous-time model for a
quantity that fluctuates around a long-run level. Its Gaussian transition
distribution makes it possible to distinguish errors due to an estimation
approximation from errors due to simulating the process itself. The model also
has an established financial interpretation: the mean-reverting Gaussian
short-rate specification in [Vasicek (1977)](#references) is an OU process.
Observation times in applications need not be equally spaced. Estimation must
then account for both the duration of each transition and the information
available over the observed time span.

Three estimation procedures illustrate different uses of these observations. Exact
conditional likelihood uses the duration of every transition. A mean-gap
pretend-fixed maximum likelihood (PFML) procedure replaces those durations by
their observed arithmetic mean. An Euler Gaussian procedure uses the actual
durations but approximates the conditional mean and variance. These choices
change the statistical objective, and an interval derived from the curvature
of a working objective need not describe its estimator's repeated-sample
uncertainty accurately.

Random sampling and its interaction with estimation are established subjects.
[Aït-Sahalia and Mykland (2003)](#references) distinguish likelihood procedures
that use different amounts of timing information, analyze their asymptotic
properties and include finite-sample simulations. The PFML terminology used
here follows that work, but the present procedure uses each dataset's realized
mean gap and jointly estimates the long-run mean. Their theoretical PFML
comparison uses the true unconditional mean gap, so the finite-sample and
nuisance-parameter conventions are not identical. Their
[2004 paper](#references) develops a general
estimating-equation treatment of randomly spaced diffusion observations. The
present experiments do not implement their integrated-out likelihood or their
formal information decomposition. Their OU specializations and nuisance
parameter conventions also differ from the joint three-parameter fit studied
here. The IID-gap assumption in the 2004 theory does not supply a
dependent-gap coverage benchmark automatically.

This study asks how observation count, nominal sampling scale and gap
variability jointly affect finite-sample estimation and interval reliability.
It compares exact, realized-mean-gap PFML and actual-gap Euler procedures on
the same simulated paths within each design cell. Bias, root mean squared
error (RMSE), interval availability and coverage are reported with Monte Carlo
uncertainty. Regular-spacing controls match the nominal design; an exploratory
extension compares clustered gaps with independent gaps having the same entire
one-gap marginal distribution. This distinction allows dependence to be
examined without changing that marginal, while retaining consequences for
realized observation spans.

The contribution is a reproducible comparison across a fixed design grid,
supported by transparent accounting for numerical outcomes and mathematical
interpretation of the working objectives. Population targets help explain
approximation bias; empirical variability and reported standard errors help
describe interval calibration. The study does not require any estimator to
dominate, irregularity to harm every parameter, or a sharp performance boundary
to emerge. Its practical classifications apply to the tested designs and
specified interval procedures. The scope is a stationary, noise-free OU process
with exogenous observations; measurement error, endogenous sampling and
real-data validation are separate questions.

## 2 Model, estimators and mathematical benchmarks

### 2.1 Model and observation conventions

Let the data-generating process satisfy

$$
dX_t=\theta_0(\mu_0-X_t)\,dt+\sigma_0\,dW_t,
\qquad \theta_0>0,\quad \sigma_0>0,
$$

where $W_t$ is standard Brownian motion. The stationary variance is
$V=\sigma_0^2/(2\theta_0)$, and simulation begins with
$X_0\sim N(\mu_0,V)$. The observation schedule is independent of the initial
state and the entire Brownian innovation stream. Write $x_i=X_{t_i}$ for
$i=0,\ldots,K$, with $t_0=0$, $n=K+1$ observations and positive observed gaps
$\Delta_i=t_i-t_{i-1}$.

For candidate parameters $\psi=(\theta,\mu,\sigma)$, the exact conditional
transition mean and variance are

$$
m_i(\psi)=\mu+(x_{i-1}-\mu)e^{-\theta\Delta_i},
\qquad
v_i(\psi)=\frac{\sigma^2}{2\theta}
             (1-e^{-2\theta\Delta_i}).
$$

Each simulated path uses this exact transition at the stored observation times
and independent standard-normal innovations. Euler discretization is never
used to generate the data. All fits condition on $x_0$ and the timestamps;
the stationary density of $x_0$ is excluded from the fitted likelihood.
Stationary initialization and a conditional fitted likelihood are therefore
distinct conventions.

Three time scales must be distinguished. Let $m$ be the nominal design gap,
$\widehat m=K^{-1}\sum_{i=1}^K\Delta_i$ its realized sample mean, and
$\bar m=E[Y]$ the population mean of the gap law actually used after any floor.
The design scale is $q=\theta_0m$. Nominal span is $Km$, while realized span
is $T=\sum_i\Delta_i$; their dimensionless versions are $Kq$ and $\theta_0T$.
The nominal and realized spans need not be equal. For example,
$n=250,q=0.1$ gives nominal dimensionless span $24.9$, because there are
$249$ transitions.

### 2.2 Exact, PFML and Euler objectives

The exact negative conditional log-likelihood is

$$
\mathcal L_X(\psi)
=\frac12\sum_{i=1}^K
\left\{\log(2\pi v_i(\psi))+
\frac{(x_i-m_i(\psi))^2}{v_i(\psi)}\right\}.
$$

The three parameters are jointly estimated. For irregular observations the
implementation profiles the drift intercept $\theta\mu$ and diffusion
variance at fixed $\theta$, then searches over $\log\theta$. The search
examines expanded endpoints, limiting objectives and sampled local minima.
It records unresolved boundary behavior rather than certifying a finite
interior optimum from solver termination alone. This numerical search is not
a proof of global optimality for every path. Regular exact fitting uses the
conditional AR(1) regression solution when it is interior.

PFML replaces every $\Delta_i$ in the exact moments by $\widehat m$. Its
conditional regression has slope $\rho=e^{-\theta\widehat m}$ and intercept
$a_P=\mu(1-\rho)$. For an admissible ordinary least-squares solution
$0<\widehat\rho<1$, let $\widehat a_P$ be the fitted intercept and
$\widehat v_P=K^{-1}\sum_i(x_i-\widehat a_P-
\widehat\rho x_{i-1})^2$ the residual variance. Then

$$
\widehat\theta_P=-\frac{\log\widehat\rho}{\widehat m},\qquad
\widehat\mu_P=\frac{\widehat a_P}{1-\widehat\rho},\qquad
\widehat\sigma_P^2=
\frac{2\widehat\theta_P\widehat v_P}{1-\widehat\rho^2}.
$$

Boundary and degenerate regressions remain explicit outcomes. For a fixed
observation vector, the mean-gap factor does not directly enter the displayed
formula for $\widehat\mu_P$. Sampling nevertheless changes its repeated-sample
behavior through the joint distribution of the observations.

Euler uses the Gaussian working moments

$$
m_{i,E}=x_{i-1}+\theta(\mu-x_{i-1})\Delta_i,
\qquad v_{i,E}=\sigma^2\Delta_i,
$$

in the same Gaussian negative-log-likelihood form. Weighted least squares
regresses $(x_i-x_{i-1})/\sqrt{\Delta_i}$ on
$\sqrt{\Delta_i}$ and $-x_{i-1}\sqrt{\Delta_i}$, yielding coefficients
$\theta\mu$ and $\theta$. The residual sum of squares divided by $K$ gives
the fitted diffusion variance. A nonpositive fitted drift rate is a boundary
outcome. A negative transition coefficient $1-\theta\Delta_i$ by itself
does not invalidate this Euler objective. Exact simulation, exact likelihood,
PFML and Euler are thus separated at the level of their statistical definitions.

### 2.3 Numerical coordinates and interval construction

Optimization and differentiation use observed-data normalization
$y=(x-a)/b$, $u=t/\widehat m$, where $a$ is the observed state mean and $b>0$
its centered root mean square. No simulation truth enters these scales. The
normalized working coordinates are
$\eta_u=(\log\theta_u,\mu_u,\log\sigma_u)$, with maps

$$
\theta_x=\theta_u/\widehat m,\qquad
\mu_x=a+b\mu_u,\qquad
\sigma_x=b\sigma_u/\sqrt{\widehat m}.
$$

For a fixed dataset these normalization constants are held fixed during
differentiation. Let $H$ be the full $3\times3$ Hessian of the summed
conditional negative log-likelihood in $\eta_u$, evaluated at the fitted
point. It is neither the scalar profile Hessian nor an optimizer's approximate
inverse Hessian. When all interval-validity checks pass, the model-based
natural-scale covariance is

$$
\widehat C_x=JH^{-1}J^\mathsf T,
\qquad J=\operatorname{diag}(\widehat\theta_x,b,\widehat\sigma_x).
$$

Because $H$ differentiates the summed objective, there is no additional
division by $K$. For each parameter $\psi_j$, the primary interval is
$\widehat\psi_j\pm z_{0.975}\sqrt{\widehat C_{x,jj}}$, with
$z_{0.975}=1.959963984540054$. Negative lower endpoints for positive parameters
are retained. Log-scale, profile, sandwich and bias-corrected intervals are
different procedures and are not substituted.

The registered numerical checks use central-difference steps
$h_j=0.001\max(1,|\eta_{u,j}|)$ and compare Hessians at $h$ and $h/2$.
The half-step Hessian must be finite, symmetric positive definite, have
condition number at most $10^{10}$, and pass relative Frobenius discrepancy
$0.005$ and whitened discrepancy $0.01$ gates. With
$H_{h/2}=BB^\mathsf T$, the latter bounds the largest absolute eigenvalue of
$B^{-1}(H_h-H_{h/2})B^{-\mathsf T}$. The standardized gradient measure
$\sqrt{g^\mathsf T H_{h/2}^{-1}g}$, with $g$ the half-step numerical gradient,
must be at most $0.01$.
An invalid point, rejected Hessian or invalid mapped variance withholds the
interval; failed matrices are not regularized to obtain an interval.

These gates assess the implemented numerical procedure. They do not establish
statistical calibration. In particular, PFML and Euler retain inverse-Hessian
working covariance under possible misspecification. A covariance calculation
appropriate to misspecification and dependence would additionally require
score variability, including serial covariance and relevant estimated nuisance
quantities.

### 2.4 Population targets and empirical calibration

For the used positive gap $Y$, define

$$
L(s)=E[e^{-sY}],\qquad
\bar m=E[Y],\qquad
H_0=E\left[\frac{1-e^{-\theta_0Y}}{Y}\right]
    =\int_0^{\theta_0}L(s)\,ds.
$$

Stationarity and exogeneity give adjacent centered-state covariance
$VL(\theta_0)$. The expected PFML objective with mean gap $\bar m$ has targets

$$
\theta_{P,*}=-\frac{\log L(\theta_0)}{\bar m},\qquad
\mu_{P,*}=\mu_0,\qquad
\sigma_{P,*}^2=2V\theta_{P,*}.
$$

The corresponding expected actual-gap Euler objective has targets

$$
\theta_{E,*}=\frac{1-L(\theta_0)}{\bar m},\qquad
\mu_{E,*}=\mu_0,\qquad
\sigma_{E,*}^2=V\left\{2H_0-
\frac{(1-L(\theta_0))^2}{\bar m}\right\}.
$$

These are expected-objective minimizers. Identifying them with estimator
probability limits also requires appropriate convergence and ergodicity
conditions, including convergence of the realized mean gap. They are not
finite-sample mean, variance or coverage predictions. Simulation errors and
coverage continue to be assessed against the true OU parameters.

For an unfloored Gamma gap with mean $m$ and CV $c>0$,
$L(s)=(1+smc^2)^{-1/c^2}$. A floor $f$ changes both the transform and the
mean. For a continuous raw gap $D$ with density $p_D$ and CDF $F_D$,

$$
L_Y(s)=e^{-sf}F_D(f)+\int_f^\infty e^{-sd}p_D(d)\,dd,
\qquad
E[Y]=fF_D(f)+E[D\,\mathbf1\{D\ge f\}].
$$

The floor creates an atom at $f$; it is not rejection truncation. The
random-gap Euler diffusion target includes $H_0$ and cannot generally be
replaced by the constant-gap expression $V(1-L(\theta_0)^2)/\bar m$.
Matched full gap marginals share the displayed expected-objective targets
under the stated stationarity and exogeneity assumptions. Serial dependence
can still change variance, realized-span distributions and coverage.

Empirical calibration compares the sample standard deviation of estimates
with the root mean square of reported standard errors on the same
interval-valid subset. The mean reported standard error is also retained
separately. This descriptive comparison neither estimates operational coverage
nor supplies a dependence-robust covariance formula. Population targets and
empirical calibration are the adopted explanatory benchmarks; a demonstrated
departure from fully matched first-order coverage theory is not established.

### 2.5 Mean precision under clustered observation gaps

For an interior PFML regression solution, let $K=N-1$. The intercept normal
equation yields the exact identity

$$
\widehat\mu_P=\frac{\sum_{i=1}^{K}X_i-\widehat\rho\sum_{i=0}^{K-1}X_i}
 {K(1-\widehat\rho)}
=\overline X_{0:K-1}+\frac{X_K-X_0}{K(1-\widehat\rho)}.
$$

Thus the estimator is an event-index average plus an endpoint correction.
Interior observations retain equal weights when many occur within a short
burst. Since nearby OU observations are highly correlated, their count does
not measure independent information about the long-run mean. Exogenous
sampling preserves each observation's stationary mean; this weighting
mechanism affects precision without requiring a mean-bias shift.

Holding a coefficient $\rho\in(0,1)$ fixed gives a linear diagnostic with
$w_0=-\rho/[K(1-\rho)]$, $w_i=1/K$ for $1\le i<K$, and
$w_K=1/[K(1-\rho)]$. Conditional on an exogenous schedule, its variance is exactly

$$
\operatorname{Var}\{\widetilde\mu(\rho)\mid t\}
=V\sum_{i,j=0}^{K}w_iw_j e^{-\theta_0|t_i-t_j|},
\qquad V=\frac{\sigma_0^2}{2\theta_0}.
$$

This expression integrates over stationary initial values, consistently with
the repeated-path experiment, although each fitted likelihood conditions on
its initial observation. An OU innovations recursion evaluates it in linear
time; dense covariance calculations provide an independent numerical check.

For the event-stationary parent schedule, define $A=\Lambda-Q$, floor $f$,
and the floor-aware joint Laplace kernel

$$
B_s=e^{-sf}\left[A^{-1}(I-e^{-Af})\Lambda
 +e^{-Af}(A+sI)^{-1}\Lambda\right].
$$

For $f=0$, this reduces to $(A+sI)^{-1}\Lambda$. Put
$\phi=\alpha B_{\theta_0}\mathbf1$. The event-index OU correlations are
$c_h=\alpha B_{\theta_0}^{h}\mathbf1$ under MMPP gaps and $c_h=\phi^h$
under independently reset parent gaps. Consequently the ordinary observation
mean has exact fixed-count variance

$$
\operatorname{Var}(\overline X)=\frac{V}{N^2}
 \left[N+2\sum_{h=1}^{N-1}(N-h)c_h\right].
$$

The weighted diagnostic uses the corresponding finite sum of its weight
cross-products. Matched one-gap marginals make $c_1$ identical but do not
make higher-lag correlations identical. The population AR(1) working variance
$V(1+\phi)/[K(1-\phi)]$ therefore need not include the actual multi-gap covariance.
These specialized calculations use standard MAP/phase-type identities [4].
They do not extend the IID assumptions in [3] without additional limit theory.

At the reference coefficient $\phi$, the actual fitted mean differs by

$$
\widehat\mu_P-\widetilde\mu(\phi)=
\frac{(X_K-X_0)(\widehat\rho-\phi)}
 {K(1-\widehat\rho)(1-\phi)}.
$$

This gives a direct empirical approximation check. The exact covariance
statements concern the linear diagnostic, while predictions for jointly fitted
PFML and studentized Wald coverage remain approximations. For random-count
fixed-horizon schedules we evaluate the conditional expression on each saved
schedule; we do not apply a fixed-count correlation formula to a stopped
sample. The parent-law $\phi$ is a reference coefficient, not an asserted
finite-window pseudo-target. Full derivations and limits appear in the
[mechanism note](../docs/week10_mechanism.md).

## 3 Experimental methods

### 3.1 Registered core designs

All studies use $(\theta_0,\mu_0,\sigma_0)=(1,0,0.5)$. E1 and E2 follow the
[27 September 2026 preregistration](../preregistration.md), adopted after the
disclosed exploratory development and before the headline observations. Each
cell has $R=2000$ independent source paths, with all three estimators fitted
to each path.

| Study | Observation design | Cells | Root seed |
| --- | --- | ---: | ---: |
| E1 | $n=1000$, $q=0.5$, $c\in\{0,0.25,0.5,1,1.5,2,3\}$ | 7 | 20261001 |
| E2 | $n\in\{250,1000,4000\}$, $q\in\{0.1,0.5,1,2\}$, $c\in\{0,0.5,1,2,3\}$ | 60 | 20261002 |

At $c=0$, gaps equal $m=q/\theta_0$. At $c>0$, independent raw gaps have
Gamma shape $c^{-2}$ and scale $mc^2$. When nominal $c\ge1.5$, each gap is
replaced by $\max(D,10^{-6}m)$ without subsequent renormalization. Raw- and
used-gap diagnostics, clipping fractions, realized moments and spans are
retained with the seeds and hashes needed to reconstruct the sources; full
gap and observation arrays are not stored in the fit partitions.
The regular reference matches $n$, truth and nominal mean gap, not realized
or exactly expected span. These comparisons consequently include changes in
span distributions and the specified floor transformation.

E1 contains 14,000 paths and 42,000 fits; E2 contains 120,000 paths and
360,000 fits. Overlapping nominal cells in the two studies have independent
seeds and remain separate analyses. No selection between them is used to
obtain a preferred classification.

### 3.2 Exploratory diagnostics and floor sensitivity

E3 follows a [dated exploratory protocol](../docs/week7_exploratory_protocol.md)
and does not replace the registered results. Its fixed selected-path
investigation examines 53 existing sources and replays 159 original fits.
The selection is diagnostic rather than a random sample for estimating the
frequency of numerical problems.

An independent follow-up uses $n=250,q=2,c\in\{0,0.5\}$ with 10,000 new paths
per cell and root seed 20261007. Its 20,000 paths and 60,000 fits retain the
original estimation and interval procedures. A separate floor study uses
$n=1000,q=0.5,c\in\{0,1.5,2,3\}$ and
$\epsilon\in\{10^{-8},10^{-6},10^{-4}\}$, with 2000 paths per variant and
root seed 20261009. Floors apply only at $c\ge1.5$.

Within each CV and replication, the three floor variants share raw Gamma
draws, the stationary initial draw and standardized OU innovations. They
have 24,000 path realizations and 72,000 fits but only 8000 independent
coupling groups. Contrasts to $\epsilon=10^{-6}$ use this pairing. The regular
variants are the same shared negative control. All headline, diagnostic and follow-up
records retain their original study identities.

### 3.3 Clustered sampling and independent full-marginal controls

E4 follows its [prospective exploratory protocol](../docs/week8_exploratory_protocol.md).
It uses $n=1000,q=0.5,R=2000$ and root seed 20261011. Three cells use a
two-state Markov-modulated Poisson process (MMPP) with intensity ratios
$r\in\{1,5,20\}$, three use independent gaps matching the respective entire
event-stationary marginal laws, and one uses regular gaps. The seven cells
contain 14,000 independent paths and 42,000 fits.

In state order calm, burst, the phase generator and stationary time
distribution are

$$
Q=m^{-1}\begin{pmatrix}-1/10&1/10\\1/2&-1/2\end{pmatrix},
\qquad \pi=(5/6,1/6).
$$

The arrival rates are $\lambda_c=6/((r+5)m)$ and
$\lambda_b=6r/((r+5)m)$, so the raw stationary time-average arrival rate is
$1/m$. With $\Lambda=\operatorname{diag}(\lambda_c,\lambda_b)$, time zero
is an event and its phase distribution is

$$
\alpha=\frac{\pi\Lambda}{\pi\Lambda\mathbf1}
       =\left(\frac5{r+5},\frac r{r+5}\right).
$$

MMPP simulation retains the arrival phase between gaps. The independent
control instead redraws a phase from $\alpha$ before every gap, while
retaining the same competing arrival/switch mechanism within each gap.
Both have raw survival function
$\Pr(D>t)=\alpha e^{(Q-\Lambda)t}\mathbf1$; their serial laws differ.
At $r=1$, both mechanisms give exponential gaps and are independent samples
from that common law, rather than regular observations.

The derived raw CV satisfies
$\mathrm{CV}^2=1+50(r-1)^2/[3(r^2+70r+25)]$.
Only $r=20$ crosses the floor trigger $\mathrm{CV}\ge1.5$. Both arms then
use $Y=\max(D,10^{-6}m)$, with no renormalization or hidden-clock advancement
through the inserted padding. The dependent arm is consequently described
as floored MMPP-derived gaps. The transformed full marginal remains matched,
while its mean need not equal $m$ exactly.

The primary contrasts compare MMPP with its corresponding independent
full-marginal control at each $r$. Cells are independent even when replication
labels coincide. This design assesses gap dependence including consequences
for realized span; it does not hold span fixed. Comparisons along the
intensity-ratio curve or against the regular reference change more than
dependence and are interpreted separately.

### 3.4 Metrics, availability and Monte Carlo uncertainty

Every planned source/estimator record is retained. Solver termination, point
validity, interval availability and pipeline errors are distinct. For one
parameter, let $S$ be the point-valid count, $N_I$ the interval-valid count,
and $C$ the number of valid intervals containing truth. Point errors
$e_j=\widehat\psi_j-\psi_0$ give

$$
\widehat{\mathrm{bias}}=S^{-1}\sum_{j=1}^{S}e_j,\qquad
\widehat{\mathrm{RMSE}}=\sqrt{S^{-1}\sum_{j=1}^{S}e_j^2},
\qquad
\widehat p_{\mathrm{op}}=C/R,\quad
\widehat p_{\mathrm{cond}}=C/N_I.
$$

Bias and RMSE condition on valid point fits. Operational coverage counts
unavailable intervals as noncoverage among all planned paths; conditional
coverage uses available intervals and is undefined when $N_I=0$. Point
availability $S/R$ and interval availability $N_I/R$ accompany these summaries.
Relative bias, RMSE and their uncertainty divide by
$|\psi_0|$ only for nonzero truth. Thus relative-to-truth metrics for $\mu_0=0$
are undefined, although an RMSE ratio to a positive control RMSE is defined.

On the point-valid subset, Monte Carlo standard errors are
$s_e/\sqrt S$ for bias and
$s_{e^2}/(2\sqrt S\,\widehat{\mathrm{RMSE}})$ for RMSE, where $s$ denotes
sample standard deviation. They require at least two valid estimates; the
RMSE approximation also requires positive RMSE. Coverage uses plug-in
binomial Monte Carlo standard errors and 95% Wilson intervals. A zero plug-in
standard error at observed coverage zero or one does not remove Wilson
uncertainty.

Within-path estimator contrasts use common-valid point subsets and paired
uncertainty; operational-coverage differences use paired indicators over all
planned paths. Independent-cell contrasts use independent-sample uncertainty
and retain each cell's own valid counts. Positive RMSE ratios use log-delta
intervals. For independent cells A and B, if their RMSE estimates are $r_A,r_B$
with Monte Carlo standard errors $u_A,u_B$, the log-ratio standard error is
$\sqrt{(u_A/r_A)^2+(u_B/r_B)^2}$; exponentiating its normal interval gives
a positive ratio interval. Coverage differences use the appropriate paired
or independent variance, rather than subtraction of Wilson endpoints.
Contrasts sharing a control are not mutually independent. All research
uncertainty statements are pointwise, without simultaneous error control.

### 3.5 Practical classification and research discipline

For $\theta$ and $\sigma$, the registered practical criteria are operational
coverage at least 80% and absolute relative bias at most 10%. A supported
breach requires the 95% Wilson upper coverage bound below 0.80, or the entire
approximate 95% relative-bias interval above 0.10 or below -0.10. A within-criteria
classification requires the Wilson lower bound at least 0.80 and the whole
bias interval inside $[-0.10,0.10]$. Other cases remain unresolved. The
overall practical classification is inapplicable to $\mu_0=0$ because its
relative-bias component is undefined.

These are broad practical criteria, not a certificate of nominal 95%
calibration. Classifications use unrounded bounds and the complete discrete
grid. They do not establish performance at untested coordinates or justify a
monotone maximum acceptable CV. Thresholds, seeds, replications, eligibility
rules and headline grids are not adjusted to obtain desired outcomes.
Exploratory analyses retain that designation, and all original failures and
unresolved cases remain available.

### 3.6 Validation and reproducibility

Scientific identities, source streams, configurations, raw outcomes and
execution environments accompany each archived study. E1/E2 and the
independent E3 follow-up address separate gap and OU streams through
`SeedSequence(root_seed, spawn_key=(4, *eight_cell_id_words, replication)).spawn(2)`.
E4 uses namespace 8 with its own full scientific cell identities. The paired
floor experiment uses its declared anchor identity to preserve coupling.
Execution order and worker identity do not determine the source streams.

The registered regular-grid checkpoint and selected bootstrap diagnostics
precede the core experiments. The bootstrap checks use two fixed source paths
and three estimators, with 500 draws per source/estimator pair. They are six
selected-path diagnostics, not a repeated-dataset bootstrap coverage study.
The published audit reports state their actual source-reconstruction and
fit-replay scope; selected endpoint-fit replay is not described as replaying
every fit in a full study.

Publication figures and tables are generated from saved results. The
[Week 9 package](../artifacts/week9_publication/README.md) binds its inputs,
claims, plotting data and recorded visual review. Its
[reproduction instructions](../docs/week9_reproduction.md) distinguish
saved-data rendering from original numerical simulation and replay. The
recorded macOS rendering environment differs from the original Linux
scientific environment; exact byte replay across different CPU dispatch,
precision or graphics environments is not assumed. The archived package,
failure dispositions and limitations remain unchanged by this manuscript.

### 3.7 Post-result span and mechanism supplement

After observing E4, we specified a separate exploratory Week 10 supplement.
The existing E4 archive and all original estimator definitions were preserved.
Saved E4 results were stratified into fixed span bins with boundaries
$0,400,450,475,500,525,550,600,\infty$. Common support required at least 30
source paths per arm in a bin. Both arms were standardized to pooled source
counts on that support. MSE was averaged before taking its square root.
Operational coverage retained unavailable intervals as noncoverage; point
metrics used valid estimates within each stratum. MC intervals condition on
observed support, counts and weights. All bins and excluded mass are reported;
this descriptive procedure does not identify a causal direct effect.

The independent fixed-horizon experiment retains the E4 parent gap laws,
truth, nominal mean gap and floor rule, with seven cells and 2,000 paths per
cell. Every source is observed at zero and $T=499.5$, with parent arrivals
strictly before $T$ retained. The origin uses event-Palm initialization.
The first crossing parent arrival is recorded for audit; a deterministic
observation at $T$ closes the window. Its terminal interval is not a new
parent draw and is not floored. Thus observed span is exactly fixed while
$N$ is random. Matching the parent complete-gap law does not assert identical
finite-window retained-gap distributions after stopping and endpoint insertion.

The full seed is 20261101 in namespace 10; development and source-only seeds
are 21261101 and 22261101. Streams are independent across cells and shared
across estimators within a source. All original numerical fitting, Hessian and
95% Wald rules remain unchanged. No coefficient from the generating truth is
used in estimation: truth enters only the external analytical benchmark.
The source, protocol and runtime were hash-bound before full execution.
Every full source was reconstructed, PFML means independently checked by OLS,
and all three endpoint fits at replications 0 and 1999 in each cell replayed.
The [supplement protocol](../docs/week10_supplement_protocol.md) and
[reproduction instructions](../docs/week10_reproduction.md) specify the complete
boundaries, provenance and validation scope.

## 4 Focused results: span controls and the PFML mean mechanism

This is a post-result exploratory extension of E4, authorized after the Week 10
methods baseline. The original fixed-N study is preserved. Two additional
analyses address span variation: descriptive reweighting of the saved original
records, and an independent fixed-horizon experiment. A separate covariance
calculation explains the mean-estimation mechanism without fitting its
parameters to the observed error curves.

### 4.1 Equal-span results

The new experiment contains seven cells, 2,000 independent sources per cell,
14,000 sources and 42,000 estimator fits. Every observed span is exactly
T=499.5. All points and intervals were available; there were zero pipeline
errors. Observation counts vary, and the same three estimators and original
full-Hessian Wald inference are used.

| r | Analysis | PFML mean RMSE ratio, MMPP/IID [95% MC interval] | Coverage, MMPP / IID | Coverage difference, percentage points [95% MC interval] |
| --- | --- | --- | --- | --- |
| 5 | Original fixed N | 1.0839 [1.0375, 1.1323] | 92.75% / 95.45% | -2.700 [-4.158, -1.242] |
| 5 | Span-standardized original | 1.0813 [1.0326, 1.1323] | 92.61% / 95.30% | -2.690 [-4.252, -1.129] |
| 5 | New fixed T | 1.0560 [1.0101, 1.1039] | 92.65% / 94.50% | -1.850 [-3.369, -0.331] |
| 20 | Original fixed N | 1.1579 [1.1073, 1.2108] | 90.70% / 94.55% | -3.850 [-5.466, -2.234] |
| 20 | Span-standardized original | 1.1450 [1.0925, 1.2001] | 90.98% / 94.34% | -3.368 [-5.048, -1.687] |
| 20 | New fixed T | 1.1564 [1.1063, 1.2087] | 90.40% / 94.40% | -4.000 [-5.638, -2.362] |

The fixed-horizon PFML mean RMSE and coverage contrasts at r=5 and r=20 have
pointwise MC intervals excluding equal performance. Thus variation in total
observed span is not necessary for the adverse PFML mean effect in these
tested designs. The r=5 point estimate is smaller than in original E4, but
these independent experiments differ in observation-count distribution and
endpoint convention. Their point estimates do not quantify a causal fraction
explained by span, and no formal attenuation claim follows from this table.

At r=1 the two parent processes both reduce to Poisson sampling. The new PFML
mean RMSE ratio is 1.0184 [0.9736,
1.0653], and its coverage difference is
-0.750 pp [-2.071,
0.571]. Both intervals include no
difference. This is a useful negative control, not an equivalence test.

Across all seven cells and all three parameters, exact-likelihood coverage
ranges from 93.80% to 96.20%.
The complete parameter tables retain PFML/Euler theta and sigma failures and
Euler mean overcoverage. Those parameter-specific findings must not be
collapsed into a claim that every estimator fails because of clustering.

### 4.2 What the saved-span diagnostic establishes

Bins and the common-support rule were specified before generating the new
analysis tables, after seeing original E4. The same pooled span-bin weights
were used for both arms. The excluded fractions and weighted mean spans are:

| r | MMPP excluded | IID excluded | Weighted MMPP span | Weighted IID span |
| --- | --- | --- | --- | --- |
| 1 | 0.25% | 0.20% | 499.853 | 499.981 |
| 5 | 9.85% | 1.65% | 499.820 | 499.798 |
| 20 | 5.50% | 0.30% | 498.471 | 498.547 |

PFML mean differences persist under this coarse standardization. These
intervals condition on the observed support, bin counts and weights; they
omit uncertainty from estimating those choices. Within-bin span variation
and differences between schedules remain. Unequal tail exclusions change the
target population. The analysis is therefore corroborating descriptive
evidence, not a matched-span causal experiment or a decomposition of the
original effect into percentages due to different mechanisms.

### 4.3 A quantitative covariance mechanism

At an interior PFML solution, with K=N-1,

$$
\widehat\mu_P=\overline X_{0:K-1}+
\frac{X_K-X_0}{K(1-\widehat\rho)}.
$$

Thus interior observations have equal event-index weights. Bursts produce
many observations of closely related OU states; those extra observations
provide less independent information about the long-run mean than their
count suggests. Exogenous sampling does not systematically select high or
low OU states, so this is a precision and interval-calibration mechanism,
without a demonstrated mean-bias shift.

For a coefficient rho held fixed, the corresponding linear estimator has
weights -rho/[K(1-rho)] at the first observation, 1/K in the interior, and
1/[K(1-rho)] at the last. Its schedule-conditional variance is exactly

$$
V\sum_{i,j} w_iw_j\exp(-\theta_0|t_i-t_j|),
\qquad V=\frac{\sigma_0^2}{2\theta_0}.
$$

In the original event-stationary fixed-count design, let B be the matrix
Laplace transform of the joint used-gap/ending-phase kernel at theta_0 and
phi=alpha B 1. The lag-h OU correlation is alpha B^h 1 for MMPP and phi^h
for IID gaps. The matched one-gap law makes lag one and the population
working variance agree; it does not equate higher-lag correlations. The
working AR(1) variance consequently misses additional serial covariance.
The floor-aware kernel and finite-count weighted sums are derived and tested
in the companion mathematical note.

| r | Design | Predicted mean SD ratio | Observed PFML RMSE ratio | Approximate MMPP coverage | Observed MMPP coverage |
| --- | --- | --- | --- | --- | --- |
| 5 | Original fixed N | 1.0776 | 1.0839 | 93.11% | 92.75% |
| 5 | New fixed T | 1.0757 | 1.0560 | 93.11% | 92.65% |
| 20 | Original fixed N | 1.1624 | 1.1579 | 90.82% | 90.70% |
| 20 | New fixed T | 1.1602 | 1.1564 | 90.89% | 90.40% |

The fixed-count prediction integrates the known generating law analytically.
For the new fixed-horizon study, the prediction averages exact conditional
variances over its saved schedules. The latter has MC uncertainty from those
schedules and uses the parent-law coefficient as a reference, without claiming
it is the exact finite-window pseudo-target. Both calculations are benchmarks
for a fixed-coefficient diagnostic; coverage additionally substitutes a
population working variance and a Gaussian approximation.

The diagnostic was also evaluated on every new OU path. Across the seven
cells, its observed MSE minus its average predicted variance is at most
0.94 estimated MC standard errors from zero. The RMS difference between
the fitted PFML mean and the fixed-coefficient diagnostic is at most
0.77% of the fitted mean's empirical SD. This directly
checks the approximation that links the covariance argument to the actual
jointly fitted estimator. It does not establish exact finite-sample coverage,
a general dependent-gap CLT, or a new robust confidence-interval method.

### 4.4 Observation-count and boundary scope

| Mechanism | r | Mean N | Observed N range |
| --- | --- | --- | --- |
| MMPP | 1 | 1000.783 | 895–1099 |
| MMPP | 5 | 1003.398 | 815–1205 |
| MMPP | 20 | 1005.499 | 660–1341 |
| IID | 1 | 1000.456 | 902–1115 |
| IID | 5 | 1000.232 | 869–1141 |
| IID | 20 | 1002.335 | 792–1250 |
| REGULAR | 0 | 1000.000 | 1000–1000 |

Both endpoints are observed. The stochastic parent process starts at an event
with its Palm phase distribution; it is not initialized at an arbitrary
calendar-time origin. Arrivals strictly before T are retained and X(T) is
appended. Accordingly, the r=1 Poisson design has expected N=2+T/m=1001,
whereas the regular design has exactly 1000 observations. Event-Palm MMPP
initialization can also affect finite-window counts near the origin.

The terminal observation interval is not a fresh parent interarrival and is
not subjected to the parent-gap floor. Stopping and endpoint insertion can
change the finite-window retained-gap distribution even when the underlying
complete-gap laws match. Random N and its distribution remain features of
this contrast. These results establish persistence at a common observed
span under the declared boundary rule, not a universal dependence-only
effect across all possible observation designs.

### 4.5 Research conclusion and remaining manuscript work

The supplemental evidence strengthens the original E4 interpretation:
clustered exogenous schedules can reduce PFML mean precision and lead its
working Wald intervals to undercover, even when every path spans the same
time interval. Higher-lag OU covariance provides a quantitative explanation
of the measured size. The matched one-gap population target for the mean
remains correct, while working uncertainty need not reflect actual spread.

All uncertainty statements are pointwise and the supplement is exploratory.
No simulation seed, estimator, interval formula or replication count was
changed to obtain these outcomes. Week 11 still needs the integrated E1–E4
Results, practical guidance, broader related work and full discussion, followed
by academic review. Week 12 remains revision, final layout and submission
preparation. This supplement does not certify journal novelty or acceptance.


## References

1. Vasicek, O. (1977). An equilibrium characterization of the term structure.
   *Journal of Financial Economics*, 5(2), 177–188.
   [Publisher DOI](https://doi.org/10.1016/0304-405X(77)90016-2).
2. Aït-Sahalia, Y., and Mykland, P. A. (2003). The effects of random and discrete
   sampling when estimating continuous-time diffusions. *Econometrica*,
   71(2), 483–549. [DOI](https://doi.org/10.1111/1468-0262.t01-1-00416);
   [author-hosted article](https://www.princeton.edu/~yacine/random.pdf).
3. Aït-Sahalia, Y., and Mykland, P. A. (2004). Estimators of diffusions with
   randomly spaced discrete observations: A general theory.
   *The Annals of Statistics*, 32(5), 2186–2222.
   [DOI](https://doi.org/10.1214/009053604000000427);
   [author-hosted article](https://www.princeton.edu/~yacine/gmmrandom.pdf).

4. Asanjarani, A., and Nazarathy, Y. (2022). Stationary Markovian arrival
   processes: Results and open problems. *ANZIAM Journal*.
   [DOI](https://doi.org/10.1017/S1446181122000013).

The [source map](week10_source_map.md) identifies the research protocols,
mathematical notes and implementation files behind these methods. The
[Week 9 literature comparison](../docs/week9_literature_comparison.md)
provides the broader comparison and its access limits for the later related-work
section.

Supplementary figures: [S1 span controls](../artifacts/week10_supplement/figures/S1_span_controls.png),
[S2 covariance mechanism](../artifacts/week10_supplement/figures/S2_mean_mechanism.png),
and [S3 full parameter/count context](../artifacts/week10_supplement/figures/S3_fixed_horizon_context.png).
Numerical Section 4 content is generated from the
[bound supplement claims](../artifacts/week10_supplement/interpretation/claims.json).
