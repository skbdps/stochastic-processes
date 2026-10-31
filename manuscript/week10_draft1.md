# Parameter estimation for the Ornstein–Uhlenbeck process under irregular sampling

Suryansh Kumar

**Draft 1 · 1 November 2026 · Week 10 methods baseline**

This first manuscript draft contains the introduction, model, estimators,
mathematical benchmarks and experimental methods. It is based on the completed
research release `4bddb527f3c8bb88bc1ebec28c8818e03dc251aa`. Results, practical
guidance, a full related-work discussion, conclusions and the final abstract
will be developed in the subsequent manuscript stage. This is an initial
working draft for author review, not a submitted or academically reviewed paper.

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
claims, plotting data and actual AI visual review. Its
[reproduction instructions](../docs/week9_reproduction.md) distinguish
saved-data rendering from original numerical simulation and replay. The
recorded macOS rendering environment differs from the original Linux
scientific environment; exact byte replay across different CPU dispatch,
precision or graphics environments is not assumed. The archived package,
failure dispositions and limitations remain unchanged by this manuscript.

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

The [source map](week10_source_map.md) identifies the research protocols,
mathematical notes and implementation files behind these methods. The
[Week 9 literature comparison](../docs/week9_literature_comparison.md)
provides the broader comparison and its access limits for the later related-work
section.
