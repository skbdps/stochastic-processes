# Week 9: population targets and empirical uncertainty calibration

Dated 31 October 2026. This note fixes the interpretation of the saved studies;
it introduces no new estimator, interval procedure or fitted experiment. Its
scope follows the population-target and empirical-calibration option in the
research proposal. The full mathematical convention remains in the
[Week 3](https://github.com/skbdps/stochastic-processes/blob/b6edb85d1e7a21ed19fb54b5c71213ed04703309/docs/week3_mathematics.md), [Week 4](https://github.com/skbdps/stochastic-processes/blob/b6edb85d1e7a21ed19fb54b5c71213ed04703309/docs/week4_mathematics.md),
[Week 5](https://github.com/skbdps/stochastic-processes/blob/b6edb85d1e7a21ed19fb54b5c71213ed04703309/docs/week5_mathematics.md), and
[prospective E4 target](https://github.com/skbdps/stochastic-processes/blob/b6edb85d1e7a21ed19fb54b5c71213ed04703309/docs/week8_population_targets.md) notes.

## Model and estimator conventions

Data follow the stationary noise-free OU model

$$
dX_t=\theta_0(\mu_0-X_t)\,dt+\sigma_0\,dW_t,
\qquad V=\sigma_0^2/(2\theta_0).
$$

The schedule is exogenous to the initial value and entire innovation stream.
Simulation begins stationary, while each fitted likelihood conditions on X0
and excludes its stationary density. The three parameters are jointly estimated.
Exact likelihood uses each stored positive time difference. PFML replaces every
gap in its exact Gaussian transition objective by the observed arithmetic mean
gap. Euler uses the actual gap in its Gaussian mean and variance.

For a fixed observation vector, PFML is an AR(1) regression with intercept.
If its slope is rho_hat and intercept a_hat, mu_hat=a_hat/(1-rho_hat).
The mean-gap scale maps rho_hat to theta_hat, but it does not mechanically enter
this mu formula. Timing affects mu's repeated-sample behavior through the
dependence of the observations. The implementation's boundary outcomes and
eligibility rules remain part of the stated procedure.

## Expected-objective targets

Let Y denote a gap after the specified floor transformation, and put

$$
\bar m=E[Y],\quad L(s)=E[e^{-sY}],\quad
H(\theta_0)=E[(1-e^{-\theta_0Y})/Y]
=\int_0^{\theta_0}L(s)\,ds.
$$

Conditional stationarity and exogeneity give a mean-zero centered state with
variance V and adjacent covariance V L(theta0). Minimizing the expected PFML
working objective yields

$$
\theta_{P,*}=-\log L(\theta_0)/\bar m,\qquad
\mu_{P,*}=\mu_0,\qquad \sigma_{P,*}^2=2V\theta_{P,*}.
$$

For the actual-gap Euler objective, the expected profiled weighted residual
variance is

$$
V\{2H(\theta_0)-2\theta(1-L(\theta_0))+\theta^2\bar m\}
+\theta^2\bar m(\mu-\mu_0)^2.
$$

Its minimizers are

$$
\theta_{E,*}=(1-L(\theta_0))/\bar m,\qquad
\mu_{E,*}=\mu_0,\qquad
\sigma_{E,*}^2=V\{2H(\theta_0)-(1-L(\theta_0))^2/\bar m\}.
$$

The H term is essential under random gaps. The constant-gap expression
V(1-L(theta0)^2)/mbar is not the random-gap Euler diffusion target. An
appropriate ergodic law is needed to identify an expected-objective minimizer
with an estimator's probability limit. These targets are not finite-n mean or
coverage predictions. All experiment bias and coverage summaries continue to
use the true OU parameters, without recentering at a working target.

For raw Gamma gaps with mean m and CV c, L(s)=(1+s*m*c^2)^(-1/c^2).
At c=0, use the separate constant-gap law. When a floor f applies, its atom and
upper tail belong in both L and E[Y]:

$$
L_Y(s)=e^{-sf}F_D(f)+\int_f^\infty e^{-sd}f_D(d)\,dd,\qquad
E[Y]=fF_D(f)+E[D\,1\{D\ge f\}].
$$

The floor is not rejection truncation and is not followed by normalization.
Raw-Gamma targets must therefore be labelled as such when compared with floored
designs. The prospective E4 note already evaluates the actual phase-type law,
including the r=20 floor, and records independent deterministic calculations.

## Dependence and covariance

E4's MMPP and independent control share the entire used one-gap law at each r.
They therefore share these expected-objective targets under the declared
stationarity and exogeneity. Their finite-sample precision and coverage can
still differ because serial gap dependence changes the joint observation law,
score covariance, span distribution and validity composition. The matched
contrast includes realized-span consequences; span is not held fixed.

The original Wald covariance is the inverse full observed Hessian in
observed-data normalized coordinates, mapped to natural units by the recorded
Jacobian. PFML/Euler covariance is model-based under a possibly misspecified
working objective. A misspecification/dependence-aware covariance generally
also needs covariance between estimating-equation contributions, including the
effect of estimated nuisance quantities. No such robust method is substituted.

Ait-Sahalia and Mykland's [2004 paper](https://www.princeton.edu/~yacine/gmmrandom.pdf)
uses IID intervals in Assumption 2 and treats specific OU nuisance conventions
in Section 4.3/Table 1. Its coefficients do not automatically apply to the
present jointly fitted mean, finite-n random-span design or dependent E4 gaps.
A dependent-gap CLT and a fully matched covariance benchmark remain
unestablished here. Consequently the current evidence does not demonstrate a
departure from a matched first-order coverage prediction.

## Saved empirical calibration

For each study/cell/method/parameter, compare the sample SD of estimates among
interval-valid paths with the RMS reported SE among those same paths. Also
retain the mean reported SE as a separate quantity. The descriptive ratio
empirical_SD/RMS_reported_SE indicates relative spread on this eligibility
subset. It is not an estimator of operational coverage, an equivalence test or
a derived long-run covariance benchmark. Do not assign it an uncertainty band
borrowed from a different statistic.

The saved summaries already contain these components. At E4 r=20, PFML mu's
empirical SD is about 0.0354051 for MMPP versus 0.0305712 for IID; RMS reported SE
is about 0.0309250 versus 0.0307701. These values are descriptive evidence for
the interval calibration discussion. Publication values must be generated from
the bound source rows rather than copied from this rounded illustration.

Read calibration beside operational coverage, conditional coverage and
availability. E2's 85 invalid points and one further interval rejection yield
86 unavailable intervals, not 171 distinct failures. E3's follow-up similarly
retains its unavailable outcomes. Strong conditional calibration cannot erase
the operational denominator.

## Permitted conclusions

The study can report population-target explanations for approximation bias,
measured finite-sample error, availability, and empirical interval calibration
within the tested designs. The exact likelihood remains correctly specified
conditional on the schedule, while finite-sample bias remains possible. E2's
short-span exact-theta breaches are bias-driven; its anticipated below-80%
coverage collapse did not occur. E4's PFML mu differences concern precision and
coverage, without a demonstrated mean-bias shift. No finding here requires a
method to win, a failure boundary to appear, or reported coverage to be tuned.

Dependence-aware sandwich/HAC inference, block bootstrap calibration,
bias-corrected procedures, new span-matched studies and general dependent-gap
limit theory would be separately specified research extensions. The current
manuscript should state the benchmark's present scope directly.
