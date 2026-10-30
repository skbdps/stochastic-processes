# Week 8 prospective population targets for the unchanged estimators

Dated 30 October 2026, after the source-only validation and before any fitted
E4 validation or full E4 experiment. This is an interpretation calculation,
not a change to the protocol, estimators, inference gates or truth-based
simulation summaries. It explains why a matched MMPP/IID comparison can have
different finite-sample behavior despite identical population targets.

## 1. Definitions and binding to the existing implementation

The frozen `ou_irregular/estimators.py` has SHA-256
`75848a02d9479991bdd7d2dcf384aa4075f6c96a6638adeb661113134a2558e9`.
PFML replaces the actual gaps by their realized arithmetic mean and fits an
AR(1) regression with an intercept. Euler uses every actual gap and fits the
weighted regression of Delta X/sqrt(Y) on sqrt(Y) and -X_previous sqrt(Y).
Normalizing Euler's units by the realized mean gap does not change that
weighted objective into a constant-gap objective. Both condition on X0.

The actual-gap Euler Gaussian objective agrees with equation (38) of
[Ait-Sahalia and Mykland (2004), Annals of Statistics 32, 2186–2222](https://www.princeton.edu/~yacine/gmmrandom.pdf),
DOI [10.1214/009053604000000427](https://doi.org/10.1214/009053604000000427).
That paper's Assumption 2 uses IID gaps. The dependent-gap conclusion below
comes from the explicit stationary expected-objective calculation here; it
does not extend the paper's IID asymptotic results without further conditions.

Let the true parameters be (theta0,mu0,sigma0), let Y be an actually used gap
after any declared floor, and define

$$
V=\frac{\sigma_0^2}{2\theta_0},\qquad
\bar m=E[Y],\qquad L(s)=E[e^{-sY}],\qquad \phi=L(\theta_0),
$$

$$
H(\theta_0)=E\left[\frac{1-e^{-\theta_0Y}}Y\right]
=\int_0^{\theta_0}L(s)\,ds.
$$

The sampling law is held fixed as candidate fitted parameters vary. The
schedule is independent of the OU initial value and entire innovation stream,
and X0 is stationary. For Z_i=X_i-mu0, conditional on the entire schedule,

$$
E[Z_{i-1}\mid\text{schedule}]=0,\qquad
E[Z_{i-1}^2\mid\text{schedule}]=V,\qquad
E[Z_{i-1}Z_i]=V\phi.
$$

These conditional stationarity statements do not assume independent
observations or independent successive gaps.

## 2. Mean-gap PFML

At limiting gap mbar, write rho=exp(-theta mbar) and
v=sigma^2(1-rho^2)/(2theta). The expected squared working residual is

$$
V(1+\rho^2-2\rho\phi)+(1-\rho)^2(\mu-\mu_0)^2.
$$

The interior minimizer is rho*=phi, mu*=mu0 and v*=V(1-phi^2). Therefore

$$
\boxed{\theta_{P,*}=-\frac{\log L(\theta_0)}{\bar m},\quad
\mu_{P,*}=\mu_0,\quad \sigma_{P,*}^2=2\theta_{P,*}V.}
$$

An appropriate law of large numbers makes the realized mean gap converge to
mbar. Its finite-sample fluctuations still affect the fitted estimator and
its uncertainty. The population formula does not replace the observed mean
by the nominal mean in any actual fit.

## 3. Euler at actual gaps

Euler profiles sigma^2 using the mean of R_i^2/Y_i, where
R_i=Delta X_i-theta Y_i(mu-X_{i-1}). Its expectation is

$$
E[R_i^2/Y_i]
=V\{2H(\theta_0)-2\theta(1-\phi)+\theta^2\bar m\}
+\theta^2\bar m(\mu-\mu_0)^2.
$$

Consequently,

$$
\boxed{\theta_{E,*}=\frac{1-L(\theta_0)}{\bar m},\quad
\mu_{E,*}=\mu_0,\quad
\sigma_{E,*}^2=V\left\{2H(\theta_0)-
\frac{(1-L(\theta_0))^2}{\bar m}\right\}.}
$$

The H term is essential for the actual-gap Euler estimator. The expression
V(1-phi^2)/mbar applies only to constant gaps; it is not its diffusion target
under random gaps. No separate E[1/Y] assumption is needed for this residual
expectation because 0<(1-exp(-theta0 Y))/Y<=theta0.

Jensen's inequality and 1-phi<=-log(phi) give
theta_E,*<=theta_P,*<=theta0. Both diffusion population targets are also at
most sigma0 under this exogenous stationary setup. These inequalities are
population statements, not restrictions imposed on individual fitted values.

## 4. Phase-type transform and the declared floor

Use the event-stationary quantities A=Lambda-Q, alpha and
S_T(t)=alpha exp(-At)1 from the [fixed mathematical appendix](week8_mathematics.md).
Integrating the density gives

$$
L_T(s)=\alpha(A+sI)^{-1}\Lambda\mathbf1
=1-s\alpha(A+sI)^{-1}\mathbf1.
$$

For Y=max(T,f), include the atom at f:

$$
L_Y(s)=e^{-sf}\left\{F_T(f)+\alpha e^{-Af}(A+sI)^{-1}\Lambda\mathbf1\right\}
=e^{-sf}\left\{1-s\alpha e^{-Af}(A+sI)^{-1}\mathbf1\right\},
$$

$$
\bar m=f+\alpha e^{-Af}A^{-1}\mathbf1.
$$

The transformed mean belongs in the target formulas. It is not reset to the
nominal raw mean. For the adopted two-state calibration, put z=sm and

$$
\tau_r=\frac35+\frac{6(r+1)}{r+5},\qquad
\kappa_r=\frac35+\frac{36r}{(r+5)^2}.
$$

Then the raw transform is

$$
L_T(s)=\frac{\kappa_r+(\tau_r-\kappa_r)z}
{z^2+\tau_r z+\kappa_r}.
$$

At m=1/2 and s=1 it equals 2/3, 66/95 and 1848/2411 for r=1,5,20.
The r=20 transformed value includes the small floor correction.

## 5. Deterministic reference values

These are large-n targets at fixed gap laws and truth (1,0,0.5), not predictions
of the exact mean estimate at n=1000. MMPP and its matched IID control share
each row's targets because their full one-gap marginals agree. Every mu target
is zero; relative mu error remains undefined.

| Gap law | E[Y] | L(1) | H(1) | PFML theta | PFML sigma | Euler theta | Euler sigma |
|---|---:|---:|---:|---:|---:|---:|---:|
| MMPP or IID, r=1 | 0.500000000000000 | 0.666666666666667 | 0.810930216216329 | 0.810930216216329 | 0.450258319250275 | 0.666666666666667 | 0.418275957086114 |
| MMPP or IID, r=5 | 0.500000000000000 | 0.694736842105263 | 0.823040480694884 | 0.728444299148231 | 0.426744741955959 | 0.610526315789474 | 0.427157724127463 |
| Floored MMPP-derived or IID, r=20 | 0.500000000000972 | 0.766486934880820 | 0.856601742685284 | 0.531875251960585 | 0.364648890016337 | 0.467026130237453 | 0.447792750919450 |
| Regular | 0.500000000000000 | 0.606530659712633 | 0.786938680574733 | 1.000000000000000 | 0.500000000000000 | 0.786938680574733 | 0.397530048810325 |

The independent calculation used 80-digit Decimal exponential-mixture
formulas, with double-precision matrix-resolvent, matrix-exponential and
quadrature agreement within 2e-14. It generated no schedules or fits. Root
independently checked both estimator objectives and the expected-residual
derivations. The r=20 floor shifts PFML theta by about +1.50e-12 and Euler
theta by +1.04e-12; that defined effect is included in the table.

## 6. What matched marginals do and do not imply

The two expected objectives depend only on one adjacent observation pair and
its gap. Under exogeneity and stationary OU initialization, that pair law is
determined by the one-gap marginal. The matched MMPP/IID designs therefore
have identical population minimizers. The r-curve changes the marginal and
its targets; it does not by itself isolate serial dependence.

Serial dependence can still change covariance between estimating-equation
contributions, realized mean gap and span variation, finite-sample bias,
boundary events and valid-fit composition, long-run variance, and coverage.
Matching targets neither proves equal variances nor justifies model-based
Wald covariance for a misspecified objective. The existing covariance method
is retained without a sandwich correction.

Interpreting population minimizers as probability limits requires the relevant
ergodic laws of large numbers; a central limit theorem needs additional
dependence and moment conditions. No general dependent-gap CLT or variance
formula is certified here. Registered simulation summaries continue to assess
the true OU parameters, with all planned outcomes and the existing eligibility
rules. No result is recentered at these pseudo-targets.
