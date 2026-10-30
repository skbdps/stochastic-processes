# Week 8 mathematical basis

Prepared with the dated [E4 protocol](week8_exploratory_protocol.md), before
its implementation and observations. Row-vector conventions and state order
(calm,burst) apply throughout. These are project derivations and implementation
requirements, not a certification of the owner's reading or theory mastery.

## Event sampling and a full-marginal independent control

The background generator has exit rates 1/(10m),1/(2m), hence stationary
probabilities pi=(5/6,1/6). Solving lambda_b=r*lambda_c and
pi_c*lambda_c+pi_b*lambda_b=1/m gives the protocol's intensities. Weighting
calendar-time state probability by the arrival intensity gives Palm
alpha=(5/(r+5),r/(r+5)). Bursts at r=20 account for 1/6 of time but 4/5 of
arrivals. Time and event stationarity answer different initialization questions.

For D0=Q-Lambda, absence of an arrival over t has transient matrix exp(D0*t).
The joint density kernel of gap duration and next arrival phase is
K(t)=exp(D0*t)*Lambda. Integrating gives P=(-D0)^(-1)*Lambda; alpha*P=alpha.
The first-arrival survival is alpha*exp(D0*t)*1. Resetting to an independent
alpha after every arrival preserves this complete one-gap law and removes
dependence between gaps. The MMPP joint gap density is
alpha*K(t1)*...*K(tK)*1, while the IID control factorizes into the product of
alpha*K(ti)*1. The matrix order cannot be reversed.

These general identities and phase-type moments are supported by
[Asanjarani and Nazarathy (2022), Eq.2.1, Proposition 2.2 and the following moment formulas](https://doi.org/10.1017/S1446181122000013).
The [PH renewal and MMPP constructions in Osogami's thesis](https://www.cs.cmu.edu/~osogami/thesis/html/node72.html)
also distinguish phase resetting, hidden switches and arrivals. The
calibration and specialized expressions below follow by substitution and
matrix integration; they are not new methodological claims.

## Closed calibration and dependence checks

Write A=-D0, d=r^2+70r+25 and beta=60r/d. Direct inversion gives

$$
A^{-1}\mathbf1=\frac{m(r+5)}{d}
\begin{pmatrix}11r+5\\r+15\end{pmatrix},\quad
P=\frac1d\begin{pmatrix}5(13r+5)&r(r+5)\\5(r+5)&r(r+65)\end{pmatrix}.
$$

The checks pi*Q=0, P*1=1, alpha*P=alpha and alpha*A^(-1)=m*pi imply
E[T]=m. Integration of t^k*exp(-A*t) gives
E[T^k]=k!*alpha*A^(-k)*1. Consequently

$$
\frac{E[T^2]}{m^2}=\frac{8(r+5)(7r+5)}{3d},\qquad
c^2=1+\frac{50(r-1)^2}{3d}.
$$

For M1=A^(-2)*Lambda, integration of the two-gap kernel gives
E[T_j*T_(j+k)]=alpha*M1*P^(k-1)*M1*1. Since P has eigenvalues 1 and beta,
its two-state decomposition simplifies this to

$$
\operatorname{Cov}(T_j,T_{j+k})=m^2\beta^k\frac{c^2-1}{2},\qquad
\rho_k=\beta^k\frac{c^2-1}{2c^2},\quad k\ge1.
$$

| r | Raw gap CV squared | Raw gap CV | Raw lag-one correlation |
|---:|---:|---:|---:|
| 1 | 1 | 1 | 0 |
| 5 | 5/3 | 1.2909944487 | 0.15 |
| 20 | 941/219 | 2.0728732843 | 0.2522527768 |

At r=1, equal arrival intensities give independent exponential gaps with mean
m, despite an evolving hidden chain. This is a Poisson limit, not a regular
grid. For a calendar-time origin the first residual has survival
pi*exp(D0*t)*1 and mean m*(1+c^2)/2; its ending phase has law pi*P, generally
not alpha. Discarding one residual is therefore not an event-stationary
initialization procedure.

## Identical floor transforms preserve the matched marginal

Let f be the applicable positive floor (zero means no transformation) and
Y=max(T,f). The transformed law has mass 1-alpha*exp(D0*f)*1 at f. Integrating
its survival, split at f, gives

$$
E[Y]=f+\alpha e^{D_0f}A^{-1}\mathbf1,
$$

$$
E[Y^2]=f^2+2\alpha e^{D_0f}(fA^{-1}+A^{-2})\mathbf1.
$$

For transformed lagged moments integrate the same joint kernel with the
transformed duration:

$$
R_f=\int_0^\infty\max(t,f)K(t)\,dt
=fP+e^{D_0f}A^{-2}\Lambda,
$$

$$
E[Y_jY_{j+k}]=\alpha R_fP^{k-1}R_f\mathbf1.
$$

For the independent control this cross-moment is instead E[Y]^2. These
formulas reduce to the raw checks when f=0. Applying the same transformation
to both mechanisms preserves their full marginal match but changes the
sampling clock; the floored sequence is not the original MMPP. No extra CTMC
evolution is inserted into floor padding. The implied transformed long-run
rate is 1/E[Y], not exactly the original 1/m.

## Conditional OU inference and interpretation

The observation mechanism is exogenous and independent of the OU process.
Given the realized times and initial observation, the exact transition
likelihood therefore remains the product of the same conditional Gaussian
transitions used by the core. The stationary initial OU density is excluded.
Gap dependence changes the sampling design, not the registered inference
definition. PFML/Euler still use model-based working covariance; matching gap
marginals is not a justification for claiming their Wald intervals are robust.

The matched designs have equal expected transformed span K*E[Y], but may have
different span variance because of gap dependence. Thus their comparison
includes the consequences of the changed random span distribution. It does
not condition on or force equal realized spans. Between-cell contrasts must
use independent-sample uncertainty; only methods on one source path are paired.
