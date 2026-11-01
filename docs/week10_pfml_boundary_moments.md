# Week 10 critical review: the fitted PFML mean has a boundary pole

Dated 1 November 2026. This is an analytical correction to the
[mechanism argument](week10_mechanism.md), with deterministic regression
checks in [test_boundary_argument.py](../tests/week10/test_boundary_argument.py).
It does not modify the estimators, the sealed experiment, or its saved results.
All derivations below are explicit consequences of the conditional Gaussian
regression and the nondegenerate OU transition density.

## 1. What the covariance argument does and does not prove

The fixed-coefficient diagnostic is a linear function of the observations.
Its schedule-conditional Gaussian law and the covariance formulas in the
mechanism note remain exact. The jointly fitted mean is a different,
nonlinear estimator:

$$
\widehat\mu_P=\overline X_{0:K-1}
 +\frac{X_K-X_0}{K(1-\widehat\rho)},\qquad K=N-1.
$$

Its denominator can approach zero while the fitted slope is strictly inside
the admissible interval. For the ideal estimator in real arithmetic, merely
requiring $0<\widehat\rho<1$ does not ensure finite moments. In particular,
the finite variance of the fixed-coefficient diagnostic cannot be transferred
to the fitted estimator by an assertion that their difference is small in
probability or small on the saved simulation paths.

## 2. A conditional proof of the pole

Fix an exogenous schedule with positive gaps and $K\ge3$. Condition on the
first $K$ observations $z=(X_0,\ldots,X_{K-1})$, and write $x=X_K$,
$\bar z=K^{-1}\sum_i z_i$, and $S=\sum_i(z_i-\bar z)^2$. For a nonconstant
prefix, the unrestricted regression slope is affine in its last response:

$$
\widehat\rho(x)=a+bx,\qquad
a=\frac{\sum_{i=0}^{K-2}(z_i-\bar z)z_{i+1}}S,
\qquad b=\frac{z_{K-1}-\bar z}S.
$$

The event $b=0$ is a hyperplane in a prefix with a nondegenerate Gaussian
density, so it has probability zero. To show that the numerator does not
cancel the pole, evaluate the slope at $x=z_0$. The response vector is then
the circular shift of the predictor vector. With $z_K=z_0$,

$$
1-\widehat\rho(z_0)
=\frac{1}{2S}\sum_{i=0}^{K-1}(z_{i+1}-z_i)^2>0.
$$

Consequently the finite root $x_*=(1-a)/b$ of $\widehat\rho(x_*)=1$ is
different from $z_0$. The exact estimator identity now gives

$$
\widehat\mu_P(x)
=\bar z-\frac{x-z_0}{Kb(x-x_*)}
=\left(\bar z-\frac1{Kb}\right)
 -\frac{x_*-z_0}{Kb(x-x_*)}.
$$

The pole coefficient is nonzero. On one side of $x_*$, a sufficiently small
neighborhood satisfies $0<\widehat\rho(x)<1$. At the pole the residual SSE
is positive for a generic prefix: zero SSE at slope one would require all
successive increments of $z$ to be equal, an arithmetic-progression subset
of Gaussian probability zero when $K\ge3$. Thus interior fits approaching
the pole need not be degenerate regressions.

Conditional on the prefix and schedule, the OU transition law for $X_K$ is
Gaussian with positive variance
$V(1-e^{-2\theta_0(t_K-t_{K-1})})$. Its density is strictly positive and
continuous at $x_*$. Therefore, on the interior side alone, the conditional
integral of $|\widehat\mu_P-\mu_0|^p$ diverges for every $p\ge1$:

$$
\int_0^\varepsilon c\,u^{-p}\,du=\infty.
$$

This proves infinite second and fourth moments for the ideal interior
estimator, including when its law is conditioned on being an admissible
interior fit. It is not necessary to observe a boundary fit for this issue
to occur. Gaussian reflection about $\mu_0$ preserves the slope and
admissibility and changes the sign of the fitted mean error. Accordingly,
the ideal error has divergent positive and negative first moments; symmetry
does not define an ordinary finite expectation or finite-sample bias.
The population objective's mean target $\mu_0$ is a separate, valid statement.

The proof applies conditional on each fixed-horizon schedule with $N\ge4$
as well as to the fixed-count design. Random counts and a deterministic
terminal observation do not remove the positive final-transition variance.

## 3. A concrete nondegenerate example and implementation scope

Take $z=(0,1,3)$ and let the terminal observation be $x$. Direct regression
gives

$$
\widehat\rho(x)=\frac{5x-7}{14},\qquad
\widehat\mu_P(x)=\frac{28-2x}{21-5x}.
$$

At $x_*=4.2$, the fitted intercept is $1.4$, the residual vector is
$(-0.4,0.6,-0.2)$, and SSE is $0.56$. For $x=4.2-\delta$ with
$\delta>0$ small,

$$
\widehat\rho=1-\frac{5\delta}{14},\qquad
\widehat\mu_P=0.4+\frac{3.92}{\delta}.
$$

The deterministic tests check the unchanged fitter for
$\delta=10^{-2},\ldots,10^{-6}$: these are point-valid interior solutions,
and the fitted mean follows the pole while the mean at a fixed coefficient
remains bounded. These inputs are checks of algebra and existing behavior,
not new research replications.

The existing point-valid rule excludes slopes outside $(0,1)$, degenerate
regressions, and nonrepresentable outputs. It does not impose a scientific
lower bound $1-\widehat\rho\ge\varepsilon$. Numerical degeneracy thresholds
also leave the well-scaled example and a neighborhood of it eligible.
The interval procedure has additional Hessian and conditioning checks.
An interval rejection does not remove an otherwise valid point estimate
from point summaries.

The divergence proof concerns the ideal real-arithmetic estimator. It does
not assert literally infinite moments of a finite-precision computer
implementation or its finite random-number support. Rounding, overflow,
and representability checks change the extreme tail in implementation-
dependent ways; they are not a declared statistical regularization and do
not supply a uniform moment bound for the ideal estimator.

On a regular mesh, the exact conditional MLE and Euler regression have the
same fitted mean as PFML whenever the common AR slope lies in $(0,1)$.
They share this mean-parameter pole on that region. The finding is not a
claim that PFML alone can have a weakly identified mean near zero reversion,
nor a proof about every parameter of every irregular-grid estimator.

## 4. The exact working Hessian is still distinct from sampling variance

At a nondegenerate interior conditional regression optimum, write
$\widehat v=\mathrm{SSE}/K$, $\bar x=K^{-1}\sum_{i=0}^{K-1}X_i$, and
$S_{xx}=\sum_{i=0}^{K-1}(X_i-\bar x)^2$. In coordinates $(a,\rho,v)$,
the observed Hessian's regression block is

$$
H_{(a,\rho)}=\frac1{\widehat v}
\begin{pmatrix}K&K\bar x\\K\bar x&\sum_iX_i^2\end{pmatrix}.
$$

The cross derivatives with $v$ vanish at OLS because both residual normal
equations are zero. Mapping $\mu=a/(1-\rho)$ through this inverse Hessian
therefore gives the exact mean entry of the working covariance:

$$
\widehat v_{\mu,\mathrm{Hess}}
=\frac{\widehat v}{(1-\widehat\rho)^2}
\left\{\frac1K+
\frac{(\widehat\mu_P-\bar x)^2}{S_{xx}}\right\}.
$$

This identity includes joint estimation of slope and residual variance.
At a stationary optimum a smooth reparameterization to OU coordinates
preserves this inverse-Hessian mean entry. The study computes it numerically
in normalized log coordinates; the deterministic tests compare that actual
implementation with the formula, including affine state and time changes.
The formula is a working-information identity, not a sampling-variance
identity under misspecification, and not a justification of the fitted
mean's moments.

If $\widehat\rho\to_p\phi<1$, $\widehat v\to_p V(1-\phi^2)$,
$S_{xx}/K\to_p V>0$, and stationary endpoints are $O_p(1)$, the additional
term involving $(\widehat\mu_P-\bar x)^2$ is $O_p(K^{-3})$ inside braces.
The leading mean covariance is consequently
$V(1+\phi)/[K(1-\phi)]$, as used in the mechanism note. This is an
in-probability statement about the reported working covariance, not
convergence of the fitted estimator's second moment.

## 5. Consequences for interpretation and Monte Carlo uncertainty

The exact discrepancy identity in the mechanism note implies
$o_p(K^{-1})$ if $\widehat\rho\to_p\phi<1$, and
$O_p(K^{-3/2})$ if the slope is also root-$K$ consistent, under the stated
endpoint condition. Those stochastic orders remain valid conditional
implications. They do not imply convergence in mean square. A suitable
uniform-integrability or tail bound would be needed to transfer moments;
the pole explains why one cannot silently assume such a bound here.

The reported empirical RMSE, empirical SD, empirical bias, and saved-path
RMS fitted-minus-fixed discrepancy remain reproducible finite-sample
descriptions. Their small observed discrepancy verifies agreement on the
saved paths; it does not establish an approximation of ideal population
MSE or existence of a fitted population variance. The historical RMSE and
log-ratio delta intervals must be interpreted with that limitation.
Ordinary bias Monte Carlo standard errors require a finite error variance;
the usual delta-method RMSE standard error requires a finite fourth error
moment as well as a positive, finite second moment. Finite values computed
from 2,000 realizations do not demonstrate these population assumptions.

Interval-coverage indicators are bounded even when fitted means have extreme
tails. With independent source replications, the existing operational and
conditional coverage summaries retain their usual binomial interpretation,
with their stated denominators and pointwise uncertainty. The covariance
identities for the fixed-coefficient diagnostic also remain valid.
Any new truncation, bound away from $\rho=1$, or alternative loss for fitted
means must be declared separately; none is retroactively inserted into the
sealed experiment by this correction.
