## 7. Fit the Gaussian AR(1) {#ar-fit}

### 7.1 The observations are fixed; the model parameters are not

You observe a numerical record $x_0,x_1,\ldots,x_n$ at equally spaced ticks. Assume the model

$$
X_i=\mu+\phi(X_{i-1}-\mu)+\varepsilon_i,
\qquad \varepsilon_i\sim N(0,\tau^2),
$$

with fresh independent Gaussian innovations, independent of the start. Here $\mu$ is the mean-reversion level, $\phi$ is the one-tick retention coefficient, and $\tau^2>0$ is the innovation variance. The Gaussian notation $N(m,v)$ uses variance in its second slot.

Unlike a calculation predicting the next value from known parameters, we now seek the unknown parameter values. The data are the observed states. The candidate parameter vector is $(\mu,\phi,\tau^2)$. A candidate's predictions, residuals, and likelihood are quantities computed from that candidate and the fixed data.

We will use a **conditional** likelihood: the starting value $x_0$ is supplied as a starting condition. There are $n+1$ recorded values but exactly $n$ observed transitions. That is why the likelihood below contains $n$ factors.

The goal is not to make every residual zero. The model allows genuine random noise, and a candidate must account for both its predicted means and its stated uncertainty. Maximum likelihood finds the candidate giving the greatest joint conditional density to the entire record.

### 7.2 Write one transition density before writing a product

Given the previous observation $x_{i-1}$, the candidate mean is

$$
m_i=\mu+\phi(x_{i-1}-\mu).
$$

All that remains random is the innovation, so the transition density evaluated at the observed next value $x_i$ is

$$
g(x_i\mid x_{i-1};\mu,\phi,\tau^2)
=\frac{1}{\sqrt{2\pi\tau^2}}
\exp\left[-\frac{(x_i-m_i)^2}{2\tau^2}\right].
$$

The residual $e_i=x_i-m_i$ is the observed value minus the candidate's one-step prediction. The word “candidate” matters: changing $\mu$ or $\phi$ changes every $m_i$, hence every residual, even though the observations stay fixed.

The Markov property follows from the recursion and fresh-noise assumption. Applying the probability chain rule and then shortening each conditional history gives

$$
L_c(\mu,\phi,\tau^2)=\prod_{i=1}^n g(x_i\mid x_{i-1};\mu,\phi,\tau^2).
$$

The states are dependent. This product is valid because it multiplies **conditional densities**, not because the observations are independent.

### 7.3 Take the logarithm carefully

For one factor, the logarithm is

$$
\log g_i=-\frac12\log(2\pi\tau^2)-\frac{e_i^2}{2\tau^2}.
$$

Adding the $n$ factors' logarithms gives

$$
\boxed{
\ell(\mu,\phi,\tau^2)
=-\frac n2\log(2\pi\tau^2)
-\frac{1}{2\tau^2}\sum_{i=1}^n[x_i-\mu-\phi(x_{i-1}-\mu)]^2.
}
$$

The same variance appears in every transition because the model assumes a common innovation distribution per tick. With irregular OU gaps, the transition variances will differ, and this particular constant-variance simplification will no longer hold.

If $\tau^2$ is fixed and positive, the first term does not depend on $\mu$ or $\phi$. Maximizing the second term is equivalent to minimizing the squared-residual sum, because its coefficient is negative. Thus a Gaussian likelihood creates a least-squares problem. We still need to solve it and then estimate the variance.

### 7.4 Why the intercept-and-slope form helps

Expand the conditional mean:

$$
\begin{aligned}
m_i&=\mu+\phi x_{i-1}-\phi\mu\\
&=\mu(1-\phi)+\phi x_{i-1}.
\end{aligned}
$$

Define

$$
a=\mu(1-\phi),\qquad b=\phi.
$$

Also name the two observed arrays

$$
u_i=x_{i-1},\qquad y_i=x_i,\qquad i=1,\ldots,n.
$$

The $u_i$ are the lagged values and the $y_i$ are their successors. The same interior observation occurs once as a response and once as a later predictor; that does not make them separate independent samples.

The prediction is now $a+bu_i$, and the residual sum of squares is

$$
S(a,b)=\sum_{i=1}^n(y_i-a-bu_i)^2.
$$

This is an ordinary intercept-and-slope least-squares objective. The transformation is reversible when $b\ne1$:

$$
\phi=b,\qquad \mu=\frac{a}{1-b}.
$$

At $\phi=1$, the original centered model becomes a random walk with no level parameter; $\mu$ disappears. The reverse formula's division by zero reflects that genuine identification problem.

### 7.5 Derive the intercept equation

Differentiate $S$ with respect to $a$. For each residual $r_i=y_i-a-bu_i$, the derivative of $r_i^2$ is $2r_i$ times the derivative of $r_i$, which is $-1$. Therefore

$$
\frac{\partial S}{\partial a}
=-2\sum_i(y_i-a-bu_i).
$$

Set the derivative to zero:

$$
\sum_i y_i-na-b\sum_i u_i=0.
$$

Define the two means separately,

$$
\bar u=\frac1n\sum_i u_i,\qquad \bar y=\frac1n\sum_i y_i.
$$

Divide the equation by $n$ and solve for $a$:

$$
\boxed{a=\bar y-b\bar u.}
$$

For any chosen slope, this intercept makes the residuals sum to zero. The fitted line passes through the point $(\bar u,\bar y)$. We have reduced a two-variable minimization to a one-variable problem in $b$.

### 7.6 Derive the slope equation and prove it minimizes the objective

Substitute $a=\bar y-b\bar u$ into the residual:

$$
\begin{aligned}
y_i-a-bu_i
&=y_i-(\bar y-b\bar u)-bu_i\\
&=(y_i-\bar y)-b(u_i-\bar u).
\end{aligned}
$$

Introduce three sums, all computed solely from the observed arrays:

$$
S_{uu}=\sum_i(u_i-\bar u)^2,
\qquad
S_{uy}=\sum_i(u_i-\bar u)(y_i-\bar y),
$$

$$
S_{yy}=\sum_i(y_i-\bar y)^2.
$$

Expanding the squared residuals gives

$$
\begin{aligned}
S(b)
&=\sum_i[(y_i-\bar y)-b(u_i-\bar u)]^2\\
&=S_{yy}-2bS_{uy}+b^2S_{uu}.
\end{aligned}
$$

Differentiate:

$$
S'(b)=-2S_{uy}+2bS_{uu}.
$$

Provided $S_{uu}>0$, setting this to zero gives

$$
\boxed{\hat b=\frac{S_{uy}}{S_{uu}}},
\qquad
\boxed{\hat a=\bar y-\hat b\bar u.}
$$

The second derivative is $S''(b)=2S_{uu}>0$, so this is the unique minimum of the profiled quadratic. Another explicit check is to complete the square:

$$
S(b)=S_{yy}-\frac{S_{uy}^2}{S_{uu}}
+S_{uu}\left(b-\frac{S_{uy}}{S_{uu}}\right)^2.
$$

The last term is nonnegative and becomes zero exactly at the derived slope. This proves global minimization over unrestricted real slopes, rather than merely locating a stationary point.

If $S_{uu}=0$, all lagged observations are equal. Then the data only identify the combination $a+bu_i$; slope and intercept cannot be separated by this conditional record. The formula's zero denominator diagnoses a lack of information, not a software inconvenience.

### 7.7 Derive the innovation-variance estimate

Write $v=\tau^2$ temporarily so there is no ambiguity about whether we differentiate in variance or standard deviation. After choosing $\hat a,\hat b$, let

$$
S_*=S(\hat a,\hat b).
$$

The log-likelihood as a function of $v>0$ is

$$
\ell(v)=-\frac n2\log(2\pi v)-\frac{S_*}{2v}.
$$

Differentiate term by term:

$$
\ell'(v)=-\frac{n}{2v}+\frac{S_*}{2v^2}
=\frac{S_*-nv}{2v^2}.
$$

If $S_*>0$, the derivative is positive when $v<S_*/n$ and negative when $v>S_*/n$. The unique maximum is therefore

$$
\boxed{\hat\tau^2=\frac{S_*}{n}.}
$$

The denominator is the number of conditional transition factors, $n$. A degrees-of-freedom correction answers a different estimation question; it is not what differentiating this likelihood produces.

If $S_*=0$, the fit interpolates the record perfectly. The likelihood then increases without bound as $v$ decreases toward zero, and no positive interior variance MLE exists. A robust implementation should report this degeneracy rather than silently return a valid-looking standard model fit.

Because minimizing $S$ was the best choice for every fixed $v>0$, combining the least-squares minimizer with $\hat v=S_*/n$ maximizes the joint conditional likelihood whenever the solution is admissible and nondegenerate.

### 7.8 Convert back to the AR parameters

The conditional estimates are

$$
\boxed{\hat\phi=\hat b,\qquad
\hat\mu=\frac{\hat a}{1-\hat b},\qquad
\hat\tau^2=\frac1n\sum_i(y_i-\hat a-\hat b u_i)^2.}
$$

The fitted mean-reversion level is not generally the average of all the recorded states. The intercept-and-slope fit is made to **paired successive observations**, with different lagged and successor means. The transformation $\hat\mu=\hat a/(1-\hat b)$ is required by the model's definition.

When $\hat b$ is close to one, a small change in the fitted intercept or slope can cause a large change in $\hat\mu$. The denominator $1-\hat b$ makes this visible. A persistent short record may estimate local movement better than a distant long-run level.

### 7.9 A complete numerical fit, including every intermediate quantity

Take the record

$$
(x_0,x_1,x_2,x_3)=(0,1,1,2).
$$

There are three transitions, so $n=3$. The lagged array is $u=(0,1,1)$ and the successor array is $y=(1,1,2)$. Their means are $\bar u=2/3$ and $\bar y=4/3$.

| $i$ | $u_i$ | $y_i$ | $u_i-\bar u$ | $y_i-\bar y$ | Product of deviations |
|---:|---:|---:|---:|---:|---:|
| 1 | 0 | 1 | $-2/3$ | $-1/3$ | $2/9$ |
| 2 | 1 | 1 | $1/3$ | $-1/3$ | $-1/9$ |
| 3 | 1 | 2 | $1/3$ | $2/3$ | $2/9$ |

Thus

$$
S_{uy}=\frac29-\frac19+\frac29=\frac13,
\qquad
S_{uu}=\frac49+\frac19+\frac19=\frac23.
$$

The slope is $(1/3)/(2/3)=1/2$, and the intercept is $4/3-(1/2)(2/3)=1$. Convert to the mean-reversion level:

$$
\hat\phi=\frac12,\qquad
\hat\mu=\frac{1}{1-1/2}=2.
$$

The fitted predictions are $(1,1.5,1.5)$. Subtract them from the observed successors $(1,1,2)$ to obtain residuals $(0,-0.5,0.5)$. Their squared sum is $0+0.25+0.25=0.5$, hence

$$
\hat\tau^2=\frac{0.5}{3}=\frac16.
$$

The result is the conditional fit $(\hat\mu,\hat\phi,\hat\tau^2)=(2,1/2,1/6)$. The observed average of all four states is $1$, not $2$. There is no contradiction: the fitted level is inferred from the relation between each state and its successor, not defined as the unweighted record average.

If we initialized a process at stationarity using these fitted parameters, its state variance would be $\hat\tau^2/(1-\hat\phi^2)=(1/6)/(3/4)=2/9$. That is a **derived model quantity**, distinct from the fitted innovation variance $1/6$.

### 7.10 Admissibility constraints can change the estimate

The unconstrained regression calculation allows any real slope. A stable causal AR(1) interpretation requires $|\phi|<1$. A fixed-gap mean-reverting scalar OU interpretation requires $0<\phi<1$.

The record $(0,2,1,3)$ gives $\bar u=1$, $\bar y=2$, $S_{uy}=-1$, and $S_{uu}=2$. Thus its fitted slope is $-1/2$. This is an admissible oscillating stable AR(1), but not an exact discretization of the positive mean-reverting scalar OU model. The observed data do not have to produce an admissible OU slope merely because we hoped to use that model.

For a closed allowed slope interval $[b_{\min},b_{\max}]$ that excludes the singular value one, the quadratic expression for $S(b)$ shows how to solve the constrained problem. The optimum is the unrestricted slope projected onto that interval; **then** the intercept and residual variance must be recomputed at that constrained slope. Changing only $\hat\phi$ while leaving the old $\hat\mu$ and variance in place generally does not give the constrained likelihood optimum.

If the allowed interval is open, such as $(0,1)$, and the best slope would occur at an excluded boundary, the supremum may not be attained. Numerical bounds such as $10^{-6}\le\phi\le1-10^{-6}$ create a different closed optimization problem. A result at such a bound should be reported as boundary-limited, not as an ordinary interior optimum.

### 7.11 Conditional and stationary full MLE are not identical calculations

The calculation above conditioned on $x_0$. A stationary Gaussian AR(1) model instead supplies the initial density

$$
X_0\sim N\left(\mu,\frac{\tau^2}{1-\phi^2}\right),\qquad |\phi|<1.
$$

Its full log-likelihood adds

$$
-\frac12\log\left(\frac{2\pi\tau^2}{1-\phi^2}\right)
-\frac{(x_0-\mu)^2(1-\phi^2)}{2\tau^2}
$$

to the conditional log-likelihood. This term depends on all three candidate parameters. It is not a constant merely because $x_0$ is observed. Its presence can change the optimum, so the conditional regression formulas should not be presented as universal formulas for stationary full MLE.

In the short worked record, the conditional fitted level is $2$ while the start is $0$. A stationary initial-density term would evaluate how plausible that start is under each candidate's stationary law. Conditional fitting deliberately does not ask that additional question.

### 7.12 Why an algebraic regression fit does not settle statistical uncertainty

The optimization objective has become a regression sum of squares, but the predictors are lagged values of a dependent process. An innovation $\varepsilon_i$ is independent of the previous state used to predict $X_i$, yet it influences $X_i$, which becomes the predictor for the next transition. Thus the full predictor array is not an externally fixed collection independent of every error.

This matters when discussing distributions of estimators, standard errors, and small-sample bias. The algebra above proves what maximizes the conditional Gaussian likelihood. It does not by itself justify every classical fixed-design regression inference formula. Nor does it make $\hat\phi$ unbiased: it is a ratio of random centered sums, and taking the expectation of a ratio is not generally the ratio of expectations.

The fitted slope also is not generally the sample Pearson correlation. If $r_{uy}$ denotes that correlation and $s_u,s_y$ are sample standard deviations calculated with the same divisor, then

$$
\hat b=r_{uy}\frac{s_y}{s_u}.
$$

The population stationary relation $\operatorname{Corr}(X_i,X_{i+1})=\phi$ does not force the two finite observed arrays to have equal standard deviations.

The exact lesson is limited but powerful: **under the stated conditional Gaussian AR(1) model, estimating the one-step mean reduces to a fully derivable least-squares problem**. Irregular observation times will preserve the likelihood principle while changing the form of that regression.

<p class="source-note">Source connection: <a href="#ref-l8">MIT Lecture 8</a> covers likelihood-based time-series estimation. The conditional normal equations, numerical examples, degeneracies, and constrained-slope calculation are derived explicitly above.</p>
