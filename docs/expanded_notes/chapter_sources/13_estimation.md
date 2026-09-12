## 13. Keep the observation gaps: exact OU estimation {#estimation}

### 13.1 What are the data, and what are the unknown parameters?

Suppose the quantity we are observing follows an OU model

$$
dX_t=\theta(\mu-X_t)dt+\sigma dB_t,
\qquad \theta>0,\quad\sigma>0.
$$

Here $\mu$ is the return level, $\theta$ the return speed, and $\sigma$ the continuous-time noise scale. The actual values of these parameters are unknown. We observe the process at times

$$
t_0<t_1<\cdots<t_n
$$

and record numerical values

$$
x_0,x_1,\ldots,x_n.
$$

The data are the **time–value pairs**, not just the values without timestamps. Define each observed time gap by

$$
\Delta_i=t_i-t_{i-1}>0,\qquad i=1,\ldots,n.
$$

There are $n+1$ recorded values but only $n$ transitions between consecutive values. The gaps are known numbers calculated from the timestamps. We are not estimating them in this model. We are estimating

$$
\eta=(\theta,\mu,\sigma).
$$

For example, one small dataset is

| Observation index | Time $t_i$ | Observed value $x_i$ | Gap from previous observation |
|---:|---:|---:|---:|
| 0 | 0.0 | 29.0 | Not applicable |
| 1 | 0.2 | 27.5 | 0.2 |
| 2 | 1.0 | 26.0 | 0.8 |
| 3 | 1.3 | 25.8 | 0.3 |

These four observations are a teaching example, not a claim that four readings provide reliable estimation of all OU parameters. Every candidate parameter triple will be evaluated on **exactly this same dataset**.

The maximization does not require us to know the true parameter before starting. As in the coin example, we ask what probability density the observed record would have under each hypothetical candidate. We then search for the candidate giving the largest likelihood, or equivalently the smallest negative log-likelihood.

### 13.2 The transition law needed for the likelihood

For a positive gap $\Delta_i$, the OU solution has the restart form

$$
X_{t_i}=\mu+e^{-\theta\Delta_i}(X_{t_{i-1}}-\mu)+\varepsilon_i,
$$

where the fresh noise term is

$$
\varepsilon_i=\sigma\int_{t_{i-1}}^{t_i}e^{-\theta(t_i-s)}\,dB_s.
$$

The weight is deterministic. A weighted Brownian integral has mean zero and variance given by the integral of the squared weight, so

$$
\begin{aligned}
\operatorname{Var}(\varepsilon_i)
&=\sigma^2\int_{t_{i-1}}^{t_i}e^{-2\theta(t_i-s)}ds\\
&=\frac{\sigma^2}{2\theta}(1-e^{-2\theta\Delta_i}).
\end{aligned}
$$

Different observation intervals use disjoint Brownian increments. Under the sampling assumptions below, their fresh noise terms are independent conditional on the recorded observation grid. Conditional on the current observed value, the next value therefore has the exact Gaussian distribution

$$
X_{t_i}\mid X_{t_{i-1}}=x_{i-1}
\sim N(m_i,v_i),
$$

with

$$
\phi_i=e^{-\theta\Delta_i},\qquad
m_i=\mu+\phi_i(x_{i-1}-\mu),\qquad
v_i=\frac{\sigma^2}{2\theta}(1-\phi_i^2).
$$

The symbols $\phi_i,m_i,v_i$ are functions calculated from a candidate parameter triple, the current observed value and the known gap. They are not extra free parameters to fit independently for every transition. All transitions share the same $\theta,\mu,\sigma$.

### 13.3 State the observation-time assumption before multiplying densities

The transition formula applies directly when observation times are fixed in advance. It also applies conditional on a realized random observation grid when the whole grid is independent of the entire OU path, including its initial value. In this second situation, first imagine drawing the observation times independently, then observing the process at those times.

Independence of the sampling design is not the same as the Markov property of the process. To see the distinction, suppose we decide to record the next observation when a continuous process first reaches a specified threshold $a$. At that observation time, its recorded value is $a$ by construction. Conditional on this sampling rule and the hitting time, the observed endpoint is not an unrestricted Gaussian variable over the real line. The fact that the underlying process is Markov does not remove the information introduced by the observation rule.

The likelihood below therefore conditions on **fixed or path-independent observation times**. State-dependent observation schedules, noisy measurements and hidden states require the appropriate observation model. Also, if the distribution of the observation grid itself depends on unknown parameters, a full joint likelihood for times and values can contain additional parameter information in the time distribution. The conditional likelihood for values given times is a specifically chosen objective, not automatically the full likelihood of every possible data-collection scheme.

For the rest of this chapter, the initial value is either fixed or independent of future Brownian increments, the observed values have no additional measurement noise, and the observation grid satisfies the stated independence condition.

### 13.4 Derive the conditional likelihood rather than multiply marginals

The Gaussian transition density evaluated at the actual next observation is

$$
g_i(x_i\mid x_{i-1};\eta,\Delta_i)
=\frac{1}{\sqrt{2\pi v_i}}
\exp\!\left[-\frac{(x_i-m_i)^2}{2v_i}\right].
$$

The residual $r_i=x_i-m_i$ is the difference between the recorded value and this candidate model's conditional mean. Different candidate parameters change both the predictions and their claimed variances.

Condition on the starting value $x_0$ and all observation times. The probability chain rule writes the joint density of the remaining observations as

$$
\begin{aligned}
L_c(\eta)
={}&p_\eta(x_1\mid x_0,\text{times})\\
&\times p_\eta(x_2\mid x_0,x_1,\text{times})\times\cdots\\
&\times p_\eta(x_n\mid x_0,\ldots,x_{n-1},\text{times}).
\end{aligned}
$$

The OU Markov property shortens each conditional distribution to dependence on the immediately preceding value and the elapsed gap. Therefore

$$
\boxed{L_c(\eta)=\prod_{i=1}^n g_i(x_i\mid x_{i-1};\eta,\Delta_i).}
$$

This is a product of **conditional** densities. The consecutive observed values are not independent. Multiplying their stationary marginal densities would describe a different statistical objective that discards the temporal dependence needed to estimate return speed properly.

Although we condition on $x_0$, we do not erase it from the calculation: it is used in $m_1$. Conditioning removes the need to assign a probability density to the supplied starting value. It does not remove the first transition from the dataset.

### 13.5 Take logs and keep the variance normalization

Take the logarithm of one transition factor:

$$
\begin{aligned}
\log g_i
&=\log\left((2\pi v_i)^{-1/2}\right)
  +\log\left(\exp[-r_i^2/(2v_i)]\right)\\
&=-\frac12\log(2\pi v_i)-\frac{r_i^2}{2v_i}.
\end{aligned}
$$

The first equality uses the rule that the logarithm of a product is the sum of logarithms. The second uses $\log(a^c)=c\log a$ for positive $a$, and $\log(e^z)=z$.

Logs turn the full product into a sum. Negating it gives the negative log-likelihood, abbreviated NLL:

$$
\boxed{\operatorname{NLL}(\eta)
=\frac12\sum_{i=1}^n\left[
\log(2\pi v_i)+\frac{(x_i-m_i)^2}{v_i}\right].}
$$

Minimizing this is equivalent to maximizing the original likelihood because the log is increasing and negation reverses the ordering.

Both terms matter. The squared residual divided by $v_i$ measures the discrepancy relative to the uncertainty claimed by the candidate. The term $\log v_i$ accounts for the density becoming more spread out when variance increases. Removing that term would reward arbitrarily large variances: all residual penalties would approach zero as $\sigma^2$ grew. The proper likelihood balances fit with the width of the prediction distribution.

The constant $n\log(2\pi)/2$ can be omitted during optimization because it does not depend on any candidate parameter. But $\log v_i$ cannot generally be omitted: $v_i$ depends on $\theta$ and $\sigma$.

### 13.6 Evaluate one candidate on a concrete dataset

Use the four readings from §13.1 and the candidate values $\theta=1$, $\mu=25$, $\sigma=2$. For the first gap $\Delta_1=0.2$,

$$
\phi_1=e^{-0.2}\approx0.81873075,
$$

$$
m_1=25+0.81873075(29-25)\approx28.27492301,
$$

and

$$
v_1=\frac{4}{2}(1-e^{-0.4})\approx0.65935991.
$$

The actual next value is 27.5, so the residual is $27.5-28.27492301\approx-0.77492301$. The first negative-log contribution is

$$
\frac12\left[\log(2\pi\times0.65935991)
+\frac{(-0.77492301)^2}{0.65935991}\right]
\approx1.16606583.
$$

Now repeat the calculation for each actual transition. Importantly, the second prediction starts from the **observed** 27.5, not from the model's previous predicted mean 28.2749. The likelihood conditions on the record that actually occurred.

| Transition | Gap | Predicted mean $m_i$ | Variance $v_i$ | Residual $r_i$ | NLL contribution |
|---|---:|---:|---:|---:|---:|
| $29.0\to27.5$ | 0.2 | 28.274923 | 0.659360 | −0.774923 | 1.166066 |
| $27.5\to26.0$ | 0.8 | 26.123322 | 1.596207 | −0.123322 | 1.157518 |
| $26.0\to25.8$ | 0.3 | 25.740818 | 0.902377 | 0.059182 | 0.869518 |

Adding the three terms gives $\operatorname{NLL}(1,25,2)\approx3.19310101$. The conditional likelihood density is consequently $e^{-3.19310101}\approx0.04104439$.

Keep $\mu=25$, $\sigma=2$ fixed and compare several candidate speeds on the same data:

| Candidate $\theta$ | Conditional NLL | Conditional likelihood density |
|---:|---:|---:|
| 0.5 | 3.978405 | 0.018715 |
| 1.0 | 3.193101 | 0.041044 |
| 2.0 | 2.461343 | 0.085320 |

Among **these three candidates with the other two parameters fixed**, $\theta=2$ has the largest likelihood. This comparison does not establish that 2 is the MLE over all positive speeds, much less the joint MLE over all three parameters. We have evaluated three points, not solved the entire optimization problem.

### 13.7 Why the parameters cannot simply be estimated separately

The conditional means depend on $\theta$ and $\mu$, while the variances depend on $\theta$ and $\sigma$. Changing the return speed therefore changes both the expected return toward the level and the amount of accumulated noise. This couples the fitting problem.

Nevertheless, we can simplify the optimization without pretending that the parameters are independent. **Profiling** means: fix one parameter temporarily, find the best values of the others at that fixed choice, and then compare those optimized scores across the first parameter.

For example, at candidate $\theta=1$, we can find the best $\mu$ and $\sigma$ conditional on that speed. At candidate $\theta=2$, we must generally find a different best $\mu$ and $\sigma$. The result is a function of $\theta$ alone, but it still incorporates the adjustment of the other parameters.

We will derive this reduction for the conditional likelihood with unrestricted real $\mu$ and positive $\sigma$. A parameter-dependent initial-density factor, extra constraints, a prior, or a different observation model may change the formulas.

### 13.8 Rewrite the likelihood at a fixed candidate speed

Fix $\theta>0$ and calculate

$$
\phi_i=e^{-\theta\Delta_i},\qquad
b_i=1-\phi_i,\qquad
c_i=\frac{1-\phi_i^2}{2\theta},\qquad
d_i=x_i-\phi_i x_{i-1}.
$$

All four quantities are now known for this fixed candidate $\theta$. Because every gap is positive and $\theta>0$, we have $0<\phi_i<1$, $b_i>0$, and $c_i>0$.

Expand the conditional mean:

$$
m_i=\mu+\phi_i(x_{i-1}-\mu)
=\phi_i x_{i-1}+(1-\phi_i)\mu.
$$

Therefore

$$
x_i-m_i=d_i-b_i\mu,
\qquad
v_i=\sigma^2c_i.
$$

The likelihood has become a weighted fitting problem for the linear expression $b_i\mu$. Define

$$
S_\theta(\mu)=\sum_{i=1}^n\frac{(d_i-b_i\mu)^2}{c_i}.
$$

This is a weighted sum of squared residuals. The weights $1/c_i$ are positive and account for different amounts of transition noise at different gaps. Substituting into NLL gives

$$
\operatorname{NLL}(\theta,\mu,\sigma)
=\frac n2\log(2\pi)
+\frac n2\log\sigma^2
+\frac12\sum_i\log c_i
+\frac{S_\theta(\mu)}{2\sigma^2}.
$$

At the fixed $\theta$, the $c_i$ are constants for the remaining optimization. The only dependence on $\mu$ is through $S_\theta(\mu)$.

### 13.9 Derive the profiled level estimate

Differentiate $S_\theta(\mu)$ with respect to $\mu$. For one term, the chain rule gives

$$
\frac{d}{d\mu}(d_i-b_i\mu)^2
=2(d_i-b_i\mu)(-b_i).
$$

Hence

$$
\begin{aligned}
S_\theta'(\mu)
&=-2\sum_i\frac{b_i(d_i-b_i\mu)}{c_i}\\
&=-2\sum_i\frac{b_i d_i}{c_i}
 +2\mu\sum_i\frac{b_i^2}{c_i}.
\end{aligned}
$$

Set this derivative equal to zero and move the first sum to the other side:

$$
\mu\sum_i\frac{b_i^2}{c_i}=\sum_i\frac{b_i d_i}{c_i}.
$$

Dividing by the positive denominator produces

$$
\boxed{\widehat\mu(\theta)
=\frac{\sum_i b_i d_i/c_i}{\sum_i b_i^2/c_i}.}
$$

This is the unique minimum in $\mu$, not merely a stationary point, because

$$
S_\theta''(\mu)=2\sum_i\frac{b_i^2}{c_i}>0.
$$

The notation $\widehat\mu(\theta)$ is important. It is a best-fitting level **as a function of the speed being tried**. It becomes part of a joint estimate only after a speed is selected. Also, this expression is not generally the arithmetic mean of the recorded $x_i$: both the preceding values and gaps affect it.

Write

$$
S(\theta)=S_\theta(\widehat\mu(\theta))
$$

for the minimized weighted squared error at that speed.

### 13.10 Derive the profiled noise variance and the remaining objective

For fixed $\theta$ and its best level, let $s=\sigma^2>0$. The terms of NLL that depend on $s$ are

$$
\frac n2\log s+\frac{S(\theta)}{2s}.
$$

Differentiate:

$$
\frac{d}{ds}\left(\frac n2\log s+\frac{S}{2s}\right)
=\frac{n}{2s}-\frac{S}{2s^2}
=\frac{ns-S}{2s^2}.
$$

When $S>0$, the derivative is negative for $s<S/n$, zero at $s=S/n$, and positive for $s>S/n$. Thus the unique minimum is

$$
\boxed{\widehat\sigma^2(\theta)=\frac{S(\theta)}{n}.}
$$

The divisor is $n$, the number of conditional transitions. This is a likelihood estimate, not an unbiasedness correction involving a regression degrees-of-freedom divisor.

At this estimate, $S/(2\widehat\sigma^2)=n/2$. Substitute back into the full expression:

$$
\boxed{
\operatorname{NLL}_{\mathrm{profile}}(\theta)
=\frac12\left\{
 n\left[\log(2\pi)+1+\log\left(\frac{S(\theta)}n\right)\right]
 +\sum_i\log c_i
\right\}.}
$$

We have reduced a three-parameter problem to a one-dimensional numerical search over $\theta>0$. Once a minimizing speed is found, plug it into the derived formulas to obtain the associated level and noise estimates.

The case $S(\theta)=0$ needs separate treatment. The fitted means then match every observed transition exactly at that candidate speed. As $\sigma^2\downarrow0$, the negative log-likelihood tends to $-\infty$, so there is no ordinary positive interior noise-variance MLE at that speed. A numerical variance floor changes the permitted model and must not be hidden as though it were the unconstrained MLE.

Profiling also does not prove that the remaining objective has one global minimum or that its minimum lies in the interior. Boundaries and weak identification can still matter. A one-dimensional optimizer is a tool for finding the remaining optimum, not a proof that every run has found the correct global answer.

### 13.11 Verify a profile calculation with the example data

At $\theta=1$, the three observed gaps give approximately

$$
b=(0.18126925,\;0.55067104,\;0.25918178),
$$

$$
c=(0.16483998,\;0.39905174,\;0.22559418),
$$

and

$$
d=(3.75680816,\;13.64345349,\;6.53872626).
$$

Substituting these into the weighted formula gives $\widehat\mu(1)\approx24.24077917$. The minimized weighted error is $S(1)\approx2.97204085$. Therefore

$$
\widehat\sigma^2(1)\approx0.99068028,
\qquad
\widehat\sigma(1)\approx0.99532923.
$$

The profiled NLL is approximately 2.13753954. This is lower than 3.19310101, the value obtained at the same speed with the arbitrarily fixed $\mu=25$ and $\sigma=2$. That is what profiling should do: it improves or preserves the score by choosing the best remaining parameters at that speed.

For comparison:

| Fixed candidate $\theta$ | Best level $\widehat\mu(\theta)$ | Best noise variance $\widehat\sigma^2(\theta)$ | Profile NLL |
|---:|---:|---:|---:|
| 0.5 | 21.819132 | 1.498571 | 3.036276 |
| 1.0 | 24.240779 | 0.990680 | 2.137540 |
| 2.0 | 25.404013 | 0.290023 | −0.177628 |

A negative NLL is not an error. This likelihood uses **densities**, which can exceed 1; a positive log-density therefore gives a negative negative-log-density. These remain only three profiled candidate evaluations, not a claim that the global joint MLE has been found.

### 13.12 Keep numerical evaluation faithful to the formula

For a small positive number $z$, $e^{-z}$ is very close to 1. Directly subtracting it from 1 can lose significant digits in floating-point arithmetic. The function `expm1(w)` evaluates $e^w-1$ accurately near zero. Thus compute

$$
b_i=-\operatorname{expm1}(-\theta\Delta_i),
\qquad
c_i=\frac{-\operatorname{expm1}(-2\theta\Delta_i)}{2\theta}.
$$

An equivalent expression is

$$
c_i=\Delta_i\frac{-\operatorname{expm1}(-2z_i)}{2z_i},
\qquad z_i=\theta\Delta_i.
$$

The ratio tends to 1 as $z_i\to0$, making the short-gap limit $c_i\approx\Delta_i$ explicit. This is numerical care for the same mathematical objective, not a new estimator.

Before evaluating a candidate, check that the data and parameters are finite, the number of values is one more than the number of gaps, every gap is positive, and $\theta,\sigma$ satisfy the parameter restrictions. If the data contain duplicated timestamps, the positive-gap transition-density formula cannot simply be applied with a zero variance.

A practical search may optimize $\alpha=\log\theta$ and recover $\theta=e^\alpha$. This enforces positivity, but does not remove the need to examine search bounds and convergence. The bounds have scientific meaning through $\theta$'s units and the observation horizon. Expanding the search range can reveal whether a reported optimum was merely sitting at an imposed limit.

A fully executable scalar implementation is supplied alongside these notes. Its functions evaluate the direct and profiled objectives; they deliberately do not label a few grid evaluations as a global optimizer. Their validation checks compare the direct objective with the profiled expression at identical parameter values.

### 13.13 Compare exact, mean-gap and Euler likelihoods

The three methods share a Gaussian-likelihood structure but use different transition means and variances. Keeping these differences explicit prevents two separate approximation errors from being confused.

**Exact gap-aware likelihood.** Use each observed $\Delta_i$ in

$$
m_i^{\mathrm{exact}}=\mu+e^{-\theta\Delta_i}(x_{i-1}-\mu),
\qquad
v_i^{\mathrm{exact}}=\frac{\sigma^2}{2\theta}(1-e^{-2\theta\Delta_i}).
$$

These are the exact endpoint transitions of the stated OU model under the observation assumptions.

**Mean-gap approximation.** Calculate the realized mean gap

$$
\overline\Delta=\frac1n\sum_i\Delta_i=\frac{t_n-t_0}{n}.
$$

Then replace every $\Delta_i$ in the exact formulas by this single number. The method uses one coefficient $e^{-\theta\overline\Delta}$ and one innovation variance for all transitions. It effectively treats the irregular record as though successive readings were evenly spaced. The earlier notes label this mean-gap approach PFML. Its relevant feature here is the replacement of the actual gaps, not the acronym.

For our small dataset, the actual gaps are $(0.2,0.8,0.3)$, while the average is $1.3/3\approx0.4333$. The second transition, which had 0.8 time units to return toward the level, would be treated as though it had only about 0.4333. The first transition would be treated as having more than twice its actual time. This changes both expected decay and noise accumulation.

**Euler pseudo-likelihood.** Keep each actual gap, but approximate the OU drift as remaining at its starting value throughout that gap:

$$
X_{t_i}\approx x_{i-1}+\theta(\mu-x_{i-1})\Delta_i
+\sigma(B_{t_i}-B_{t_{i-1}}).
$$

Because the Brownian increment is $N(0,\Delta_i)$, the approximation gives

$$
m_i^{\mathrm{Euler}}=x_{i-1}+\theta(\mu-x_{i-1})\Delta_i,
\qquad
v_i^{\mathrm{Euler}}=\sigma^2\Delta_i.
$$

The mean-gap method changes the **recorded clock**; the Euler method approximates the **dynamics over each gap**. They are not the same approximation. When the gaps are all equal, the exact and mean-gap objectives agree identically, but Euler generally still differs from the exact transition.

Euler accuracy is controlled by the scale of $\theta\Delta_i$, not by the mean gap alone. Some large gaps can remain problematic even when the average is small. At a fixed Euler step $h$, the discrete autoregressive coefficient is $1-\theta h$. It is negative when $\theta h>1$, and has magnitude at least 1 when $\theta h\geq2$. The exact coefficient $e^{-\theta h}$ remains in $(0,1)$ for every positive $h$. This algebra shows how a coarse Euler scheme can introduce oscillating or unstable behavior absent from the exact mean-reverting transition.

### 13.14 Why averaging the gaps can systematically distort return speed

We can explain one important effect with a population calculation. This is an additional large-sample illustration with stronger assumptions, not a universal finite-sample claim.

Assume a stationary OU process with true parameters $\theta,\mu,\sigma$ and stationary variance $q=\sigma^2/(2\theta)$. Let successive gaps $D_i$ be independent, identically distributed positive variables, independent of the entire process, with finite mean $d=E[D_i]$. Also assume that the relevant sample means and cross-products converge to their population values as the observation horizon grows; this is the long-run averaging condition needed for the regression-limit argument.

Let $Y_i=X_{t_i}-\mu$. Stationarity and independent sampling give $E[Y_i]=0$ and $E[Y_i^2]=q$. Conditional on a particular gap $D_i$, the transition relation is

$$
Y_i=e^{-\theta D_i}Y_{i-1}+\varepsilon_i.
$$

The fresh noise has conditional mean zero. The independently chosen gap is independent of the preceding sampled value. Therefore

$$
\begin{aligned}
E[Y_{i-1}Y_i]
&=E[e^{-\theta D_i}Y_{i-1}^2]+E[Y_{i-1}\varepsilon_i]\\
&=E[e^{-\theta D_i}]\,E[Y_{i-1}^2]+0\\
&=qE[e^{-\theta D_i}].
\end{aligned}
$$

A constant-coefficient regression of $Y_i$ on $Y_{i-1}$ has population slope equal to covariance divided by predecessor variance. Its target is consequently

$$
\phi_* = \frac{qE[e^{-\theta D_i}]}{q}=E[e^{-\theta D_i}].
$$

If the mean-gap estimator interprets this slope as $e^{-\theta_*d}$, then

$$
\boxed{\theta_*=-\frac{\log E[e^{-\theta D_i}]}{d}.}
$$

The asterisk denotes the long-run target of this misspecified constant-gap fitting procedure. It is not the true speed by definition, and it is not the value of every finite-data estimate.

### 13.15 Derive the direction of that distortion

The function $g(x)=e^{-\theta x}$ is strictly convex for $\theta>0$, because

$$
g''(x)=\theta^2e^{-\theta x}>0.
$$

A differentiable convex function lies above its tangent at the mean $d$:

$$
g(x)\geq g(d)+g'(d)(x-d).
$$

Substitute $x=D_i$ and take expectations. The linear term disappears because $E[D_i-d]=0$. This gives Jensen's inequality in the form needed here:

$$
E[e^{-\theta D_i}]\geq e^{-\theta d}.
$$

The logarithm preserves the inequality, while multiplying by the negative factor $-1/d$ reverses it:

$$
\log E[e^{-\theta D_i}]\geq-\theta d
\quad\Longrightarrow\quad
\boxed{\theta_*\leq\theta.}
$$

When the gap distribution is nondegenerate, strict convexity makes the inequality strict. Shorter and longer gaps do not average out linearly inside an exponential. Their average surviving fraction is larger than the surviving fraction at the average gap. The constant-gap interpretation consequently infers a slower return speed.

For a numerical closed-form example, let $D$ have an exponential distribution with mean $d$, so its density is $f_D(u)=d^{-1}e^{-u/d}$ for $u>0$. Then

$$
\begin{aligned}
E[e^{-\theta D}]
&=\int_0^\infty e^{-\theta u}\frac1d e^{-u/d}\,du\\
&=\frac1d\int_0^\infty e^{-(\theta+1/d)u}\,du\\
&=\frac1d\cdot\frac1{\theta+1/d}\\
&=\frac1{1+\theta d}.
\end{aligned}
$$

Thus

$$
\theta_* =\frac{\log(1+\theta d)}{d}.
$$

At true speed $\theta=1$ and mean gap $d=0.5$, this gives $\theta_*=\log(1.5)/0.5\approx0.81093022$. It predicts a long-run distortion in this specified setting; it does not claim that every fitted value must be below 1 or equal to 0.81093.

The implied noise scale can also shift. The population residual variance of the best constant-coefficient regression is

$$
\begin{aligned}
\tau_*^2
&=E[(Y_i-\phi_*Y_{i-1})^2]\\
&=q-2\phi_*(q\phi_*)+\phi_*^2q\\
&=q(1-\phi_*^2).
\end{aligned}
$$

Mapping the fitted constant-gap model back to OU parameters gives

$$
\sigma_*^2=\frac{2\theta_*\tau_*^2}{1-\phi_*^2}
=2\theta_*q
=\sigma^2\frac{\theta_*}{\theta}.
$$

Thus the same limiting fit can preserve the stationary marginal variance $q$ while distorting the speed and continuous-time noise scale. Looking only at whether the fitted stationary histogram seems reasonable would not reveal all of the temporal misspecification.

### 13.16 Exact transition formulas do not guarantee unbiased estimates

An exact likelihood uses the correct transition density for the stated model and observation scheme. That is a statement about how the objective is specified. It is not the same as saying its maximizer has no finite-sample bias or perfectly calibrated uncertainty.

For an estimator $\widehat\theta$ computed from a random dataset, its bias is

$$
\operatorname{Bias}(\widehat\theta)=E[\widehat\theta]-\theta.
$$

Its mean squared error is $E[(\widehat\theta-\theta)^2]$, and its root mean squared error, RMSE, is the square root of that quantity. Write $\widehat\theta-\theta=(\widehat\theta-E[\widehat\theta])+(E[\widehat\theta]-\theta)$ and expand the square. The cross term has expectation zero, giving

$$
E[(\widehat\theta-\theta)^2]
=\operatorname{Var}(\widehat\theta)+\operatorname{Bias}(\widehat\theta)^2.
$$

This decomposition distinguishes systematic error from dispersion across datasets. A method can have low bias but high variability, or vice versa.

For an interval procedure that produces $[L(D),U(D)]$ from dataset $D$, coverage is the repeated-sampling probability

$$
P_\theta\big(L(D)\leq\theta\leq U(D)\big).
$$

A procedure labelled “95%” is intended to have coverage near 0.95 under its stated conditions. Using the exact Gaussian transition alone does not establish that a particular asymptotic interval achieves this target on a short, persistent sample. Simulation studies of bias, RMSE and coverage address different questions from verifying the transition formula.

The number of observations and the total observation horizon are also different pieces of information. Recording many closely spaced values can create a large $n$ while adjacent values remain strongly dependent: the stationary correlation over a gap is $e^{-\theta\Delta}$. More timestamps do not turn a short observed time span into the same experiment as many independent long-separated observations. Precise consistency and uncertainty claims require their own asymptotic assumptions; they are not inferred merely from the size of $n$.

### 13.17 Full initial-state likelihood and exact simulation

If we choose to model $X_0$ by the stationary law rather than condition on it, its density contributes another factor. With $q=\sigma^2/(2\theta)$,

$$
L_{\mathrm{full}}(\eta)
=\frac1{\sqrt{2\pi q}}
  \exp\!\left[-\frac{(x_0-\mu)^2}{2q}\right]L_c(\eta).
$$

The corresponding additional NLL term is

$$
\frac12\left[\log(2\pi q)+\frac{(x_0-\mu)^2}{q}\right].
$$

This depends on all three parameters through $\mu$ and $q$. It cannot be discarded as a constant merely because $x_0$ is an observed fixed number. Conditional and full likelihoods are both legitimate when matched to their modeling choices, but they are not numerically identical objectives. The profiling formulas derived above should not be reused without checking the additional term.

For simulation at the observation times, the exact transition provides a simpler alternative to hidden small Euler steps. Generate independent standard normal numbers $z_i\sim N(0,1)$ and set

$$
x_i=\mu+e^{-\theta\Delta_i}(x_{i-1}-\mu)
+\sqrt{\frac{\sigma^2}{2\theta}(1-e^{-2\theta\Delta_i})}\,z_i.
$$

Because shifting and scaling a standard normal gives the required conditional Gaussian law, this simulates the exact joint distribution on the selected grid. It does not reconstruct every unobserved point of a continuous sample path. Joining simulated endpoints with straight lines supplies a visual interpolation, not an exact conditional Brownian or OU bridge between them.

At the true parameter values, the standardized innovations $(X_{t_i}-m_i)/\sqrt{v_i}$ are independent standard normals conditional on the independent observation grid. Residuals formed after estimating parameters from the same data are not automatically an exactly independent standard-normal sample. The shared parameter fitting introduces additional dependence and finite-sample effects.

### 13.18 Reconstruct the entire estimation logic

A useful final check is to reproduce the argument as a sequence of mathematical decisions, not as a memorized objective. First identify the observations and their timestamps. Specify the OU model and how observation times were selected. Choose whether to condition on the starting value or model its distribution. For each candidate parameter triple, calculate each transition's mean and variance using its actual gap. Use the probability chain rule and the Markov property to form the product of conditional densities. Take logs, keep all parameter-dependent normalization terms, and optimize over the stated parameter space.

The profiled version makes one additional calculation: at each candidate speed, analytically minimize the weighted residual sum over the level and then fit the noise variance. It does not change the underlying conditional likelihood. Comparing it with mean-gap and Euler procedures then becomes a controlled comparison of precisely identified assumptions and approximations.

::: {.worked}
**Final checks to answer without looking at a formula card.**

If every gap equals $h$, the exact method and mean-gap method use identical means and variances for every candidate parameter, so their conditional likelihoods coincide. Euler need not coincide, because $1-\theta h$ is generally not $e^{-\theta h}$ and $\sigma^2h$ is generally not $q(1-e^{-2\theta h})$.

If all observed values and gaps are unchanged but the candidate $\theta$ changes, the data have not changed. The model's forecast decay and transition variances have changed, so the likelihood can change. This is exactly the same candidate-parameter reasoning as the coin example, now applied to dependent continuous observations.

If a routine returns a parameter at a search bound, that is a numerical outcome to investigate and report. It is not evidence that the unconstrained likelihood has an interior maximum at that value. If a simulation verifies a transition variance, that verifies the model calculation being checked; it does not by itself prove that an estimator is unbiased or that confidence intervals have the advertised coverage.
:::

<p class="source-note">Source connection: the OU transition follows from <a href="#ref-ss">Särkkä and Solin, Example 6.2</a>. Random sampling and estimation are studied in <a href="#ref-am">Aït-Sahalia and Mykland (2003)</a>. The numerical examples, conditional profiling calculation and explicit population mean-gap illustration are derived here under the assumptions stated in the text; they are not presented as new research estimators or as completed finite-sample research results.</p>
