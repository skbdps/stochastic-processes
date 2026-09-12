## 8. Refine the clock {#continuous}

### 8.1 Step count is not elapsed time

A random walk $W_k=S_1+\cdots+S_k$ is indexed by the number of steps. The symbol $k=100$ means one hundred updates; it says nothing about whether those updates fill one second or one century.

To connect a numerical model to physical time, specify a tick duration $\delta>0$. At grid time $t=k\delta$, the number of completed ticks is $k=t/\delta$. If the steps have mean $m_\delta$ and variance $v_\delta$, independence gives

$$
E[W_k]=km_\delta=\frac{t}{\delta}m_\delta,
\qquad
\operatorname{Var}(W_k)=kv_\delta=\frac{t}{\delta}v_\delta.
$$

The subscripts on $m_\delta$ and $v_\delta$ allow the step distribution to change when we refine the tick size. That change is not a nuisance; it is necessary if different numerical resolutions are meant to approximate the same continuous-time model.

A genuinely discrete model can have a physically meaningful tick duration, such as one transaction or one generation. There is nothing wrong with that. Here the goal is different: we want to make the numerical tick arbitrarily small without changing the intended mean and variance over a fixed physical duration.

### 8.2 Derive the drift scaling

Suppose the desired deterministic drift rate is $a$ units of state per unit of time. Over elapsed time $t$, the expected accumulated change should be $at$.

Our grid model gives expected change $(t/\delta)m_\delta$. Equating it with $at$ and canceling $t$ gives

$$
\frac{m_\delta}{\delta}=a
\quad\Longrightarrow\quad
\boxed{m_\delta=a\delta.}
$$

Halving the tick duration doubles the number of ticks within a fixed interval, so the mean change per tick must halve. This is the ordinary rate-times-duration relationship from elementary calculus.

For example, with drift rate $a=3$ and duration $t=2$, two ticks of length one each contribute mean three, while two hundred ticks of length $0.01$ each contribute mean $0.03$. Both give total expected change six.

### 8.3 Derive the noise scaling separately

Suppose the desired variance over elapsed time $t$ is $\sigma^2t$, where $\sigma>0$ is a noise scale. The grid model gives variance $(t/\delta)v_\delta$. Equating the two yields

$$
\frac{t}{\delta}v_\delta=\sigma^2t
\quad\Longrightarrow\quad
\boxed{v_\delta=\sigma^2\delta.}
$$

Variance is proportional to duration. Standard deviation is its square root, so the random step's standard deviation must be $\sigma\sqrt\delta$, not $\sigma\delta$.

To implement both mean and noise scaling, choose i.i.d. variables $\xi_i$ with mean zero and variance one and define a tick's change by

$$
\boxed{S_{i,\delta}=a\delta+\sigma\sqrt\delta\,\xi_i.}
$$

The first term is deterministic. The second has mean zero and variance $\sigma^2\delta$, since multiplying by $\sigma\sqrt\delta$ multiplies variance by its square. At $t=k\delta$, the accumulated process has mean $at$ and variance $\sigma^2t$.

For unit noise scale and one unit of elapsed time:

| Tick length $\delta$ | Number of ticks | Step standard deviation $\sqrt\delta$ | Step variance | Total variance |
|---:|---:|---:|---:|---:|
| $1$ | $1$ | $1$ | $1$ | $1$ |
| $1/4$ | $4$ | $1/2$ | $1/4$ | $1$ |
| $1/100$ | $100$ | $1/10$ | $1/100$ | $1$ |

The distribution shapes may still differ at finite resolution, but their intended total variance is calibrated consistently.

### 8.4 Why the square root is not an arbitrary convention

Suppose instead that random step amplitude were proportional to $\delta^\alpha$, so a centered step were $\sigma\delta^\alpha\xi_i$. Its variance would be $\sigma^2\delta^{2\alpha}$. Summing $t/\delta$ independent steps would produce

$$
\operatorname{Var}(W_{t/\delta})
=\sigma^2t\delta^{2\alpha-1}.
$$

For a finite positive variance limit in this finite-variance diffusion scaling, the exponent must be zero: $2\alpha-1=0$, hence $\alpha=1/2$.

If $\alpha=1$, noise is scaled like ordinary drift. Then the total variance is $\sigma^2t\delta$, which tends to zero: the random fluctuations disappear at fixed times in mean square. If $\alpha=0$, the step amplitude is not reduced at all, and the total variance grows like $1/\delta$, rather than staying on the intended finite diffusion scale.

This variance calculation selects a regime; it is not a theorem that every possible stochastic scaling limit must be Brownian. Heavy-tailed steps, dependence between steps, or other scaling choices can lead to other limits. Our construction assumes independent finite-variance steps and aims at a continuous Gaussian-noise model.

The dimensions also agree. A drift rate has units “state divided by time,” so $a\delta$ has state units. The noise scale has units “state divided by square root of time,” so $\sigma\sqrt\delta$ also has state units. Adding the two is dimensionally consistent.

### 8.5 Construct a whole interpolated path

Take independent fair steps $\xi_i\in\{-1,1\}$ and let $W_k=\sum_{i=1}^k\xi_i$. Choose resolution $n$ ticks per unit time, so $\delta=1/n$. Define a process on the grid by

$$
Z_n(k/n)=\frac{W_k}{\sqrt n}.
$$

Equivalently, each tick has an effective step $\xi_i/\sqrt n$ of variance $1/n$. The normalization can be placed in each step or outside the accumulated sum; those constructions are identical.

To obtain a value between grid times, join successive grid positions by a straight line. If $t=(k+r)/n$ for $0\le r<1$, the linear interpolation is

$$
\boxed{Z_n(t)=\frac{W_k+r\xi_{k+1}}{\sqrt n}.}
$$

At $r=0$ this gives the left endpoint $W_k/\sqrt n$. At $r=1$ it gives the right endpoint $W_{k+1}/\sqrt n$. The interpolated sample path is continuous for each finite $n$, but it is piecewise linear and not yet a Brownian path.

### 8.6 Check finite-resolution properties rather than assign limiting ones early

At a grid point $t=k/n$, the variance is $k/n=t$. Between grid points, independence of $W_k$ and $\xi_{k+1}$ gives

$$
\operatorname{Var}(Z_n(t))
=\frac{\operatorname{Var}(W_k)+r^2\operatorname{Var}(\xi_{k+1})}{n}
=\frac{k+r^2}{n}.
$$

The target variance $t$ equals $(k+r)/n$, not $(k+r^2)/n$. Their difference is

$$
t-\operatorname{Var}(Z_n(t))=\frac{r-r^2}{n},
$$

which is at most $1/(4n)$ because $r-r^2=r(1-r)$ has maximum $1/4$. Thus the discrepancy vanishes with refinement, but it is not zero at every finite resolution.

For $n=4$ and $t=1/8$, we have $k=0$ and $r=1/2$. The interpolated variance is $(1/4)/4=1/16$, while $t=1/8$. This concrete discrepancy prevents us from claiming that the interpolated walk already has exact Brownian variance at all times.

Likewise, two disjoint subintervals lying inside one grid cell share the same random slope $\xi_{k+1}\sqrt n$. Their increments are not independent. Increments over disjoint grid-aligned blocks are independent because they use disjoint steps. The distinction between grid points and interpolated points must be respected before taking a limit.

### 8.7 What the ordinary central limit theorem establishes

The central limit theorem for i.i.d. zero-mean unit-variance steps states

$$
\frac{W_k}{\sqrt k}\xrightarrow{d}N(0,1).
$$

At time one, $Z_n(1)=W_n/\sqrt n$, so the theorem directly gives convergence to $N(0,1)$.

For fixed $t>0$, let $k=\lfloor nt\rfloor$, the greatest integer not exceeding $nt$. Then

$$
\frac{W_k}{\sqrt n}
=\frac{W_k}{\sqrt k}\sqrt{\frac kn}.
$$

The first factor approaches a standard normal distribution; the deterministic second factor approaches $\sqrt t$. The resulting distribution approaches $N(0,t)$. The interpolation contribution has magnitude at most $1/\sqrt n$ for our $\pm1$ steps, so it disappears as $n$ grows. Thus the interpolated $Z_n(t)$ has the same fixed-time limit.

This establishes a **one-time distributional limit**. It does not mean a particular walk outcome becomes a Gaussian number, nor does it by itself prove convergence of the entire path. A probability distribution is a description across repeated outcomes, not a shape drawn by one trajectory.

### 8.8 Why convergence of paths is a stronger statement

To describe a process, we need joint behavior at different times. For ordered times $t_1<t_2<\cdots<t_m$, disjoint increments of the scaled walk come from separate blocks of independent steps. Their variances approach the interval lengths. This suggests jointly independent Gaussian increments in the limit, and therefore the desired joint distributions of positions obtained by adding those increments.

Even matching all finite collections of times is not, on its own, a complete proof about continuous functions under a path topology. One must also control the paths' oscillations so that probability does not concentrate on increasingly irregular behavior invisible to a finite observation set.

The functional central limit theorem, also called Donsker's theorem, supplies the stronger result: on a fixed bounded time interval, these normalized, linearly interpolated i.i.d. finite-variance walks converge in distribution as random continuous functions to Brownian motion. The usual topology measures the maximum difference between functions over the whole interval.

This theorem is stated here, not proved in its full functional-analytic generality. What we have derived is the normalization, the finite-time limits, and the reason a further path-level theorem is required. None of the later OU calculations will pretend that a one-variable CLT proved an entire process construction.

### 8.9 Understand a process as a random function

For one experiment outcome $\omega$, a stochastic process assigns a trajectory $t\mapsto X_t(\omega)$. There are three distinct operations.

If you fix the time $t$ and let the experiment vary, $X_t$ is a random variable. If you fix the outcome $\omega$ and vary time, $X_t(\omega)$ is an ordinary function of $t$. If you fix both, $X_t(\omega)$ is one number.

For a three-step fair walk, the outcome could be the step sequence $(1,-1,1)$. Its path is the function on times $0,1,2,3$ taking values $0,1,0,1$. The outcome is not the same object as its value at time three, which is just the number one. The random variable $W_3$ is the rule that reads the final value from whichever outcome is realized.

One can choose the sample space to be a space of paths itself. Then an event is a set of paths, such as “paths ending above zero” or “paths that stay nonnegative throughout the interval.” Assigning probabilities to such sets is the path-space viewpoint.

### 8.10 Marginal distributions do not determine a process

A marginal distribution describes the state at one time without describing its relationship to states at other times. Two continuous processes can have the same marginal distribution at every time and still have very different dependence.

Let $Z,Z_1,Z_2$ be independent standard normal variables. Define

$$
U_t=Z,
\qquad
V_t=Z_1\cos t+Z_2\sin t.
$$

For every fixed $t$, $U_t\sim N(0,1)$. The variable $V_t$ is a linear combination of independent Gaussians, so it is Gaussian with mean zero and variance $\cos^2t+\sin^2t=1$. Thus $V_t\sim N(0,1)$ too.

But $U_t$ is a constant path: once its value is drawn, it never changes. The second process generally varies with time. Their covariances expose the difference:

$$
\operatorname{Cov}(U_s,U_t)=1,
$$

whereas

$$
\begin{aligned}
\operatorname{Cov}(V_s,V_t)
&=\cos s\cos t+\sin s\sin t\\
&=\cos(t-s).
\end{aligned}
$$

Therefore giving a distribution such as $X_t\sim N(0,t)$ for every $t$ would not, by itself, specify Brownian motion. We also need the joint relationship between increments and path continuity.

### 8.11 Numerical ticks and observation gaps are not interchangeable

The symbol $\delta$ in this chapter is a numerical resolution used to build an approximation. Later, $\Delta_i=t_i-t_{i-1}$ will denote an actual gap between recorded observations. A fine simulation grid and a sparse observation grid serve different purposes.

For a model with an exact transition distribution, such as OU, one can simulate directly between the required observation times. There is no need to insert a tiny Euler grid merely to obtain correct endpoint observations. The continuum model describes the process between recordings; its transition law already summarizes the effect of the unobserved interval.

::: {.worked}
**Worked check.** You want drift rate $a=2$, noise scale $\sigma=3$, and total time $T=1$. With tick size $\delta=1/100$, each step should be $0.02+0.3\xi_i$, where $E[\xi_i]=0$ and $\operatorname{Var}(\xi_i)=1$. One hundred independent steps have mean $2$ and variance $100(0.3^2)=9$. Using $0.03\xi_i$ instead would give variance $0.09$, one hundred times too small. That error comes from scaling noise by $\delta$ instead of $\sqrt\delta$.
:::

<p class="source-note">Source connection: the path-space viewpoint and Brownian scaling construction are discussed in <a href="#ref-l17">MIT Lecture 17, sections 1–2</a>. The normalization and finite-interpolation corrections are calculated explicitly here. The functional limit theorem is identified as an additional theorem rather than silently inferred from a scalar CLT.</p>
