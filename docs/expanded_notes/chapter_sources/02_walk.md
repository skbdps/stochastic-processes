## 2. Random kicks, dependent positions {#walk}

### 2.1 A random variable gives one value; a process gives values over time

Imagine a running balance. At the end of each day, a random adjustment is added to it. An adjustment can be positive or negative. We want to distinguish two questions: what is the adjustment on a particular day, and what is the accumulated balance after several days?

A **random variable** is a rule that assigns a numerical value to the outcome of an experiment. A **stochastic process** is a collection of such random variables indexed by time. The notation $W_0,W_1,W_2,\ldots$ describes the uncertain balance at successive times. A realized list, such as $0,1,0,1,2$, is one sample path of that process.

For a simple example, toss a coin independently at every step. A head adds one; a tail subtracts one. Define the random step $S_i$ by

$$
P(S_i=1)=p,\qquad P(S_i=-1)=1-p.
$$

The parameter $p$ is the same for all steps. Independence means that any finite collection of steps has a joint probability equal to the product of their individual probabilities. Identical distribution means each step uses the same two-point probability law. These are assumptions of the model; they are not automatic properties of quantities called “noise.”

Start at the fixed position $W_0=0$. After one step $W_1=S_1$; after two, $W_2=S_1+S_2$; after three, $W_3=S_1+S_2+S_3$. In general,

$$
\boxed{W_k=\sum_{i=1}^k S_i},
\qquad
W_{k+1}=W_k+S_{k+1}.
$$

This process is a random walk. The index $k$ counts steps; it does not yet specify whether a step lasts a second or a year.

### 2.2 A complete small experiment makes the two layers visible

For three fair coin tosses, each of the eight possible step sequences has probability $1/8$. Here are all of them.

| Steps $(S_1,S_2,S_3)$ | Positions $(W_0,W_1,W_2,W_3)$ |
|---|---|
| $(1,1,1)$ | $(0,1,2,3)$ |
| $(1,1,-1)$ | $(0,1,2,1)$ |
| $(1,-1,1)$ | $(0,1,0,1)$ |
| $(1,-1,-1)$ | $(0,1,0,-1)$ |
| $(-1,1,1)$ | $(0,-1,0,1)$ |
| $(-1,1,-1)$ | $(0,-1,0,-1)$ |
| $(-1,-1,1)$ | $(0,-1,-2,-1)$ |
| $(-1,-1,-1)$ | $(0,-1,-2,-3)$ |

The third step is $1$ in four rows and $-1$ in four rows. Even after restricting attention to rows with the first two steps fixed, the third step is still equally likely to have either sign. That is the steps' independence.

The positions tell a different story. If you know $W_2=2$, the last position can only be $1$ or $3$. Without that information, it could also be $-1$ or $-3$. Knowledge of one position changes the distribution of a later position. Thus the positions are dependent even though the steps are independent.

There is no contradiction because $W_2$ and $W_3$ contain shared randomness:

$$
W_2=S_1+S_2,
\qquad
W_3=S_1+S_2+S_3.
$$

Both include $S_1$ and $S_2$. Only $S_3$ is new. This distinction between independent new disturbances and dependent accumulated states is the central structural idea of these notes.

### 2.3 Expectation: derive the center of the distribution

For a discrete random variable $U$, its expectation is the probability-weighted average of its possible values:

$$
E[U]=\sum_u uP(U=u),
$$

provided the sum is well-defined. For one step,

$$
\begin{aligned}
E[S_i]&=1\cdot p+(-1)\cdot(1-p)\\
&=p-1+p\\
&=2p-1.
\end{aligned}
$$

A fair walk has $p=1/2$, hence zero step mean. A walk with $p>1/2$ has positive step mean and tends to move upward on average.

Expectation is linear. To see why for two discrete variables, expand using their joint probabilities:

$$
\begin{aligned}
E[U+V]
&=\sum_{u,v}(u+v)P(U=u,V=v)\\
&=\sum_u u\sum_vP(U=u,V=v)
 +\sum_v v\sum_uP(U=u,V=v)\\
&=E[U]+E[V].
\end{aligned}
$$

Summing the joint probabilities over the unwanted variable gives the corresponding marginal probability. Notice that we did **not** replace a joint probability by a product; independence was not used. The same linearity extends to finite sums and constant multiples whenever the expectations exist.

Apply it to the walk:

$$
E[W_k]=E[S_1]+\cdots+E[S_k]=k(2p-1).
$$

More generally, if the i.i.d. steps have any distribution with finite mean $m_S$, then $E[W_k]=km_S$. The mean of an accumulated total is the accumulated mean of its contributions.

### 2.4 Variance: why independence matters here

Expectation locates a distribution's center. Variance measures its spread around that center:

$$
\operatorname{Var}(U)=E[(U-E[U])^2].
$$

Expand the square and use linearity:

$$
\begin{aligned}
\operatorname{Var}(U)
&=E[U^2-2U E[U]+(E[U])^2]\\
&=E[U^2]-2(E[U])^2+(E[U])^2\\
&=E[U^2]-(E[U])^2.
\end{aligned}
$$

For a $\pm1$ step, $S_i^2=1$ regardless of which outcome occurs, so $E[S_i^2]=1$. Therefore

$$
\begin{aligned}
\operatorname{Var}(S_i)
&=1-(2p-1)^2\\
&=1-(4p^2-4p+1)\\
&=4p(1-p).
\end{aligned}
$$

For $p=1/2$, the variance is $1$. At $p=0$ or $p=1$, there is no randomness in a step, and the variance is zero, as it should be.

For two variables, expanding the centered square gives

$$
\operatorname{Var}(U+V)
=\operatorname{Var}(U)+\operatorname{Var}(V)+2\operatorname{Cov}(U,V),
$$

where

$$
\operatorname{Cov}(U,V)
=E[(U-E[U])(V-E[V])]
=E[UV]-E[U]E[V]
$$

is their covariance. A covariance term records whether the deviations of the two variables tend to reinforce or offset one another.

If $U$ and $V$ are independent with finite second moments, their joint probability factors, giving $E[UV]=E[U]E[V]$. Hence their covariance is zero. Repeating the expansion for many variables gives

$$
\operatorname{Var}\left(\sum_{i=1}^k S_i\right)
=\sum_{i=1}^k\operatorname{Var}(S_i)
+2\sum_{i<j}\operatorname{Cov}(S_i,S_j).
$$

For independent steps, every off-diagonal covariance vanishes. If each step has variance $v_S$, we obtain

$$
\boxed{\operatorname{Var}(W_k)=kv_S.}
$$

Thus a $\pm1$ walk has $\operatorname{Var}(W_k)=4kp(1-p)$. The fair case has $\operatorname{Var}(W_k)=k$.

The standard deviation is the square root of the variance. A fair walk therefore has standard deviation $\sqrt{k}$. After four times as many steps, its standard deviation is twice as large, not four times as large. A zero expected position does not mean that a typical path stays at zero: positive and negative positions balance in the average while their spread increases.

### 2.5 Derive dependence quantitatively through covariance

Suppose $j\le k$. Separate the first $j$ steps from the later ones:

$$
W_k=W_j+R,
\qquad
R=S_{j+1}+\cdots+S_k.
$$

The variable $W_j$ depends only on the first block of steps, while $R$ depends only on the second block. Functions of disjoint sets of independent random variables are independent. Thus $\operatorname{Cov}(W_j,R)=0$.

Covariance is linear in each argument when the other is held fixed. This follows by distributing products in its centered definition. Consequently,

$$
\begin{aligned}
\operatorname{Cov}(W_j,W_k)
&=\operatorname{Cov}(W_j,W_j+R)\\
&=\operatorname{Cov}(W_j,W_j)+\operatorname{Cov}(W_j,R)\\
&=\operatorname{Var}(W_j)+0\\
&=jv_S.
\end{aligned}
$$

Interchanging $j$ and $k$ if necessary gives the symmetric formula

$$
\boxed{\operatorname{Cov}(W_j,W_k)=v_S\min(j,k).}
$$

The number of shared steps is $\min(j,k)$. That is precisely how much covariance the two accumulated positions share.

When $v_S>0$ and $1\le j\le k$, the correlation divides covariance by the product of standard deviations:

$$
\operatorname{Corr}(W_j,W_k)
=\frac{jv_S}{\sqrt{jv_S}\sqrt{kv_S}}
=\sqrt{\frac jk}.
$$

For example, the correlation between positions after nine and ten steps is $\sqrt{9/10}$, which is close to one. Most of their randomness is shared. The correlation between positions after one and one hundred steps is only $1/10$, because the hundredth position contains much more fresh randomness beyond the first step.

### 2.6 An increment is a change over an interval

An increment of the walk from step $a$ to step $b$, with $a<b$, is

$$
W_b-W_a=S_{a+1}+\cdots+S_b.
$$

The cancellation of the first $a$ steps explains why increments are more convenient than positions for describing fresh randomness. Consider

$$
U=W_5-W_2=S_3+S_4+S_5,
$$

and

$$
V=W_8-W_5=S_6+S_7+S_8.
$$

These increments use disjoint blocks of independent steps, so they are independent. The common time label $5$ is only a shared endpoint; no step is counted twice.

By contrast,

$$
V'=W_6-W_4=S_5+S_6
$$

shares $S_5$ with $U$. Expanding covariance term by term leaves only $\operatorname{Cov}(S_5,S_5)=v_S$, so $\operatorname{Cov}(U,V')=v_S$. Simply calling two quantities “increments” does not make them independent. Their intervals must use disjoint steps.

**Independent increments** means that increments over disjoint time intervals are mutually independent. **Stationary increments** means that the distribution of an increment depends only on the interval length. For i.i.d. steps,

$$
W_{k+h}-W_k=S_{k+1}+\cdots+S_{k+h}
\overset{d}=S_1+\cdots+S_h=W_h,
$$

where $\overset d=$ means equality in distribution, not equality of realized values. Every block of $h$ steps has the same joint law. That establishes stationary increments.

Do not confuse stationary increments with a stationary process. A stationary process has time-shift-invariant joint distributions, and in particular time-invariant one-time distributions. This walk's variance is $kv_S$, which changes with $k$ whenever $v_S>0$. Its increments can be stationary while its position process is not.

### 2.7 Derive the exact endpoint distribution

Let $K_k$ be the number of upward steps among the first $k$ steps. Then there are $k-K_k$ downward steps. The net position is

$$
W_k=K_k-(k-K_k)=2K_k-k.
$$

Therefore ending at $m$ requires

$$
K_k=\frac{k+m}{2}.
$$

This immediately gives two restrictions: $|m|\le k$, and $k+m$ must be even. A walk of four $\pm1$ steps cannot end at position $1$, because its possible endpoints are $-4,-2,0,2,4$.

If $r=(k+m)/2$ is an integer between zero and $k$, there are $\binom{k}{r}$ choices of which steps go up. Each corresponding sequence has probability $p^r(1-p)^{k-r}$. Thus

$$
\boxed{P(W_k=m)=\binom{k}{(k+m)/2}
 p^{(k+m)/2}(1-p)^{(k-m)/2}.}
$$

For all other $m$, the probability is zero. For a fair four-step walk, the probabilities at $-4,-2,0,2,4$ are $1,4,6,4,1$, respectively, divided by $16$. Central positions can be reached in more different orders than extreme positions. That counting explains the shape in this particular model.

### 2.8 Where the Gaussian approximation enters—and what it does not prove

The central limit theorem says that for i.i.d. steps with finite mean $m_S$ and finite positive variance $v_S$,

$$
\frac{W_k-km_S}{\sqrt{kv_S}}
\ \xrightarrow{d}\ N(0,1)
\qquad\text{as }k\to\infty.
$$

The centering subtracts the expected accumulated position. The denominator divides by its standard deviation, making the normalized variable have mean zero and variance one. Convergence in distribution here means that probabilities below fixed thresholds converge to those of the standard normal distribution at its continuity points. It does not mean that a particular random path becomes a bell curve.

The central limit theorem itself is an imported theorem, not proved by counting the outcomes of one finite walk. We can nevertheless see exactly how to use it. Undoing the centering and scaling suggests

$$
W_k\ \text{is approximately }N(km_S,kv_S)
$$

for large $k$. A finite $\pm1$ walk is still discrete, so the continuous normal distribution is an approximation for probabilities over ranges, not an exact identity of point masses. For a lattice calculation, a continuity correction may improve an approximation; the exact binomial formula remains available in this simple example.

The average step is $\overline S_k=W_k/k$. Dividing a variable by $k$ divides its variance by $k^2$, so

$$
E[\overline S_k]=m_S,
\qquad
\operatorname{Var}(\overline S_k)=v_S/k.
$$

There is no conflict between an accumulated sum spreading out and an average concentrating. They are different variables. The sum has variance proportional to $k$; the average has variance proportional to $1/k$.

### 2.9 The next-position distribution uses only the current position

The recursion $W_{k+1}=W_k+S_{k+1}$ gives a useful prediction. If $W_k=w$, then the next position is $w+1$ with probability $p$ and $w-1$ with probability $1-p$. Knowing earlier positions adds nothing once the current position is supplied, because the next step is independent of all earlier steps.

This is the Markov property in its simplest form. It does not say that $W_{k+1}$ is independent of $W_k$. Its distribution is centered around the supplied $w$, so it clearly depends on the current state. It says that a whole history can be replaced by its last position for the purpose of predicting the next position.

::: {.worked}
**Worked exercise: steps of unequal magnitude.** Suppose $P(S=2)=1/4$ and $P(S=-1)=3/4$. The mean is $2/4-3/4=-1/4$. The second moment is $4/4+3/4=7/4$. Therefore the variance is $7/4-(-1/4)^2=27/16$. For one hundred independent such steps, the mean position is $-25$ and its variance is $100(27/16)=168.75$. A normal approximation is $N(-25,168.75)$. For positions at steps 40 and 100, the covariance is $40(27/16)=67.5$. All of these calculations concern different features of the same accumulated process.
:::

We have built a dependent process without making its new disturbances dependent. We have also learned to separate endpoint values, interval changes, and entire paths. The next topic uses the path viewpoint to ask which boundary the process reaches first.

<p class="source-note">Source connection: random walks, their increment properties, and their relation to Markov chains appear in <a href="#ref-l5">MIT Lecture 5, sections 1–3</a>. The moment, covariance, and counting calculations are derived here from the step model.</p>
