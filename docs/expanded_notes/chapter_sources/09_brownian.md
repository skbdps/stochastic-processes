## 9. Read Brownian motion correctly {#brownian}

### 9.1 Define the process before using its notation

Brownian motion is a mathematical model for continuously accumulated random fluctuations. The word “continuous” refers to its sample paths, not to smoothness or differentiability. We define the process through its probability properties rather than by pretending there is a smallest next instant of time.

A **standard Brownian motion** $(B_t)_{t\ge0}$ has the following properties. It starts at $B_0=0$. Its paths are continuous with probability one. Its increments over disjoint time intervals are mutually independent. Finally, for every $0\le s<t$,

$$
\boxed{B_t-B_s\sim N(0,t-s).}
$$

The notation $N(m,v)$ means Gaussian with mean $m$ and variance $v$. Thus the increment's standard deviation is $\sqrt{t-s}$. The increment distribution depends only on elapsed time, not on when the interval starts. This is stationary increments.

“Standard” fixes the variance rate to one. Multiplying by a constant $\sigma$ gives increments of variance $\sigma^2(t-s)$. Adding a deterministic drift $at$ changes the increment mean to $a(t-s)$.

Existence of a process satisfying the definition is a theorem. The normalized random-walk construction motivates it and, with a functional limit theorem, provides one route to it. The defining properties will be the starting assumptions for the derivations in this chapter.

### 9.2 The distribution at one time

Because $B_0=0$, the state at time $t$ is itself the increment over $[0,t]$:

$$
B_t=B_t-B_0\sim N(0,t).
$$

Therefore

$$
E[B_t]=0,\qquad \operatorname{Var}(B_t)=t,
\qquad E[B_t^2]=t.
$$

The last equality follows from $\operatorname{Var}(U)=E[U^2]-(E[U])^2$ and the zero mean. It is not a new independent assumption.

At time four, the distribution is $N(0,4)$ and its standard deviation is two. To calculate a probability, standardize:

$$
P(B_4>2)=P\left(\frac{B_4}{2}>1\right)=P(Z>1),
\qquad Z\sim N(0,1).
$$

If $\Phi(z)=P(Z\le z)$ denotes the standard normal distribution function, this is $1-\Phi(1)$, approximately $0.159$. The process mean staying zero does not imply that a realized path stays close to zero; its one-time standard deviation grows as $\sqrt t$.

An interval such as $[-1.96\sqrt t,1.96\sqrt t]$ contains approximately $95\%$ of the distribution at a **specified single time** $t$. It is not automatically a band with $95\%$ probability of containing the entire path over many times. A pathwise event places multiple simultaneous restrictions and is a different probability question.

### 9.3 Independent increments do not mean independent positions

Take $0<s<t$ and split the later position into two pieces:

$$
B_t=B_s+(B_t-B_s).
$$

The first piece is the increment over $[0,s]$; the second is the increment over $[s,t]$. These pieces are independent. Nevertheless, $B_s$ and $B_t$ are not independent because the second position includes the first piece.

Covariance quantifies the shared contribution. By its centered definition, covariance is linear in either argument and $\operatorname{Cov}(U,U)=\operatorname{Var}(U)$. Therefore

$$
\begin{aligned}
\operatorname{Cov}(B_s,B_t)
&=\operatorname{Cov}(B_s,B_s+(B_t-B_s))\\
&=\operatorname{Var}(B_s)+\operatorname{Cov}(B_s,B_t-B_s)\\
&=s+0\\
&=s.
\end{aligned}
$$

The zero cross term uses independence of the increments. Reversing the roles of $s$ and $t$ when necessary gives

$$
\boxed{\operatorname{Cov}(B_s,B_t)=\min(s,t).}
$$

For positive $s\le t$, the correlation is

$$
\operatorname{Corr}(B_s,B_t)
=\frac{s}{\sqrt s\sqrt t}=\sqrt{\frac st}.
$$

At $s=0$, the initial value has variance zero, so a correlation coefficient involving $B_0$ is undefined; the covariance formula still gives zero. The distinction between covariance and correlation includes their denominators and domains, not just their names.

### 9.4 Joint Gaussianity: why means and covariances become enough

Choose times $0=t_0<t_1<\cdots<t_m$. Define increments $U_j=B_{t_j}-B_{t_{j-1}}$. By the Brownian definition, these are independent centered Gaussian variables with variances $t_j-t_{j-1}$. The positions are their partial sums:

$$
B_{t_1}=U_1,\quad B_{t_2}=U_1+U_2,\quad\ldots,\quad B_{t_m}=U_1+\cdots+U_m.
$$

Any linear combination of the positions is therefore a linear combination of independent Gaussian increments, and is Gaussian. This is the defining property of a jointly Gaussian vector. For such vectors, the mean vector and covariance matrix determine the distribution.

For example,

$$
\begin{pmatrix}B_1\\B_3\\B_4\end{pmatrix}
\quad\text{has mean }\begin{pmatrix}0\\0\\0\end{pmatrix}
\quad\text{and covariance }
\begin{pmatrix}1&1&1\\1&3&3\\1&3&4\end{pmatrix}.
$$

The off-diagonal entries are the lengths of shared history. The diagonal entries are the variances at their respective times.

The qualification “jointly Gaussian” cannot be dropped. Knowing that each separate state is Gaussian does not by itself imply that every linear combination is Gaussian or that covariances determine all dependence. Here joint Gaussianity was derived from the independent Gaussian increments.

### 9.5 Derive the conditional transition law and the Markov property

Suppose the value at time $s$ is known to be $b$. The decomposition $B_t=B_s+(B_t-B_s)$ becomes a fixed starting value plus a fresh Gaussian increment:

$$
\boxed{B_t\mid B_s=b\sim N(b,t-s).}
$$

For example, before any observation, $B_5\sim N(0,5)$. After observing $B_2=3$, the remaining uncertainty is the increment over three time units, so $B_5\mid B_2=3\sim N(3,3)$. Conditioning changes both the mean and the variance.

More generally, the future increments are independent of the information generated by all Brownian values up to time $s$. That accumulated information is often written $\mathcal F_s$ and called the Brownian filtration at $s$. A filtration is an increasing family of information sets: more time reveals more observations.

Once $B_s$ is supplied, the older Brownian history cannot improve the distributional prediction of $B_t$. This establishes the Markov property. It is a consequence of the increment construction, not a separate appeal to independence of the positions.

### 9.6 Covariance of arbitrary increments is overlap length

Let $U=B_b-B_a$ and $V=B_d-B_c$, with $a<b$ and $c<d$. Expand covariance using the position covariance already derived:

$$
\begin{aligned}
\operatorname{Cov}(U,V)
&=\operatorname{Cov}(B_b,B_d)-\operatorname{Cov}(B_b,B_c)\\
&\quad-\operatorname{Cov}(B_a,B_d)+\operatorname{Cov}(B_a,B_c)\\
&=\min(b,d)-\min(b,c)-\min(a,d)+\min(a,c).
\end{aligned}
$$

This expression equals the length of the overlap of the intervals $[a,b]$ and $[c,d]$. If the intervals are disjoint apart from an endpoint, the covariance is zero, and Brownian independent increments give independence. For overlapping intervals, the common Brownian fluctuation contributes to both changes.

For $U=B_3-B_1$ and $V=B_4-B_2$, the overlap is $[2,3]$, of length one. Substitution confirms it:

$$
\operatorname{Cov}(U,V)=3-2-1+1=1.
$$

Each increment has variance two, so their correlation is $1/2$. This example distinguishes “an increment” from “an independent increment relative to another specified interval.”

### 9.7 Why Brownian motion is not stationary

A stationary process has joint distributions unchanged when all time arguments are shifted together. In particular, its one-time distribution cannot depend on absolute time. Brownian motion fails that condition because $\operatorname{Var}(B_t)=t$.

Its increments are stationary because $B_{s+h}-B_s\sim N(0,h)$ for every starting time $s$. The quantity being shifted in this statement is an **interval change**, not an absolute position. These two notions of stationarity concern different objects.

A related but different symmetry is Brownian scaling. For a positive constant $c$, define $\widetilde B_t=B_{ct}/\sqrt c$. Its increments satisfy

$$
\widetilde B_t-\widetilde B_s
=\frac{B_{ct}-B_{cs}}{\sqrt c}
\sim N\left(0,\frac{c(t-s)}c\right)=N(0,t-s).
$$

They remain independent over disjoint intervals, the paths remain continuous, and the start remains zero. Hence $\widetilde B$ is also standard Brownian motion. This is equality of process laws, not an assertion that $B_{ct}/\sqrt c$ and $B_t$ are equal on each particular path.

### 9.8 Why there is no ordinary Brownian velocity

For an ordinary differentiable function $x(t)$, the difference quotient $[x(t+h)-x(t)]/h$ approaches a finite derivative as $h\to0$. Brownian increments behave differently.

For a fixed deterministic time $t$ and $h>0$, define

$$
Q_h=\frac{B_{t+h}-B_t}{h}.
$$

The numerator is $N(0,h)$. Dividing by $h$ multiplies variance by $1/h^2$, giving

$$
\boxed{Q_h\sim N(0,1/h).}
$$

Its typical magnitude is $1/\sqrt h$, rather than a finite slope. But variance divergence alone is not a general proof that random variables cannot converge. We can make a more precise argument using the actual Gaussian law.

For any fixed finite $M>0$,

$$
P(|Q_h|\le M)=P(|Z|\le M\sqrt h),\qquad Z\sim N(0,1).
$$

The standard normal density is at most $1/\sqrt{2\pi}$. The interval $[-M\sqrt h,M\sqrt h]$ has length $2M\sqrt h$, so

$$
P(|Q_h|\le M)\le\frac{2M\sqrt h}{\sqrt{2\pi}}\to0.
$$

If a finite derivative existed at this fixed $t$, then along the sequence $h=1/n$ the quotients would eventually be bounded by some integer $M$. For fixed integers $M,K$, the event that $|Q_{1/n}|\le M$ for every $n\ge K$ is contained in each individual event $|Q_{1/n}|\le M$. Its probability is therefore no larger than arbitrarily small such probabilities, and hence is zero. Taking the countable union over $M$ and $K$ shows that a finite derivative at this fixed time has probability zero.

The stronger theorem is that a Brownian path is almost surely nowhere differentiable: with probability one there is no time with a finite ordinary derivative. The fixed-time calculation alone does not prove that stronger statement, because an uncountable union of probability-zero events need not have probability zero. The stronger theorem is a genuine additional regularity result. We will use the fact, not pretend that the fixed-time argument proved more than it did.

### 9.9 Squared increments do not disappear

A second calculation explains why ordinary differential manipulations can fail. Partition $[0,T]$ into deterministic intervals with endpoints $0=t_0<t_1<\cdots<t_m=T$. Let $h_j=t_{j+1}-t_j$ and $\Delta B_j=B_{t_{j+1}}-B_{t_j}$. Consider

$$
Q=\sum_{j=0}^{m-1}(\Delta B_j)^2.
$$

Each increment is $N(0,h_j)$, so $E[(\Delta B_j)^2]=h_j$. Therefore

$$
E[Q]=\sum_j h_j=T.
$$

To see whether $Q$ concentrates near $T$, calculate its variance. If $Z\sim N(0,1)$, then $E[Z^4]=3$. This follows by integration by parts: writing the standard normal density as $\varphi$ with $\varphi'(z)=-z\varphi(z)$,

$$
\int z^4\varphi(z)\,dz
=-\int z^3\varphi'(z)\,dz
=3\int z^2\varphi(z)\,dz=3,
$$

where the integrals run over the real line and the Gaussian tails make the boundary term vanish. Thus $E[(\Delta B_j)^4]=3h_j^2$, and

$$
\operatorname{Var}((\Delta B_j)^2)=3h_j^2-h_j^2=2h_j^2.
$$

The squared increments are independent because the underlying increments are independent. Hence

$$
\operatorname{Var}(Q)=2\sum_jh_j^2
\le2\left(\max_jh_j\right)\sum_jh_j
=2T\max_jh_j.
$$

As the largest interval length tends to zero,

$$
\boxed{E[(Q-T)^2]\to0.}
$$

This is convergence in **mean square**: the expected squared error tends to zero. In particular it implies convergence in probability, because $P(|Q-T|>\epsilon)\le E[(Q-T)^2]/\epsilon^2$. The limiting value $T$ is Brownian motion's quadratic variation along these refining deterministic partitions.

For a smooth function with bounded derivative $|x'|\le C$ on $[0,T]$, each increment satisfies $|\Delta x_j|\le Ch_j$. Thus $\sum_j(\Delta x_j)^2\le C^2\sum_jh_j^2\to0$. Brownian motion instead leaves a nonzero accumulated squared-increment term. This difference is the source of the extra second-order terms in Itô calculus.

### 9.10 The meaning of the informal symbol $(dB)^2=dt$

The expression $(dB_t)^2=dt$ is often used as a mnemonic for quadratic variation. It is not a pathwise equality saying the square of each small Brownian increment equals its interval length. A single increment is random, and its square fluctuates.

The rigorous calculation above says that the **sum** of squared increments over a refining partition approaches total elapsed time in a specified mode of convergence. That is the fact a valid stochastic-calculus argument must ultimately use. We will not treat $dB$ as an ordinary infinitesimal whose algebra is self-explanatory.

### 9.11 Exact simulation at a finite set of times

For deterministic times $0=t_0<t_1<\cdots<t_n$, draw independent standard normal variables $Z_1,\ldots,Z_n$ and set

$$
B_{t_i}=B_{t_{i-1}}+\sqrt{t_i-t_{i-1}}\,Z_i.
$$

Each increment has the required Gaussian law and is independent of the others. Thus the resulting finite vector has the exact Brownian distribution at those times. A line drawn between the simulated points is only a visual interpolation; it does not recreate the unobserved Brownian path inside each interval.

::: {.worked}
**Worked check.** Compute $\operatorname{Cov}(B_2,B_5)$ and the law of $B_5$ after observing $B_2=-1$. The covariance is the shared time length two. The future increment is independent $N(0,3)$, so the conditional law is $N(-1,3)$. The unconditional law is $N(0,5)$. Now compare increments $B_2-B_0$ and $B_5-B_2$: they are independent, even though the positions $B_2$ and $B_5$ are not. Every answer follows from identifying what randomness is shared and what randomness is new.
:::

<p class="source-note">Source connection: Brownian motion and its regularity are discussed in <a href="#ref-l17">MIT Lecture 17</a>; quadratic variation motivates <a href="#ref-l18">Lecture 18</a>. Covariance, conditional prediction, the fixed-time derivative argument, and the mean-square quadratic-variation calculation are derived above.</p>
