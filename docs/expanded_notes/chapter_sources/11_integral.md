## 11. Weight the noise: the deterministic Wiener integral {#integral}

### 11.1 What object are we trying to define?

A Brownian motion $B_t$ starts at zero and has independent Gaussian increments: over a time interval of length $h$, the increment is distributed as $N(0,h)$. Here and throughout these notes, the second argument of $N(m,v)$ is the **variance**. An increment is the change $B_{t+h}-B_t$, not the value $B_{t+h}$ by itself.

Suppose disturbances arriving at different times do not have equal influence on a final outcome. Perhaps an early disturbance has mostly worn off, whereas a recent disturbance is still strongly felt. We need an accumulated-noise model that can assign different weights to different time intervals. The notation for that object is

$$
I(f)=\int_a^b f(s)\,dB_s.
$$

The function $f$ specifies the weights. In this chapter it is **deterministic**: once we choose the function and the interval, its values are fixed, not selected by looking at the realized Brownian path. The resulting integral is nevertheless random because the Brownian increments are random.

For comparison, $\int_a^b f(s)\,ds$ weights small amounts of **time**. The integral $\int_a^b f(s)\,dB_s$ weights small **changes in Brownian motion**. Replacing $dB_s$ by $ds$ changes the object and generally produces the wrong answer.

We will not define this integral as an antiderivative of a nonexistent ordinary Brownian derivative. We will first define it for weights that are constant on a finite collection of intervals. Then we will explain precisely how to take limits to obtain more general deterministic weights.

### 11.2 Start with weights that are constant on each interval

Choose a partition

$$
a=u_0<u_1<\cdots<u_m=b.
$$

A partition is simply a list of increasing times dividing the interval into smaller intervals. Let $f$ equal the constant $c_j$ on $(u_j,u_{j+1}]$. Such an $f$ is called a **step function**. It can be written compactly as

$$
f(s)=\sum_{j=0}^{m-1}c_j\mathbf 1_{(u_j,u_{j+1}]}(s),
$$

where the indicator $\mathbf 1_A(s)$ is 1 when $s$ belongs to the set $A$, and 0 otherwise. The indicator notation is only a way to express “use weight $c_j$ on the $j$th interval.”

For this function, define the integral by the finite sum

$$
\boxed{I(f)=\sum_{j=0}^{m-1}c_j(B_{u_{j+1}}-B_{u_j}).}
$$

There is nothing limiting or informal in this definition. Every Brownian value exists, so every difference and the finite weighted sum exist.

For example, take $a=0$, $b=3$, and use weight 2 for the first unit of time and weight $-1$ for the next two units. Then

$$
f(s)=\begin{cases}2,&0<s\leq1,\\-1,&1<s\leq3,\end{cases}
\qquad
I(f)=2(B_1-B_0)-(B_3-B_1).
$$

Because $B_0=0$, this equals $2B_1-(B_3-B_1)$. A negative weight reverses the sign of a disturbance; it does not make its variance negative.

### 11.3 Derive the mean and variance of the finite sum

Write $\Delta B_j=B_{u_{j+1}}-B_{u_j}$ and $h_j=u_{j+1}-u_j$. Brownian motion gives $E[\Delta B_j]=0$, $E[(\Delta B_j)^2]=h_j$, and independence between increments from distinct intervals.

Linearity of expectation therefore gives

$$
E[I(f)]=\sum_j c_j E[\Delta B_j]=0.
$$

Since the mean is zero, the variance equals the second moment. Expand the square before taking its expectation:

$$
I(f)^2
=\sum_j c_j^2(\Delta B_j)^2
+2\sum_{j<k}c_jc_k\Delta B_j\Delta B_k.
$$

For $j\ne k$, independence permits multiplication of the expectations:

$$
E[\Delta B_j\Delta B_k]
=E[\Delta B_j]E[\Delta B_k]=0.
$$

Consequently every cross term disappears, while each diagonal term contributes its interval length:

$$
\begin{aligned}
\operatorname{Var}(I(f))
&=\sum_j c_j^2 E[(\Delta B_j)^2]\\
&=\sum_j c_j^2h_j\\
&=\int_a^b f(s)^2\,ds.
\end{aligned}
$$

The last equality is an ordinary integral of a step function: its value is the sum of rectangle heights $c_j^2$ multiplied by their widths $h_j$.

For our two-weight example, the variance is $2^2(1)+(-1)^2(2)=6$. The two coefficients must be **squared individually**. Neither $\int f(s)\,ds$ nor $(\int f(s)\,ds)^2$ is the variance. In fact, this example has $\int_0^3 f(s)\,ds=2-2=0$ even though the random integral has variance 6. Opposite weights cancel in the ordinary integral; independent random disturbances do not cancel deterministically.

Each summand is Gaussian, and a sum of independent Gaussian variables is Gaussian. Therefore the finite-sum definition gives the exact law

$$
I(f)\sim N\!\left(0,\int_a^b f(s)^2\,ds\right).
$$

This is not an appeal to a central limit theorem. Even two weighted independent Brownian increments already give an exactly Gaussian sum.

### 11.4 Comparing two integrals explains how limits will work

Suppose $f$ and $g$ are step functions. Their breakpoints may differ, but we can combine both sets of breakpoints into a common partition. On each resulting interval, write their weights as $c_j$ and $d_j$.

The finite sums show directly that

$$
I(\alpha f+\beta g)=\alpha I(f)+\beta I(g)
$$

for deterministic numbers $\alpha,\beta$. This property is called **linearity**. It follows because each Brownian increment is multiplied by $\alpha c_j+\beta d_j$, and distributing that multiplication separates the sums.

They also give a covariance identity. Both integrals have mean zero, so

$$
\begin{aligned}
\operatorname{Cov}(I(f),I(g))
&=E\!\left[\left(\sum_j c_j\Delta B_j\right)
             \left(\sum_k d_k\Delta B_k\right)\right]\\
&=\sum_j\sum_k c_jd_kE[\Delta B_j\Delta B_k]\\
&=\sum_j c_jd_jh_j\\
&=\int_a^b f(s)g(s)\,ds.
\end{aligned}
$$

As before, different increments contribute zero cross expectation. Only the same-interval pairs survive.

Apply the variance formula to $f-g$. By linearity, its integral is $I(f)-I(g)$, and hence

$$
\boxed{E[(I(f)-I(g))^2]=\int_a^b(f(s)-g(s))^2\,ds.}
$$

This is the crucial bridge. It says that **a small integrated squared error in the weights produces an equally small expected squared error in the random integrals**. That is the control needed to define integrals by approximation.

### 11.5 What “square-integrable” and “mean-square convergence” mean

A deterministic function $f$ is square-integrable on $[a,b]$ when

$$
\int_a^b|f(s)|^2\,ds<\infty.
$$

The collection of such functions is conventionally called $L^2([a,b])$. The notation $L^2$ is not another stochastic model; it describes a class of functions with a finite squared-integral size. Functions that differ only on a set of time length zero represent the same element for this purpose.

For the smooth functions used later, such as exponentials on finite intervals, square-integrability is immediate: a bounded function on a finite interval has a finite integral of its square. More singular functions can also qualify. For example, $f(s)=s^{-1/4}$ on $(0,1]$ has $\int_0^1f(s)^2ds=\int_0^1s^{-1/2}ds=2$. In contrast, $s^{-1/2}$ on $(0,1]$ does not qualify because its square has divergent integral $\int_0^1s^{-1}ds$.

A sequence of random variables $U_m$ converges to $U$ **in mean square** when

$$
E[(U_m-U)^2]\longrightarrow0.
$$

This is a statement about the average squared discrepancy across random outcomes. It does not mean every approximation equals the limit, nor does it automatically assert convergence along every individual path.

Mean-square convergence implies convergence in probability. Indeed, for every $\varepsilon>0$,

$$
P(|U_m-U|>\varepsilon)
\leq\frac{E[(U_m-U)^2]}{\varepsilon^2}\longrightarrow0.
$$

The inequality follows because on the event $|U_m-U|>\varepsilon$, the squared error is at least $\varepsilon^2$. This elementary argument also explains why expected squared error is a useful way to control approximation.

### 11.6 Construct the integral for a general deterministic weight

Take a deterministic square-integrable $f$. A standard approximation theorem says that step functions $f_m$ can be chosen with

$$
\int_a^b(f_m(s)-f(s))^2\,ds\longrightarrow0.
$$

For continuous $f$, taking its values on increasingly fine partitions gives an intuitive construction. The approximation theorem extends the statement to all square-integrable $f$; we use that functional-analysis result rather than prove it here.

We already know how to define each $I(f_m)$. The difference identity gives

$$
E[(I(f_m)-I(f_k))^2]=\int_a^b(f_m-f_k)^2\,ds\longrightarrow0
\quad\text{as }m,k\to\infty.
$$

Thus the random integrals form a **Cauchy sequence in mean square**: late members become arbitrarily close to one another in expected squared error. Completeness of the space of square-integrable random variables guarantees a mean-square limit. Completeness is the general mathematical theorem that “Cauchy sequences in this space have a limit in the same space.”

We **define** $I(f)$ to be this limit. There is still a question: would a different sequence of step approximations produce another integral? Let $g_m$ also approximate $f$ in squared-integral norm. Then

$$
E[(I(f_m)-I(g_m))^2]=\int_a^b(f_m-g_m)^2\,ds\longrightarrow0.
$$

Therefore their limits agree in mean square and hence almost surely. The integral does not depend on which valid approximation sequence we chose. This uniqueness is essential: otherwise the symbol $\int f\,dB$ would not identify a well-defined random variable.

Linearity extends to the limits. The mean remains zero because

$$
|E[I(f_m)]-E[I(f)]|
\leq E[|I(f_m)-I(f)|]
\leq\sqrt{E[(I(f_m)-I(f))^2]}\longrightarrow0.
$$

The second inequality is the Cauchy–Schwarz inequality applied to the random error and the constant 1. The finite-sum second-moment identity extends as well:

$$
\boxed{E\!\left[\left(\int_a^bf(s)\,dB_s\right)^2\right]
=\int_a^bf(s)^2\,ds.}
$$

This is the **Itô isometry**, here in its deterministic-integrand form. An isometry preserves a notion of size or distance. The preceding difference identity explains the name more fully: squared-integral distance between deterministic functions becomes mean-square distance between their stochastic integrals.

### 11.7 Why the limiting integral is still Gaussian

Knowing the mean and variance of a random variable does not generally determine its distribution. We therefore need to justify the additional Gaussian conclusion, rather than infer it just from the isometry.

Each approximating integral has distribution $N(0,v_m)$, where $v_m=\int f_m^2$. Squared-integral convergence of $f_m$ to $f$ implies $v_m\to v=\int f^2$. For example, using Cauchy–Schwarz,

$$
\left|\int(f_m^2-f^2)\right|
=\left|\int(f_m-f)(f_m+f)\right|
\leq\left(\int(f_m-f)^2\right)^{1/2}
     \left(\int(f_m+f)^2\right)^{1/2}\longrightarrow0.
$$

The first factor tends to zero and the second stays bounded. Thus the Gaussian distributions of the approximations tend to $N(0,v)$. At the same time, mean-square convergence makes $I(f_m)$ converge in probability, and therefore in distribution, to $I(f)$. Uniqueness of the limiting distribution gives

$$
\boxed{\int_a^bf(s)\,dB_s\sim N\!\left(0,\int_a^bf(s)^2\,ds\right).}
$$

Here we use the standard facts that probability convergence implies distribution convergence and that a distributional limit is unique. The argument is specific: **Gaussian approximations whose variances converge**, not an assertion that every limit of arbitrary weighted noise must be Gaussian. When $v=0$, the integral is zero almost surely, conventionally a degenerate Gaussian.

The covariance identity also extends:

$$
\operatorname{Cov}\!\left(\int f\,dB,\int g\,dB\right)=\int fg.
$$

Moreover, any finite collection of these deterministic-weight integrals is jointly Gaussian. To see why, any linear combination is the integral of the corresponding linear combination of deterministic weights, so it is Gaussian. This is one characterization of a jointly Gaussian vector.

For jointly Gaussian variables, zero covariance implies independence. Consequently, deterministic-weight integrals over disjoint intervals are independent: their weights have disjoint supports, so $fg=0$ except possibly at endpoints, and the covariance integral is zero. An integral over $(s,t]$ is also independent of the Brownian information up to time $s$. Its approximating sums use only future independent increments. If an initial value $X_0$ is independent of future Brownian noise, adding that initial value to the past information preserves this conclusion.

### 11.8 A special integration-by-parts formula requires no Brownian derivative

For a deterministic continuously differentiable $f$, the Wiener integral can also be written using an ordinary time integral of the Brownian path. This is especially useful for solving OU without invoking the full stochastic product rule.

Begin with the exact finite-partition identity

$$
\begin{aligned}
f(u_{j+1})B_{u_{j+1}}-f(u_j)B_{u_j}
={}&f(u_j)(B_{u_{j+1}}-B_{u_j})\\
&+B_{u_{j+1}}(f(u_{j+1})-f(u_j)).
\end{aligned}
$$

Check it by expanding the right side: the two terms involving $f(u_j)B_{u_{j+1}}$ cancel. Sum over all intervals. The left side telescopes, leaving only endpoint products. Rearranging gives

$$
\begin{aligned}
\sum_j f(u_j)(B_{u_{j+1}}-B_{u_j})
={}&f(b)B_b-f(a)B_a\\
&-\sum_j B_{u_{j+1}}(f(u_{j+1})-f(u_j)).
\end{aligned}
$$

As the partition becomes fine, the first sum on the left converges in mean square to the Wiener integral, because the left-step approximations converge to the continuous $f$. On the right, use

$$
f(u_{j+1})-f(u_j)=\int_{u_j}^{u_{j+1}}f'(s)\,ds.
$$

Continuity of the Brownian path ensures that replacing $B_{u_{j+1}}$ by $B_s$ within these small intervals changes the sum by a quantity tending to zero. More explicitly, its absolute error is at most the maximum within-interval oscillation of $B$ multiplied by $\int_a^b|f'(s)|ds$. Uniform continuity on the finite interval makes that maximum oscillation tend to zero almost surely.

Thus the right side has a pathwise limit. Since the left side also has a mean-square limit, both limits agree almost surely. We obtain

$$
\boxed{\int_a^bf(s)\,dB_s
=f(b)B_b-f(a)B_a-\int_a^bf'(s)B_s\,ds.}
$$

The derivative in this formula is $f'$, not $B'$. Brownian motion appears inside an ordinary integral as a continuous function. This is why the calculation is legal despite Brownian paths having no ordinary derivative.

### 11.9 Two smooth examples, worked all the way through

First use $f(s)=s$ on $[0,t]$. The Gaussian law and the isometry give

$$
E\!\left[\int_0^t s\,dB_s\right]=0,
\qquad
\operatorname{Var}\!\left(\int_0^t s\,dB_s\right)
=\int_0^t s^2ds
=\left[\frac{s^3}{3}\right]_0^t=\frac{t^3}{3}.
$$

Therefore $\int_0^t s\,dB_s\sim N(0,t^3/3)$. Independently, the integration-by-parts formula, with $f'(s)=1$, gives

$$
\int_0^t s\,dB_s=tB_t-\int_0^t B_s\,ds.
$$

These are two descriptions of the same random variable: one gives its distribution, and one expresses it in terms of a realized path. The ordinary integral on the right remains random because its integrand is the random path.

Now fix a terminal time $t$ and a positive number $\theta$, and use $f(s)=e^{-\theta(t-s)}$. Even though the symbol $t$ appears in the weight, the weight is deterministic once the terminal time is fixed. Its derivative with respect to the integration variable $s$ is

$$
f'(s)=\theta e^{-\theta(t-s)}.
$$

The sign is positive: $-\theta(t-s)=-\theta t+\theta s$. Integration by parts gives

$$
\boxed{\int_0^t e^{-\theta(t-s)}\,dB_s
=B_t-\theta\int_0^t e^{-\theta(t-s)}B_s\,ds.}
$$

The variance follows by squaring the weight, substituting $r=t-s$, and integrating:

$$
\begin{aligned}
\int_0^t e^{-2\theta(t-s)}ds
&=\int_0^t e^{-2\theta r}dr\\
&=\left[-\frac{e^{-2\theta r}}{2\theta}\right]_0^t\\
&=\frac{1-e^{-2\theta t}}{2\theta}.
\end{aligned}
$$

Both the pathwise identity and this variance calculation will be used in the OU solution.

### 11.10 What fails when the weight itself is random?

The deterministic restriction is substantive. A random weight may depend on the same noise that appears in the increment sum, so the finite sum is no longer a fixed linear combination of independent Gaussian variables. A general adapted Itô integral can still exist and satisfy a suitable isometry, but its distribution need not be Gaussian.

A concrete illustration is the random integrand $B_s$ itself. On a partition, the elementary identity

$$
B_{u_{j+1}}^2-B_{u_j}^2
=2B_{u_j}(B_{u_{j+1}}-B_{u_j})
 +(B_{u_{j+1}}-B_{u_j})^2
$$

gives, after summing and using $B_0=0$,

$$
\sum_jB_{u_j}\Delta B_j
=\frac12\left(B_t^2-\sum_j(\Delta B_j)^2\right).
$$

For Brownian motion the sum of squared increments tends to $t$ in mean square as a deterministic partition's mesh tends to zero. The reason is that its mean is $\sum_jh_j=t$ and its variance is $2\sum_jh_j^2\leq2t\max_jh_j\to0$. Thus the left-endpoint sums here converge to

$$
\int_0^tB_s\,dB_s=\frac12(B_t^2-t).
$$

For $t>0$, this is not a Gaussian random variable: it is nonconstant but can never be less than $-t/2$, whereas any nondegenerate Gaussian has possible values arbitrarily far in both directions. The example is not a construction of every random-integrand integral. It is a direct counterexample to the idea that “integrating against Brownian motion always produces a Gaussian.”

### 11.11 A complete covariance example

::: {.worked}
Let $f=2$ on $(0,1]$ and $f=-1$ on $(1,3]$. Let $g=1$ on $(0,2]$ and $g=0$ on $(2,3]$. Define $U=\int_0^3f\,dB$ and $V=\int_0^3g\,dB$.

Use the common partition $0,1,2,3$. The weights are $(2,-1,-1)$ for $U$ and $(1,1,0)$ for $V$, on intervals of unit length. Therefore

$$
\operatorname{Var}(U)=4+1+1=6,
\qquad
\operatorname{Var}(V)=1+1+0=2,
$$

and

$$
\operatorname{Cov}(U,V)=2(1)+(-1)(1)+(-1)(0)=1.
$$

Thus $\operatorname{Corr}(U,V)=1/\sqrt{12}$. The pair is jointly Gaussian, but not independent because the covariance is nonzero. Also $U+V=\int(f+g)dB$ has weights $(3,0,-1)$, hence variance $9+0+1=10$. This agrees with the general variance formula $6+2+2(1)=10$.
:::

The central lesson is now precise. Deterministic weights turn Brownian increments into a Gaussian accumulated noise whose variance is the ordinary integral of the **squared** weights. We know how the object is constructed, why the construction is independent of the approximations, and why its law has the stated form.

<p class="source-note">Source connection: deterministic Wiener integration and the Itô isometry are treated in <a href="#ref-l18">MIT Lecture 18</a> and <a href="#ref-ss">Särkkä and Solin, §4.1</a>. The finite-sum proofs, approximation explanation and deterministic integration-by-parts derivation are given explicitly here. The general density of step functions and completeness of the square-integrable spaces are identified as standard theorems, not presented as proofs completed in these notes.</p>
