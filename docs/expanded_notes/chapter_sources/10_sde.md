## 10. What an SDE actually says {#sde}

### 10.1 There are two different unknowns in this project

A stochastic differential equation, or SDE, specifies how a random process evolves. For example,

$$
dX_t=\theta(\mu-X_t)\,dt+\sigma\,dB_t
$$

is the Ornstein–Uhlenbeck model. The state $X_t$ is a random quantity at time $t$. The parameters $\theta>0$, $\mu\in\mathbb R$, and $\sigma>0$ control the model. The process $B_t$ is standard Brownian motion: it starts at zero, has continuous paths, and has independent Gaussian increments of variance equal to elapsed time.

**Solving the SDE** means expressing the process $X$ in terms of a given starting value, given parameters, and the driving Brownian noise. **Estimating the parameters** means using observations of $X$ to infer the unknown values of $\theta,\mu,\sigma$. These tasks are connected, but they are not the same operation. We usually solve or characterize the model first so that we can later calculate its likelihood.

Before doing either task, we need to give the displayed differential equation a precise meaning. In particular, $dB_t$ cannot be treated as an ordinary derivative times $dt$, because Brownian paths do not possess an ordinary finite derivative.

### 10.2 Start with the meaning of an ordinary differential equation

Consider the ordinary equation

$$
x'(t)=g(t),\qquad x(0)=x_0,
$$

with continuous $g$. Its integral form, by the fundamental theorem of calculus, is

$$
x(t)=x_0+\int_0^t g(s)\,ds.
$$

The integral accumulates changes over time. On a partition $0=t_0<\cdots<t_m=t$, the sum

$$
\sum_{j=0}^{m-1}g(t_j)(t_{j+1}-t_j)
$$

approximates the accumulated change. Each summand is a rate multiplied by a duration. As the largest interval length tends to zero, the sum tends to the ordinary integral.

The shorthand $dx=g(t)dt$ is compatible with this accumulation interpretation. It is not a license to regard $dx$ and $dt$ as arbitrary independent algebraic symbols. The underlying derivative or integral statement determines what operations are valid.

For example, $dx=3dt$ with $x(0)=1$ means $x(t)=1+3t$. The equation $dx=2t\,dt$ means $x(t)=x_0+t^2$. Both are solved by integrating known functions of time.

### 10.3 A random path can still have an ordinary time integral

Suppose $X_s$ is a continuous random process. For one realized outcome $\omega$, the function $s\mapsto X_s(\omega)$ is a continuous ordinary function. Therefore

$$
\int_0^t X_s(\omega)\,ds
$$

exists as an ordinary time integral on a bounded interval. The result is random across different outcomes, but the integration variable is still ordinary time.

Thus the OU drift integral

$$
\theta\int_0^t(\mu-X_s)\,ds
$$

does not require differentiating $X$. A continuous path can be integrated even when it cannot be differentiated. This distinction is important: roughness of the path does not make every integral involving it mysterious.

The new issue is an integral **against Brownian increments**, written $\int b_s\,dB_s$. Its increments come from the random integrator $B$, rather than from time intervals $ds$. We will define that object explicitly in Chapter 11, beginning with constant and step-function weights.

### 10.4 Constant Brownian weight: the sum telescopes exactly

For any partition of $[0,t]$,

$$
\begin{aligned}
\sum_{j=0}^{m-1}(B_{t_{j+1}}-B_{t_j})
&=(B_{t_1}-B_0)+(B_{t_2}-B_{t_1})+\cdots+(B_t-B_{t_{m-1}})\\
&=B_t-B_0\\
&=B_t.
\end{aligned}
$$

Every interior Brownian value appears once with a plus sign and once with a minus sign. The cancellation is exact at every partition; no derivative is involved.

Multiplying by a constant $\sigma$ gives

$$
\int_0^t\sigma\,dB_s=\sigma(B_t-B_0)=\sigma B_t.
$$

For this constant integrand, the stochastic integral is already completely understandable as an accumulated increment sum. We do not need a general Itô integral construction merely to read this case.

Contrast this with $\int_0^t\sigma B_s\,ds$. That is an ordinary time integral of the Brownian path and is not generally equal to $\sigma B_t$. The placement of $B_s$ and the differential matters: weighting elapsed time by a Brownian value is different from accumulating Brownian increments.

### 10.5 The precise meaning of additive-noise SDEs

For a constant noise scale $\sigma$, the equation

$$
dX_t=a(X_t,t)\,dt+\sigma\,dB_t,\qquad X_0=x_0,
$$

means

$$
\boxed{X_t=x_0+\int_0^t a(X_s,s)\,ds+\sigma B_t.}
$$

A solution is a process that satisfies this integral equation, with a continuous version and the specified initial value. Its value at time $t$ uses the initial condition and noise up to time $t$, not future Brownian increments.

A process that uses only information available up to its current time is called **adapted** to that information flow. Formally, the information available at time $t$ is represented by a filtration $\mathcal F_t$, an increasing family of collections of events. For our linear model, one can think of $\mathcal F_t$ as information from $X_0$ and all Brownian values observed by time $t$. No advanced filtration calculation is needed to understand the examples here.

The integral equation is the definition behind the differential shorthand. We should not divide the displayed SDE by $dt$ and interpret $dB_t/dt$ as a finite ordinary velocity. Over a short interval of length $h$, the Brownian increment has standard deviation $\sqrt h$. Dividing by $h$ produces standard deviation $1/\sqrt h$, not a finite derivative.

### 10.6 Solve constant drift plus constant noise

Consider

$$
dX_t=a\,dt+\sigma\,dB_t,\qquad X_0=x_0,
$$

where $a$ and $\sigma$ are fixed constants. Its integral form is

$$
X_t=x_0+\int_0^t a\,ds+\sigma B_t=x_0+at+\sigma B_t.
$$

The right side now contains only known constants and the given Brownian process. That is an explicit solution.

Because $B_t\sim N(0,t)$, scaling and shifting give

$$
\boxed{X_t\sim N(x_0+at,\sigma^2t).}
$$

For a later time $s+h$, subtract the solution at $s$:

$$
X_{s+h}-X_s=ah+\sigma(B_{s+h}-B_s).
$$

The future Brownian increment is independent of the past, so

$$
X_{s+h}\mid X_s=x\sim N(x+ah,\sigma^2h).
$$

Take $x_0=1$, $a=3$, and $\sigma=2$. Then $X_t=1+3t+2B_t$ and $X_t\sim N(1+3t,4t)$. At $t=2$, the mean is seven and the variance is eight. The noise **standard deviation** is $2\sqrt2$, not eight and not four.

If the drift is a deterministic time function instead, such as $a(s)=2s$, the same reading gives $X_t=x_0+t^2+\sigma B_t$. Its mean is $x_0+t^2$ and its variance remains $\sigma^2t$. Randomness does not interfere with the ordinary integration of a known drift function.

### 10.7 Read the OU equation term by term

For OU, the drift is $a(x,t)=\theta(\mu-x)$ and the noise scale is constant. Its integral equation is

$$
\boxed{X_t=X_0+\theta\int_0^t(\mu-X_s)\,ds+\sigma B_t.}
$$

The starting value contributes $X_0$. At each time $s$, the discrepancy from the target is $\mu-X_s$. The parameter $\theta$ converts that discrepancy into a rate of expected correction. Integrating the rate over time accumulates the drift contribution. The last term accumulates the independent Brownian disturbances.

If the state is above $\mu$, the drift coefficient is negative; if below $\mu$, it is positive. But the noise can move the state in either direction, so this does not require every observed increment to point toward the target.

The equation is not yet an explicit solution because $X_s$ appears inside the integral on the right. It is a self-consistent relationship that the whole process must satisfy. Merely giving it an integral meaning does not remove that unknown. Chapter 12 will solve this dependence using an ordinary integrating-factor argument after separating the additive Brownian term.

### 10.8 Exact increments and frozen-drift approximations are different

Subtract the OU integral equation at time $t$ from the equation at time $t+h$:

$$
\boxed{
X_{t+h}-X_t
=\theta\int_t^{t+h}(\mu-X_s)\,ds
+\sigma(B_{t+h}-B_t).
}
$$

This is an exact equality. The drift depends on the path throughout the interval, not only on its starting value.

A simple approximation replaces the evolving $X_s$ inside the drift integral by the starting state $X_t$:

$$
X_{t+h}-X_t
\approx\theta(\mu-X_t)h+\sigma(B_{t+h}-B_t).
$$

The approximation is called an Euler–Maruyama step. The exact difference between the true drift integral and its frozen version is

$$
\theta\int_t^{t+h}(X_t-X_s)\,ds.
$$

This is generally nonzero at a finite $h$. Under continuity it becomes small relative to $h$ along a fixed path as $h\to0$, because

$$
\left|\theta\int_t^{t+h}(X_t-X_s)ds\right|
\le\theta h\sup_{t\le s\le t+h}|X_t-X_s|.
$$

The supremum tends to zero by continuity. This local observation explains the approximation; establishing global error rates for a numerical scheme is a further question.

An Euler **approximation process** on grid times can be defined exactly by

$$
X^{E}_{i+1}=X^{E}_i+\theta(\mu-X^{E}_i)h_i+\sigma\sqrt{h_i}\,Z_i,
$$

where $Z_i$ are independent standard normals. The equality is exact as the definition of the numerical recursion, not as a claim that its finite-gap distribution equals the true OU transition.

### 10.9 Why a short-time drift does not mean the noise is negligible

Over a small duration $h$, a fixed drift rate contributes order $h$, while Brownian noise has standard deviation proportional to $\sqrt h$. Since $\sqrt h$ is larger than $h$ when $0<h<1$, random fluctuations can dominate individual very short increments even though the drift accumulates a systematic effect over time.

This explains why seeing one short movement away from the mean is not evidence against mean reversion. The drift describes a conditional expected tendency embedded in noisier individual changes. Estimation uses many transitions and their time gaps to separate these effects.

The dimensions reinforce the distinction. The return speed $\theta$ has units of inverse time, while $\sigma$ has units of state divided by square root of time. The quantity $\theta(\mu-X_t)h$ and the quantity $\sigma\sqrt h Z$ both have units of the state, although they scale differently with duration.

### 10.10 State-dependent diffusion needs an additional integral definition

A general Itô SDE is written

$$
dX_t=a(X_t,t)\,dt+b(X_t,t)\,dB_t.
$$

Its meaning is

$$
X_t=X_0+\int_0^t a(X_s,s)\,ds+\int_0^t b(X_s,s)\,dB_s.
$$

Now the Brownian weight can change with the random state. The weighted increment sum no longer telescopes. An Itô integral uses suitable nonanticipating approximations, typically left-endpoint values, and a probabilistic limit. Different integral conventions can differ when the integrand is random and varies with the driving noise.

For the OU model, the original diffusion coefficient is constant. Its explicit solution will involve nonconstant but **deterministic** exponential weights. That restricted stochastic-integral theory is enough for our derivation and will be developed in Chapter 11. We do not need to assume general nonlinear Itô calculus just to solve this linear additive-noise model.

### 10.11 Meaning, existence, and uniqueness are separate claims

Writing an integral equation specifies what a solution would have to satisfy. It does not automatically prove that a solution exists or is unique for every choice of coefficients.

One useful sufficient condition in general SDE theory is a Lipschitz condition: the coefficients do not change faster than a fixed multiple of the change in state. For the OU drift,

$$
|\theta(\mu-x)-\theta(\mu-y)|=\theta|x-y|.
$$

The constant diffusion has zero difference between states. The coefficients also have at most linear growth. These properties place OU within standard existence-and-uniqueness results, under the usual initial-data conditions. More concretely, Chapter 12 constructs the solution and proves uniqueness directly for this model, so its solvability will not rest on an unexplained general theorem.

::: {.worked}
**Worked check.** For $dX_t=(2t+1)dt+3dB_t$ and $X_0=4$, the integral form is $X_t=4+\int_0^t(2s+1)ds+3B_t=4+t^2+t+3B_t$. Therefore $X_t\sim N(4+t^2+t,9t)$. Replacing the drift by $2-X_t$ would prevent this direct substitution because the unknown process would then appear inside the time integral. The problem would still be meaningful, but it would require solving a self-consistent equation.
:::

<p class="source-note">Source connection: the integral interpretation of SDEs is discussed in <a href="#ref-l21">MIT Lecture 21</a> and <a href="#ref-ss">Särkkä and Solin, Chapter 4</a>. The examples and exact-versus-Euler distinction are developed here from the integral equations.</p>
