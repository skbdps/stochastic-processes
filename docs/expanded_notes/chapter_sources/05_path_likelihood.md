## 5. The likelihood of a dependent record {#path-likelihood}

### 5.1 State what is observed and what is unknown

Suppose you observe a process at $n+1$ successive recording times and obtain the numerical record

$$
D=(x_0,x_1,\ldots,x_n).
$$

The corresponding uncertain values before observation are $X_0,X_1,\ldots,X_n$. A parameter vector $\eta$ specifies the model's transition probabilities or densities. The actual parameter vector is unknown; for likelihood evaluation, we insert candidate values of $\eta$ and calculate how well each candidate assigns probability or density to this **fixed record**.

For an independent coin sample, the likelihood is a product of individual probabilities. For a process, the values are typically dependent. Today's position helps predict tomorrow's position. Multiplying the individual marginal probabilities $P_\eta(X_i=x_i)$ would ignore that structure.

We need a factorization that remains valid for dependent observations. The tool is the probability chain rule. The Markov property will simplify the chain rule, but the two ideas must be introduced separately.

### 5.2 Derive the product rule from the definition of conditional probability

For events $A$ and $B$ with $P(A)>0$, conditional probability is

$$
P(B\mid A)=\frac{P(A\cap B)}{P(A)}.
$$

Multiplying by $P(A)$ gives

$$
P(A\cap B)=P(A)P(B\mid A).
$$

This is an identity. It does not require independence. If $A$ and $B$ happen to be independent, then $P(B\mid A)=P(B)$ and the formula reduces to $P(A)P(B)$. Without independence, the conditional factor must be retained.

The same identity works inside a world where another event $C$ is already supplied:

$$
P(A\cap B\mid C)=P(A\mid C)P(B\mid A,C).
$$

You can verify it by inserting the definitions and canceling denominators. Supplying a model parameter does not change this algebra; every probability is simply evaluated under the same candidate probability law $P_\eta$.

### 5.3 Build the chain rule for three and then many observations

Begin with three observations. Apply the product rule to the block $(X_0,X_1)$ and the final value $X_2$:

$$
\begin{aligned}
&P_\eta(X_0=x_0,X_1=x_1,X_2=x_2)\\
&\quad=P_\eta(X_0=x_0,X_1=x_1)
 P_\eta(X_2=x_2\mid X_0=x_0,X_1=x_1).
\end{aligned}
$$

Apply it again to the first block:

$$
\begin{aligned}
&P_\eta(X_0=x_0,X_1=x_1,X_2=x_2)\\
&=P_\eta(X_0=x_0)
P_\eta(X_1=x_1\mid X_0=x_0)\\
&\qquad\times P_\eta(X_2=x_2\mid X_0=x_0,X_1=x_1).
\end{aligned}
$$

Each new factor asks for the next observation conditional on everything preceding it. Continuing in this way gives

$$
\boxed{
P_\eta(X_0=x_0,\ldots,X_n=x_n)
=P_\eta(X_0=x_0)
\prod_{i=1}^n P_\eta(X_i=x_i\mid X_0=x_0,\ldots,X_{i-1}=x_{i-1}).
}
$$

For compactness we will sometimes write this as

$$
p_\eta(x_0,\ldots,x_n)
=p_\eta(x_0)\prod_{i=1}^n p_\eta(x_i\mid x_0,\ldots,x_{i-1}).
$$

The lowercase $p_\eta$ here denotes a mass function, not the coin's single scalar heads probability. Its arguments identify what distribution it describes.

The important fact is that the chain rule is always available when the relevant conditional laws are defined. It does not make a process Markov. It simply writes a joint probability as a sequence of conditional probabilities.

### 5.4 Now use the Markov assumption—and only now

A first-order Markov model says that the conditional law of the next state given its entire past equals its conditional law given only its immediately preceding state:

$$
p_\eta(x_i\mid x_0,\ldots,x_{i-1})
=p_\eta(x_i\mid x_{i-1}).
$$

Substitute that equality into each conditional factor of the chain rule:

$$
\boxed{
p_\eta(x_0,\ldots,x_n)
=p_\eta(x_0)\prod_{i=1}^n p_\eta(x_i\mid x_{i-1}).
}
$$

The observations have not become independent. Each factor still depends on the previous observed state. What has disappeared is the need to carry every older observation as a conditioning argument.

This separates the two mathematical ingredients cleanly. **The chain rule justifies multiplication. The Markov property shortens the history in each conditional.** Neither can be replaced by the vague explanation that “the observations are independent.”

### 5.5 Work through one fully specified discrete record

Use a two-state homogeneous chain with

$$
p_{11}=0.7,\quad p_{12}=0.3,\quad p_{21}=0.4,\quad p_{22}=0.6,
$$

where $p_{ij}$ is the probability of moving from state $i$ to state $j$. Suppose the observed record is $(1,1,2,1)$ and the initial distribution assigns probability $1/2$ to state $1$.

The record contains the three observed transitions $1\to1$, $1\to2$, and $2\to1$. Its full probability is

$$
\begin{aligned}
P(X_0=1,X_1=1,X_2=2,X_3=1)
&=P(X_0=1)p_{11}p_{12}p_{21}\\
&=\tfrac12(0.7)(0.3)(0.4)\\
&=0.042.
\end{aligned}
$$

If the starting state is supplied as a condition rather than modeled as random, then

$$
P(X_1=1,X_2=2,X_3=1\mid X_0=1)
=(0.7)(0.3)(0.4)=0.084.
$$

The first calculation asks for the probability of the start **and** the subsequent record. The second asks for the probability of the subsequent record **given** that start. They are different questions, so different numerical answers are appropriate.

### 5.6 Endpoint probability versus observed-record probability

Suppose instead that you observe only $X_0=1$ and $X_2=1$, with the middle state unobserved. There are two possible routes: $(1,1,1)$ and $(1,2,1)$. The conditional endpoint probability is

$$
P(X_2=1\mid X_0=1)=p_{11}^2+p_{12}p_{21}=0.49+0.12=0.61.
$$

If the full route $(1,1,1)$ was observed, its conditional probability would be only $p_{11}^2=0.49$. We would not add the other route because it did not occur in the recorded data.

The rule is: **sum over unobserved alternatives; retain observed intermediate values in the product**. Matrix powers perform these route sums for finite homogeneous chains. A continuous-state transition over a long interval likewise already accounts for all unobserved behavior within that interval; one should not invent unobserved intermediate observations and multiply their densities into a likelihood.

A useful related question is reconstruction of the missing middle state. Take a weather chain with $S\to S=0.8$, $S\to R=0.2$, $R\to S=0.4$, and $R\to R=0.6$. Given $S$ at time zero and $R$ at time two, the chance of rain at the middle time is

$$
\begin{aligned}
P(R_1\mid S_0,R_2)
&=\frac{P(R_1,R_2\mid S_0)}{P(R_2\mid S_0)}\\
&=\frac{0.2(0.6)}{0.8(0.2)+0.2(0.6)}\\
&=\frac{0.12}{0.28}=\frac37.
\end{aligned}
$$

The numerator is the one route consistent with a rainy middle state. The denominator is the sum of all routes consistent with the observations. Every factor has a specific conditioning interpretation.

### 5.7 Move from masses to densities without changing the structure

For continuous observations, the exact-value probability $P(X_i=x_i)$ is zero. The probability of a small interval is obtained by integrating a density over that interval. Therefore likelihood uses a **joint density** evaluated at the recorded values, not the probability of the exact-value event.

The density version of conditional probability is

$$
f_\eta(x,y)=f_\eta(x)f_\eta(y\mid x),
$$

where a conditional density can be obtained by the ratio $f_\eta(x,y)/f_\eta(x)$ wherever the denominator is positive. Repeating the identity gives the same chain factorization. For a Markov process with initial density $h_\eta$ and transition density $g_\eta$,

$$
\boxed{
f_\eta(x_0,\ldots,x_n)
=h_\eta(x_0)\prod_{i=1}^n g_\eta(x_i\mid x_{i-1}).
}
$$

For each fixed current state $x$, the function $y\mapsto g_\eta(y\mid x)$ is a full probability density over the next state: it is nonnegative and integrates to one. At the observed next value $y$, its height supplies the transition's likelihood contribution.

If the transition depends on an elapsed time $\Delta_i$, write $g_\eta(x_i\mid x_{i-1};\Delta_i)$. The time gap is observed information controlling the transition, not an extra unknown parameter unless a separate sampling model is being estimated. We will work through this explicitly for OU observations in Chapter 13.

### 5.8 A continuous example with every likelihood factor visible

Consider the recursion

$$
X_i=25+0.5(X_{i-1}-25)+\varepsilon_i,
\qquad \varepsilon_i\sim N(0,4),
$$

where the innovations $\varepsilon_i$ are mutually independent and independent of the starting value. The notation $N(m,v)$ means a Gaussian distribution with mean $m$ and variance $v$.

Given $X_{i-1}=x$, the only remaining random quantity is $\varepsilon_i$. Therefore

$$
X_i\mid X_{i-1}=x\sim N\bigl(25+0.5(x-25),4\bigr).
$$

The corresponding transition density is

$$
g(y\mid x)=\frac{1}{\sqrt{8\pi}}
\exp\left[-\frac{(y-25-0.5(x-25))^2}{8}\right].
$$

Now observe $(x_0,x_1,x_2)=(29,27,25.5)$ and condition on $x_0=29$. For the first transition, the predicted mean is $25+0.5(29-25)=27$, so the observed error is zero. The factor is $g(27\mid29)=1/\sqrt{8\pi}$.

For the second transition, the **observed** current value is $27$, not an unobserved prediction. The predicted next mean is $25+0.5(27-25)=26$. The actual next observation is $25.5$, so the error is $-0.5$ and its square is $0.25$. Thus

$$
g(25.5\mid27)=\frac{1}{\sqrt{8\pi}}\exp\left[-\frac{0.25}{8}\right].
$$

Multiply the two factors:

$$
\boxed{L=\frac{1}{8\pi}e^{-1/32}.}
$$

This is a joint **conditional density value**, not a dimensionless probability of those exact real numbers. In this calculation the parameters were fixed at $(25,0.5,4)$ simply to evaluate a candidate. For estimation, replace them by $(\mu,\phi,\tau^2)$, keep $(29,27,25.5)$ fixed, and maximize the resulting function of the candidate parameters.

### 5.9 Conditional likelihood does not discard the starting value

The conditional likelihood is

$$
L_c(\eta)=\prod_{i=1}^n g_\eta(x_i\mid x_{i-1}).
$$

The full likelihood is

$$
L_f(\eta)=h_\eta(x_0)L_c(\eta).
$$

Notice that $x_0$ is still used in the first conditional transition of $L_c$. Conditioning removes the **density factor for how the start was generated**, not the start's numerical role in predicting $x_1$.

There are several legitimate setups. If the initial value was fixed by experimental design, no initial random density is needed. If it was drawn from a known distribution independent of $\eta$, its density factor is constant with respect to the optimization and does not alter the maximizing parameter. If it was drawn from a stationary distribution depending on $\eta$, that factor generally changes the full likelihood and may change the estimate.

The fact that $x_0$ is observed does not make $h_\eta(x_0)$ constant in $\eta$. All observations in a likelihood are observed and fixed; their probabilities or densities still vary across candidate models. This is exactly the data-versus-parameter distinction from the coin example.

“Exact transition likelihood” and “full likelihood” also answer different questions. A conditional likelihood can use exact transitions. A full likelihood additionally specifies the initial law. Neither word should silently stand in for the other.

### 5.10 Derive a two-state chain MLE from transition counts

Let $a=p_{11}$ and $b=p_{21}$ be unknown. Then $p_{12}=1-a$ and $p_{22}=1-b$. The unknown parameter vector is $(a,b)$; the observed data are a state sequence. Condition on its start.

Use the concrete record

$$
1,1,2,1,2,2,1,1,1,2,1.
$$

It contains ten transitions. Counting them gives $n_{11}=3$, $n_{12}=3$, $n_{21}=3$, and $n_{22}=1$. Each $1\to1$ contributes $a$, each $1\to2$ contributes $1-a$, and so on. Consequently,

$$
L_c(a,b)=a^3(1-a)^3 b^3(1-b).
$$

Taking logs produces

$$
\ell(a,b)=3\log a+3\log(1-a)+3\log b+\log(1-b).
$$

The terms involving $a$ and $b$ are separate. Differentiating with respect to $a$ while holding $b$ fixed gives

$$
\frac3a-\frac3{1-a}=0
\quad\Longrightarrow\quad \hat a=\frac36=\frac12.
$$

Differentiating in $b$ gives

$$
\frac3b-\frac1{1-b}=0
\quad\Longrightarrow\quad \hat b=\frac34.
$$

Both second derivatives are negative in the interior, so the solutions maximize their respective blocks. The estimate of an outgoing transition probability is the number of observed departures of that type divided by the total number of observed departures from the same state.

For a general finite chain, this yields

$$
\hat p_{ij}=\frac{n_{ij}}{\sum_k n_{ik}},
$$

provided the denominator is positive. If no departure from state $i$ is observed, the conditional likelihood contains no factor involving its outgoing distribution, so that distribution is not identified by these data. It is incorrect to manufacture an estimate by dividing zero by zero.

The transition counts summarize the dependence of this conditional likelihood on the parameters. This is the relevant meaning of a sufficient statistic here. It does not mean that the observed transitions are i.i.d., or that any arbitrary array of counts can come from one path. For example, in a two-state path the difference $n_{12}-n_{21}$ can only be $-1$, $0$, or $1$, depending on the starting and ending states. Our concrete record obeys that constraint.

### 5.11 The logarithm changes computation, not the probability model

For the general Markov likelihood,

$$
\ell_f(\eta)=\log h_\eta(x_0)+\sum_{i=1}^n\log g_\eta(x_i\mid x_{i-1}).
$$

Taking logs converts a product of many factors into a sum. It does not remove dependence, justify a Markov assumption, or change the initial-condition decision. Those were modeling choices made before the logarithm was taken.

The advantage is that each recorded transition contributes one evaluable term. A long record does not force the conditioning history of each term to grow, provided the process is truly Markov in the recorded state. If the state is incomplete, such as a noisy measurement of an unobserved Markov state, the recorded series need not inherit this simple factorization using the latent state's transition density. A different observation model may be needed.

The Gaussian AR(1) model supplies a concrete transition density for which this likelihood can be derived and optimized by hand. That is the subject of the next two chapters.

<p class="source-note">Source connection: <a href="#ref-l5">MIT Lecture 5</a> supplies the Markov framework, and <a href="#ref-l8">MIT Lecture 8</a> discusses time-series likelihoods. Every factorization and numerical example above is developed from the stated conditional laws.</p>
