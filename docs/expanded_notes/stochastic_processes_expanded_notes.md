## Before the first chapter {#orientation}

A stochastic process is a model for quantities that evolve with randomness over time. This book develops the ideas needed to understand one particular example: a process that tends to return toward a level while receiving continuing random disturbances. The eventual goal is not just to write down that model, but to understand how observations at unequal time gaps can be used to estimate its parameters.

The route contains thirteen chapters. It begins with the meaning of a likelihood, constructs dependent observations from independent random steps, explains how a Markov state summarizes predictive information, and introduces the Gaussian AR(1). It then develops the continuous-time side: scaling the clock, Brownian motion, the integral meaning of an SDE, deterministic Wiener integration, and the OU solution. The last chapter brings these pieces together in a fully specified estimation problem.

The exposition distinguishes four kinds of statements. A **definition** gives a precise meaning to a word or symbol. A **model assumption** specifies the mechanism we choose to study. A **derivation** obtains a consequence of those assumptions. An **approximation** replaces an exact object by a simpler one and must be identified as such. Confusing these categories is one reason stochastic-process formulas can feel disconnected. For example, the OU equation is a model specification; its Gaussian transition is a derived consequence; the Euler transition is an approximation; the MLE is an estimation rule applied to observed data.

Each chapter begins with the object and question it studies. Symbols are introduced where they are used, and the needed probability identities are restated or derived locally. The chapters use elementary algebra and single-variable calculus. Calculus operations are written out in the derivations; expectations, covariance, conditional distributions and matrix operations are introduced through the particular calculations that use them. The text does not require a prior course in stochastic processes or an earlier conversation to supply missing definitions. For the general continuous-time foundations, some major theorems are stated rather than proved: these include the central limit theorem, existence of Brownian motion, its nowhere-differentiability theorem, and completeness results used in defining the Wiener integral. The text identifies those dependencies rather than presenting a heuristic as a proof. The main model calculations—likelihoods, recursions, moments, transitions and profile estimates—are worked through explicitly.

The full derivations are part of the visible chapter text. The navigation can take you directly to a topic, and each chapter also has its own list of subsections. Worked checks contain both a question and its resolution; covering the resolution and reproducing the steps is more useful than memorizing the final equation. The reference card near the end is for returning to a result after understanding its derivation, not for replacing the chapters.

### A notation convention used throughout

Capital letters such as $X_t$ usually mean random variables before their values are observed. Lowercase letters such as $x_i$ mean recorded numerical values. A hat, as in $\widehat\theta$, indicates an estimate. During likelihood evaluation, an unhatted parameter such as $\theta$ is a candidate value being tried; the true unknown parameter can be denoted $\theta_0$ when that distinction matters.

$N(m,v)$ means a Gaussian distribution with mean $m$ and **variance** $v$. Its standard deviation is $\sqrt v$. The Brownian motion symbol is $B_t$; $W_k$ is reserved for a discrete random walk. A numerical discretization tick is $\delta$, while a gap between actual observations is $\Delta_i$. These two durations need not be equal.

For the OU model, $\theta$ is the positive return speed, $\mu$ is the return level, and $\sigma$ is the positive continuous-time noise scale. The full parameter vector is $\eta=(\theta,\mu,\sigma)$. For AR(1), $\phi$ is the retained fraction of the previous deviation, and $\tau^2$ is the innovation variance. The notation is repeated when the models are introduced, so the table is not a prerequisite lookup exercise.


## 1. What a likelihood measures {#likelihood}

### 1.1 The question comes before the formula

You have a coin, and you do not know its probability of producing heads. You toss it ten times and record

$$
H,H,T,H,T,H,H,H,T,H.
$$

There are seven heads and three tails. Those results have already happened. You can count them, write them down, or store them in a file. They are the **data**.

The coin also has a property that you have not directly observed: its underlying probability of heads. Call that probability $p$. A coin with $p=0.5$ is fair; a coin with $p=0.8$ favors heads. Thus $p$ measures the coin's bias, rather than merely labeling it “fair” or “unfair.” The probability of tails is $1-p$, so we do not need a second independent parameter for tails.

It is important not to identify $p$ with the observed fraction of heads before doing any reasoning. The number $7/10$ describes this particular record. The number $p$ describes the mechanism that produced the record. A fair coin can produce seven heads in ten tosses, and a biased coin can produce five. Randomness prevents one short record from revealing the true probability with certainty.

Our task is therefore: **use the observed record to estimate the unobserved probability $p$**. Maximum likelihood is one rule for doing that. Before applying it, we need a probability model telling us what records different coins would tend to produce.

### 1.2 Data, random variables, parameters, and assumptions

Represent toss $i$ by a random variable $Y_i$:

$$
Y_i=\begin{cases}1,&\text{if toss }i\text{ is heads},\\0,&\text{if toss }i\text{ is tails}.\end{cases}
$$

Before a toss, $Y_i$ is uncertain. After observing it, its recorded value is a number $y_i$, either zero or one. Capital letters will usually denote random variables; lowercase letters will denote their observed values.

The model makes two assumptions. First, each toss has the same probability of heads, $P_p(Y_i=1)=p$. Second, the tosses are independent under that model: learning the outcomes of some tosses does not change the probabilities of the others when $p$ is held fixed. The subscript in $P_p$ labels the probability law belonging to a candidate coin with parameter $p$.

The first assumption is **identical distribution**. The second is **independence**. Together they are often abbreviated “independent and identically distributed,” or i.i.d. The abbreviation names two separate requirements, not one vague assertion of randomness.

| Object | Value or definition | Role |
|---|---|---|
| $n$ | $10$ | Number of tosses in the observed experiment |
| $(y_1,\ldots,y_{10})$ | $(1,1,0,1,0,1,1,1,0,1)$ | Ordered observed data |
| $k=\sum_i y_i$ | $7$ | Observed number of heads |
| $p$ | Any candidate in $[0,1]$ | Parameter we vary while comparing models |
| $p_0$ | Unknown | The actual probability of heads, assuming the model is correct |
| $\hat p$ | To be calculated | Estimate selected from the observed data |

We know the **form** of the probability model without knowing its true parameter value. This is no more contradictory than knowing that a line has the form $y=ax+b$ without yet knowing its slope and intercept.

### 1.3 How can we calculate a probability “given $p$” when $p$ is unknown?

The expression $P_p(\text{data})$ does not require us to know the true $p_0$. It asks a hypothetical question:

> Suppose the coin had this candidate value of $p$. What probability would that model assign to the record we observed?

We can answer that question at $p=0.3$, then at $p=0.5$, then at $p=0.7$, without knowing which candidate, if any, is the truth. The calculation is conditional on a proposed model, not on knowledge that the proposed model is correct.

A simple algebra analogy is useful. You can write $f(a)=a^2-4a+5$ and find which $a$ minimizes it even though you did not know that minimizing value in advance. The unknown is the input you are searching over. Knowing how a function depends on its input is different from already knowing which input is best.

In maximum likelihood, the analogous function is the **likelihood**:

$$
L(p;D)=P_p(D).
$$

Here $D$ denotes the fixed observed data. Writing $L(p;D)$ emphasizes that $p$ is the variable of the function and $D$ is what we have already observed. You will also see $P(D\mid p)$. In this frequentist setting, that notation is shorthand for “the probability of $D$ under the model labeled by $p$.” It does not mean that we have placed a probability distribution on $p$ itself.

### 1.4 Derive the likelihood of the ordered record

Under a candidate coin, a head contributes a factor $p$ and a tail contributes a factor $1-p$. Independence lets us multiply the probabilities of the ten individual outcomes. For our specific sequence,

$$
\begin{aligned}
L_{\mathrm{sequence}}(p)
&=p\cdot p\cdot(1-p)\cdot p\cdot(1-p)\\
&\qquad\cdot p\cdot p\cdot p\cdot(1-p)\cdot p\\
&=p^7(1-p)^3.
\end{aligned}
$$

The exponents are not mysterious model parameters. They are counts extracted from the data. Seven observed heads give seven factors of $p$; three observed tails give three factors of $1-p$.

For a general binary record of length $n$, each factor can be written as $p^{y_i}(1-p)^{1-y_i}$. Check the two cases: when $y_i=1$, the factor is $p$; when $y_i=0$, it is $1-p$. Thus

$$
\begin{aligned}
L_{\mathrm{sequence}}(p)&=\prod_{i=1}^n p^{y_i}(1-p)^{1-y_i}\\
&=p^{\sum_i y_i}(1-p)^{\sum_i(1-y_i)}\\
&=p^k(1-p)^{n-k}.
\end{aligned}
$$

The product symbol $\prod$ means “multiply all the indicated factors,” just as $\sum$ means “add all the indicated terms.” We have written out the multiplication before compressing it into notation.

### 1.5 Why does another version contain $\binom{10}{7}$?

The event “this exact ordered sequence occurs” is different from the event “exactly seven heads occur, in any order.” To calculate the second probability, we must add the probabilities of all sequences with seven heads.

How many such sequences exist? Choose which seven of the ten positions hold heads. The number of choices is

$$
\binom{10}{7}=\frac{10!}{7!\,3!}=120.
$$

The factorial $m!$ means $m(m-1)\cdots 1$. Every one of these 120 sequences has probability $p^7(1-p)^3$. They cannot occur simultaneously, so the probability of their union is the sum:

$$
L_{\mathrm{count}}(p)=P_p(K=7)=120p^7(1-p)^3,
$$

where $K=Y_1+\cdots+Y_{10}$ is the random count before the experiment. This is a binomial probability.

The sequence likelihood and the count likelihood have different heights. However, the multiplier 120 does not depend on $p$. Multiplying all candidates' scores by the same positive number does not change which candidate ranks highest. Therefore both versions produce the same maximum-likelihood estimate of $p$.

Do not add the binomial factor merely because you see repeated observations in any likelihood. It appears here because we have changed the event from one sequence to a collection of sequences. In a dependent time series, changing the order can change the likelihood itself.

### 1.6 What “argmax” asks us to return

For the count likelihood, candidate comparisons give approximately:

| Candidate $p$ | $P_p(K=7)$ |
|---:|---:|
| $0.1$ | $0.000008748$ |
| $0.3$ | $0.009001692$ |
| $0.5$ | $0.117187500$ |
| $0.7$ | $0.266827932$ |
| $0.9$ | $0.057395628$ |

Each row concerns the **same event**, $K=7$. Only the candidate coin changes. Among these candidates, $0.7$ gives the largest score. To establish the maximum over the whole interval $[0,1]$, rather than merely over this small table, we will use calculus.

The notation distinguishes the answer we want from the score it achieves:

$$
\max_{p\in[0,1]}L(p)
\quad\text{is the greatest likelihood value},
$$

whereas

$$
\operatorname*{argmax}_{p\in[0,1]}L(p)
\quad\text{is the input where that value occurs}.
$$

For this example, the first is about $0.2668$ and the second is $0.7$. We are estimating a probability-of-heads parameter, so the second is the useful answer. In general there can be several maximizers or no attained maximum. The notation does not guarantee uniqueness or existence. We will check both in this example.

### 1.7 Why logarithms preserve the estimate

The natural logarithm is strictly increasing on positive numbers. If $L(p_1)>L(p_2)>0$, then $\log L(p_1)>\log L(p_2)$. Thus the log preserves the ranking of candidate parameters. It changes the heights, but not the location of the highest point.

The logarithm also transforms products into sums:

$$
\log(ab)=\log a+\log b,
\qquad \log(a^r)=r\log a.
$$

Applying these identities gives the log-likelihood

$$
\ell(p)=\log 120+7\log p+3\log(1-p),\qquad 0<p<1.
$$

The term $\log 120$ can be omitted while locating the maximizer because its value is the same for every $p$. Its derivative is zero. Terms depending on the parameter cannot be discarded this way.

There is also a numerical reason for working with logs. Multiplying many small probabilities can produce a number too close to zero for a computer to represent. Adding their logarithms avoids that particular failure. This argument does not imply that every likelihood factor must be below one: continuous **density** values, introduced below, can exceed one.

### 1.8 Derive the maximizing value, one differentiation at a time

A derivative measures how the score changes when we change the candidate parameter by a small amount. A positive derivative means the score is increasing locally; a negative derivative means it is decreasing. At a smooth peak inside the permitted interval, the first-order rate of change must be zero. That gives a candidate location for the maximum, which we will then check.

Differentiate the log-likelihood with respect to the candidate $p$. The derivative of $\log p$ is $1/p$. For $\log(1-p)$, the chain rule supplies an extra factor of $-1$:

$$
\frac{d}{dp}\log(1-p)=\frac{1}{1-p}\frac{d}{dp}(1-p)=-\frac{1}{1-p}.
$$

Consequently,

$$
\ell'(p)=\frac{7}{p}-\frac{3}{1-p}.
$$

At an interior maximum, a differentiable function must have derivative zero. Solving that equation gives

$$
\begin{aligned}
\frac{7}{p}-\frac{3}{1-p}&=0,\\
\frac{7}{p}&=\frac{3}{1-p},\\
7(1-p)&=3p,\\
7-7p&=3p,\\
7&=10p,\\
p&=0.7.
\end{aligned}
$$

The cross-multiplication is valid because $p$ and $1-p$ are positive in the interior of the parameter interval. It has located a stationary point; we still need to show it is the maximum.

Differentiate again:

$$
\ell''(p)=-\frac{7}{p^2}-\frac{3}{(1-p)^2}<0,\qquad 0<p<1.
$$

A strictly negative second derivative means the log-likelihood is strictly concave: its slope decreases continuously, so it can turn from rising to falling at most once. Our stationary point is therefore its unique interior maximum. At $p=0$, observing heads is impossible, and at $p=1$, observing tails is impossible. Both endpoint likelihoods are zero. The interior likelihood is positive. Thus no endpoint defeats the interior solution:

$$
\boxed{\hat p_{\mathrm{MLE}}=0.7.}
$$

For $0<k<n$, exactly the same differentiation gives

$$
\begin{aligned}
\frac{k}{p}-\frac{n-k}{1-p}&=0,\\
k(1-p)&=(n-k)p,\\
\hat p&=k/n.
\end{aligned}
$$

When $k=0$, the likelihood $(1-p)^n$ decreases with $p$, so the maximum is at $p=0$. When $k=n$, the likelihood $p^n$ increases, so the maximum is at $p=1$. Those are boundary solutions. If the parameter space were restricted to the open interval $(0,1)$, these two records would have a supremum but no attained MLE in that space. The parameter space matters.

### 1.9 Why software often minimizes instead

Many optimization routines are written as minimizers. Define the negative log-likelihood $J(p)=-\ell(p)$. Every peak of $\ell$ becomes a valley of $J$ at the **same horizontal location**. Therefore

$$
\boxed{\operatorname*{argmax}_{p}L(p)=\operatorname*{argmax}_{p}\ell(p)=\operatorname*{argmin}_{p}[-\ell(p)].}
$$

The minus sign acts on the function, not on the answer. The estimate remains $+0.7$; it does not become $-0.7$.

In a model with several unknowns, the same rule searches over a vector. For an OU process we will use $\eta=(\theta,\mu,\sigma)$, with $\theta>0$, $\mu\in\mathbb R$, and $\sigma>0$. The notation $\operatorname{argmax}_\eta$ then means searching over all three components simultaneously, subject to those restrictions. The basic idea is unchanged: keep the observed data fixed and compare candidate parameter values.

### 1.10 What changes when the data are continuous?

Suppose a measuring instrument records a real-valued quantity. A Gaussian model has density

$$
f(y;m,v)=\frac{1}{\sqrt{2\pi v}}\exp\left[-\frac{(y-m)^2}{2v}\right],\qquad v>0.
$$

Here $m$ is the mean and $v$ is the variance. The notation $N(m,v)$ uses **variance**, not standard deviation, in its second slot. The standard deviation is $\sqrt v$.

For a variable with a density, a particular exact value has probability zero:

$$
P(Y=y)=\int_y^y f(u;m,v)\,du=0.
$$

This is not a statement that an observation cannot occur. It means that probability is assigned to intervals by integrating the density, rather than by assigning positive mass to each point. For a small interval of width $h$ around an observed $y$, and a density continuous there,

$$
P(y\le Y\le y+h)\approx f(y;m,v)h.
$$

To compare candidates using the same observed interval, its width $h$ is common to all candidates. The density value therefore supplies the relevant relative score as the interval becomes small. For independent observations $y_1,\ldots,y_n$, the continuous likelihood is

$$
L(m,v)=\prod_{i=1}^n f(y_i;m,v).
$$

A density value can exceed one because probability is **area**, not height. For example, a uniform density on an interval of length $0.1$ has height $10$ and total area $10\times0.1=1$. A likelihood made from densities is not itself a point probability.

The Gaussian log-likelihood illustrates two effects that later reappear in every OU transition:

$$
\ell(m,v)=-\frac n2\log(2\pi v)-\frac{1}{2v}\sum_{i=1}^n(y_i-m)^2.
$$

The squared-error term penalizes predictions far from the observations. The $\log v$ term penalizes spreading the density too broadly. Omitting it would let us make the squared-error penalty arbitrarily small simply by sending $v$ to infinity. The variance-normalization term is part of the statistical model, not cosmetic algebra.

For fixed $v$, differentiating in $m$ gives

$$
\frac{\partial\ell}{\partial m}=\frac1v\sum_i(y_i-m)=0
\quad\Longrightarrow\quad \hat m=\frac1n\sum_i y_i.
$$

Writing $S=\sum_i(y_i-\hat m)^2$ and differentiating in $v$ gives

$$
\frac{\partial\ell}{\partial v}=-\frac{n}{2v}+\frac{S}{2v^2}=0
\quad\Longrightarrow\quad \hat v=S/n,
$$

provided $S>0$. This derivation explains why a Gaussian variance MLE has denominator $n$, not the $n-1$ associated with a different goal, unbiased variance estimation. If all the observations are identical, $S=0$ and this unconstrained Gaussian likelihood can increase without bound as $v$ approaches zero. Again, writing an optimization problem does not guarantee a regular interior solution.

### 1.11 What the estimate does—and does not—tell us

The statement $\hat p=0.7$ means that the candidate coin with heads probability $0.7$ maximizes the likelihood of our observed record within the assumed model. It does not mean that the probability the true parameter equals $0.7$ is $70\%$. It does not guarantee seven heads in the next ten tosses. It does not, by itself, constitute a hypothesis test rejecting fairness.

An **estimator** is the rule $\hat p(D)=K/n$ before the data are observed. Because $K$ varies across experiments, the estimator is a random variable across experiments. An **estimate** is its realized value, here $0.7$. This distinction is how we later discuss bias and uncertainty: repeat the experiment conceptually and examine how the estimates vary.

For the independent coin model, $E[K]=np_0$ and $\operatorname{Var}(K)=np_0(1-p_0)$. Dividing by $n$ gives

$$
E[\hat p]=p_0,\qquad \operatorname{Var}(\hat p)=\frac{p_0(1-p_0)}{n}.
$$

Thus this particular MLE is unbiased. That conclusion comes from calculating its repeated-sampling expectation, not from its being an MLE. Other maximum-likelihood estimators can be biased. Chapter 2 derives the expectation and variance rules used in this calculation.

::: {.worked}
**Worked check: identify every object before differentiating.** A coin is tossed twenty times and produces twelve heads. The observed count is $k=12$, the fixed sample size is $n=20$, and the unknown parameter is the true heads probability $p_0$. The candidate likelihood is $L(p)=\binom{20}{12}p^{12}(1-p)^8$. Taking logs and differentiating gives $12/p-8/(1-p)=0$, hence $12-12p=8p$, and therefore $\hat p=0.6$. The maximizing parameter is $0.6$; the maximum likelihood value is whatever number results from inserting $0.6$ into $L$. They answer different questions.
:::

The idea to carry forward is precise: **a likelihood is an evaluable function of hypothetical parameter values, even while the actual parameter remains unknown**. For a time series, our next challenge is not this optimization principle. It is constructing the correct likelihood when observations depend on one another.

<p class="source-note">Source connection: see <a href="#ref-psu">Penn State STAT 504, §1.5</a> for the Bernoulli/binomial likelihood.  the likelihood viewpoint and Gaussian time-series likelihood are covered in <a href="#ref-l8">MIT Lecture 8</a>. The coin calculations and the differentiation above are worked explicitly here.</p>


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


## 3. Think one step ahead {#first-step}

### 3.1 State the experiment and the event separately

A gambler begins with $m$ units of fortune, where $m$ is an integer between $0$ and a target $N$. On each round the fortune increases by one with probability $p$ or decreases by one with probability $q=1-p$. The rounds are independent and use the same probabilities. The game stops as soon as the fortune reaches either $0$ or $N$.

For this chapter, the unknown is **not** the coin parameter $p$. Treat $p$ and $N$ as specified features of the game. The unknown we want to calculate is a probability of an event: reaching $N$ before reaching $0$.

Let $X_k$ denote the fortune after $k$ rounds, with $X_0=m$. Define the stopping time

$$
T=\inf\{k\ge0:X_k\in\{0,N\}\}.
$$

The notation $\inf$ means the earliest member of the indicated set; if the boundary is never reached, take $T=\infty$. Calling $T$ a stopping time means that, after watching the first $k$ rounds, we can tell whether the game has stopped by round $k$. No knowledge of future rounds is required to make that decision.

Now define

$$
f(m)=P(X_T=N\mid X_0=m).
$$

The input of $f$ is a **starting fortune**, not a time. The output is a **success probability**, not an expected fortune. Distinguishing those two roles prevents several tempting but incorrect recursions.

The boundaries require no calculation:

$$
f(0)=0,\qquad f(N)=1.
$$

Starting at ruin means success has already failed. Starting at the target means success has already occurred.

### 3.2 Why the game really does end

Before setting success plus failure equal to one, it is worth checking that the game cannot avoid both boundaries forever with positive probability. Assume $0<p<1$. From any interior fortune, a sequence of $N$ upward steps would reach the target within at most $N$ rounds. The probability of that specified sequence is $p^N>0$.

Thus, conditional on any history that has not yet stopped, the probability of stopping within the next $N$ rounds is at least $p^N$. Consequently, the probability of surviving $r$ complete blocks of $N$ rounds is at most

$$
P(T>rN)\le(1-p^N)^r.
$$

The right side tends to zero as $r$ increases. Therefore $P(T<\infty)=1$: the game ends almost surely. “Almost surely” means with probability one, not that we have found a fixed deterministic upper bound on its duration.

The same geometric bound also ensures a finite expected duration. In the special cases $p=1$ and $p=0$, the game ends deterministically at the upper and lower boundary, respectively. We can therefore discuss success and ruin as complementary outcomes.

### 3.3 Derive the first-step equation from total probability

At an interior fortune $m$, the next round has only two possible outcomes. Split the success event according to which occurs first. The law of total probability gives

$$
\begin{aligned}
f(m)
&=P(\text{success}\mid X_0=m,\text{first step up})\,p\\
&\quad+P(\text{success}\mid X_0=m,\text{first step down})\,q.
\end{aligned}
$$

After an upward first step, the fortune is $m+1$. All remaining steps are independent of the first step and have the same probabilities as before. Therefore the remaining experiment has exactly the law of a fresh game started at $m+1$. Its success probability is $f(m+1)$. The downward case gives $f(m-1)$.

Substitution yields

$$
\boxed{f(m)=p f(m+1)+q f(m-1)},
\qquad 1\le m\le N-1.
$$

This equation conditions on the **first step out of the current state**. We did not assume that eventual success is independent of that step; indeed its probability depends on the step. We used total probability to average the two conditional success probabilities.

Why is there no time argument, such as $f(k,m)$? Because the future game has the same rules regardless of the round number. If the upward probability changed with time, the success probability from fortune $m$ could depend on the current round, and the time-free recursion would not generally be valid.

### 3.4 Solve the fair case without guessing the answer

For a fair game, $p=q=1/2$. The recursion becomes

$$
f(m)=\frac12f(m+1)+\frac12f(m-1).
$$

Multiply by two and rearrange:

$$
2f(m)=f(m+1)+f(m-1),
$$

so

$$
f(m+1)-f(m)=f(m)-f(m-1).
$$

Define a successive difference

$$
d_m=f(m)-f(m-1).
$$

The rearranged equation states $d_{m+1}=d_m$. Every successive difference is therefore the same number, say $c$. Starting from $f(0)=0$,

$$
f(1)=c,\quad f(2)=2c,\quad\ldots,\quad f(m)=mc.
$$

The other boundary says $1=f(N)=Nc$, hence $c=1/N$. We have obtained, rather than assumed,

$$
\boxed{f(m)=\frac mN.}
$$

At fortune $30$ with boundaries $0$ and $100$, success has probability $0.3$ and ruin has probability $0.7$. At $m=0$ the formula gives zero; at $m=N$ it gives one; increasing the initial fortune increases success probability. These are elementary but useful checks on a formula's direction.

A function equal to the weighted average of its neighboring values is called harmonic for the corresponding walk. In this one-dimensional fair example, the harmonic equation reduces to equal successive differences. We do not need more advanced harmonic-function theory to solve it.

### 3.5 Why a fair game can still have an unequal chance of success

A fair step has expected gain zero. That is a statement about the expected change on **one round**. It does not say that every pair of terminal outcomes has equal probability.

Starting with $30$ and stopping at $0$ or $100$, failure loses $30$ units, while success gains $70$. The success probability $0.3$ balances those unequal outcomes:

$$
0.3(70)+0.7(-30)=21-21=0.
$$

Equivalently, the expected terminal fortune is $0.3(100)+0.7(0)=30$. This is consistent with fairness. It is a check on the first-step calculation, not a replacement for the derivation or an unrestricted use of a stopping theorem.

### 3.6 Derive the biased-game formula

For $0<p<1$, return to

$$
f(m)=p f(m+1)+q f(m-1),\qquad p+q=1.
$$

Replace the coefficient of $f(m)$ by $p+q$ and collect terms:

$$
p[f(m+1)-f(m)]=q[f(m)-f(m-1)].
$$

Using the same differences $d_m=f(m)-f(m-1)$, we get

$$
d_{m+1}=\frac qp d_m.
$$

Let $r=q/p$. Successive differences are now geometric:

$$
d_1=c,\quad d_2=cr,\quad d_3=cr^2,\quad d_j=cr^{j-1}.
$$

Sum the differences from fortune zero to fortune $m$:

$$
f(m)-f(0)=\sum_{j=1}^m d_j
=c(1+r+\cdots+r^{m-1}).
$$

For $r\ne1$, the finite geometric sum follows by subtracting $r$ times the sum from itself:

$$
(1-r)(1+r+\cdots+r^{m-1})=1-r^m.
$$

Thus $f(m)=c(1-r^m)/(1-r)$. The boundary $f(N)=1$ gives $c=(1-r)/(1-r^N)$. Substituting $c$ eliminates the unknown difference:

$$
\boxed{f(m)=\frac{1-(q/p)^m}{1-(q/p)^N}},
\qquad p\ne\frac12.
$$

When $p=1/2$, $r=1$ and the displayed expression has the indeterminate form $0/0$. That does not mean the success probability is undefined; it means the geometric-sum shortcut divided by $1-r$. Use the preceding sum directly, where all the differences are equal, to recover $m/N$.

For a numerical example, take $p=2/3$, $q=1/3$, $N=3$, and $m=1$. Then $q/p=1/2$ and

$$
f(1)=\frac{1-1/2}{1-1/8}=\frac{1/2}{7/8}=\frac47.
$$

At $m=2$, the probability is $(1-1/4)/(1-1/8)=6/7$. Check the first-step equation at $m=1$: $(2/3)(6/7)+(1/3)(0)=4/7$, as required.

### 3.7 Backward questions and forward questions use different equations

It is easy to reverse the neighboring coefficients because a different, equally legitimate question uses arrivals rather than departures.

Our success function $f(m)$ asks: starting at $m$, what happens eventually? Its recursion follows the next departure:

$$
f(m)=p f(m+1)+q f(m-1).
$$

Now consider an unrestricted walk and let $g_k(m)=P(X_k=m)$ be the probability of occupying $m$ at time $k$. To be at $m$ next round, the walk must arrive from $m-1$ by an upward step or from $m+1$ by a downward step:

$$
g_{k+1}(m)=p g_k(m-1)+q g_k(m+1).
$$

This is a forward distribution update. The coefficient $p$ multiplies the lower neighboring position because upward moves arrive from below. In the success equation, $p$ multiplies $f(m+1)$ because an upward first move takes the current fortune to the higher neighbor.

The fair case hides the distinction because both coefficients are $1/2$. A biased example exposes it. Before writing a recursion, state whether your function measures a future outcome **from** a current state or current occupancy **arriving at** a state.

### 3.8 Use the same method to derive an expected duration

First-step analysis can calculate expectations as well as probabilities. Let

$$
u(m)=E[T\mid X_0=m]
$$

be the expected number of rounds until absorption. At a boundary, $u(0)=u(N)=0$: the game has already ended. In the interior, the first round consumes one unit of time. After that round, the expected remaining duration is $u(m+1)$ or $u(m-1)$. For the fair game,

$$
u(m)=1+\frac12u(m+1)+\frac12u(m-1).
$$

Rearranging gives

$$
u(m+1)-2u(m)+u(m-1)=-2.
$$

The left side is a second difference. A linear function has second difference zero. A quadratic $m^2$ has second difference two, because

$$
(m+1)^2-2m^2+(m-1)^2=2.
$$

Therefore $-m^2$ supplies the required $-2$. Add a general linear part, which does not change the second difference:

$$
u(m)=-m^2+am+b.
$$

The boundary $u(0)=0$ gives $b=0$. The boundary $u(N)=0$ gives $-N^2+aN=0$, hence $a=N$. Thus

$$
\boxed{u(m)=m(N-m).}
$$

For $m=30$ and $N=100$, the expected duration is $2100$ rounds. Individual games can end much sooner or much later. The expectation is an average over possible games, not a promised stopping time.

Why does verifying this candidate settle the solution? The difference of any two solutions satisfies the homogeneous equation $h(m)=[h(m+1)+h(m-1)]/2$ with zero boundary values. We already solved that equation: the differences are constant and the zero boundaries force the whole function to be zero. Thus the candidate is unique. The finite-expectation argument above ensures the actual expected duration is among these solutions.

### 3.9 Shifted boundaries and a full check

Suppose the boundaries are $-A$ and $B$, where $A,B>0$ are integers, and the current fortune is $m$. Add $A$ to every fortune. The lower boundary becomes zero, the upper boundary becomes $A+B$, and the starting point becomes $m+A$.

For a fair game,

$$
P(\text{hit }B\text{ before }-A\mid X_0=m)
=\frac{m+A}{A+B}.
$$

From zero, success is $A/(A+B)$ and ruin is $B/(A+B)$. The expected duration is

$$
(m+A)\bigl[(A+B)-(m+A)\bigr]=(m+A)(B-m).
$$

::: {.worked}
**Worked check.** Begin at zero, stop at $-2$ or $3$, and use fair steps. The success probability is $2/5$. The expected number of rounds is $2\cdot3=6$. Starting instead at $m=1$, the success probability becomes $3/5$ and the expected duration becomes $3\cdot2=6$. The equal durations do not imply equal success probabilities; one quantity measures time, the other measures the terminal event.
:::

The method is now reusable. Define the quantity you want as a function of the present state, supply the boundary values, condition on the next move, and use the same function for the remaining task. A Markov model is the general setting in which this “present state is enough” step can be formalized.

<p class="source-note">Source connection: the fair boundary-hitting problem and first-step recursion appear in <a href="#ref-l5">MIT Lecture 5, section 2</a>. The biased-case and duration calculations above are direct extensions derived in full.</p>


## 4. What the present remembers {#markov}

### 4.1 Begin with a prediction problem

Suppose a machine is inspected once per day. At each inspection it is either working, labeled state $1$, or broken, labeled state $2$. Let $X_n$ denote its state on day $n$. A possible observed record is $1,1,2,1,2$.

To predict tomorrow, you could ask for the conditional probability

$$
P(X_{n+1}=j\mid X_0=x_0,\ldots,X_n=x_n).
$$

The vertical bar means that the observations to its right are supplied information. The expression asks how likely tomorrow's state $j$ is after seeing the entire recorded history. A general stochastic process may genuinely require much of that history.

A first-order Markov model makes a simplifying assertion: once today's state is known, the older states provide no additional information for tomorrow's distribution. Formally,

$$
\boxed{
P(X_{n+1}=j\mid X_0=x_0,\ldots,X_n=x_n)
=P(X_{n+1}=j\mid X_n=x_n).
}
$$

For discrete states, this equality is required for histories with positive probability. In a continuous-state formulation it is an equality of appropriate conditional distributions, up to the usual probability-zero qualifications.

“First-order” means one current state is sufficient; it does not mean the process is independent across time. If tomorrow's probabilities change with today's state, the consecutive states are dependent. Markovness concerns whether adding **older** history improves the prediction after the present is already supplied.

### 4.2 Conditional independence is different from independence

Consider this transition rule:

$$
P(X_{n+1}=1\mid X_n=1)=0.7,
\qquad
P(X_{n+1}=1\mid X_n=2)=0.4.
$$

The next state depends on the current one because the probabilities $0.7$ and $0.4$ differ. A Markov assumption adds that yesterday's state does not alter either probability once today's state is specified. For example,

$$
P(X_{n+1}=1\mid X_n=1,X_{n-1}=2)=0.7
$$

and

$$
P(X_{n+1}=1\mid X_n=1,X_{n-1}=1)=0.7.
$$

If both histories are possible, these equalities say that today's state summarizes their predictive effect. They do not say yesterday has no relation to tomorrow. Without conditioning on today, yesterday changes the distribution of today, which then changes the distribution of tomorrow.

This distinction is called **conditional independence**: the next state and the earlier past are independent conditional on the present state. The condition is essential. Removing it changes the assertion.

### 4.3 Choosing the state is part of choosing the model

Suppose a traffic signal stays green for two ticks, yellow for two ticks, and red for two ticks. Seeing only the color “green” does not reveal whether it will remain green or change next tick. You also need to know which of the two green ticks is current.

A rigorous example can randomize the initial phase uniformly over the six phases. At a given time, the history “red, green” means this is the first green tick, so the next color is green. The history “green, green” means this is the second green tick, so the next color is yellow. Both histories have the same current color but different next-color predictions. The color process is not first-order Markov under this initialization.

Enlarge the state to include the phase: $(\text{color},\text{first or second tick})$. The next enlarged state is now determined by the current enlarged state. The enlarged process is Markov.

The same issue arises in numerical models. If

$$
X_{n+1}=X_n+0.8(X_n-X_{n-1})+\varepsilon_{n+1},
$$

with fresh independent noise, then today's level alone is generally insufficient. Today's level and yesterday's level specify the momentum term. The pair $Z_n=(X_n,X_{n-1})$ is a natural Markov state because

$$
Z_{n+1}=(X_n+0.8(X_n-X_{n-1})+\varepsilon_{n+1},\ X_n)
$$

depends only on $Z_n$ and the fresh noise. Calling a process Markov is therefore a statement about a specified state representation, not about a system independently of how it is recorded.

### 4.4 The transition probabilities are conditional distributions

For a time-homogeneous finite-state chain, write

$$
p_{ij}=P(X_{n+1}=j\mid X_n=i).
$$

The first subscript is the departure state; the second is the destination. **Time-homogeneous** means that this number does not change with the day $n$. A process can be Markov without being homogeneous: it may use a different transition rule on different days.

For the two-state example, complete the missing probabilities by making each outgoing distribution sum to one:

| Current state | Next is $1$ | Next is $2$ |
|---|---:|---:|
| $1$ | $p_{11}=0.7$ | $p_{12}=0.3$ |
| $2$ | $p_{21}=0.4$ | $p_{22}=0.6$ |

Each row here is a conditional probability distribution. Therefore $p_{i1}+p_{i2}=1$ for each current state $i$.

We will retain the **column matrix convention** used in MIT Lecture 5. Store the outgoing distribution from state $i$ in column $i$:

$$
A_{ji}=p_{ij},\qquad
A=\begin{pmatrix}0.7&0.4\\0.3&0.6\end{pmatrix}.
$$

The first column describes departures from state $1$; the second describes departures from state $2$. Each column sums to one. Many books instead store $P_{ij}=p_{ij}$ in a row-stochastic matrix $P$. That matrix is $A^\mathsf T$, the transpose. Neither convention is more correct, but mixing the probability-vector convention from one with the matrix from the other produces wrong calculations.

### 4.5 Derive the distribution update from total probability

Let

$$
r_n=\begin{pmatrix}P(X_n=1)\\P(X_n=2)\end{pmatrix}
$$

be the column vector describing uncertainty about today's state. This vector is not itself the machine's state. A realized state is either $1$ or $2$; the vector records probabilities of those alternatives.

For tomorrow to be in state $1$, today's state must be either $1$ or $2$. The alternatives are mutually exclusive and exhaustive, so total probability gives

$$
\begin{aligned}
P(X_{n+1}=1)
&=P(X_{n+1}=1\mid X_n=1)P(X_n=1)\\
&\quad+P(X_{n+1}=1\mid X_n=2)P(X_n=2)\\
&=0.7P(X_n=1)+0.4P(X_n=2).
\end{aligned}
$$

The analogous calculation for state $2$ gives $0.3P(X_n=1)+0.6P(X_n=2)$. Together these are exactly the matrix-vector product

$$
\boxed{r_{n+1}=Ar_n.}
$$

For $r_0=(0.5,0.5)^\mathsf T$, one step gives $r_1=(0.55,0.45)^\mathsf T$. The entries remain nonnegative, and their sum is one. Algebraically, the column sums of $A$ ensure total probability mass is preserved.

Repeating the same update gives $r_2=A(Ar_0)=A^2r_0$, and by induction,

$$
\boxed{r_n=A^nr_0.}
$$

The homogeneity assumption is what allows the same $A$ to be used at every step. With matrices $A_0,A_1,\ldots$ changing over time, the update would instead be $r_n=A_{n-1}\cdots A_1A_0r_0$.

### 4.6 Why matrix multiplication sums over possible routes

Suppose the starting state is definitely $1$. To be in state $1$ two ticks later, the path must be either $1\to1\to1$ or $1\to2\to1$.

The first route has probability

$$
P(X_1=1,X_2=1\mid X_0=1)
=P(X_1=1\mid X_0=1)P(X_2=1\mid X_1=1,X_0=1)
=0.7\cdot0.7.
$$

The multiplication comes from the probability product rule. The Markov property removes $X_0$ from the second conditional, and homogeneity identifies the resulting factor as $p_{11}$. We did not assert that the consecutive states are independent.

The second route has probability $0.3\cdot0.4$. Add the mutually exclusive routes:

$$
P(X_2=1\mid X_0=1)=0.7^2+0.3\cdot0.4=0.61.
$$

For arbitrary departure $i$ and destination $j$,

$$
P(X_2=j\mid X_0=i)=\sum_k p_{ik}p_{kj}.
$$

The intermediate state $k$ is unobserved, so we sum it out. In our convention,

$$
(A^2)_{ji}=\sum_k A_{jk}A_{ki}
=\sum_k p_{kj}p_{ik},
$$

which is the same route sum. Explicit multiplication yields

$$
A^2=
\begin{pmatrix}
0.7^2+0.4(0.3)&0.7(0.4)+0.4(0.6)\\
0.3(0.7)+0.6(0.3)&0.3(0.4)+0.6^2
\end{pmatrix}
=
\begin{pmatrix}0.61&0.52\\0.39&0.48\end{pmatrix}.
$$

Thus column $i$ of $A^n$ is the distribution after $n$ ticks conditional on starting in $i$. It is not the probability of one specific $n$-step path. The likelihood of an observed path keeps its intermediate states, a distinction developed in Chapter 5.

### 4.7 A stationary distribution is an unchanged probability distribution

A distribution $\pi$ is stationary when feeding it through one transition produces the same distribution:

$$
\boxed{A\pi=\pi.}
$$

This does not mean that each realized machine stays in its current state. Individual paths continue moving. It means that the fractions of probability assigned to the possible states stay the same.

Write $\pi=(\pi_1,\pi_2)^\mathsf T$, with $\pi_1+\pi_2=1$. The first component equation is

$$
0.7\pi_1+0.4\pi_2=\pi_1.
$$

Move $0.7\pi_1$ to the right:

$$
0.4\pi_2=0.3\pi_1.
$$

Substitute $\pi_2=1-\pi_1$:

$$
0.4(1-\pi_1)=0.3\pi_1
\quad\Longrightarrow\quad
0.4=0.7\pi_1
\quad\Longrightarrow\quad
\pi_1=4/7.
$$

Therefore $\pi_2=3/7$. Substituting into both component equations confirms the result. The equality $0.3\pi_1=0.4\pi_2$ can also be read as balancing the flow from $1$ to $2$ against the flow from $2$ to $1$.

The equation $A\pi=1\cdot\pi$ makes $\pi$ an eigenvector with eigenvalue one. An eigenvector is a nonzero vector whose direction is unchanged by a matrix multiplication; the multiplier is its eigenvalue. Probability normalization selects the scaling for this particular eigenvector.

### 4.8 Derive convergence for this example, rather than confuse it with stationarity

Finding a fixed point does not yet show that every initial distribution approaches it. We can prove convergence directly for our two-state chain.

Let $a_n=P(X_n=1)$, so $P(X_n=2)=1-a_n$. The update becomes

$$
a_{n+1}=0.7a_n+0.4(1-a_n)=0.4+0.3a_n.
$$

The stationary value $a_*=4/7$ satisfies $a_*=0.4+0.3a_*$. Subtract this fixed-point equation from the update:

$$
a_{n+1}-a_*=0.3(a_n-a_*).
$$

Apply the relation repeatedly:

$$
\boxed{a_n=\frac47+0.3^n\left(a_0-\frac47\right).}
$$

Because $0.3^n\to0$, the initial probability $a_0$ eventually stops mattering. Starting definitely in state $1$ gives $a_0=1$ and $a_2=4/7+0.09(3/7)=0.61$, agreeing with the route calculation. Starting definitely in state $2$ gives $a_0=0$ and $a_2=4/7-0.09(4/7)=0.52$.

If $a_0=4/7$ initially, the error term is zero at every time: that is stationary initialization. If $a_0\ne4/7$, the distribution changes and only converges toward stationarity. A process started at a fixed state is not automatically stationary just because it has a stationary distribution.

### 4.9 General conditions, explained through failures

For a finite-state chain, **irreducible** means that every state can reach every other state along some sequence of transitions with positive probability. A chain with two disconnected groups fails this condition. Its initial group cannot be forgotten.

A state has a period determined by the common divisor of its possible return times. A chain is **aperiodic** when that period is one. The practical issue is whether a rigid cycle prevents convergence. For example,

$$
A=\begin{pmatrix}0&1\\1&0\end{pmatrix}
$$

has the stationary distribution $(1/2,1/2)^\mathsf T$, but a definite start alternates forever between $(1,0)^\mathsf T$ and $(0,1)^\mathsf T$. The chain is irreducible but periodic; the time-specific distribution has no limit.

For finite chains, irreducibility ensures a unique stationary distribution, and adding aperiodicity gives convergence to it from every initial distribution. This general result is a theorem. Our earlier recurrence supplies a complete proof for the running two-state example, not a proof for arbitrary matrices. Having every transition probability strictly positive is an easy sufficient condition for both properties: every state communicates directly, and self-transitions remove a rigid return cycle. It is sufficient, not necessary.

At the other extreme, let every column of $A$ already equal the same probability vector $\pi$. Then $Ar=\pi$ for every initial probability vector $r$. Convergence occurs after one step, exactly. It is therefore false that finite matrix powers can never equal their stationary limit.

### 4.10 Three words that should not be used interchangeably

**Markov** describes how much information is needed for conditional prediction. **Time-homogeneous** describes whether the transition rule changes with absolute time. **Stationary** describes whether time shifts change the process's joint distributions.

For a homogeneous chain initialized in an invariant distribution, the joint probability of a consecutive block is

$$
P(X_n=i_0,\ldots,X_{n+r}=i_r)
=\pi_{i_0}p_{i_0i_1}\cdots p_{i_{r-1}i_r}.
$$

There is no dependence on $n$, so shifting the block in time does not change its distribution. For nonconsecutive times, the same argument uses the appropriate powers of $A$ between observations. This explains why homogeneous dynamics plus stationary initialization produce a stationary process, rather than treating the phrase as an unexplained label.

::: {.worked}
**Worked check.** A chain has $p_{11}=0.9$, $p_{12}=0.1$, $p_{21}=0.2$, and $p_{22}=0.8$. In column convention, $A=\begin{pmatrix}0.9&0.2\\0.1&0.8\end{pmatrix}$. Stationarity requires $0.1\pi_1=0.2\pi_2$ and total mass one, so $\pi=(2/3,1/3)^\mathsf T$. Writing $a_n=P(X_n=1)$ gives $a_{n+1}=0.2+0.7a_n$, hence $a_n=2/3+0.7^n(a_0-2/3)$. This proves convergence for this example and shows explicitly how its starting distribution is forgotten.
:::

<p class="source-note">Source connection: <a href="#ref-markov">Berkeley Prob 140</a> states the general finite-chain convergence theorem.  the column convention, route sums, invariant distributions, and positive-entry convergence theorem are in <a href="#ref-l5">MIT Lecture 5, section 3</a>. The explicit recurrences here prove the stated numerical examples without requiring spectral theory.</p>


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


## 6. Give the process a home {#ar}

### 6.1 What a random walk cannot express

A random walk evolves as $W_{i+1}=W_i+S_{i+1}$, where $S_{i+1}$ is a fresh independent step. If the step mean is zero, knowing $W_i=w$ gives expected next position $w$. The process has no preferred level. An unusually high current position does not, by itself, create a downward expected correction.

Now imagine a quantity that tends to return toward a normal level, while still receiving unpredictable disturbances. To model that behavior, specify a target level $\mu$ and retain only a fraction $\phi$ of the current deviation from it. Then add new noise:

$$
\boxed{X_i=\mu+\phi(X_{i-1}-\mu)+\varepsilon_i.}
$$

For this chapter, initially treat the parameters as known. We want to understand the model's behavior. Estimating the parameters from observations is a separate problem, treated in Chapter 7.

The quantity $X_i$ is the state at tick $i$. The parameter $\mu$ is the level toward which the conditional mean tends. The coefficient $\phi$ is the fraction of a deviation retained for one tick. The random innovation $\varepsilon_i$ is the fresh unpredictable contribution at that tick. We assume

$$
\varepsilon_i\sim N(0,\tau^2),\qquad \tau^2>0,
$$

with all innovations mutually independent and independent of the initial value $X_0$. The notation $N(m,v)$ means Gaussian with mean $m$ and variance $v$. Thus $\tau$ is the innovation standard deviation, and $\tau^2$ is its variance.

“Autoregressive of order one,” abbreviated AR(1), means that the next value is expressed using one lag of the same process. Rewriting the equation gives the familiar intercept-and-slope form

$$
X_i=\underbrace{\mu(1-\phi)}_{\text{intercept}}+\phi X_{i-1}+\varepsilon_i.
$$

### 6.2 Conditional prediction: what stays fixed and what remains random

Suppose the current value is observed as $X_{i-1}=x$. Under this condition, the quantity $\mu+\phi(x-\mu)$ is a fixed number. The only remaining uncertainty in the next value is the innovation. Adding a constant to a Gaussian shifts its mean and leaves its variance unchanged. Therefore

$$
\boxed{X_i\mid X_{i-1}=x\sim N\bigl(\mu+\phi(x-\mu),\tau^2\bigr).}
$$

Conditional expectation averages over the remaining randomness while holding the supplied current value fixed. Since $E[\varepsilon_i]=0$,

$$
E[X_i\mid X_{i-1}=x]=\mu+\phi(x-\mu).
$$

Conditional variance measures the uncertainty around that prediction. Because the predictable part is constant under the condition,

$$
\operatorname{Var}(X_i\mid X_{i-1}=x)=\tau^2.
$$

Take $\mu=25$, $\phi=1/2$, and $\tau^2=4$. If the current observation is $29$, then the deviation is $29-25=4$. Half remains in the forecast, giving predicted next value $25+2=27$. The distribution is $N(27,4)$, with standard deviation $2$.

If the realized innovation is $-1$, the next value is $26$. If it is $3$, the next value is $30$. Both are compatible with the same mean-reverting model. In fact,

$$
P(X_i>29\mid X_{i-1}=29)
=P\left(Z>\frac{29-27}{2}\right)=P(Z>1),
$$

where $Z\sim N(0,1)$. This probability is about $0.159$. Mean reversion does not prohibit an upward move from an already high value. It describes the **average conditional move**, not every realized move.

### 6.3 Derive the expected pull toward the target

Subtract the supplied current value $x$ from the conditional mean:

$$
\begin{aligned}
E[X_i-X_{i-1}\mid X_{i-1}=x]
&=\mu+\phi(x-\mu)-x\\
&=\mu-\phi\mu+\phi x-x\\
&=(1-\phi)\mu-(1-\phi)x\\
&=-(1-\phi)(x-\mu).
\end{aligned}
$$

When $0<\phi<1$, the factor $1-\phi$ is positive. If $x>\mu$, the expected change is negative; if $x<\mu$, it is positive. If $x=\mu$, there is no expected drift, although noise still produces movement.

For $x\ne\mu$ and $0<\phi<1$, the conditional mean lies strictly between $x$ and $\mu$. A smaller positive $\phi$ means less of the deviation survives one tick, hence faster mean reversion on that tick scale.

The parameter range matters. At $\phi=0$, the next value is $\mu+\varepsilon_i$ and forgets the current value completely. For $-1<\phi<0$, the conditional deviation changes sign at each tick while its magnitude shrinks: the model oscillates around $\mu$. At $\phi=1$, the level parameter cancels and the model becomes $X_i=X_{i-1}+\varepsilon_i$, a Gaussian random walk. At $|\phi|>1$, the deterministic effect of a deviation grows rather than decays.

The stable causal AR(1) model with fresh independent innovations uses $|\phi|<1$. The mean-reverting scalar OU model sampled at a positive fixed time gap corresponds specifically to $0<\phi<1$; it does not produce a negative coefficient of this kind.

### 6.4 Why the model is Markov

Given the entire history through tick $i-1$, the new innovation is still independent Gaussian noise. The recursion only uses that history through $X_{i-1}$. Therefore the conditional distribution of $X_i$ given $X_0,\ldots,X_{i-1}$ is the same Gaussian as the one given $X_{i-1}$ alone.

This proves the first-order Markov property for the specified model. The argument depends on the fresh-noise assumption. Merely writing a formula with a symbol $\varepsilon_i$ would not prove Markovness if the noise retained extra dependence on older states.

The innovations are independent; the states are not. An innovation at an early time is carried forward through the recursion and affects several later states. To see exactly how, we now solve the recursion.

### 6.5 Unroll the equation instead of skipping to a formula

Center the process at the target by defining $Y_i=X_i-\mu$. Then

$$
Y_i=\phi Y_{i-1}+\varepsilon_i.
$$

Start from a fixed $X_0=x_0$, hence $Y_0=x_0-\mu$. The first three expansions are

$$
Y_1=\phi Y_0+\varepsilon_1,
$$

$$
\begin{aligned}
Y_2&=\phi Y_1+\varepsilon_2\\
&=\phi(\phi Y_0+\varepsilon_1)+\varepsilon_2\\
&=\phi^2Y_0+\phi\varepsilon_1+\varepsilon_2,
\end{aligned}
$$

and

$$
Y_3=\phi^3Y_0+\phi^2\varepsilon_1+\phi\varepsilon_2+\varepsilon_3.
$$

Every existing contribution is multiplied by $\phi$ each time it advances one more tick. The newest innovation has multiplier one; an innovation from $r$ ticks ago has multiplier $\phi^r$.

The resulting general expression is

$$
\boxed{X_k=\mu+\phi^k(x_0-\mu)+\sum_{j=1}^k\phi^{k-j}\varepsilon_j.}
$$

It can be checked by induction: multiply the centered expression at tick $k$ by $\phi$, append $\varepsilon_{k+1}$, and the formula at tick $k+1$ follows. Thus the expansion is a solution of the recursion, not just a visual pattern.

### 6.6 Derive the finite-time mean

The starting deviation is a fixed number, while each innovation has mean zero. Linearity of expectation gives

$$
\begin{aligned}
E[X_k\mid X_0=x_0]
&=\mu+\phi^k(x_0-\mu)+\sum_{j=1}^k\phi^{k-j}E[\varepsilon_j]\\
&=\mu+\phi^k(x_0-\mu).
\end{aligned}
$$

For $|\phi|<1$, the factor $\phi^k$ tends to zero. Hence the mean approaches $\mu$ as the number of ticks grows. That is why $\mu$ is called the long-run mean or mean-reversion level.

It is not necessarily the mean at every finite time. Starting at $29$ with $\mu=25$ and $\phi=1/2$, the successive means are $27$, $26$, $25.5$, and so on. Their destination is $25$, but their finite-time values are not all $25$.

### 6.7 Derive the finite-time variance and the geometric sum

For a fixed start, only the innovation sum contributes variance. If a variable is multiplied by $a$, its variance is multiplied by $a^2$, since

$$
E[(aU-E[aU])^2]=a^2E[(U-E[U])^2].
$$

Distinct innovations are independent, so the covariance terms vanish. Therefore

$$
\begin{aligned}
\operatorname{Var}(X_k\mid X_0=x_0)
&=\sum_{j=1}^k\operatorname{Var}(\phi^{k-j}\varepsilon_j)\\
&=\tau^2\sum_{j=1}^k\phi^{2(k-j)}\\
&=\tau^2(1+\phi^2+\phi^4+\cdots+\phi^{2(k-1)}).
\end{aligned}
$$

Let $G_k=1+r+\cdots+r^{k-1}$. Then $rG_k=r+r^2+\cdots+r^k$. Subtracting gives $(1-r)G_k=1-r^k$, so $G_k=(1-r^k)/(1-r)$ for $r\ne1$. With $r=\phi^2$,

$$
\boxed{\operatorname{Var}(X_k\mid X_0=x_0)
=\tau^2\frac{1-\phi^{2k}}{1-\phi^2}},\qquad |\phi|<1.
$$

The square in $1-\phi^2$ is essential. It appears because **variance squares the innovation weights**. Using $1-\phi$ would be summing the wrong powers.

For the numerical model $(\mu,\phi,\tau^2)=(25,1/2,4)$, the first variance is $4$. The second is $4(1+1/4)=5$. The third is $4(1+1/4+1/16)=5.25$. The variance is increasing from a known initial value, but it is approaching a finite limit rather than increasing without bound as in a random walk.

### 6.8 Why the finite-time law is exactly Gaussian

The expansion is a constant plus a finite weighted sum of independent Gaussian variables. Such a sum is Gaussian. This is an exact closure property, not a central-limit approximation.

One way to justify it is through exponential moments. If $U\sim N(m,v)$, completing the square in its density integral gives

$$
E[e^{rU}]=e^{rm+r^2v/2}.
$$

To see where that Gaussian exponential moment comes from, insert the Gaussian density into the definition of expectation. For $v>0$,

$$
E[e^{rU}]=\int_{-\infty}^{\infty}
 e^{ru}\frac{1}{\sqrt{2\pi v}}e^{-(u-m)^2/(2v)}du.
$$

Complete the square in the exponent:

$$
ru-\frac{(u-m)^2}{2v}
=-\frac{(u-m-rv)^2}{2v}+rm+\frac{r^2v}{2}.
$$

The terms $rm+r^2v/2$ do not involve $u$, so their exponential can be taken outside the integral. The remaining integrand is the density of $N(m+rv,v)$, whose integral is 1. This proves the displayed expression for $E[e^{rU}]$. It is a way of describing a distribution through a function of the auxiliary number $r$, not an extra model parameter being estimated.

For independent $U_j$, the exponential moment of $\sum_j a_jU_j$ factors into $\prod_jE[e^{ra_jU_j}]$. Multiplying the expressions above gives

$$
\exp\left[r\sum_j a_jm_j+\frac{r^2}{2}\sum_j a_j^2v_j\right],
$$

which is the exponential-moment function of a Gaussian with those mean and variance parameters. The standard uniqueness result for such moment-generating functions identifies the distribution. The essential mechanism is independence plus the quadratic form of the Gaussian exponential moment.

Thus, for a fixed initial value,

$$
X_k\mid X_0=x_0
\sim N\left(\mu+\phi^k(x_0-\mu),\ \tau^2\frac{1-\phi^{2k}}{1-\phi^2}\right).
$$

If the initial value is random and non-Gaussian, an unconditional finite-time state need not be Gaussian. Its distribution contains the transformed initial randomness as well as Gaussian innovations. The conditioning or initialization assumption must accompany the Gaussian claim.

### 6.9 Stationary variance as a balance of old and new uncertainty

For $|\phi|<1$, the finite-time variance tends to

$$
q=\frac{\tau^2}{1-\phi^2}.
$$

Here $q$ denotes a **stationary state variance**, not an innovation variance. The distinction is physical: $\tau^2$ is new uncertainty injected on one tick; $q$ includes the surviving contributions of all past innovations.

We can derive the same formula by asking which variance would remain unchanged under a transition. Suppose the current state has variance $q$ and is independent of the next innovation. The next state has variance $\phi^2q+\tau^2$. For it to remain $q$, we require

$$
q=\phi^2q+\tau^2
\quad\Longrightarrow\quad
(1-\phi^2)q=\tau^2
\quad\Longrightarrow\quad
q=\frac{\tau^2}{1-\phi^2}.
$$

If the current state is Gaussian with mean $\mu$ and variance $q$, the next is Gaussian too. Its mean remains $\mu$ and its variance remains $q$, so the entire one-time distribution is preserved:

$$
\boxed{\pi=N\left(\mu,\frac{\tau^2}{1-\phi^2}\right).}
$$

Initialize $X_0$ from this distribution independently of all future innovations. Because the transition rule is the same at every tick, every finite block of states then has a distribution unchanged by a common shift in time. This is a stationary process.

Starting from a fixed $x_0$ produces convergence toward this distribution, not stationarity from the first tick. Also, convergence of distributions does not mean individual paths eventually stop at $\mu$. Fresh noise continues to enter forever.

### 6.10 Derive autocovariance and autocorrelation

In stationarity, each state has variance $q$. To compare states $k$ ticks apart, restart the expansion at tick $i$:

$$
X_{i+k}-\mu=\phi^k(X_i-\mu)
+\sum_{j=1}^k\phi^{k-j}\varepsilon_{i+j}.
$$

The future innovations are independent of $X_i$. Covary both sides with $X_i$:

$$
\begin{aligned}
\operatorname{Cov}(X_i,X_{i+k})
&=\phi^k\operatorname{Var}(X_i)
+\sum_{j=1}^k\phi^{k-j}\operatorname{Cov}(X_i,\varepsilon_{i+j})\\
&=\phi^kq+0.
\end{aligned}
$$

Autocovariance means covariance of the process with itself at a different time. The stationary lag-$k$ autocovariance is therefore $q\phi^k$ for $k\ge0$. Correlation divides by the two standard deviations, both $\sqrt q$:

$$
\boxed{\rho(k)=\frac{q\phi^k}{\sqrt q\sqrt q}=\phi^k.}
$$

The coefficient governing the fraction of a deviation retained each tick also governs how correlation decays across ticks. For $\phi=1/2$, lag-one correlation is $1/2$, lag-two correlation is $1/4$, and lag-three correlation is $1/8$.

With a fixed start instead, the same covariance argument yields $\operatorname{Cov}(X_i,X_{i+k})=\phi^k\operatorname{Var}(X_i)$, but the variances at $i$ and $i+k$ are not generally equal. One cannot automatically divide by $q$ and declare the correlation to be $\phi^k$ before stationarity is established.

### 6.11 What to distinguish when reading an AR(1) formula

For a known current observation, the next-step mean is $\mu+\phi(x-\mu)$ and the next-step variance is $\tau^2$. For a fixed start several ticks earlier, accumulated innovation noise gives the larger finite-time variance derived above. Under stationary initialization, the unconditional variance is $q=\tau^2/(1-\phi^2)$.

These are three different conditioning setups. Calling all three quantities “the variance” without saying what is conditioned on makes correct formulas appear contradictory.

::: {.worked}
**Worked check.** In the model $\mu=25$, $\phi=1/2$, $\tau^2=4$, with fixed $X_0=29$, calculate the law of $X_2$. Expanding gives $X_2=25+(1/2)^2(4)+(1/2)\varepsilon_1+\varepsilon_2$. Its mean is $26$, and its variance is $(1/2)^2(4)+4=5$. Hence $X_2\sim N(26,5)$. The stationary law would instead be $N(25,16/3)$. Knowing $X_1=27$ makes the next-step conditional law $X_2\mid X_1=27\sim N(26,4)$. The same mean $26$ in two of these calculations does not make their variances equal: conditioning on $X_1$ removes uncertainty about the first innovation's contribution.
:::

<p class="source-note">Source connection: the AR(1) model and stationarity framework are in <a href="#ref-l8">MIT Lecture 8</a>. The formulas here are derived from the recursion, including the required squared coefficient in the stationary-variance denominator.</p>


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


## 12. Solve the Ornstein–Uhlenbeck model {#ou}

### 12.1 State the model and separate the two possible tasks

The Ornstein–Uhlenbeck process, or OU process, is defined here by

$$
dX_t=\theta(\mu-X_t)\,dt+\sigma\,dB_t,
\qquad \theta>0,\quad\sigma>0,\quad\mu\in\mathbb R.
$$

The parameter $\mu$ is the level toward which the process tends to return. The positive parameter $\theta$ determines the speed of this return. The positive parameter $\sigma$ scales the Brownian disturbances. Standard Brownian motion has independent increments $B_{t+h}-B_t\sim N(0,h)$, so the noise accumulated over a short interval has standard deviation proportional to $\sqrt h$.

There are two different tasks involving these parameters. **Solving the model** means finding a process $X_t$ for specified parameter values, an initial value and Brownian noise. **Estimating the model** means inferring parameter values from observed process values. This chapter solves the model first. The resulting transition distribution will make estimation possible in Chapter 13.

The differential equation is shorthand for the integral equation

$$
\boxed{X_t=X_0+\theta\int_0^t(\mu-X_s)\,ds+\sigma B_t.}
$$

Initially take $X_0=x_0$ fixed. We will later allow a random initial value independent of future Brownian increments. The ordinary time integral on the right is well-defined along a continuous path. The challenge is that it contains the unknown $X_s$ itself. We want to remove that self-reference.

### 12.2 First solve the same return mechanism without noise

Ignore the noise temporarily and consider the ordinary equation

$$
x'(t)=\theta(\mu-x(t)),\qquad x(0)=x_0.
$$

Shift the variable by its target: let $y(t)=x(t)-\mu$. Because $\mu$ is a constant, $y'(t)=x'(t)$, and the equation becomes

$$
y'(t)=-\theta y(t).
$$

To solve it, multiply by $e^{\theta t}$. Why that factor? The ordinary product rule gives

$$
\begin{aligned}
\frac{d}{dt}\big(e^{\theta t}y(t)\big)
&=\theta e^{\theta t}y(t)+e^{\theta t}y'(t)\\
&=\theta e^{\theta t}y(t)-\theta e^{\theta t}y(t)\\
&=0.
\end{aligned}
$$

The exponential was chosen so that its derivative generates exactly the term needed to cancel the original decay. A function with derivative zero is constant, so $e^{\theta t}y(t)=y(0)=x_0-\mu$. Hence

$$
x(t)=\mu+e^{-\theta t}(x_0-\mu).
$$

The deviation from the target is multiplied by $e^{-\theta t}$. It becomes smaller as time passes because $\theta>0$. We expect this deterministic decay to remain in the stochastic solution, but the noise creates a new issue: $X_t$ is not ordinarily differentiable. We must not simply pretend that the same derivative calculation applies to the entire noisy process without justification.

### 12.3 Remove the rough noise before using ordinary differentiation

Center the stochastic process in the same way:

$$
Y_t=X_t-\mu.
$$

The integral equation becomes

$$
Y_t=Y_0-\theta\int_0^tY_s\,ds+\sigma B_t.
$$

Now introduce a second process

$$
V_t=Y_t-\sigma B_t.
$$

This subtraction is useful because it removes the explicitly rough Brownian contribution. Substituting $Y_s=V_s+\sigma B_s$ into the integral equation gives

$$
V_t=Y_0-\theta\int_0^t\big(V_s+\sigma B_s\big)\,ds.
$$

For any fixed continuous Brownian path, the right side is an ordinary integral of a continuous function. Therefore $V_t$ is ordinarily differentiable along that path, even though $Y_t$ need not be. Differentiating this ordinary integral equation is legitimate and yields

$$
V'(t)=-\theta V_t-\theta\sigma B_t,
\qquad V_0=Y_0.
$$

Notice what has **not** happened: we have not differentiated $B_t$. It is merely a continuous forcing term on the right of an ordinary differential equation. An equation can be driven by a nondifferentiable continuous input and still have a differentiable accumulated response.

This converts the stochastic problem into an ordinary linear equation with a random but fixed-in-the-calculation input path. We can now apply the integrating-factor method without treating Brownian motion as differentiable.

### 12.4 Solve the transformed ordinary equation step by step

Multiply the equation for $V$ by $e^{\theta t}$ and apply the ordinary product rule:

$$
\begin{aligned}
\frac{d}{dt}\big(e^{\theta t}V_t\big)
&=e^{\theta t}\big(V'(t)+\theta V_t\big)\\
&=-\theta\sigma e^{\theta t}B_t.
\end{aligned}
$$

Integrate from 0 to $t$:

$$
e^{\theta t}V_t-V_0
=-\theta\sigma\int_0^t e^{\theta s}B_s\,ds.
$$

Since $V_0=Y_0$, dividing by $e^{\theta t}$ gives

$$
V_t=e^{-\theta t}Y_0
-\theta\sigma\int_0^t e^{-\theta(t-s)}B_s\,ds.
$$

Finally recover $Y_t=V_t+\sigma B_t$ and then $X_t=\mu+Y_t$:

$$
\boxed{
X_t=\mu+e^{-\theta t}(X_0-\mu)
+\sigma\left[B_t-\theta\int_0^t e^{-\theta(t-s)}B_s\,ds\right].}
$$

Everything in this expression is already defined using Brownian values and an ordinary time integral. It is an explicit solution: the unknown process no longer appears inside its own integral. This formula is also a direct construction of a continuous process, because $B$ and the ordinary integral term are continuous.

### 12.5 Express the solution as a weighted stochastic integral

The usual OU solution is written in a more interpretable weighted-noise form. For a deterministic continuously differentiable $f$, integration by parts against Brownian motion gives

$$
\int_0^tf(s)\,dB_s=f(t)B_t-f(0)B_0-\int_0^tf'(s)B_s\,ds.
$$

This identity does not require $B'$: it follows from a finite telescoping product identity and continuity, as derived in Chapter 11. In the present calculation take $f(s)=e^{-\theta(t-s)}$, treating the terminal time $t$ as fixed. Then $f(t)=1$, $B_0=0$, and $f'(s)=\theta e^{-\theta(t-s)}$. Hence

$$
B_t-\theta\int_0^t e^{-\theta(t-s)}B_s\,ds
=\int_0^t e^{-\theta(t-s)}\,dB_s.
$$

Substitute this into the explicit pathwise solution:

$$
\boxed{
X_t=\mu+e^{-\theta t}(X_0-\mu)
+\sigma\int_0^t e^{-\theta(t-s)}\,dB_s.}
$$

The first term is the target level. The second is the surviving part of the initial deviation. The third is all the later noise, with each disturbance damped according to its age. A disturbance arriving near $s=t$ has weight close to 1. A disturbance arriving long before $t$ has a smaller weight $e^{-\theta(t-s)}$.

This is the continuous-time version of an AR(1) expansion in which old innovations carry powers of an autoregressive coefficient. The exponential weight is not guessed solely from analogy: we derived it from the equation.

### 12.6 Verify that the construction is a solution and that it is unique

The ordinary equation for $V$ supplies a short verification. Define $V$ by the explicit expression in §12.4. Ordinary differentiation shows it satisfies

$$
V'(t)=-\theta(V_t+\sigma B_t),\qquad V_0=X_0-\mu.
$$

Integrating this equation and defining $X_t=\mu+V_t+\sigma B_t$ gives

$$
\begin{aligned}
X_t
&=\mu+V_0-\theta\int_0^t(V_s+\sigma B_s)ds+\sigma B_t\\
&=X_0-\theta\int_0^t(X_s-\mu)ds+\sigma B_t\\
&=X_0+\theta\int_0^t(\mu-X_s)ds+\sigma B_t.
\end{aligned}
$$

That is exactly the defining integral equation. No approximation has been used.

For uniqueness, suppose $X$ and $\widetilde X$ are continuous solutions with the **same** initial value and the **same** Brownian path. Their difference $D_t=X_t-\widetilde X_t$ satisfies

$$
D_t=-\theta\int_0^tD_s\,ds.
$$

The Brownian terms cancel. This difference is therefore differentiable, with $D'(t)=-\theta D_t$ and $D_0=0$. The ordinary solution from §12.2 is $D_t=e^{-\theta t}D_0=0$. Thus the solutions agree at all times along the continuous path.

Different Brownian paths produce different OU paths; uniqueness does not say every realization looks the same. It says that once the initial value and driving noise are fixed, the model does not leave multiple possible responses.

### 12.7 Derive the distribution from a fixed starting value

Let $X_0=x_0$ be fixed. A deterministic weighted Brownian integral is a centered Gaussian, with variance equal to the ordinary integral of the squared weight. Consequently, the random term

$$
Z_t=\sigma\int_0^t e^{-\theta(t-s)}\,dB_s
$$

has mean zero. Its variance is

$$
\begin{aligned}
\operatorname{Var}(Z_t)
&=\sigma^2\int_0^t e^{-2\theta(t-s)}\,ds\\
&=\sigma^2\int_0^t e^{-2\theta r}\,dr
\qquad(r=t-s)\\
&=\sigma^2\left[-\frac{e^{-2\theta r}}{2\theta}\right]_0^t\\
&=\frac{\sigma^2}{2\theta}(1-e^{-2\theta t}).
\end{aligned}
$$

The factor $2\theta$ appears because the exponential weight was **squared** in the variance calculation. Taking the expectation of the solution removes only its centered random term:

$$
E[X_t\mid X_0=x_0]=\mu+e^{-\theta t}(x_0-\mu).
$$

Adding a deterministic number to a Gaussian shifts its mean without changing its variance. Therefore

$$
\boxed{
X_t\mid X_0=x_0
\sim N\!\left(
\mu+e^{-\theta t}(x_0-\mu),
\frac{\sigma^2}{2\theta}(1-e^{-2\theta t})
\right).}
$$

At $t=0$, the conditional variance is zero because the starting value is known. For any $t>0$, the variance is positive. Mean reversion controls the accumulation of noise; it does not eliminate noise.

### 12.8 Derive the transition from an arbitrary observed time

To use data, we need the distribution of $X_t$ given an earlier value $X_s$, not only a distribution from time zero. Split the noise integral in the solution at time $s<t$:

$$
\int_0^t e^{-\theta(t-u)}dB_u
=\int_0^s e^{-\theta(t-u)}dB_u
 +\int_s^t e^{-\theta(t-u)}dB_u.
$$

For $u\leq s$, the exponential factors as

$$
e^{-\theta(t-u)}=e^{-\theta(t-s)}e^{-\theta(s-u)}.
$$

Therefore the portion containing the initial value and all noise before $s$ is

$$
\begin{aligned}
&e^{-\theta t}(X_0-\mu)
 +\sigma\int_0^s e^{-\theta(t-u)}dB_u\\
&\qquad=e^{-\theta(t-s)}
\left[e^{-\theta s}(X_0-\mu)
 +\sigma\int_0^s e^{-\theta(s-u)}dB_u\right]\\
&\qquad=e^{-\theta(t-s)}(X_s-\mu).
\end{aligned}
$$

Substitution yields the exact restart relation

$$
\boxed{X_t=\mu+e^{-\theta(t-s)}(X_s-\mu)
+\sigma\int_s^t e^{-\theta(t-u)}dB_u.}
$$

The final integral depends only on Brownian increments after $s$. Its deterministic weight makes it Gaussian, and independent increments make it independent of the information available by $s$, including an independent initial value. Once $X_s$ is supplied, the earlier path adds no further information to this conditional distribution. This is the **Markov property**, now justified by the solution rather than assumed without explanation.

Write $\Delta=t-s>0$ and suppose $X_s=x$. The conditional distribution is

$$
\boxed{
X_{s+\Delta}\mid X_s=x
\sim N\big(m_\Delta(x),v_\Delta\big),}
$$

where

$$
m_\Delta(x)=\mu+e^{-\theta\Delta}(x-\mu),
\qquad
v_\Delta=\frac{\sigma^2}{2\theta}(1-e^{-2\theta\Delta}).
$$

The transition depends on the elapsed duration $\Delta$, not the calendar time $s$. That is time-homogeneity in continuous time. It does not require every pair of recorded observations to have the same gap.

### 12.9 Identify the stationary distribution and explain initialization

Define

$$
q=\frac{\sigma^2}{2\theta}.
$$

This is a convenient name for a variance, not an extra independent parameter. For a fixed starting value, as $t\to\infty$, the mean tends to $\mu$ and the variance tends to $q$. The Gaussian distributions therefore approach $N(\mu,q)$.

To check that this is an invariant distribution, suppose the current value is already distributed as $N(\mu,q)$ and is independent of future noise. Set $\phi=e^{-\theta\Delta}$. The transition can be expressed as

$$
X_{s+\Delta}=\mu+\phi(X_s-\mu)+\varepsilon,
\qquad \varepsilon\sim N(0,q(1-\phi^2)),
$$

with $\varepsilon$ independent of $X_s$. The next value is Gaussian, its mean is $\mu$, and its variance is

$$
\phi^2q+q(1-\phi^2)=q.
$$

Thus the transition preserves $N(\mu,q)$ for every positive duration. Starting with

$$
X_0\sim N(\mu,q)
$$

independently of future Brownian increments gives a stationary OU process. Its joint laws are invariant under time shifts because the initial marginal is invariant and all subsequent conditional laws depend only on elapsed gaps.

Starting at a fixed $x_0$ gives a different situation. Its variance begins at zero and increases toward $q$, so that process is not stationary from the start. Also, convergence of distributions to a stationary law does not say an individual path settles permanently at $\mu$. In stationarity there is still positive variance and continuing random movement.

More generally, for an independent random initial value with finite mean $m_0$ and variance $v_0$, the solution gives

$$
E[X_t]=\mu+e^{-\theta t}(m_0-\mu),
$$

and, because the initial value and future noise are independent,

$$
\operatorname{Var}(X_t)
=e^{-2\theta t}v_0+q(1-e^{-2\theta t}).
$$

Setting $m_0=\mu$ and $v_0=q$ makes these moments time-invariant. To conclude Gaussian stationarity from the initialization argument above, take the **Gaussian** initial distribution, not merely an arbitrary initial distribution with those two moments.

### 12.10 Derive covariance and correlation, including the fixed-start case

For $s\leq t$, use the restart relation and take covariance with $X_s$:

$$
\begin{aligned}
\operatorname{Cov}(X_s,X_t)
&=\operatorname{Cov}\big(X_s,\mu+e^{-\theta(t-s)}(X_s-\mu)
     +\text{future noise}\big)\\
&=e^{-\theta(t-s)}\operatorname{Var}(X_s).
\end{aligned}
$$

Constants contribute zero covariance, and the future noise is independent of $X_s$. In stationarity, $\operatorname{Var}(X_s)=q$, so

$$
\boxed{\operatorname{Cov}(X_s,X_t)=q e^{-\theta|t-s|}.}
$$

Both standard deviations are $\sqrt q$. Dividing covariance by their product gives

$$
\boxed{\operatorname{Corr}(X_s,X_t)=e^{-\theta|t-s|}.}
$$

The same exponential controls how much a current deviation survives in a future forecast and how strongly two stationary values are correlated.

For a fixed initial value, use $\operatorname{Var}(X_s)=q(1-e^{-2\theta s})$ instead. For $s\leq t$,

$$
\begin{aligned}
\operatorname{Cov}(X_s,X_t)
&=q e^{-\theta(t-s)}(1-e^{-2\theta s})\\
&=q\big(e^{-\theta(t-s)}-e^{-\theta(t+s)}\big).
\end{aligned}
$$

The symmetric form is $q(e^{-\theta|t-s|}-e^{-\theta(t+s)})$. It depends on the two times separately, not just their difference. This is a direct mathematical indication that the fixed-start process is not stationary. At time zero the value has zero variance, so its correlation with a later value is undefined even though its covariance is zero.

### 12.11 Check short gaps, long gaps, half-life and units

For small $\Delta$, use the first terms of the exponential expansion:

$$
e^{-\theta\Delta}=1-\theta\Delta+O(\Delta^2).
$$

The notation $O(\Delta^2)$ means that for fixed parameters, the magnitude of the omitted remainder is bounded by a constant times $\Delta^2$ for sufficiently small $\Delta$. It identifies the order of the approximation error. Substituting into the conditional mean gives

$$
m_\Delta(x)=x+\theta(\mu-x)\Delta+O(\Delta^2).
$$

Similarly,

$$
\begin{aligned}
v_\Delta
&=\frac{\sigma^2}{2\theta}
   \left(2\theta\Delta-2\theta^2\Delta^2+O(\Delta^3)\right)\\
&=\sigma^2\Delta-\sigma^2\theta\Delta^2+O(\Delta^3).
\end{aligned}
$$

The Euler approximation keeps only the leading mean change and leading variance: mean $x+\theta(\mu-x)\Delta$ and variance $\sigma^2\Delta$. These are approximations to the exact transition, not its definition.

As $\Delta\to\infty$, the conditional mean tends to $\mu$ and the conditional variance tends to $q$. The process eventually loses predictive memory of the supplied current value, but not its randomness.

The **half-life** of a conditional deviation is the duration $h$ solving $e^{-\theta h}=1/2$. Taking logs gives $-\theta h=-\log 2$, hence

$$
h=\frac{\log 2}{\theta}.
$$

This measures halving of the forecast deviation, not a guaranteed time for a random sample path to reach a particular level.

If time is measured in days and $X$ in temperature units, $\theta$ has units of $1/\text{day}$ and $\sigma$ has units of temperature divided by $\sqrt{\text{day}}$. Then $\theta\Delta$ is dimensionless, as required inside an exponential, and $\sigma^2/(2\theta)$ has units of squared temperature, as required for a variance.

Another useful limit holds $\mu,\sigma$ fixed and lets $\theta\downarrow0$. Then $m_\Delta(x)\to x$ and $v_\Delta\to\sigma^2\Delta$: the model approaches scaled Brownian motion with no mean-reverting drift. The stationary variance diverges; at exactly $\theta=0$ there is no finite invariant Gaussian variance. The limiting behavior is consistent with losing the restoring force.

### 12.12 Recover AR(1) at observation times

For a fixed observation gap $\Delta$, write $X_i=X_{i\Delta}$. The restart relation becomes

$$
X_i=\mu+\phi(X_{i-1}-\mu)+\varepsilon_i,
$$

with

$$
\phi=e^{-\theta\Delta},
\qquad
\varepsilon_i\sim N(0,\tau^2),
\qquad
\tau^2=q(1-\phi^2).
$$

The innovations are independent because they are deterministic-weight Brownian integrals over disjoint observation intervals. Thus a fixed-gap sampled OU process is an exact Gaussian AR(1), not just an approximation to one.

Conversely, given a fixed-gap AR(1) with $0<\phi<1$ and $\tau^2>0$, take logs to recover

$$
\theta=-\frac{\log\phi}{\Delta}.
$$

Rearranging $\tau^2=[\sigma^2/(2\theta)](1-\phi^2)$ then gives

$$
\sigma^2=\frac{2\theta\tau^2}{1-\phi^2}.
$$

The restriction $0<\phi<1$ matters. A stable AR(1) may have a negative coefficient, but no positive real $\theta$ and positive gap produce a negative $e^{-\theta\Delta}$. Such an AR fit is not a mean-reverting scalar OU model of the form used here.

With irregular gaps $\Delta_i$, the same derivation gives $\phi_i=e^{-\theta\Delta_i}$ and $\tau_i^2=q(1-\phi_i^2)$. The coefficients now change with the recorded gaps. This does not mean the underlying continuous-time OU mechanism changed; it means we observed that mechanism over different durations.

### 12.13 Two shorter transitions must agree with one longer transition

This consistency check is useful both conceptually and in code. Define $\phi(h)=e^{-\theta h}$ and $v(h)=q(1-\phi(h)^2)$. Starting from $x$, transition for duration $h$ and then duration $k$:

$$
X_{s+h}=\mu+\phi(h)(x-\mu)+\varepsilon_1,
$$

$$
\begin{aligned}
X_{s+h+k}
&=\mu+\phi(k)(X_{s+h}-\mu)+\varepsilon_2\\
&=\mu+\phi(k)\phi(h)(x-\mu)+\phi(k)\varepsilon_1+\varepsilon_2.
\end{aligned}
$$

The two innovations are independent, with variances $v(h)$ and $v(k)$. The retained deviation is $\phi(k)\phi(h)=e^{-\theta(h+k)}=\phi(h+k)$. The combined variance is

$$
\begin{aligned}
\phi(k)^2v(h)+v(k)
&=q\left[\phi(k)^2(1-\phi(h)^2)+1-\phi(k)^2\right]\\
&=q\left[1-\phi(k)^2\phi(h)^2\right]\\
&=v(h+k).
\end{aligned}
$$

Thus the two-stage Gaussian transition has exactly the same mean and variance, and hence the same Gaussian distribution, as the single transition of duration $h+k$. This is an explicit instance of the transition-composition property, also called the semigroup property.

### 12.14 A numerical transition with every quantity identified

::: {.worked}
Suppose the current observed value is $x=29$, and use the candidate model parameters $\theta=1$, $\mu=25$, $\sigma=2$. For a gap $\Delta=0.5$, the surviving fraction is

$$
\phi=e^{-0.5}\approx0.60653066.
$$

The current deviation from the target is $29-25=4$, so the predicted next mean is

$$
m=25+(0.60653066)(4)\approx27.42612264.
$$

The stationary variance corresponding to these parameters is $q=2^2/(2\cdot1)=2$. The conditional variance over this half-unit gap is smaller:

$$
v=2(1-e^{-1})\approx1.26424112.
$$

Its standard deviation is $\sqrt v\approx1.12438$. Therefore the next value, conditional on the current reading and these parameters, has distribution approximately $N(27.4261,1.2642)$. The prediction 27.4261 is not a guaranteed next value.

| Gap $\Delta$ | Retained fraction | Conditional mean | Conditional variance |
|---:|---:|---:|---:|
| 0.10 | 0.9048 | 28.6193 | 0.3625 |
| 0.50 | 0.6065 | 27.4261 | 1.2642 |
| 1.00 | 0.3679 | 26.4715 | 1.7293 |
| 3.00 | 0.0498 | 25.1991 | 1.9950 |

The same starting value has different predictions and uncertainties at different horizons. This is why timestamps carry information about the return speed. The next chapter keeps that information inside every likelihood factor.
:::

<p class="source-note">Source connection: the OU solution is treated in <a href="#ref-l21">MIT Lecture 21, §1.2</a> and <a href="#ref-ss">Särkkä and Solin, Examples 4.5 and 6.2</a>. Here the solution is obtained through an ordinary equation for the noise-subtracted process, then verified directly; its distributions, stationarity and covariance are derived from the resulting formula.</p>


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


## Reference card: a formula and the conditions that make it true {#reference-card}

This card is a return point after working through the derivations. In particular, distinguish the variance of fresh noise from a state variance, a conditional law from a stationary marginal law, and a model parameter from its estimate.

| Topic | Result | Conditions and interpretation |
|---|---|---|
| Coin MLE | $\widehat p=k/n$ | Independent tosses with a common probability; $n>0$; the data are fixed during maximization. Boundary cases use the permitted parameter space. |
| Random walk | $E[W_k]=km_S$, $\operatorname{Var}(W_k)=kv_S$ | Independent identically distributed finite-variance steps and fixed start zero. Positions are not independent merely because steps are. |
| Fair gambler's ruin | $f(m)=m/N$, $E_m[T]=m(N-m)$ | Unit up/down steps, fair probabilities, absorbing boundaries $0,N$. |
| Markov propagation | $r_{n+1}=Ar_n$ | Here $A_{ji}=p_{ij}$: outgoing distributions are columns. Homogeneity gives $r_n=A^nr_0$. |
| Stationary distribution | $A\pi=\pi$ | Fixed-point equation for a probability vector; convergence from other initial distributions is a separate claim. |
| Markov likelihood | $L_c=\prod_i p_\eta(x_i\mid x_{i-1})$ | Supplied starting value; correct observed-state Markov model. Include gaps when transitions depend on them. |
| Stable AR(1) | $q=\tau^2/(1-\phi^2)$ | $|\phi|<1$; stationary state variance, not one-step innovation variance. |
| Stationary AR correlation | $\rho(k)=\phi^k$ | Stationary initialization and independent innovations; fixed-start correlations can differ. |
| Gaussian AR fit | $\widehat b=S_{uy}/S_{uu}$, $\widehat a=\bar y-\widehat b\bar u$ | Unconstrained conditional regression, $S_{uu}>0$; check whether the solution is in the required parameter domain. |
| Diffusion scaling | noise per tick $=\sigma\sqrt\delta\,\xi_i$ | Centered unit-variance independent kicks; finite-variance diffusion scaling, not a universal description of all stochastic limits. |
| Brownian covariance | $\operatorname{Cov}(B_s,B_t)=\min(s,t)$ | Standard Brownian motion. This is covariance of positions; disjoint increments have zero covariance. |
| Deterministic Wiener integral | $\int f\,dB\sim N(0,\int f^2)$ | Deterministic square-integrable weights. General adapted random weights need not produce Gaussian integrals. |
| OU conditional mean | $m_\Delta(x)=\mu+e^{-\theta\Delta}(x-\mu)$ | $\theta>0$ and a positive fixed or exogenously selected gap. |
| OU conditional variance | $v_\Delta=\frac{\sigma^2}{2\theta}(1-e^{-2\theta\Delta})$ | Fresh accumulated noise over the supplied gap; exact for the stated OU model. |
| Stationary OU | $q=\sigma^2/(2\theta)$, $\rho(h)=e^{-\theta|h|}$ | Initialize with the invariant Gaussian independently of future noise for stationarity from time zero. |
| AR-to-OU conversion | $\theta=-\log\phi/\Delta$, $\sigma^2=2\theta\tau^2/(1-\phi^2)$ | Fixed positive gap and $0<\phi<1$. A negative stable AR coefficient is not an embedding in this scalar OU model. |
| Exact conditional NLL | $\frac12\sum_i[\log(2\pi v_i)+(x_i-m_i)^2/v_i]$ | Use all actual gaps and keep parameter-dependent normalization. Exact model specification does not imply unbiased finite-sample estimation. |

### What is proved here, and what is used as a theorem? {#proof-status}

The book explicitly derives the coin maximization, random-walk moments and covariance, first-step boundary probabilities and fair expected duration, two-state stationary distribution and convergence, Markov likelihood factorization, Gaussian AR regression estimates, time-scaling calculations, Brownian covariance and quadratic variation along deterministic partitions, deterministic-integral mean and variance, the OU solution and uniqueness, its transition and stationary laws, and conditional OU profile formulas. Each derivation has its modeling and initialization conditions attached.

Some general mathematical results are used with their role identified: the central limit theorem and functional central limit theorem; existence of Brownian motion; the global nowhere-differentiability theorem; approximation by step functions and completeness in square-integrable spaces; standard uniqueness facts for distributional limits or moment-generating functions; and the general finite irreducible aperiodic Markov-chain convergence theorem. For several of these, the notes also prove a narrower calculation illustrating the result—for example the explicit two-state convergence formula and a fixed-time Brownian derivative argument. A narrower calculation is not presented as a proof of the entire general theorem.

## Sources and further reading {#sources}

These notes are an independently written teaching exposition with original worked calculations. The following primary teaching and research sources provide the standard definitions, results and wider context. Links require internet access; no source PDFs, copied textbook figures, external rendering libraries or font files are bundled with this document.

### Likelihood and discrete processes

<p id="ref-psu"><strong>Penn State, STAT 504, Lesson 1, §1.5: Maximum Likelihood Estimation.</strong> The Bernoulli and binomial discussion separates the observed sample from the unknown parameter and compares likelihoods for the ordered sample and its count. <a href="https://online.stat.psu.edu/stat504/Lesson01">Open the official lesson.</a></p>

<p id="ref-l5"><strong>MIT 18.S096, Fall 2013, Lecture 5: Stochastic Processes I.</strong> Random walks, gambler's ruin, first-step analysis and finite Markov chains. The matrix convention in these notes follows that lecture's column orientation. <a href="https://ocw.mit.edu/courses/18-s096-topics-in-mathematics-with-applications-in-finance-fall-2013/f5784e4facf3de690210d17c97358eba_MIT18_S096F13_lecnote5.pdf">Open the lecture PDF.</a></p>

<p id="ref-markov"><strong>Berkeley Prob 140 textbook: Long Run Behavior.</strong> A supplementary course reference for the finite irreducible aperiodic chain convergence theorem and its stationary-distribution interpretation. That text uses its own matrix convention; transpose consistently when comparing it with the column convention here. <a href="https://data140.org/fa18/textbook/chapters/Chapter_10/03_Long_Run_Behavior">Open the course chapter.</a></p>

<p id="ref-l8"><strong>MIT 18.S096, Fall 2013, Lecture 8: Time Series Analysis I.</strong> AR(1) and likelihood-based estimation. The stationary variance in this book is derived as τ²/(1−φ²); the squared coefficient follows from the variance of a scaled variable. <a href="https://ocw.mit.edu/courses/18-s096-topics-in-mathematics-with-applications-in-finance-fall-2013/1926c83ecd7ea700f7cb63914c6d7c0f_MIT18_S096F13_lecnote8.pdf">Open the lecture PDF.</a></p>

### Brownian motion, stochastic integration and OU

<p id="ref-l17"><strong>MIT 18.S096, Fall 2013, Lecture 17: Stochastic Processes II.</strong> Path space, Brownian motion and path regularity. The normalized walk used here is Wₖ/√n at time k/n, and the finite-resolution interpolation is explicitly distinguished from the Brownian limit. <a href="https://ocw.mit.edu/courses/18-s096-topics-in-mathematics-with-applications-in-finance-fall-2013/3b97c6b0c282dd9dc024c4c7ffe3fba8_MIT18_S096F13_lecnote17.pdf">Open the lecture PDF.</a></p>

<p id="ref-l18"><strong>MIT 18.S096, Fall 2013, Lecture 18: Itô Calculus.</strong> Motivation for stochastic integration, deterministic and adapted integrands, and the Itô isometry. The weighted-sum construction in Chapter 11 explains the deterministic case before discussing why random integrands are different. <a href="https://ocw.mit.edu/courses/18-s096-topics-in-mathematics-with-applications-in-finance-fall-2013/ef2c66c8079ba656210ad1fd4a5e2fa8_MIT18_S096F13_lecnote18.pdf">Open the lecture PDF.</a></p>

<p id="ref-l21"><strong>MIT 18.S096, Fall 2013, Lecture 21: Stochastic Differential Equations.</strong> Integral interpretation, existence-and-uniqueness conditions, and the OU solution in §1.2. Chapter 12 here gives a direct pathwise derivation by subtracting the Brownian term and solving the resulting ordinary equation. <a href="https://ocw.mit.edu/courses/18-s096-topics-in-mathematics-with-applications-in-finance-fall-2013/a671d0d2abe626c371fd5850ad670397_MIT18_S096F13_lecnote21.pdf">Open the lecture PDF.</a></p>

<p id="ref-ss"><strong>Simo Särkkä and Arno Solin, <em>Applied Stochastic Differential Equations</em>, Cambridge University Press, 2019.</strong> The author-hosted book provides a broader treatment. Relevant places include §4.1, Example 4.5 (OU solution), Example 6.2 (exact OU discretization), and §6.5 (steady-state behavior). <a href="https://users.aalto.fi/~ssarkka/pub/sde_book.pdf">Open the author-hosted PDF.</a></p>

### Irregular observation times and estimation

<p id="ref-am"><strong>Yacine Aït-Sahalia and Per A. Mykland, “The Effects of Random and Discrete Sampling When Estimating Continuous-Time Diffusions,” <em>Econometrica</em> 71 (2003), 483–549.</strong> The NBER working-paper version is Technical Working Paper 0276 (2002). This is a primary research reference for separating discrete-sampling effects from random-spacing effects and the consequences of ignoring random sampling. The elementary population illustration in Chapter 13 states its own assumptions and is not a replacement for that paper's general asymptotic theory. <a href="https://www.nber.org/papers/t0276">Open the NBER paper page and publication record.</a></p>

## Reproducible likelihood example {#code-example}

The following Python program reproduces the numerical candidate and profile tables in Chapter 13 using only the standard library. Save it as `ou_likelihood_examples.py` and run it with Python 3.10 or later. It evaluates likelihoods; it deliberately does not call its three candidate speeds a global MLE search. The observation assumptions are part of the mathematics, not conditions that this code can infer merely from the numerical arrays.


```python
"""Reproduce the OU likelihood calculations in the expanded teaching notes.

Python 3.10+; standard library only. No external packages are needed.

The observation times must be fixed or independent of the entire OU path.
These are conditional likelihoods: the initial value is supplied, not modeled
by a parameter-dependent stationary density. Functions evaluate objectives;
this module does not claim to find a global maximum-likelihood estimate.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from collections.abc import Sequence


@dataclass(frozen=True)
class ProfileFit:
    theta: float
    mu: float
    sigma2: float
    nll: float


def validate_data(
    values: Sequence[float], gaps: Sequence[float]
) -> tuple[tuple[float, ...], tuple[float, ...]]:
    """Require n positive gaps and n+1 finite, unnoised observations."""
    x, dt = tuple(map(float, values)), tuple(map(float, gaps))
    if not dt or len(x) != len(dt) + 1:
        raise ValueError("Need n >= 1 gaps and exactly n+1 observations")
    if not all(math.isfinite(v) for v in x):
        raise ValueError("Observations must be finite")
    if not all(math.isfinite(h) and h > 0 for h in dt):
        raise ValueError("Gaps must be finite and strictly positive")
    return x, dt


def transition_coefficients(theta: float, gap: float) -> tuple[float, float, float]:
    """Return phi, b=1-phi, c such that the transition variance is sigma^2*c."""
    if not math.isfinite(theta) or theta <= 0:
        raise ValueError("theta must be finite and positive")
    if not math.isfinite(gap) or gap <= 0:
        raise ValueError("gap must be finite and positive")
    z = theta * gap
    if not math.isfinite(z) or z == 0:
        raise ArithmeticError("theta * gap is outside the supported numeric range")
    phi = math.exp(-z)
    b = -math.expm1(-z)
    # expm1 avoids cancellation at short gaps. The long-gap expression
    # avoids overflow in 2*z or 2*theta when evaluating a saturated transition.
    c = 0.5 / theta if z > 350 else gap * (-math.expm1(-2 * z) / (2 * z))
    if not math.isfinite(c) or c <= 0:
        raise ArithmeticError("Nonpositive or non-finite transition coefficient")
    return phi, b, c


def conditional_nll(
    theta: float, mu: float, sigma: float,
    values: Sequence[float], gaps: Sequence[float]
) -> float:
    """Evaluate the exact conditional OU negative log-likelihood.

    The recorded predecessor, not the previous predicted mean, is used in
    each transition. All parameter-dependent variance terms are retained.
    """
    x, dt = validate_data(values, gaps)
    if not math.isfinite(mu) or not math.isfinite(sigma) or sigma <= 0:
        raise ValueError("mu must be finite; sigma must be finite and positive")
    terms = []
    for previous, observed, h in zip(x[:-1], x[1:], dt):
        phi, b, c = transition_coefficients(theta, h)
        mean = phi * previous + b * mu
        variance = (sigma * sigma) * c
        if not math.isfinite(variance) or variance <= 0:
            raise ArithmeticError("Transition variance is outside numeric range")
        standardized = (observed - mean) / math.sqrt(variance)
        term = 0.5 * (math.log(2 * math.pi) + math.log(variance)
                      + standardized * standardized)
        if not math.isfinite(term):
            raise ArithmeticError("Non-finite likelihood contribution")
        terms.append(term)
    return math.fsum(terms)


def profile_at_theta(
    theta: float, values: Sequence[float], gaps: Sequence[float]
) -> ProfileFit:
    """Analytically fit free mu and positive sigma^2 at a FIXED theta.

    Return the one-dimensional profile objective plus the associated level
    and variance. A zero residual sum is degenerate, not an interior MLE.
    """
    x, dt = validate_data(values, gaps)
    coefficients = [transition_coefficients(theta, h) for h in dt]
    b = [item[1] for item in coefficients]
    c = [item[2] for item in coefficients]
    d = [next_value - item[0] * previous
         for previous, next_value, item in zip(x[:-1], x[1:], coefficients)]
    denominator = math.fsum(bi * bi / ci for bi, ci in zip(b, c))
    numerator = math.fsum(bi * di / ci for bi, di, ci in zip(b, d, c))
    if not math.isfinite(denominator) or denominator <= 0:
        raise ArithmeticError("Degenerate profile weights")
    mu = numerator / denominator
    squared_error = math.fsum((di - bi * mu) ** 2 / ci
                             for bi, di, ci in zip(b, d, c))
    sigma2 = squared_error / len(dt)
    if not math.isfinite(mu) or not math.isfinite(sigma2) or sigma2 <= 0:
        raise ValueError("Degenerate or non-finite profiled fit")
    nll = 0.5 * (len(dt) * (math.log(2 * math.pi) + 1 + math.log(sigma2))
                 + math.fsum(math.log(ci) for ci in c))
    if not math.isfinite(nll):
        raise ArithmeticError("Non-finite profile objective")
    return ProfileFit(theta, mu, sigma2, nll)


def main() -> None:
    values = (29.0, 27.5, 26.0, 25.8)
    times = (0.0, 0.2, 1.0, 1.3)
    gaps = tuple(right - left for left, right in zip(times[:-1], times[1:]))
    print("Candidate evaluations, NOT a global MLE search")
    print("theta   fixed(mu=25,sigma=2) NLL     profiled mu    profiled sigma^2    profile NLL")
    for theta in (0.5, 1.0, 2.0):
        nll = conditional_nll(theta, 25.0, 2.0, values, gaps)
        fit = profile_at_theta(theta, values, gaps)
        print(f"{theta:4.1f}    {nll:24.8f}  {fit.mu:14.8f}  {fit.sigma2:18.8f}  {fit.nll:13.8f}")
        # The direct and profiled formulas must describe the same objective
        # when evaluated at identical parameters and data.
        direct = conditional_nll(theta, fit.mu, math.sqrt(fit.sigma2), values, gaps)
        if not math.isclose(direct, fit.nll, rel_tol=1e-10, abs_tol=1e-10):
            raise AssertionError("Profile and direct objective disagree")


if __name__ == "__main__":
    main()
```
