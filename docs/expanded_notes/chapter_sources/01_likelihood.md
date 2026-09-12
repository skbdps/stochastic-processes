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
