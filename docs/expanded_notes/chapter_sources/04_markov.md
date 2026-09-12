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
