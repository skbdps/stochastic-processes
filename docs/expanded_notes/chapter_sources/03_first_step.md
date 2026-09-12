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
