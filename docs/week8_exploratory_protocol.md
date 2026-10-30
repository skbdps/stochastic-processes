# Week 8 E4: event-stationary clustered sampling and full-marginal controls

Dated 30 October 2026, after completed Week 7 publication
`d8e921df5a395cd17b62be51c6b8b5fabd7f9221` and the
[core-completion gate](week8_gate.md). This protocol is fixed before Week 8
implementation, development observations or full E4 observations. It is an
exploratory extension, not an amendment to the registered E1/E2 experiments.

## 1. Fixed question and complete design

Compare exogenous clustered observation gaps with independent gaps having
the **same entire event-stationary marginal distribution**. A same-CV Gamma
control is insufficient for this comparison. The intensity-ratio curve alone
changes both marginal gaps and dependence and cannot isolate the latter.

Retain the original extension settings: truth (theta,mu,sigma)=(1,0,0.5),
q=theta*m=0.5, n=1000 observations and K=n-1=999 transitions. Use calm and
burst mean sojourns 10m and 2m, respectively, and intensity ratios r=1,5,20.
This resolves the original sojourn-order convention explicitly. Use exactly
R=2000 independent source paths per cell, matching the E1/E2 per-cell scale.

The **seven fixed cells** are:

1. Event-stationary MMPP, r=1,5,20: three cells.
2. IID full-Palm-marginal gaps, r=1,5,20: three cells.
3. One independent regular reference with gap m: one cell.

Total: 14,000 independent source paths and 42,000 fitted records. All three
estimators share each source path. All seven source cells use distinct
scientific identities and independent streams; matching replication numbers
does not pair different cells. The two r=1 cells are independent draws from
the same exponential-gap law, not three copies of the regular control.

## 2. Model calibration and initialization

Use state order (calm,burst), row probability vectors and

$$
Q=m^{-1}\begin{pmatrix}-1/10&1/10\\1/2&-1/2\end{pmatrix},
\qquad \pi=(5/6,1/6).
$$

Let lambda_b=r*lambda_c and calibrate the raw time-average arrival rate to 1/m:

$$
\lambda_c=\frac{6}{(r+5)m},\quad
\lambda_b=\frac{6r}{(r+5)m},\quad
\Lambda=\operatorname{diag}(\lambda_c,\lambda_b),\quad
\alpha=\frac{\pi\Lambda}{\pi\Lambda\mathbf1}
=\left(\frac5{r+5},\frac r{r+5}\right).
$$

Time zero is an event. Initialize the MMPP phase from alpha and retain the
arrival phase between gaps. This is event-stationary initialization with no
burn-in. Initializing from pi at an arbitrary calendar time instead produces
a forward residual waiting time; discarding one such gap does not in general
make the resulting phase distribution alpha. That is not this design.

While in phase j, draw the next arrival-or-switch wait as exponential with
rate lambda_j+q_j. The event is an arrival with probability
lambda_j/(lambda_j+q_j); an arrival retains j. Otherwise switch phase and
continue accumulating the gap. Generate the joint gap/ending-phase outcome;
do not sample a marginal gap and a phase transition independently.

With D0=Q-Lambda and A=-D0, the whole-gap law is

$$
\Pr(T>t)=\alpha e^{D_0t}\mathbf1,\qquad
P=A^{-1}\Lambda,\qquad
E[T^k]=k!\alpha A^{-k}\mathbf1.
$$

The IID control redraws an independent phase from alpha at **every** gap,
then runs the same competing arrival/switch mechanism to absorption. Hidden
switches within a gap remain allowed. Choosing alpha followed by a single
Exp(lambda_j) draw would give the wrong marginal. The regular reference uses
constant gap m. These standard MAP/phase-type identities are supported by
[Asanjarani and Nazarathy, Eq.2.1 and Proposition 2.2](https://doi.org/10.1017/S1446181122000013)
and the [PH-renewal/MMPP construction in Osogami's thesis](https://www.cs.cmu.edu/~osogami/thesis/html/node72.html).
The specialized derivations are in [Week 8 mathematics](week8_mathematics.md).

## 3. Floor, OU path and source addressing

Retain a nominal-CV trigger and no renormalization. Here nominal CV means the
**derived raw event-gap CV**, not a Gamma input parameter:

$$
\mathrm{CV}^2=1+\frac{50(r-1)^2}{3(r^2+70r+25)}.
$$

For r=1,5,20, CV is 1, sqrt(5/3), sqrt(941/219), respectively. Apply
Y=max(T,1e-6*m) only when this raw CV is at least 1.5: only r=20 in this
design. Apply the same transformation to its MMPP and IID control. Do not
advance the hidden CTMC through inserted floor padding. The transformed
sequence is therefore described as **floored MMPP-derived gaps**, not the
original MMPP on an unchanged clock. Its full one-gap marginal still matches
the transformed IID control, but its mean need not equal m exactly. Retain
raw/transformed moments, minimum gaps, clip counts and realized spans.

Generate exact OU transitions at the stored transformed observation times,
with an independent stationary initial OU draw. The sampling mechanism is
independent of the entire OU initial/innovation stream. Conditional fitted
likelihood excludes the stationary initial density. Preserve the original
estimators, full Hessian/Wald procedure, failure eligibility rules and
thresholds without any source-specific numerical tuning.

Full root seed is **20261011**. Scientific addresses are
`SeedSequence(seed, spawn_key=(8,*eight_cell_id_words,rep)).spawn(2)`:
child 0 supplies schedule/phase randomness, child 1 stationary X0 and OU
innovations. No worker, partition, execution order or requested total R enters
the cell address. The cell's SHA-256 identity binds the full flat scientific
configuration: source mechanism, ratio, nominal design/truth, calibration,
sojourns, floor policy version and namespace 8. The regular cell uses r=0 as
an explicit inapplicable sentinel; no MMPP formula is evaluated at that
sentinel. The old registered Gamma source/address policy is not used for E4.

## 4. Fixed development validation before full execution

Development fit runs use the complete seven-cell design, R=4 and root seed
21261011 (full seed +1,000,000). They are validation records and cannot supply
E4 results. Source-only analytical validation uses root seed **22261011**,
256 independent length-999 gap sequences for each of the six stochastic
mechanism/ratio cells, with the same namespace/address construction. It
generates no fitted OU experiments and remains separate evidence.

Before full E4 observations:

- Verify exact matrix/closed-form identities, stationary probabilities,
  calibrated means, CVs, raw and transformed moment formulas, and the r=1
  exponential limit against independent calculations.
- Check deterministic seed replay, source independence, complete identity
  binding, correct event-phase retention versus IID reset, exact regular
  gaps, and the declared posthoc floor transform on constructed inputs.
- For each stochastic validation cell, record per-sequence raw and transformed
  first and second moments, adjacent-gap cross-products, and empirical CDF
  values at t/m in {0.1,0.5,1,2,5}. Moment/CDF averages use all 999 gaps;
  adjacent products use the 998 adjacent pairs. Empirical CDF uses the
  indicator gap<=t. Compare their across-sequence means with analytical
  expectations, using six MCSEs plus 1e-10 absolute tolerance. Across the
  256 sequences, MCSE is the sample standard deviation (ddof=1) divided by
  sqrt(256).
  MCSE is calculated across the 256 independent sequences, never by treating
  correlated MMPP gaps as independent observations. Save all per-sequence
  statistics and the comparison outcomes.
- Compare exactly the first phase of each of the 256 independent sequences
  with alpha, for each mechanism/ratio, using six theoretical binomial
  standard errors with N=256. Do not pool dependent MMPP event phases into
  an independent-binomial check. Record clipping frequency, but do not require
  an empirical rare-event normal interval to contain its tiny theoretical
  probability. Floor correctness also has deterministic tests and exact
  transformed-law calculations.
- Preserve every failed validation and its original inputs. Investigate any
  discrepancy; do not change seeds, R or tolerances to obtain acceptance.

The six-MCSE diagnostic is a prespecified development check, not a research
confidence claim. Tests and validation support only their stated scope.
Publish the reviewed complete implementation freeze before full R=2000 data.

## 5. Estimands, contrasts and interpretation

Retain all planned records, including finite candidates and failures. Separate
optimizer success, point validity, interval availability and pipeline errors.
Bias/RMSE condition on valid points; operational coverage counts unavailable
intervals as noncoverage in all R planned paths; conditional coverage uses
valid intervals. Use the existing uncertainty formulas and practical criteria.
Mu relative metrics remain undefined. All inference is pointwise; PFML/Euler
Wald covariance remains model-based, not sandwich corrected.

Report complete per-cell summaries and paired within-path method comparisons.
For each r/method/parameter, contrast MMPP minus the independently generated
matched IID control, using independent-sample Monte Carlo uncertainty and
independent RMSE-ratio delta calculations. Report each stochastic cell against
the single independent regular reference separately. Contrasts sharing that
regular sample are not mutually independent. Preserve each component's valid
point count; independent cells have no common-valid paired subset.

The full-marginal MMPP/IID comparison assesses consequences of gap dependence
within this fixed design, including its effect on the distribution of realized
span. It does not hold realized span fixed. No direct dependence-only claim
may be made from an r-curve or a regular comparison alone. Null, mixed or
adverse outcomes receive the same complete reporting as favorable outcomes.

## 6. Audit, figures and completion

Use a new source-specific audit. Check full planned keys, actual E4 metadata,
exact schema including empty failure exports, source addresses, policy and
dependency hashes, raw-to-summary arithmetic and independent contrasts.
Reconstruct every one of the 14,000 source realizations and its observation/time
hashes. Replay complete original fits at fixed replications 0 and R-1 in each
cell (14 sources / 42 fits). Preserve every exception; pipeline errors block
completion. Keep E1/E2/E3 archives and source files unchanged.

Render mechanism/ratio coverage, availability, matched-control contrasts and
source-spacing diagnostics from saved tables. Bind their inputs and inspect
every PNG before guarded sealing. Preserve renderer state and invocation logs,
then publish a completion report, mathematical interpretation, concrete
reproduction commands and limitations. Do not infer completion from code or
configuration alone. No professor contact, outside approval, novelty claim,
personal reading certificate or final-manuscript claim is part of this stage.
