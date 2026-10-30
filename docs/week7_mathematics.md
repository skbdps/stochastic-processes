# Week 7 mathematical diagnostics

This note documents the selected-path component of the exploratory protocol
dated 29 October 2026 UTC. Its formulas explain diagnostics of the unchanged
Week 3/4 estimators and inference procedure. They do not define replacement
estimates, covariance matrices, intervals, headline classifications or a new
registered method. The reconstructed implementation is reviewed and tested
against its current bytes; historical tests of the lost implementation do not
certify this version.

## Conditional likelihood and fixed-theta profiling

Write the observed, normalized states as \(y_0,\ldots,y_K\), with observed
positive gaps \(d_1,\ldots,d_K\); \(K=n-1\). At a fixed positive \(\theta\),
define

\[
 \phi_i=e^{-\theta d_i},\qquad b_i=1-\phi_i,\qquad
 a_i=\frac{1-e^{-2\theta d_i}}{2\theta},\qquad
 z_i=y_i-\phi_i y_{i-1}=\Delta y_i+b_i y_{i-1}.
\]

The implementation evaluates the two subtractions from one with `expm1` and
uses the last expression for \(z_i\) to avoid cancellation when gaps are small.
Conditional on the initial observation, \(z_i\) has mean \(b_i\mu\) and
variance \(\sigma^2a_i\). No density for \(y_0\) is added. With
\(A=\operatorname{diag}(a_i)\), the conditional negative log likelihood is

\[
 \ell(\theta,\mu,\sigma)=\frac12\left[
 K\log(2\pi)+2K\log\sigma+\sum_i\log a_i+
 \frac{(z-b\mu)^T A^{-1}(z-b\mu)}{\sigma^2}\right].
\]

For fixed \(\theta\), its nuisance optimum is

\[
 \widehat\mu(\theta)=\frac{b^T A^{-1}z}{b^T A^{-1}b},\qquad
 \widehat\sigma^2(\theta)=\frac{(z-b\widehat\mu)^T A^{-1}
          (z-b\widehat\mu)}{K}.
\]

The likelihood variance divisor is **K**, not n and not a degrees-of-freedom
adjusted divisor. When the residual variance is positive, the profiled NLL is

\[
 p(\log\theta)=\frac12\left[K\{\log(2\pi)+1+
             \log\widehat\sigma^2(\theta)\}+\sum_i\log a_i\right].
\]

Degenerate normal equations, zero residual variances and nonrepresentable
values are recorded as numerical/degenerate outcomes. No truth value,
regularization or arbitrary residual-variance floor is substituted.

The independent full NLL uses the same transition identity, with residual
\(e_i=\Delta y_i+b_i(y_{i-1}-\mu)\), log variance
\(v_i^{\log}=2\eta_3+\log a_i\), and
\(\eta=(\log\theta,\mu,\log\sigma)\):

\[
 \ell(\eta)=\frac12\sum_i\{\log(2\pi)+v_i^{\log}
                       +e_i^2e^{-v_i^{\log}}\}.
\]

This is independently coded rather than calling the historical exact objective.
The historical objective remains the objective for the registered float64
Hessian ladder and the original-fit replay.

## Two likelihood limits

For the lower limit, \(\alpha=\theta\mu\) can remain finite while
\(\theta\downarrow0\) and \(\mu\) diverges. The limiting increments have
mean \(\alpha d_i\) and variance \(\sigma^2d_i\). Profiling gives

\[
 \widehat\alpha=\frac{\sum_i\Delta y_i}{\sum_i d_i},\qquad
 \widehat\sigma_B^2=\frac1K\sum_i
              \frac{(\Delta y_i-\widehat\alpha d_i)^2}{d_i},
\]

\[
 \ell_B=\frac12\left[K\{\log(2\pi)+1+\log\widehat\sigma_B^2\}
                       +\sum_i\log d_i\right].
\]

This is the **Brownian-with-free-drift** limit. Holding a finite \(\mu\)
before taking the limit would impose zero drift and answer a different
likelihood question.

For \(\theta\to\infty\), the positive-gap conditional transition correlations
vanish. Let \(v=\sigma^2/(2\theta)\) remain finite. The limiting likelihood
is that of K independent normal responses \(y_1,\ldots,y_K\), with

\[
 \widehat\mu_I=K^{-1}\sum_i y_i,\qquad
 \widehat v_I=K^{-1}\sum_i(y_i-\widehat\mu_I)^2,\qquad
 \ell_I=\tfrac K2\{\log(2\pi)+1+\log\widehat v_I\}.
\]

This closure likelihood does not supply a finite OU estimate: both theta and
sigma diverge along a fixed positive stationary variance. Neither closure is a
new point fit or a basis for a replacement interval.

## AR coefficient and scope of the boundary proof

At one constant gap, conditional OU transitions are an intercept-and-slope
Gaussian regression, with \(\phi=e^{-\theta d}\in(0,1)\). Profiling the
intercept and variance leaves a positive multiple of the logarithm of

\[
 \operatorname{SSE}(\phi)=\sum_i(y_i-\bar y_+-
                                  \phi(y_{i-1}-\bar y_-))^2.
\]

If predictor variation is positive, this is a convex quadratic, minimized over
the real line at the ordinary least-squares slope. For a nonpositive slope,
the infimum over \(0<\phi<1\) is approached at \(\phi=0\); the corresponding
upper-theta likelihood supremum has no finite interior estimate. A slope at
least one instead gives the lower-theta boundary. A zero regression residual
variance is a separate degeneracy, not a finite-variance boundary estimate.

The constant-gap model is exactly the PFML working model at the **realized
mean gap**, for any observed spacing pattern. The proof also applies to exact
likelihood when all observed gaps are equal. It does **not** establish an exact
likelihood boundary from an OLS slope under irregular gaps. The diagnostic
records exact gap equality separately from the historical numerical
near-regular shortcut, whose tolerance is
\(128K\epsilon_{64}\) for relative deviations from the observed mean gap.
Calling a sample near regular is a numerical recognition rule, not a proof
that unequal gaps are mathematically equal. The actual exact objective always
uses its observed gaps.

## Fixed finite search and numerical resolution

The profile scan follows the dated rule

\[
 L=\max[-40,\min\{-30,\log(10^{-10})-\log(\max_i d_i)\}],
\]

\[
 U=\min[40,\max\{12,\log(50)-\log(\min_i d_i)\}].
\]

Conversions to float64 round outward within the fixed caps. The saved record
states whether \(e^L\max d_i\le10^{-10}\) and
\(e^U\min d_i\ge50\) are attained. Cap-limited tails stay unresolved. A
linspace is refined enough that **every actual grid spacing is at most .05**.
Every finite interior sampled local-minimum bracket, including ties and
plateaus, is polished by bounded scalar minimization with `xatol=1e-12` and
`maxiter=500`.

Five genuine one-dimensional Nelder–Mead searches receive starting points
equal to the original valid exact fit's normalized log theta plus
\(\{-2,-1,0,1,2\}\), clipped to [L,U]. Each starts from an explicit two-point
simplex of width .1, reversing direction at the upper bound. The fixed options
are `xatol=fatol=1e-10`, `maxiter=2000`, `maxfev=4000`, `adaptive=False`.
Every start, simplex, setting, convergence outcome and returned candidate is
retained, including unsuccessful searches. A solver exception retains its
error and the existing bracket center or initial point explicitly as such;
that point is not described as a returned optimizer solution. The remaining
fixed brackets and starts are still attempted. Optimizer callbacks are float64
because the optimizer is float64; every returned candidate is independently
reevaluated in longdouble for comparisons.

The original theta with profiled nuisances, the original **full** three-parameter
objective, grid extrema, both endpoints and both limiting objectives are kept
distinct. A lower profiled value at the original theta can be caused by nuisance
profiling and must not be reported as a different theta optimum.

The diagnostic resolution is
\(10^{-7}+10^{-9}|\ell_{\rm best}|\). Extended-precision objective differences
are saved as decimal strings as well as plotting floats. Competitive limiting
likelihoods and cap-limited tails remain explicitly unresolved. Finite searches
do not prove global optimality, and a search improvement by itself does not
establish a verified software defect. This resolution rule does not replace
any registered fit or inference threshold.

## Full Hessian, profile curvature and observed-data representations

Use the unchanged `ObservationCoordinates` transformation, with recorded state
offset a, state scale b and mean-gap scale m estimated only from observations:
\(y=(x-a)/b\), \(u=t/m\). In the working parameter order
\(\eta=(\log\theta_u,\mu_u,\log\sigma_u)\), the complete conditional
NLL Hessian has dimension three. Its four diagnostic step vectors are

\[
 h_j=c\max(1,|\eta_j|),\qquad
 c\in\{.001,.0005,.00025,.000125\}.
\]

The float64 ladder delegates to the historical central-difference stencil.
The additional exact ladder computes the same full NLL and stencil independently
in longdouble, at the original fitted float64 eta and the same float64 step
vectors lifted to longdouble. All four attempts survive even if derivatives
are nonfinite, an eigenvalue is nonpositive or a score cannot be calculated.

For two adjacent completed matrices \(H_c,H_f\), the relative Frobenius
comparison is \(\|H_c-H_f\|_F/\|H_f\|_F\). If the fine matrix is positive
definite and \(H_f=CC^T\), the directional comparison is

\[
 \max|\lambda\{C^{-1}(H_c-H_f)C^{-T}\}|.
\]

This is the spectral radius of a symmetric whitened discrepancy, equivalently
the largest absolute generalized eigenvalue of \((H_c-H_f,H_f)\). It is not a
Frobenius norm of the whitened matrix or a comparison of diagonal entries.
The score is \(\sqrt{g^T H^{-1}g}\), reported only for a positive definite
matrix. No inverse regularization is permitted. Eigenvalues, condition,
gradient, Hessian and missing/failed diagnostics are all retained. Converting
the **completed** longdouble matrix and gradient to float64 for linear algebra
does not imply that their stencil arithmetic used float64.

At a regular nuisance optimum, partition the full Hessian between log theta
and nuisance parameters \(\nu=(\mu,\log\sigma)\). Differentiating the nuisance
score equation gives

\[
 p''(\log\theta)=H_{\theta\theta}
       -H_{\theta\nu}H_{\nu\nu}^{-1}H_{\nu\theta}.
\]

Thus profile curvature is a Schur complement; it is not the theta diagonal of
the full Hessian. Each completed matrix records `theta_hessian_diagonal` and
`theta_schur_complement_at_stored_point`; despite the short field names, both
refer to the first **working coordinate log(theta_u)**, not natural theta_x.
The latter is an algebraic matrix quantity at the original stored point. It
equals actual profile curvature only when nuisances are stationary on that
same data representation and their Hessian block is invertible. The recorded
gradient permits this distinction to be assessed. Singular nuisance blocks
leave the complement unavailable. The registered covariance uses the full matrix and its
existing gates, followed by the natural-scale Jacobian
\(\operatorname{diag}(\theta_x,b,\sigma_x)\). No diagnostic ladder matrix is
used to construct an alternative interval.

Two floating-point representations are purposefully distinguished:

1. The independent full-NLL Hessian ladder lifts the exact float64 normalized
   states and `diff(normalized_sample.times)` used by registered inference.
   Substituting original gaps divided by m would change some rounded numbers.
2. The profile scan converts the original float64 states and times to
   longdouble first, then applies the same recorded a,b,m in longdouble. Its
   gap calculation is `diff(original_times_longdouble)/m`.

Maximum state/gap differences and the full NLL at the original fitted eta under
both representations are recorded, together with the registered float64 NLL.
The density Jacobian is \(\ell_x=\ell_y+K\log b\). Representation rounding
must not be mistaken for a scientific change in likelihood or for evidence that
one covariance procedure has been substituted for another. The platform's
actual fraction bits, significand bits, decimal precision and storage width
are serialized; 16 storage bytes alone do not imply 128-bit precision.

## Replay, failure retention and evidence boundaries

The fixed selection binds both its literal file bytes and semantic content.
All three original fit records are retained for every selected source,
including invalid candidates and unavailable intervals. Recreated times and
states have exact SHA-256 matches. Every original fit field is replayed: exact
keys, discrete values, hashes and missing values; floats use only the fixed
rtol \(2\times10^{-10}\) and atol \(2\times10^{-12}\). Embedded diagnostic
JSON is recursively compared. A failed replay prevents diagnostic analysis and
survives as an immutable failed-source record.

Each completed or failed source is an immutable JSON file. A root manifest
binds its digest, explicit dependency files, protocol, mathematical note,
configuration files, numerical environment and recorded Git freeze. Frozen
execution verifies that the recorded commit actually contains every recorded
dependency byte, before observations or output creation. Unfrozen development
execution rejects headline inputs. Closed diagnostic outputs cannot be resumed.
For incomplete outputs, committed failures are skipped exactly like successes.
An uncommitted next-source file can be adopted only after a deterministic replay
produces exactly the same bytes; different prior evidence is never overwritten.

Before and after each root-manifest replacement, and before source/audit
commits, a byte guard compares the root with its last known state. Unexpected
state is preserved together with the expected bytes in an adjacent incident
directory. The failure handler does not overwrite that observed state. A bound
incident also blocks subsequent automatic resumption even if someone restores
an earlier root file. This guards against observed state drift; an advisory
lock and file comparisons cannot prove the absence of every noncooperating
external filesystem writer or identify the cause of a previous incident.

The selected-source audit verifies all planned source/estimator keys, hashes,
retained originals, complete original-fit replay and immutable readback. Its
scope explicitly says that it does **not** independently recalculate every
diagnostic mathematical quantity. Independent formula tests and review cover
the implementation; later reporting must distinguish those validation layers.
These selected diagnostics are exploratory evidence, not a new Monte Carlo
coverage experiment or a correction to E1/E2.

## Paired floor source: actual cells and shared random addresses

The dated floor-sensitivity protocol fixes twelve actual cells with
\(n=1000\), \(q=\theta\bar d=.5\), truth
\((\theta,\mu,\sigma)=(1,0,.5)\), nominal CV in
\(\{0,1.5,2,3\}\), and epsilon in \(\{10^{-8},10^{-6},10^{-4}\}\).
The physical gap floor and the stored scientific cell ID always refer to the
**actual** epsilon. For each nominal CV, a distinct canonical cell at
epsilon \(10^{-6}\) supplies only the seed address, called the anchor.

Let \(c\) denote the actual cell and \(a(c)\) its same-CV anchor. Convert
the anchor's 64-character hexadecimal cell hash into the eight specified
32-bit words \(w_1,\ldots,w_8\). At replication r, the source forms

\[
 S=\operatorname{SeedSequence}(s,
           \operatorname{spawn\_key}=(4,w_1,\ldots,w_8,r)),
 \qquad (S_D,S_X)=S.\operatorname{spawn}(2).
\]

The adapter composes the unchanged `OuSampleSource`. It passes every actual
physical field unchanged, replacing only the internal mapping's `cell_id`
with the anchor ID for this address construction. This internal mapping is a
source-address override, not a newly declared canonical scientific cell. It
is not exported as the stored actual cell. The physical mean gap remains
\(\bar d=q/\theta=.5\); there is no realized-span matching.

For positive nominal CV \(c_v\), each epsilon variant uses the same raw draws

\[
 D_i\sim\operatorname{Gamma}(c_v^{-2},\bar d\,c_v^2),
 \qquad d_i^{(\epsilon)}=\max(D_i,\epsilon\bar d),
 \qquad i=1,\ldots,K,
\]

with \(K=999\). All positive CVs in this particular study meet the registered
floor threshold 1.5. No renormalization follows the floor. The same inherited
implementation constructs cumulative timestamps and records clip counts,
realized moments, spans and minimum timestamp gaps, including the effects of
floating-point timestamp accumulation. For CV=0, the gaps are regular and no
floor is applied, so all three variants have exactly the same times.

The second child generator is restarted at the same address for every epsilon
of a given (CV,r). It supplies a shared stationary initial draw and shared
standard normal innovations. In mathematical form,

\[
 X_0=\mu+\frac{\sigma}{\sqrt{2\theta}}Z_0,
\]

\[
 X_i^{(\epsilon)}=\mu+
  (X_{i-1}^{(\epsilon)}-\mu)e^{-\theta d_i^{(\epsilon)}}+
  \sqrt{\frac{\sigma^2}{2\theta}
      (1-e^{-2\theta d_i^{(\epsilon)}})}Z_i.
\]

The actual recursion uses the unchanged simulator and the differences of its
actual timestamps. The shared normals are independent of the shared raw Gamma
draws because gaps and states use separate seed children. Different CVs use
different anchor IDs and hence independent scientific addresses.

This is a common-random-number coupling of exact OU transition simulations.
It does **not** construct one continuous OU trajectory observed at three
different sets of times. If the floor changes gaps, it can change the total
span and subsequent states. Two epsilon variants can nevertheless coincide
on a finite path when no raw gap is affected by their different floors. For
CV=0 all three variants necessarily coincide, including their numerical fits:
that negative control is one shared sample, not three independent samples.

The adapter requires an exact canonical `Cell.to_dict()` mapping. Field names,
Python scalar types, finiteness, the derived mean gap, the supplied actual ID
and membership in the fixed twelve-cell set are checked before delegation.
An invalid ID is rejected; it is not replaced with an invented valid one.
Input mappings are copied, not mutated. The run design separately controls
the root seed and replication count: the exploratory design uses seed
20261009 and R=2000, while development uses seed+1,000,000 and a small fixed R.
Changing R does not change any retained source address.

Returned diagnostics preserve every inherited source diagnostic and add
`actual_cell_id`, `seed_anchor_cell_id`, `seed_anchor_floor_fraction` and
`source_coupling_policy`. The policy identifier is
`week7_paired_floor_v1_shared_gamma_stationary_x0_and_standard_normal_innovations`.
There is no diagnostic `cell_id` that could overwrite the actual row ID during
composition with the unchanged fitting engine. Collisions with scientific
row or coupling-identity fields are rejected. The actual row also retains its
own `floor_fraction`; source hashes and spawn keys retain their original
meanings.

Each realized variant is fitted and normalized from its own observed states
and times by the unchanged inference services. Scales are not borrowed from
another epsilon and simulation truth is not used for numerical normalization.
The future paired comparison and source audit must recognize this explicit
anchor policy; these records must not be labelled as independent registered
source-address experiments. This source module supplies no new estimator,
interval, contrast, or claim of a completed sensitivity experiment.

## Paired epsilon comparisons and their Monte Carlo uncertainty

The comparison module accepts a complete fixed floor design, or a complete
three-floor subset for each included CV when reporting in bounded batches.
It verifies the entire Cartesian set of actual cell, replication and estimator
keys before computing anything. Every included CV must contain epsilon
\(10^{-8},10^{-6},10^{-4}\). A missing cell, missing method, missing
replication or duplicate record blocks the paired table; no group is silently
dropped because only the groups present were inspected.

Every row must retain its canonical actual cell, registered-floor anchor,
coupling policy, root seed, expected source/gap/path spawn keys, and well-formed
source hashes. The expected address uses the same-CV registered-floor hash,
the given design's seed and that row's replication. Three estimators within
one actual cell/replication must share the same source diagnostics. Different
epsilon variants at positive CV are **not** required to have identical times
or observations: their common randomness is the shared Gamma draws and normal
innovations, and their actual floor transformations can differ.

For CV=0, all variants must instead retain identical source arrays and complete
numerical records after excluding the actual cell ID, actual-ID diagnostic and
actual floor fraction. This includes invalid numerical candidates, NLLs,
statuses and interval evidence, not only successful estimates. Inconsistent
control evidence blocks the table rather than being reported as an apparent
nonzero control effect. Raw source failures lacking verified coupling evidence
likewise block a certified paired table; their underlying raw evidence is not
deleted. Point/interval statistical failures with valid source metadata remain
ordinary retained outcomes. An inference pipeline error may retain a valid
point fit, but its interval must be unavailable; a nonempty pipeline error
combined with an interval marked valid is inconsistent and blocks the table.

For one positive or zero CV, one estimator and one parameter with truth
\(\psi\), let A be a nonreference epsilon and B be epsilon \(10^{-6}\).
Direction is always **A minus B** for differences and **A over B** for ratios,
serialized as `epsilon_minus_registered_floor` and
`epsilon_over_registered_floor`. There are no reference-to-itself rows.
The full design supplies
\(4\times2\times3\times3=72\) comparison rows; a single-CV three-floor
subset supplies 18.

### Point metrics on one common-valid subset

Let \(\mathcal S\) be the intersection of the two variants' valid point-fit
replications, with size S. Recompute **both** variants' errors and metrics on
this same subset:

\[
 e_{A,r}=\widehat\psi_{A,r}-\psi,\qquad
 e_{B,r}=\widehat\psi_{B,r}-\psi,\qquad r\in\mathcal S,
\]

\[
 B_j=\overline e_j,\qquad M_j=\overline{e_j^2},\qquad
 R_j=\sqrt{M_j},\qquad j\in\{A,B\}.
\]

The output retains each variant's own point-valid and interval-valid counts,
the common point-valid and interval-valid counts, and the two common-subset
biases, MSEs and RMSEs. Metrics from the two different own-valid subsets are
not subtracted and described as paired metrics.

For a paired vector \(v_r\), define
\(s(v)^2=\sum_{r\in\mathcal S}(v_r-\bar v)^2/(S-1)\). The sample
Monte Carlo uncertainties are

\[
 \operatorname{MCSE}(B_A-B_B)=\frac{s(e_A-e_B)}{\sqrt S},\qquad
 \operatorname{MCSE}(M_A-M_B)=\frac{s(e_A^2-e_B^2)}{\sqrt S}.
\]

For positive RMSEs, the delta-method influence vector for their difference is

\[
 v_r=\frac{e_{A,r}^2}{2R_A}-\frac{e_{B,r}^2}{2R_B},\qquad
 \operatorname{MCSE}(R_A-R_B)\simeq\frac{s(v)}{\sqrt S}.
\]

For the log RMSE ratio, the corresponding vector is

\[
 w_r=\frac12\left(\frac{e_{A,r}^2}{M_A}-\frac{e_{B,r}^2}{M_B}\right),
 \qquad s_{\log}=\frac{s(w)}{\sqrt S}.
\]

Subtracting each population functional's constant influence term would not
change these sample variances. Computing each paired influence vector before
its sample variance retains the cross-variant covariance. For example, the
RMSE-difference variance contains
\(-\operatorname{Cov}(e_A^2,e_B^2)/(2R_A R_B S)\). Adding independent
variance estimates would omit that term and be incorrect for this coupling.
The common-subset squared-error sample covariance is also saved.

For \(\rho=R_A/R_B\), report

\[
 \operatorname{MCSE}(\rho)\simeq\rho s_{\log},\qquad
 [\rho e^{-z s_{\log}},\rho e^{z s_{\log}}],
 \quad z=1.959963984540054.
\]

Difference intervals use the corresponding approximate Monte Carlo interval
estimate \(\pm z\operatorname{MCSE}\). These are simulation-summary
uncertainties, not replacements for the original fitted parameter intervals.
Zero observed paired variance in a finite positive-CV sample does not prove
zero population variance or population equality between floor variants.
For S=0, point metrics are unavailable. For S<2, sample Monte Carlo uncertainty
is unavailable. A zero required RMSE makes the square-root/log delta derivatives
singular: a zero numerator with a positive denominator still has ratio zero,
but its delta uncertainty is undefined; a zero denominator has no defined ratio.
No fictitious zero standard error is supplied in these cases. Relative biases
and RMSEs divide by \(|\psi|\), and remain undefined for zero-truth mu.
The absolute mu metrics and the ratio of its two absolute RMSEs remain usable
when their own denominators permit.

### Operational coverage on all planned pairs

For each of the R planned replications, set \(I_{j,r}=1\) only when variant
j has a valid interval that covers truth; an unavailable interval contributes
zero. The paired difference is

\[
 D_r=I_{A,r}-I_{B,r},\qquad
 \widehat\Delta_C=R^{-1}\sum_r D_r,\qquad
 \operatorname{MCSE}(\widehat\Delta_C)=s(D)/\sqrt R.
\]

The output retains all four counts: both cover \(n_{11}\), only A covers
\(n_{10}\), only B covers \(n_{01}\), and neither covers \(n_{00}\).
They sum to R. For R>1, the equivalent paired variance is

\[
 \operatorname{MCSE}^2=
 \frac{n_{10}+n_{01}-(n_{10}-n_{01})^2/R}{R(R-1)}.
\]

It is a sample-variance calculation using the discordant pairs, not the sum of
two independent binomial variances. Each variant's total covered and
interval-valid counts are reported separately. Conditioning coverage on a
common-valid point subset would change the operational estimand and is not
done. The coverage difference's Monte Carlo interval is an unbounded normal
approximation and can extend outside [-1,1]; its endpoints are not clipped.
Coverage MCSE is unavailable for R<2, including a CV=0 singleton. For the
identical CV=0 control, bias, MSE and coverage differences are zero, and their
paired uncertainties are zero when the respective sample size is at least
two. Its RMSE ratio is one when its common RMSE is positive; the RMSE-difference
and ratio delta uncertainties additionally require positive common RMSE and
S>=2. It is never treated as two independent measurements.

The implementation composes the existing global record validator, semantic
summarizer, coverage indicator and pure Monte Carlo helpers. Overflow and
nonfinite arithmetic on otherwise valid quantities block table production;
such paths are not filtered away. The exact 71-column output contract is the
module's literal `COLUMNS` tuple, including directions, eligibility counts,
common metrics, paired uncertainties and coverage contingencies. This table
does not replace the separately retained per-variant summaries or failed-fit
records, and does not itself establish completion of the sensitivity study.
