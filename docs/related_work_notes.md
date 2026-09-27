# Related work — "does / doesn't" notes (Tier A)

**Aït-Sahalia & Mykland (2003).** *Does:* defines three likelihood estimators under random sampling — FIML (uses states
and intervals), IOML (integrates the intervals out), PFML (pretends the intervals are fixed at their mean) — and derives
the asymptotic cost of discreteness, the cost of randomness, and the cost of ignoring randomness, with simulations for the
OU example. *Doesn't:* report finite-sample coverage, vary the gap coefficient of variation, or give a practitioner-facing
boundary. **Relation to us:** our naive estimator *is* their PFML; our nominal equidistant comparisons are inspired by their information-loss distinctions,
not identical to their FIML/IOML/PFML decomposition. The proposed coverage study would examine finite-sample behavior.
See [Week 3 source verification and conventions](week3_mathematics.md#2-primary-literature-reading-notes-and-limits).

**Aït-Sahalia & Mykland (2004).** *Does:* a general estimating-equation theory under random sampling via a generalized
infinitesimal generator, with exact MLE and Euler schemes as examples. *Doesn't:* any coverage or simulation-grid study.

**Tang & Chen (2009); Yu (2012).** The project's earlier notes associated these works with
finite-span mean-reversion bias. Week 3 verified Tang–Chen's publisher-supplied abstract, but its full text
was inaccessible, so its exact assumptions and coefficients are not certified here. Finite-span effects motivate
future coverage checks; they do not prove a specific coverage breach. Nominal equidistant controls do not
equalize every realized span, and high-CV flooring also changes expected span. See the
[access limit and corrected interpretation](week3_mathematics.md#tang-and-chen-2009-verification-status).

**Fleming et al. (2017, 2019) — ctmm.** *Do:* compare OU/OUF estimators (ML, REML, pHREML) and check empirical 95% CI
coverage via simulation. *Don't:* vary the spacing distribution — sampling is regular, only the record length varies.

**Holý & Tomanová (2018/2025).** *Do:* OU estimation with irregular (Poisson) spacing under microstructure noise, via a
noisy-ARMA(1,1) representation; bias/MSE and trading profit. *Don't:* coverage, a CV dial, or a naive baseline.

**Astronomy DRW literature (Kozłowski 2017; Kelly 2009; Burke et al. 2021).** *Does:* Monte Carlo over cadence and
baseline for OU parameters; establishes the "baseline ≥ 10 τ" rule of thumb. *Doesn't:* compare estimators or report
coverage.
