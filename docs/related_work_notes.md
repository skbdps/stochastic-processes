# Related work — "does / doesn't" notes (Tier A)

**Aït-Sahalia & Mykland (2003).** *Does:* defines three likelihood estimators under random sampling — FIML (uses states
and intervals), IOML (integrates the intervals out), PFML (pretends the intervals are fixed at their mean) — and derives
the asymptotic cost of discreteness, the cost of randomness, and the cost of ignoring randomness, with simulations for the
OU example. *Doesn't:* report finite-sample coverage, vary the gap coefficient of variation, or give a practitioner-facing
boundary. **Relation to us:** our naive estimator *is* their PFML; our twin normalization is the empirical mirror of their
cost decomposition; our coverage map is the finite-sample counterpart of their asymptotic PFML bias result.

**Aït-Sahalia & Mykland (2004).** *Does:* a general estimating-equation theory under random sampling via a generalized
infinitesimal generator, with exact MLE and Euler schemes as examples. *Doesn't:* any coverage or simulation-grid study.

**Tang & Chen (2009); Yu (2012).** *Do:* characterize the finite-span (order 1/T) bias of the mean-reversion MLE and its
curvature for slow reversion. *Don't:* treat irregular sampling. **Relation:** explains why even the exact MLE breaches
coverage at low-span grid corners; the equidistant twin absorbs this effect by construction.

**Fleming et al. (2017, 2019) — ctmm.** *Do:* compare OU/OUF estimators (ML, REML, pHREML) and check empirical 95% CI
coverage via simulation. *Don't:* vary the spacing distribution — sampling is regular, only the record length varies.

**Holý & Tomanová (2018/2025).** *Do:* OU estimation with irregular (Poisson) spacing under microstructure noise, via a
noisy-ARMA(1,1) representation; bias/MSE and trading profit. *Don't:* coverage, a CV dial, or a naive baseline.

**Astronomy DRW literature (Kozłowski 2017; Kelly 2009; Burke et al. 2021).** *Does:* Monte Carlo over cadence and
baseline for OU parameters; establishes the "baseline ≥ 10 τ" rule of thumb. *Doesn't:* compare estimators or report
coverage.
