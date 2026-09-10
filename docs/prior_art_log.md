# Prior-art log (kill-criterion search, Aug 22 2026)

**Pre-registered kill criterion.** The project pivots if a published paper compares ≥2 of {exact MLE, naive/PFML,
Euler contrast} across ≥2 spacing distributions **with CI coverage** as a reported metric.

**Verdict: does not fire.** No paper satisfies all three axes ([E] ≥2 of our estimators, [S] spacing distribution varied,
[C] empirical CI coverage). Two papers satisfy two axes each; their intersection — coverage across gap-variability with a
mean-gap-naive baseline — is open.

| Paper | Axes | Does | Does not |
|---|---|---|---|
| Aït-Sahalia & Mykland (2003), *Econometrica* 71:483–549 | E + S | FIML / IOML / **PFML** (our naive estimator, verbatim: Δ̄ in place of the actual Δ); fixed vs random sampling; asymptotic bias/variance expansions with confirmatory simulations | finite-sample CI coverage; a CV dial; an operating boundary |
| Aït-Sahalia & Mykland (2004), *Ann. Statist.* 32:2186–2222 | E(partial) + S | general theory for randomly spaced observations; exact MLE and Euler contrasts as worked examples | coverage (purely asymptotic) |
| Fleming et al. (2017) *MEE* 8:571–579; Fleming et al. (2019) *MEE* 10:1679–1689 (ctmm) | E + C | multiple OU estimators (ML/REML/pHREML) with empirical 95% CI coverage | irregular spacing — only duration/rate is varied |
| Holý & Tomanová (2018 arXiv:1811.09312 → 2025 *Ann. Oper. Res.*) | E(partial) | noise-robust OU estimation for ultra-high-frequency data; Poisson spacing | coverage; CV dial; naive baseline (their axis is microstructure noise) |
| Kozłowski (2017) *A&A* 597:A128 and the DRW literature (Kelly 2009; Burke et al. 2021) | S | Monte Carlo over cadence and baseline for OU = damped random walk; baseline ≥ 10τ rule | coverage; estimator comparison |

**Supporting results used in the design.** Tang & Chen (2009, *J. Econometrics* 149:65–81): the mean-reversion MLE has
bias of order 1/T (T = data span), not 1/n. Yu (2012, *J. Econometrics* 169:114–122): refinement for slow mean reversion.
Florens-Zmirou (1989): the Euler contrast. Kessler (1997): higher-order Gaussian contrasts.

**Commitment.** Re-run a targeted 2024–2026 preprint scan immediately before submission.
