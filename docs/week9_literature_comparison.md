# Week 9 literature comparison and bounded contribution

Checked 31 October 2026. This is a targeted primary-source review for the
publication evidence package, not an exhaustive novelty certificate. It adds
current qualifications to the historical [prior-art log](prior_art_log.md)
without rewriting that dated record. A fresh search immediately before
submission remains necessary.

## Search record

The following queries were executed through web search on the check date.
Search hits were screened against primary publisher, author-hosted or arXiv
records. Dates shown below come from those records rather than search-engine
publication-age labels.

| Query or route | Scope of examination |
| --- | --- |
| `Ornstein Uhlenbeck irregular sampling confidence interval coverage PFML Euler 2024 2025 2026` | Recent interval-calibration and estimator comparisons |
| `Ornstein Uhlenbeck irregular observations maximum likelihood simulation coverage 2024 2025 2026` | Recent irregular-observation methods and simulations |
| `Ornstein Uhlenbeck sampling coverage 2024 2025 2026` | Broader recent screening; unrelated applications excluded by model/inference target |
| Exact-title searches for Tang--Chen, Yu, Fleming 2017/2019, Kelly, Kozlowski and Burke | Resolve versions, primary access and actual estimands |
| Author-hosted Ait-Sahalia--Mykland reprints | Verify sampling/nuisance conventions and simulation content |
| arXiv version records and HTML for Holý--Tomanová | Check the 19 October 2025 v4 rather than treating the 2018 deposit as current |
| Primary arXiv records of recent screened candidates | Check current titles, submission/version dates and model differences |

The search found a directly relevant October 2026 interval-coverage preprint.
It is included below. No result of this limited search establishes absence of
an overlapping study; incomplete access is retained as a limitation.

## Checked comparison

| Source and checked material | Relevant contribution and applicability limit |
| --- | --- |
| Ait-Sahalia and Mykland (2003), [author-hosted published article](https://www.princeton.edu/~yacine/random.pdf), *Econometrica* 71, 483--549; sampling/likelihood distinctions, OU section, Figures 6--8 and Table V | Establishes FIML/IOML/PFML distinctions and asymptotic bias/variance analysis with finite-sample simulations. PFML attribution belongs here. The OU specialization uses a fixed zero mean and specified nuisance conventions. Our nominal regular contrasts do not implement its formal information decomposition or IOML. Table V reports simulated distributional moments; describing this paper as containing no simulations would be incorrect. |
| Ait-Sahalia and Mykland (2004), [author-hosted published article](https://www.princeton.edu/~yacine/gmmrandom.pdf), *Annals of Statistics* 32, 2186--2222; Assumption 2, Section 4.3/Table 1 and Euler objective | Provides estimating-equation asymptotics for randomly spaced diffusions. Assumption 2 has IID intervals; the OU examples have particular known-nuisance settings and a zero mean. These coefficients cannot be transplanted into jointly estimated theta/mu/sigma with dependent gaps. The present expected-objective derivation supplies a narrower interpretation without claiming the dependent-gap limit theory. |
| Tang and Chen (2009), [publisher record](https://www.sciencedirect.com/science/article/pii/S030440760800208X), DOI [10.1016/j.jeconom.2008.11.001](https://doi.org/10.1016/j.jeconom.2008.11.001) | The publisher abstract describes Vasicek/CIR bias and variance expansions and bootstrap bias correction. Direct publisher full-text access returned HTTP 403 in this review. A 2008 working-paper mirror surfaced in search, but equivalence to the published text was not verified. Full-text working-depth reading and correction coefficients remain uncertified; no correction is adopted from this access-limited review. |
| Yu (2012), [publisher record](https://www.sciencedirect.com/science/article/pii/S030440761200005X), DOI [10.1016/j.jeconom.2012.01.004](https://doi.org/10.1016/j.jeconom.2012.01.004); publisher-supplied abstract/introduction visible in search | Analyzes mean-reversion bias for discretely observed OU with a known long-run mean and discusses finite-span behavior. Direct full-text fetch returned HTTP 403. Known-mean coefficients are not the benchmark for the present jointly estimated mean. A bias-order result does not imply a particular Wald-coverage threshold. |
| Fleming and Calabrese (2017), [publisher full text](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12673), *Methods in Ecology and Evolution* 8, 571--579; simulation methods and Figure 4 | The interval-coverage target is home-range **area**, not the three OU parameters. OU/OUF movement simulations are sampled every three hours with varying duration; interval tests use nominal 50% coverage for the true 95% area. This is relevant small-effective-sample inference, but calling it our three-parameter nominal-95% Wald experiment would conflate targets and levels. |
| Fleming et al. (2019), [publisher full text](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.13270), *Methods in Ecology and Evolution* 10, 1679--1689; Section 3.1 and Figures 1--3 | Compares ML, REML and perturbative/hybrid corrections using two-dimensional OU/OUF movement models and an effective-sample-size perspective. Checked simulation figures use daily/hourly observations and include 99% Monte Carlo bands for mean area or RMSE. Those bands are not themselves empirical coverage of our theta/mu/sigma intervals. The paper also discusses home-range probability coverage; its inferential objects must be named precisely. |
| Kelly, Bechtold and Siemiginowska (2009), [author preprint record](https://arxiv.org/abs/0903.5315), abstract/version record | Uses continuous-time stochastic modeling and Bayesian estimation of quasar variability time scale/amplitude. It supplies application context for irregular light curves. This review checked the abstract/record, not an exhaustive theorem or interval-coverage inventory. It therefore makes no blanket claim that the paper contains no coverage assessment. |
| Kozlowski (2017), [author preprint record](https://arxiv.org/abs/1611.08248), published DOI [10.1051/0004-6361/201629890](https://doi.org/10.1051/0004-6361/201629890); abstract and primary search excerpt | Studies DRW parameter recovery versus survey duration, cadence, noise and other observing features. The reported ten-times-timescale guidance belongs to that fitting/application setup. Publisher HTML/PDF fetches returned HTTP 403; it is not imported as a universal OU interval guarantee or our numerical boundary. |
| Burke et al. (2021), [author preprint record](https://arxiv.org/abs/2108.05389), abstract/version record | Measures an optical variability timescale and its relationship with black-hole mass in an AGN sample. This is physical application and timescale context. Abstract-level screening does not establish that all simulation or coverage questions were absent from the full article or supplement. |
| Holý and Tomanová, [arXiv v4 record](https://arxiv.org/abs/1811.09312v4) and [full HTML](https://arxiv.org/html/1811.09312v4), 19 October 2025; Section 2.4 and likelihood/ARMA definitions | Evaluates noise-sensitive and noise-robust OU estimators under Poisson observation times, with exact simulated transitions. The checked simulation comparison emphasizes mean absolute estimation errors and includes a microstructure-noise model. This differs from our noise-free exact/PFML/Euler coverage grid. Parameter notation also differs: their tau is mean-reversion rate, while our theta is the rate. |
| Li, Kenney and Gu (2026), [arXiv v1](https://arxiv.org/abs/2610.25266v1) and [full HTML](https://arxiv.org/html/2610.25266v1), submitted 18 October 2026; Sections 2--3 and Table 1 | Studies an OU state-space model with measurement error. Table 1 compares asymptotic and profile-likelihood interval coverage for log(theta), using 200 datasets at several sample sizes. Its model supports unequal intervals, while the stated simulation design in Section 2.3 uses unit time steps. This is substantive nearby coverage work; it does not establish the same PFML/Euler by gap-law experiment. General claims that OU finite-sample interval coverage is unstudied are untenable. |

## Additional 2024--2026 screening

These primary records were screened at abstract/version-record depth. They are
not a substitute for full-text comparisons if the eventual manuscript claims
extend into their domains.

| Primary record | Observed scope and disposition |
| --- | --- |
| Ben Alaya, Friesen and Kremer, [arXiv:2404.05554v3](https://arxiv.org/abs/2404.05554v3), revised 27 September 2025 | Ergodic Volterra OU drift inference, including continuous and discrete observation regimes. Different process/limit-theory setting; retained as adjacent theory. |
| Han and Zhang, [arXiv:2403.18255](https://arxiv.org/abs/2403.18255), submitted 27 March 2024 | Multi-regime threshold OU estimation and asymptotic results. Different drift model; retained as adjacent inference. |
| Kawamo and Masuda, [arXiv:2601.19071](https://arxiv.org/abs/2601.19071), submitted 27 January 2026 | Skewed stable Levy-driven OU likelihood/quasi-likelihood theory. Non-Gaussian driving noise and different asymptotic regime; outside the present Gaussian grid. |
| Ning et al., [arXiv:2610.27806v1](https://arxiv.org/abs/2610.27806v1), submitted 20 October 2026 | Recent AGN DRW variability application and physical scaling relations. Retained as current application context; no coverage-grid conclusion is drawn from abstract screening. |

## Contribution statement and claim disposition

The supported description is a reproducible finite-sample comparison of exact,
realized-mean-gap PFML and actual-gap Euler estimation for a stationary,
noise-free OU process with exogenous observations. It reports parameter-specific
error, interval availability and operational/conditional coverage with Monte
Carlo uncertainty across a fixed n/q/CV grid. Exploratory diagnostics retain
boundary, curvature and floor outcomes, and E4 compares dependent gaps with
independent controls sharing the complete one-gap marginal.

This description does not claim a new estimator, first discovery of
random-sampling effects, universal performance boundary, or established
degradation beyond a matched first-order covariance theory. Whether the full
combination is sufficiently distinct for a particular publication requires
academic assessment and further source checking. The checked recent preprint
must appear in manuscript positioning; it cannot be omitted merely because it
postdates the historical search.

The historical kill-criterion verdict is a dated search outcome. This review
does not identify an exact duplicate of the complete three-estimator by gap-law
design among the checked sources, but its scope and access limits do not
justify an exhaustive absence claim. A close overlap found later would require
revising contribution language before submission, rather than generating
outcome-driven experiments to manufacture novelty.
