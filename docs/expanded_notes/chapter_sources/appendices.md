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
