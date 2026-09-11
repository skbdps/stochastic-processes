# Revision notes: review and validation

Reviewed on 11 September 2026.

## Read the notes

Open [`../index.html`](../index.html) from a checkout, or open
[`../stochastic_processes_revision.html`](../stochastic_processes_revision.html)
directly in a browser. The original learning notes are preserved unchanged.
The standalone HTML delivered in the conversation has identical mathematics;
its supplementary repository links use GitHub URLs instead of local paths.

Sections 01–10 revise the original foundations through SDE notation.
Sections 11–13 are explicitly optional continuation into Wiener integration,
OU transitions and irregular-grid likelihoods. There are 20 expandable
self-checks. This documentation does not certify personal learning milestones
or change the computational research results.

## Mathematical/source review

The explanations and worked derivations were reviewed against the linked primary
sources and by direct calculation. The source map in the HTML points to:

- MIT 18.S096 Lecture 5, sections 2–3: random walks, first-step analysis,
  the column-oriented transition matrix and the positive-entry convergence theorem.
- Lecture 8, slides 17–18 and 28 onward: AR(1), stationary moments and likelihoods.
- Lecture 17, sections 1–2: path space, Brownian increments, the scaling construction,
  differentiability and quadratic variation.
- Lecture 18: stochastic integration and the Itô isometry.
- Lecture 21, section 1.2: the integrating-factor OU solution.
- Särkkä & Solin (2019), section 4.1, Example 4.5 (printed p. 50 / PDF p. 58),
  Example 6.2 (printed p. 81 / PDF p. 89), and section 6.5.

Two printed-source corrections were checked: Lecture 8 slide 17 omits the square
on phi in the stationary-variance denominator; Lecture 17 p. 2 omits the spatial
1/sqrt(n) normalization in its walk construction. The companion derives the
correct formulas rather than copying those displays. Some PDF rendering requests
were unavailable; their parsed source text was used alongside the successfully
rendered source pages. No third-party PDFs or copied figures are bundled.

Other explicit qualifications include stationary initialization versus convergence,
fixed-time CLT versus path convergence, grid versus off-grid interpolation,
conditional versus full likelihood, admissible AR-to-OU parameter conversion,
deterministic-integrand Gaussianity, and independence of the observation grid
from the path. The profiled likelihood and exponential-gap Jensen illustration
are derived with their assumptions stated. This is a mathematical/editorial review,
not a peer review or a formal proof-verification system.

## Executed checks

`python docs/check_revision_examples.py` passed **20 test methods** locally.
They cover worked arithmetic, gambler's-ruin recursions, Markov examples,
AR(1) moments and fitting, Brownian overlap, the OU display table,
transition composition, stationarity, small gaps, parameter conversion,
conditional profiling and invalid/degenerate inputs. These tests validate the
examples and identities, not statistical bias or interval coverage claims.

`CHROMIUM_EXECUTABLE=/usr/bin/chromium python docs/check_revision_browser.py`
passed locally in Chromium 144.0.7559.96. Checks covered:

- All 20 self-check disclosures, unique IDs, internal anchors and local link targets.
- Expand/collapse controls and native keyboard slider operation at four gap values.
- No page-level horizontal overflow at widths 320, 375, 768, 1024 and 1440 pixels.
- Print expansion/restoration and an A4 PDF render.
- No-JavaScript reading and native disclosure behavior; no captured JavaScript errors.

Screenshots and a print render were generated. Desktop and mobile screenshots
were inspected. This is not a full accessibility certification or an exhaustive
cross-browser/device audit. Long equations and tables may scroll within their
own containers on narrow screens.

The repository workflow repeats the mathematical and Chromium tests, and retains
source, JSON report, screenshots and the print render as a private Actions artifact
for 14 days. The first CI attempt failed because its test treated a range control
as a text input; the test was corrected to use keyboard interaction. Check the
latest run rather than that superseded attempt.

Reviewed repository HTML SHA-256:
`ae3bfea689bf6b134b43e611a3ed8c3deabb3f370008e35974e82380c52d6998`

## Publication boundary

The repository is private and reported `has_pages: false` at review time.
Merging these files into `main` publishes them **inside the existing repository**;
it does not create a publicly hosted GitHub Pages site or change repository
visibility. The HTML itself is self-contained and works without a web server.
External references require network access. Supplementary private-repository links
also require the reader to have GitHub access.
