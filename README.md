# Ornstein–Uhlenbeck estimation under irregular sampling

A computational study of how observation spacing affects exact maximum-likelihood,
mean-gap PFML, and Euler estimates of the Ornstein–Uhlenbeck model

$$
dX_t = \theta(\mu-X_t)\,dt + \sigma\,dW_t.
$$

**Current stage: Week 2 computational milestone.** The written OU derivations,
simulator, estimation demonstration, tests and reproducible results are available.
The project is proceeding independently; professor outreach is deferred.
See [progress and validation](docs/week2_progress.md) for the revised scope and
[OU derivations](docs/ou_derivations.md) for the mathematical reference. The owner's
optional closed-notes theory exercise is not certified by the code or these notes.

The next research stages will measure finite-sample bias, RMSE and confidence-interval
coverage across a wider spacing grid, with equidistant comparisons. The current
notebook establishes point-estimation behavior at selected settings; it does not
establish coverage certification or an operating-region boundary. The mean-gap
estimator is the PFML estimator attributed to Aït-Sahalia & Mykland (2003), not a
new estimator introduced by this project.

## Reproduce the Week 2 milestone

Use Python 3.12 and run from a checkout of this repository:

```bash
git clone https://github.com/skbdps/stochastic-processes.git
cd stochastic-processes
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python run_all.py
```

The repository is private, so cloning requires access to the owner's GitHub account.
On Windows activate with `.venv\Scripts\activate` instead. The direct dependencies
are pinned in `requirements.txt`; `requirements-lock.txt` records their complete
dependency closure for the validated Linux/Python 3.12 environment.

`run_all.py` runs pytest and executes every notebook code cell with the same Python
interpreter, regenerating the executed notebook, all four replication CSVs and the
teaser PNG. It uses an in-process IPython shell with captured text and image outputs;
no separate Jupyter server or local TCP connection is required. This runner targets
this notebook's synchronous Python cells. You can also open the notebook in a
notebook editor using a kernel from the same environment.

The seed is `20260902`, with independent child streams for each experiment. Every
reported fit must be finite and converged; a failed fit stops the run. Seeds and
version pins support reproducibility, but floating-point details and optimizer
statuses can still vary across platforms. Execution time excludes dependency
installation; recorded measurements are in the progress document.

## Contents

| Path | Purpose |
| --- | --- |
| `ou_irregular/ou.py` | Exact simulation; equidistant, exponential and Gamma gaps; exact/PFML/Euler likelihoods; multistart fitting; analytic checks |
| `notebooks/01_exact_mle_demo.ipynb` | Four Week 2 experiments and a clearly marked Euler preview, with executed outputs |
| `tests/test_ou.py` | Mathematical checks and optimizer regression tests |
| `results/*.csv` | One row per replication, including optimizer diagnostics |
| `figures/teaser.png` | Exact-versus-PFML demonstration |
| `docs/ou_derivations.md` | Wiener integral, OU solution, transition law, stationarity, likelihoods and optional exercises |
| `docs/week2_progress.md` | Current scope, validation evidence and next steps |
| `docs/prior_art_log.md`, `docs/related_work_notes.md` | Historical prior-art notes supplied with the project |
| `OU_master_plan_v2.pdf` | Original master plan; the progress document records the outreach deferral |
| `stochastic_processes_for_OU_notes.html` | Original learning notes through P2, supplemented by the derivations document |

## Teaser: exact MLE versus PFML

![Exact and PFML sampling distributions](figures/teaser.png)

Exponential observation gaps, CV = 1, mean gap = 0.5, 1000 observations and 500
replications; true theta = 1 and sigma = 0.5. PFML uses the realized average gap.
The exact estimator uses each timestamp and has small bias at this setting.

| Quantity | Exact MLE mean | PFML mean | PFML large-sample limit |
| --- | ---: | ---: | ---: |
| theta | 1.0095 | 0.8183 | 0.8109 |
| sigma | 0.5013 | 0.4515 | 0.4503 |

The exact-MLE relative biases are about +0.95% and +0.26%; PFML's are about
-18.17% and -9.71%, respectively. This result does not imply that exact MLE is
unaffected by irregular sampling: finite-sample bias and uncertainty may depend
on spacing. All table values regenerate from `results/cell4_teaser_naive_vs_exact.csv`.

## Next stage

Week 3 builds the configurable experiment harness with per-replication seed/config
records, high-CV gap-floor logging, three estimators and plotting. Week 4 handles
Fisher information, Wald intervals, certification and preregistration before any
headline experiments. The Gamma helper and Euler preview here are prototypes,
not evidence that those later milestones are finished.

Author: Suryansh Kumar.
