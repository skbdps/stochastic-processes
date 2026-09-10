# Ornstein–Uhlenbeck estimation under irregular sampling

A computational study of how observation spacing affects exact maximum-likelihood,
mean-gap PFML, and Euler estimates of the Ornstein–Uhlenbeck model

$$
dX_t = \theta(\mu-X_t)\,dt + \sigma\,dW_t.
$$

**Current stage: Week 3 experiment harness.** One YAML now generates schedules,
exact OU paths, three estimator fits, raw replication records, metrics and figures.
The project is proceeding independently; professor outreach remains deferred.
See [Week 3 completion and run guide](docs/week3_progress.md) and
[mathematical requirements and primary sources](docs/week3_mathematics.md).
The [Week 2 record](docs/week2_progress.md) and [OU derivations](docs/ou_derivations.md)
remain available. Written notes do not certify the owner's personal theory exercise.

The next research stages will measure finite-sample bias, RMSE and confidence-interval
coverage across a wider spacing grid, with equidistant comparisons. The current
notebook establishes point-estimation behavior at selected settings; it does not
establish coverage certification or an operating-region boundary. The mean-gap
estimator is the PFML estimator attributed to Aït-Sahalia & Mykland (2003), not a
new estimator introduced by this project.

## Run the Week 3 harness

After installing the environment below:

```bash
python -m ou_irregular.runner run --config configs/week3_smoke.yaml --output artifacts/my_week3_run
python -m ou_irregular.runner replot --run-dir artifacts/my_week3_run
```

The supplied smoke grid is 3 coarseness values × 3 CV values × 20 replications ×
3 estimators = **540 fits**. Its 250 observations per path have 249 gaps. Existing
outputs are protected; use a fresh directory or `--overwrite` for a known OU run.
Committed [smoke artifacts](artifacts/week3_smoke/) include estimates, spacing
and optimizer diagnostics, summaries, configurations, source hashes and three PNGs.

The harness checks saved data against the planned design before replotting and
can replay one replication by its recorded cell ID. Its coverage fields are
explicitly uncomputed: Hessians, Wald intervals, certification and preregistration
remain Week 4. The tiny smoke sample cannot establish a research threshold.

## Environment and Week 2 reproduction

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
| `ou_irregular/config.py`, `configs/week3_smoke.yaml` | Validated, serialized experiment definition |
| `ou_irregular/spacing.py`, `ou_irregular/estimators.py` | Gap flooring and strict exact/PFML/Euler fitting with retry diagnostics |
| `ou_irregular/metrics.py`, `ou_irregular/plots.py`, `ou_irregular/runner.py` | Bias/RMSE and MCSE, nominal equidistant comparisons, run/replot/replay |
| `docs/week3_progress.md`, `docs/week3_mathematics.md` | Requirements, completed work, formulas, source checks and limitations |
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

Week 3 now provides the configurable harness and an exploratory 3×3 run. Week 4 handles
Fisher information, Wald intervals, certification and preregistration before any
headline experiments. Successful smoke execution does not certify those later statistical milestones.

Author: Suryansh Kumar.
