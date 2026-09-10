# Week 2 progress and validation

Historical Week 2 snapshot. For the current harness milestone, see
[Week 3 progress](week3_progress.md).

## Current scope

Updated 10 September 2026. The owner has chosen to proceed independently and defer
the outreach proposal and professor email until later. This decision supersedes
the outreach dependency in the original master plan; the historical PDF is retained.

The completed milestone is the **Week 2 computational and written-reference
package**. Writing derivation notes does not establish the owner's personal mastery.
The optional blank-paper exercise remains unverified and does not block the
independent code work or this publication.

| Item | Status | Evidence |
| --- | --- | --- |
| Prior-art notes | Present, historical | `prior_art_log.md`, `related_work_notes.md`; no new literature scan is claimed |
| Written OU derivations | Complete | `ou_derivations.md`, including assumptions and exact/PFML/Euler likelihoods |
| Personal theory exit exercise | Unverified | Optional self-check in the derivations document; no completion is claimed |
| Simulator and four core notebook experiments | Complete | Executed notebook with assertions |
| Euler demonstration | Preview only | Cell 5, outside the Week 2 exit requirements |
| Optimizer status handling | Corrected | Converged finite starts preferred only within an explicit numerical objective tolerance |
| Reproduction and test evidence | Recorded below | Same-interpreter runner, version pins, regenerated CSVs and PNG |
| Professor proposal/email | Deferred | Owner decision; no communication sent |
| Week 3 harness and Week 4 certification | Not started | Future work; no coverage or certification claim |

## Changes from the supplied ZIP

- `fit()` preserves optimizer diagnostics, validates and honors the requested number
  of starts, and raises an error if every candidate is nonfinite or invalid.
- A converged candidate is selected over a nominally better failed candidate only
  when the negative log-likelihood difference is no greater than
  `1e-8 + 1e-10 * abs(best_nll)`. A materially better failed candidate remains failed.
- The notebook asserts finite, successful fitting for every reported replication,
  stores diagnostics separately for the exact/PFML pair, and never silently drops
  failed replications.
- Scientific statements are limited to the configurations actually tested. Removed
  claims that exact MLE is indifferent to irregularity or that the demonstration
  establishes an isolated order-1/T law.
- `run_all.py` executes the notebook with the active interpreter in an in-process
  IPython shell. It captures rich outputs, avoids a separate TCP kernel and returns
  a failure when tests or a notebook cell fail.
- Dependencies are pinned and the original plan and HTML notes retain their paths.
  The new derivation reference supplements the historical learning notes.

## Validation

Validated on Linux with Python 3.12.14. Direct dependencies and their resolved
versions are in `requirements.txt` and `requirements-lock.txt`.

- All **21 pytest tests passed**, including optimizer failure/tie-selection cases.
- `python run_all.py` passed in **35.0 seconds** in the working checkout.
- The same command passed in **34.0 seconds** in a fresh local Git clone of the
  publication files, using the installed pinned environment. This verifies clean
  checkout execution; dependency download/install time is excluded.
- All six notebook code cells executed, including the optional Euler preview;
  no error outputs remain, and the teaser is embedded in the executed notebook.
- All **1600 Monte Carlo fits** in the four result files succeeded: 200 equidistant
  exact fits, 200 exponential-gap exact fits, 500 exact/PFML pairs and 200 Euler fits.
  The separate equidistant single-fit cross-check also succeeded.
- All four CSV files and the exported PNG were byte-identical between the two runs.
  Notebook timing text may differ; cross-platform bitwise identity is not promised.
- The existing master plan and original theory HTML are unchanged.

The regenerated teaser means are theta: exact 1.0095, PFML 0.8183; sigma: exact
0.5013, PFML 0.4515. These remain consistent with the original numerical finding,
with explicit convergence handling in the revised pipeline.

## Interpretation and limits

The demonstration uses a correct conditional likelihood given the initial value,
although simulation initializes from the stationary distribution. Conditional and
stationary-initial likelihoods are different finite-sample estimators; this is
intentional and documented in the derivation reference.

For random observation times, the gaps are generated independently of the OU path.
This does not cover state-dependent/event-triggered sampling. The experiment seed
is fixed, with one child stream per cell. Per-replication seed/config records,
high-CV clipping policy, input hardening, twin-normalized metrics and interval
coverage belong to the next harness stages.

For mu = 0, relative bias is undefined, so the notebook reports its mean and Monte
Carlo standard error rather than dividing by zero. The pilot sample sizes of
200/500 replications do not satisfy the later >=2000-replication headline design.

## Next steps

1. Build `spacing.py`, `estimators.py`, `metrics.py`, `runner.py` and `plots.py` around
   a serialized experiment configuration and independent replication seeds.
2. Specify and report the high-CV gap floor and clipping fraction.
3. Add a small end-to-end smoke grid using exact, PFML and Euler estimators.
4. Add uncertainty estimation, certification and frozen preregistration in Week 4,
   before headline runs. Keep outreach deferred until the owner chooses to resume it.
