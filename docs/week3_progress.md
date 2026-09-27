# Week 3: completed harness, mathematical work, and validation

**Historical record of the 10 September implementation.** The
[27 September validation update](week3_validation_update.md) supersedes the
fitting, point-validity, provenance and plotting behavior described below.
The original smoke artifacts and their stated counts remain preserved;
the corrected same-input run is stored separately. See also the
[versioned plan amendment](plan_amendment_week3_validation.md).

Updated 10 September 2026. The **Week 3 programming exit is complete**: one YAML
runs the exact simulator, all three estimators, saved replication records, metrics
and figures. Mathematical reference notes and source checks are documented in
[week3_mathematics.md](week3_mathematics.md), with explanatory comments and
docstrings in each module. Professor outreach remains deferred.

One reading item remains explicitly qualified: Aït-Sahalia–Mykland 2003/2004
primary texts were checked at the relevant working sections, but Tang–Chen 2009's
full text was inaccessible. Its publisher-supplied abstract and related primary
research were checked; a full theorem-level reading is **not** claimed. The
implemented likelihoods and metrics do not depend on an unverified bias coefficient.

## Requirement-by-requirement completion

| Mathematical requirement | Programming delivered | Status |
| --- | --- | --- |
| Exact OU transition and stationary initialization; conditional likelihood given x0 | Reuses exact simulator in `ou.py`; strict new likelihoods in `estimators.py` | Complete |
| Gamma shape/scale from mean and CV; exponential and regular cases | `spacing.py`: CV=0 regular, CV=1 exponential, other CV Gamma | Complete |
| High-CV floor, clipped fraction, distinction between nominal/expected/realized spacing | Literal max(D, epsilon*mean), analytic expected clipped fraction/mean, raw and realized diagnostics | Complete |
| Exact, mean-gap PFML and Euler objectives in log(theta),mu,log(sigma) | Shared interface, stable variance arithmetic, validated data, 3-start optimization and deterministic 5-start retry | Complete |
| Bias, relative bias, RMSE, relative RMSE and Monte Carlo SE | `metrics.py`, including bias and RMSE MCSE and relative MCSE | Complete |
| Zero truth and failed-fit denominators | Missing relative mu metrics at mu=0; all planned fits retained; error metrics explicitly conditional on success | Complete |
| Equidistant design comparisons | Same nominal mean gap, n and truth; excess bias, RMSE ratios; missing/ambiguous controls checked | Complete |
| Reproducible independent replications | Stable scientific cell ID, per-replication SeedSequence, separate gap/path children, saved keys and replay | Complete |
| One configuration to end-to-end outputs | Strict safe YAML, `config.py`, `runner.py`, `plots.py` | Complete |
| 3×3 smoke grid | Nine cells, 20 replications, all three estimators, 540 fits | Complete |
| AS&M 2003/2004 working-section reading | Primary links, conventions and corrections in mathematics document | Complete |
| Tang–Chen full working-depth reading | Access limit documented; no coefficients copied from memory | Full-text reading outstanding |
| Coverage, Hessian/Wald validation, certification and preregistration | Coverage fields explicitly uncomputed; runner accepts smoke phase only | Week 4, not claimed complete |

The floor fraction in the smoke YAML is **exploratory**, not a frozen
preregistration. Mathematical reference derivations for profiling, floored second
moments and PFML population limits are documented as reference work; they are not
presented as implemented estimator or correction routines.

## Reproduce and inspect

Install the pinned requirements with Python 3.12, then from the repository root:

```bash
python -m pytest -q
python -m ou_irregular.runner run --config configs/week3_smoke.yaml --output artifacts/my_week3_run
python -m ou_irregular.runner replot --run-dir artifacts/my_week3_run
```

For one replication, copy its full `cell_id` from `cells.json` or `replications.csv`:

```bash
python -m ou_irregular.runner replay --run-dir artifacts/my_week3_run --cell-id FULL_CELL_ID --rep 0
```

Replay emits CSV on stdout and does not consume earlier replications. Replotting
uses saved raw estimates, checks their digest and complete planned key set, and
refreshes derived-output hashes and provenance. To regenerate an existing,
identified OU output directory, add `--overwrite` to `run`; unrelated nonempty
directories are rejected. No output directory tree is deleted.

The supplied grid uses theta=1, mu=0, sigma=0.5; theta*nominal mean gap in
{0.1, 0.5, 1}; nominal CV in {0, 1, 2}; **250 observations / 249 transitions**;
20 replications; seed 20260910; epsilon=1e-6 for CV>=1.5. The smoke plotting
interface accepts at most three coarseness values and three sample sizes; these
limits are checked before simulations start. All three estimators and a CV=0
control are required by the configuration schema.

## Output files and conventions

The checked-in example is [artifacts/week3_smoke](../artifacts/week3_smoke/).

| File | Meaning |
| --- | --- |
| `config.yaml`, `resolved_config.json` | Normalized configuration, including truth, grid, seeds, floor and optimizer settings |
| `cells.json` | Every scientific cell and its stable SHA-256 identifier |
| `replications.csv` | One row per cell, replication and estimator, including every failed fit |
| `spacing.csv` | One row per cell/replication: actual gaps' diagnostics, recorded once per simulated path |
| `spacing_summary.csv` | Per-cell means of clipping, mean-gap, CV and span diagnostics |
| `failures.csv` | Failed-estimation rows; a header-only file means none failed |
| `summary.csv` | One row per cell/estimator/parameter: errors, MCSEs, denominators and regular-grid comparisons |
| `run_metadata.json` | Status, planned/recorded counts, environment versions, code hashes, output hashes and runtime |
| `figures/*.png` | Bias curves, RMSE-ratio curves and optimizer-failure rates |

`n` always means observation count in configuration/output. `theta_mean_gap` and
`cv` are nominal design coordinates, not measured values. `realized_span` is
measured from the actual timestamps. `expected_mean_gap_after_floor` identifies
the small expectation shift introduced by flooring; the generator does not
renormalize the gaps.

The estimators use the same simulated observations within each replication.
Different cells, including regular-grid controls, have independent streams.
Adding grid cells or replications does not change existing cell/replication draws.
Cell IDs depend on data-generating inputs, not estimator order or optimizer starts.

The default retry requests 3 initial starts and, only if needed, 5 retry starts
(8 total requested starts because the first 3 are repeated). A successful retry
cannot replace a materially better failed optimum. `attempted_starts` covers both
calls; `n_starts`, `n_finite`, `n_successful` and `selected_start` describe the
selected call. Objective tolerance/gap cover both calls after a retry.

`relative_bias` divides signed bias by the magnitude of truth. At zero truth its
value is missing, as are the relative RMSE/MCSE values. RMSE MCSE uses a delta-method
approximation and is missing with fewer than two valid fits or zero RMSE.
Coverage values remain missing with `ci_status=not_computed_week3`; success rate
is never substituted for coverage.

## Validation evidence

- **141 tests passed** on the completed code in 8.17 seconds.
- The final supplied smoke run completed in **24.232 seconds**: 9 cells, 180
  simulated paths, **540 planned and recorded fits, zero final failures**.
- All 81 parameter summary rows have a matching regular-grid control. Coverage
  is uncomputed throughout; relative mu quantities are missing as defined.
- Regression tests cover closed-form fits, independently calculated Gaussian
  likelihoods, high-CV clipping and expected-floor means, bad input rejection,
  optimizer retries and exceptions, simulation-failure retention, saved-data
  integrity, per-replication replay, and rebuilding figures from raw rows.
- Hand-computable metric tests check error/MCSE denominators and missing-value
  behavior. Two design types (regular and irregular) verify affine state and
  time-unit transformations for all three likelihoods and the exact simulator.
- All three generated figures were visually reviewed. Each states that it is
  exploratory, conditions point metrics on success, and has no computed coverage.
- Week 2's simulator/API and executed notebook remain unchanged. The full suite
  includes the existing Week 2 mathematical and optimizer regression tests.

Times exclude package installation and vary with hardware. Same-environment
replay is tested; numerical identity across all platforms is not promised.
The source/output hashes refer to the actual files used in the recorded run.

## Interpretation changes and Week 4 handoff

The original plan overstated two conclusions. Our regular-grid comparisons are
not AS&M's formal FIML/IOML/PFML decomposition, and a bias expansion does not prove
that exact-MLE coverage falls below 80% at a specific small-span corner. Those
coverage outcomes remain hypotheses to measure. The new mathematics notes and
updated related-work notes state these qualifications.

Next: implement observed-Hessian diagnostics and interval construction, define
invalid-interval/failure policy, run the restricted certification cell, and freeze
preregistration before E1/E2 headline experiments. Revisit nominal versus
expected-span-matched controls and floor sensitivity before freezing the design.
Acquire the Tang–Chen full text before adopting any attributed bias correction.
