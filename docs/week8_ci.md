# Validation commands and archived CI incidents

The Week 8 release adds no active GitHub CI workflow to main. Numerical audits, exact replay assertions, actual figure review, guarded seals and complete accounting remain required research checks. Their reproducible commands and saved outcomes are retained; a CI badge is not their substitute.

The original frozen source passed the pre-full-run gate. The final fresh native invocation 41215391840 retains 2,054 passes, six design-inapplicable skips and one strict historical input-hash failure. Four complete groups (2,057 cases) passed on that host. The entire four-case historical group subsequently passed unchanged under the recorded Intel SDE instruction runtime in run 41219755639. A separate composite record reconciles every one of the 2,061 identities, with 2,055 passes and six skips across the selected complete groups. It does not relabel the original failed invocation.

The retained assembly binds both complete original artifacts, commands, package/thread settings, CPU or effective dispatch, all node IDs, JUnit, console output and source control commits. The original AVX2-only host failure, two unavailable-hardware probes and unsuccessful emulator-launch attempts remain distinct from the successful strict replay. No source, assertion, seed or tolerance was changed to obtain a preferred scientific result. Historical byte replay is environment-specific; matching version numbers alone does not guarantee it.

For saved readback from a clean full-history clone using Linux/x86, CPython 3.12.14, the package lock and recorded thread settings:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
unset MKL_NUM_THREADS VECLIB_MAXIMUM_THREADS
python -m pip install -r requirements-lock.txt
python -m scripts.week8_release_verify --output /tmp/fresh-week8-readback.json
```

For another deliberately intended full regression, prepare the historical optional-package interpreter as described in [reproduction](week8_reproduction.md), inspect CPU compatibility, and run:

```sh
OU_LEGACY_PYTHON=/path/to/historical/bin/python \
  python scripts/week8_release_regression.py /tmp/fresh-week8-regression
```

The output must be fresh and outside the checkout. This runs all unchanged tests, including prescribed small development fixtures, and retains every outcome. It starts no full E4 experiment. A historical hash mismatch must remain a failure; do not substitute a skip, weaken its comparison, or retry random hosts to obtain a badge. The successful emulator's precise launcher, official-kit digest, runtime smoke and unchanged test invocation are retained with the composite evidence. No emulator is selected automatically by the scientific code.

Earlier CI and diagnostic branches are preserved as historical provenance. Full raw GitHub job-server logs remain in the original local evidence, with focused public failure records carrying hashes and source links. Scientific CLI logs and the necessary test/runtime records are retained with the release.
