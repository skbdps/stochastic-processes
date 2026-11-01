# Reproducing and verifying the published E4 evidence

The full experiment already exists. Begin with the recovered original and
saved-evidence verification; a new full simulation is unnecessary for reading,
auditing, plotting or checking publication. The
[completion record](week8_completion.md) identifies original versus later
evidence, and the [interpretation](week8_results_interpretation.md) links the
saved tables used for numerical claims.

## Full-history checkout and runtime

Use a full-history clone, preserving file bytes. The lightweight handoff
snapshot lacks historical archives and commit objects required by provenance
checks. Use Linux x86_64, CPython 3.12.14 and the pinned lock file for the
recorded E4 environment. The original study used NumPy 2.3.5, SciPy 1.17.0,
pandas 2.2.3 and matplotlib 3.10.8 on an AMD EPYC 7763 host, with AVX2 present
and AVX-512 unavailable. Its invocation retains CPU flags, NumPy runtime and
thread settings. Versions alone do not guarantee exact source-array replay
on another CPU or operating system.

```bash
git clone -c core.autocrlf=false https://github.com/skbdps/stochastic-processes.git
cd stochastic-processes
git fetch --all --tags
python3.12 -m venv /tmp/ou-e4-current
source /tmp/ou-e4-current/bin/activate
python -m pip install -r requirements-lock.txt
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
unset MKL_NUM_THREADS VECLIB_MAXIMUM_THREADS
export MPLCONFIGDIR=/tmp/ou-e4-reproduction-mpl
```

Check that `python --version` is exactly 3.12.14. Do not edit the dependency
lock, frozen source or protocol files to make an old archive verify. The
original scientific commit is
`0ff57919f5d8345f52d8689732ff94448ca3ca1d`. The release adds publication
evidence and documentation while preserving its declared scientific bytes;
the saved manifest checks those bytes. For a scientific-freeze investigation,
make a separate detached worktree at that commit and leave the release
checkout intact.

The handoff's original ZIP has SHA-256
`e9a780206925c4085ff4bb47cf988c7ccec66ce40db038f3fdef82703fe52513`.
When using it, verify that digest before extracting into a fresh external
directory, reject absolute or escaping member paths, and hash every member.
Preserve the ZIP and its untouched extraction. Copy `week8_e4/` together with
the sibling `week8_e4-figure-renderer/` when making a working copy; renderer
state is execution evidence. A new seal changes only the working copy and
must remain distinguishable from the original pending-review artifact.

## Saved readback, independent audit and external figures

From the release checkout and matching environment:

```bash
python -m ou_irregular.week8.runner verify --run-dir artifacts/week8_e4 --require-completed > /tmp/e4-saved-verification.json
python -m ou_irregular.week8.runner audit --run-dir artifacts/week8_e4 > /tmp/e4-repeat-audit.json
python -m ou_irregular.week8.runner replot --run-dir artifacts/week8_e4 --output /tmp/e4-external-replot
```

Use fresh external filenames/directories for another invocation. `verify`
generates neither schedules nor fits. The full read-only `audit` reconstructs
all 14,000 sources, checks exact array/trace hashes and all tables, and replays
42 endpoint fits. `audit --no-fit-replay` still reconstructs every source but
does not supply a full-completion audit. Neither command writes a replacement
saved audit. A discrepancy is retained and investigated rather than relaxed.
Incompatible macOS/ARM or CPU-dispatch results must not be relabeled as exact
original replay.

`replot` consumes saved audited tables, creates all 18 PNG/PDF pairs externally
and retains a separate sibling renderer-state directory. It samples no paths
and refits no model. Compare every output to the original figures and verify
the source archive's complete hash map before and after:

```python
from pathlib import Path
import hashlib
from ou_irregular.week8 import runner

archive = Path("artifacts/week8_e4")
def hashes(root):
    return {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(root.rglob("*")) if p.is_file()}

before = hashes(archive)
runner.verify_run(archive, require_completed=True)
runner.replot(archive, Path("/tmp/e4-external-replot-check"))
assert hashes(archive) == before
original = {name: value for name, value in hashes(archive / "figures").items()
            if name.endswith((".png", ".pdf"))}
recreated = {name: value for name, value in
             hashes(Path("/tmp/e4-external-replot-check/figures")).items()
             if name.endswith((".png", ".pdf"))}
assert len(original) == len(recreated) == 36
assert original == recreated
```

Exact figure bytes require the recorded renderer versions, settings and font
environment. Never overwrite old figures or renderer state to hide drift.
The release's actual external reproduction record is
**[all 36 external PNG/PDF files byte-identical; source archive unchanged](../artifacts/week8_evidence/completion/replot_verification.json)**.

## Actual review and guarded sealing

The released archive retains its actual full-study review. The review identifies
all 18 current PNG filenames exactly once, individual hashes, current
figure-manifest digest, reviewer/time and factual inspection notes. The
record identifies the inspection scope and provenance. Neither a synthetic
unit-test review nor the R=4 review can close the full-study archive.

For a fresh working copy of the original pending-review artifact, inspect its
current images and provide an actual review before invoking the existing
guarded seal:

```bash
python -m ou_irregular.week8.runner verify --run-dir /tmp/e4-working/week8_e4
python -m ou_irregular.week8.runner seal --run-dir /tmp/e4-working/week8_e4 --review /tmp/e4-actual-review.json
python -m ou_irregular.week8.runner verify --run-dir /tmp/e4-working/week8_e4 --require-completed
```

This is a lifecycle operation on a verified copy, not a new experiment. The
runner checks the bound numerical audit, figures and artifact maps before the
completed transition. Later external mutation invalidates the seal.

## Regression and strict historical replay

The pre-full-run gate accounted for 2,061 tests: 2,055 passed and six
design-inapplicable skips. Final release validation is recorded separately:
**[composite gate: 2,061 identities, 2,055 passes, six design-inapplicable skips](../artifacts/week8_evidence/completion/regression_assembly.json)**. Use the exact recorded workflow/commands,
not a single combined pytest process with an arbitrary latest environment.
The successful development workflow is retained at
`validation/week8-runner-20261031` commit
`afcec2b8c82e6f4e1ba7fd58fbde56582a57e1ee`, in
`.github/workflows/week8-runner.yml`.

Its five disjoint groups cover core, current Week 7, Week 8, historical archive
fixtures and unchanged strict Week 3 isolation. The historical fixture
environment uses the same lock with `ipykernel`, `ipython` and `nbformat`
omitted, reproducing the original recorded optional-package presence.
Install that separate environment without editing the committed lock, then
run the release launcher from a clean full-history checkout:

```bash
python3.12 -m venv /tmp/ou-e4-legacy
python - <<'PY'
from pathlib import Path
omitted = {"ipykernel", "ipython", "nbformat"}
lines = [line for line in Path("requirements-lock.txt").read_text().splitlines()
         if line.split("==")[0].lower() not in omitted]
Path("/tmp/ou-e4-legacy-requirements.txt").write_text("\n".join(lines) + "\n")
PY
/tmp/ou-e4-legacy/bin/python -m pip install -r /tmp/ou-e4-legacy-requirements.txt
OU_LEGACY_PYTHON=/tmp/ou-e4-legacy/bin/python python scripts/week8_release_regression.py /tmp/e4-fresh-release-evidence
```

The last destination must not already exist and must be outside the checkout.
The current interpreter is the activated `/tmp/ou-e4-current` environment
installed above. The launcher checks both runtimes against saved metadata,
captures CPU/dependency state, retains collection, per-group logs, outcomes,
JUnit and failures, and reconciles every complete test identity exactly once.
Its small prescribed validation fixtures do not start a new full E4 study.
The strict historical test runs in a fresh interpreter:

```bash
python -m pytest -q -ra tests/week4/test_week3_isolation.py --junitxml=/tmp/e4-strict-replay.xml
```

The other current groups use a real-file, spawn-safe entry point that imports
`ou_irregular.runner` and then restores the recorded unset MKL/VECLIB variables
before calling pytest. This matters because importing the historical runner
has global thread-environment side effects. It changes no assertion or test.
The retained workflow supplies the exact launcher, environment comparisons,
test group boundaries, complete collection/JUnit identity accounting and
failure-log retention. Six floor-only checks on the non-paired follow-up
fixture remain inapplicable; none waives strict replay or a scientific gate.

CPU-dispatch diagnostics showed that disabling AVX-512 could change historical
source hashes on the tested host. The cause of the earlier combined-process
failure remains unproven because its CPU state was not captured. Preserve
CPU/runtime captures and exact mismatches; do not weaken comparisons or change
old manifests. Historical Week 4's recursive-source auditor must run at its
recorded revision; new explicit-file dependency checks permit unrelated later
additions.

## Optional independent new execution

The following is the exact original scientific CLI shape, documented for a
separately intended new reproduction. It is not part of recovering the saved
experiment and must use a fresh destination, committed dependencies, the
published scientific freeze and the recorded environment:

```bash
python -m ou_irregular.week8.runner run --workers 4 --partition-size 250 --output /tmp/e4-new-independent-reproduction
```

It executes the fixed 14,000 paths / 42,000 fits with seed 20261011; it has no
scientific overrides. A new output remains a new execution, with its own audit,
actual figure inspection and seal; never relabel it as the original artifact.
The complete R=4 development validation uses `run --validation` and its
separate fixed seed 21261011. Source-only validation uses seed 22261011 and
supplies no fitted E4 results.
