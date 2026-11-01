# Reproducing the Week 10 explanatory supplement

This supplement is separate from original E4 and the Week 9 figure freeze.
Read [the protocol](week10_supplement_protocol.md),
[mechanism derivation](week10_mechanism.md) and
[span-standardization method](week10_span_method.md) first. Saved results live
under `artifacts/week10_supplement/`; normal verification does not rerun a
simulation. Original E4 remains under `artifacts/week8_e4/`.

## Recorded environment

The new study uses CPython 3.12.10 on macOS ARM64, NumPy 2.3.5, SciPy 1.17.0,
pandas 2.2.3, Matplotlib 3.10.8 and PyYAML 6.0.3. Install the exact package
closure from [requirements-week10-lock.txt](../requirements-week10-lock.txt)
in a separate Python 3.12 environment. Its environment, source hashes,
uncommitted-file status, CPU information and numerical runtime are retained
in the package evidence. The scientific code was frozen by its actual bytes
before the full run; a later Git commit is not described as its execution
commit. This is a new experiment's native environment, not a claim to replay
the original Linux E4 bit for bit on macOS.

Use the following environment prefix for every command below, substituting
your Week 10 interpreter path for `python`:

```sh
env -u MKL_NUM_THREADS -u VECLIB_MAXIMUM_THREADS \
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLCONFIGDIR=/tmp/ou-week10-mpl \
  python -m ou_irregular.week10.cli verify \
  --package artifacts/week10_supplement
```

This checks the entire completed package inventory, rederives the fitted-
parameter/source/mechanism summary tables, and rederives all six span-analysis
exports from original saved E4 records. It checks source bytes and the recorded
Python/package/thread environment. Exact replay across different hardware or
numerical libraries is not promised. A mismatch is a diagnostic to retain,
not permission to rewrite the sealed evidence or weaken a comparison.

## Reconstructing sources and selected fits

```sh
python -m ou_irregular.week10.cli audit \
  --run-dir artifacts/week10_supplement/horizon_full
```

With the same environment prefix, this reconstructs all 14,000 schedules and
OU paths, verifies their hashes and scalar source diagnostics, independently
checks PFML means by ordinary least squares, and replays all three original
fits at replications 0 and 1999 in each cell. It reuses the existing successful
audit record only if the recalculation agrees; it never refreshes its timestamp
to claim a new original review. No full simulation is needed for this check.

The mechanism's per-schedule predictions can be rederived separately:

```python
from ou_irregular.week10.design import full_design
from ou_irregular.week10.reporting import verify_tables
verify_tables(full_design(), "artifacts/week10_supplement/horizon_full",
              "artifacts/week10_supplement/analysis", "artifacts/week8_e4",
              replay_mechanism=True)
```

This generates the same saved-source realizations to evaluate covariance
diagnostics; it performs no new fits and does not add research replications.

## Regenerating tables and figures outside the completed package

Always use fresh output locations. The writers refuse an existing destination.

```sh
python -m ou_irregular.week10.cli span --output /tmp/week10-rebuild/span
python -m ou_irregular.week10.cli analyze \
  --run-dir artifacts/week10_supplement/horizon_full \
  --output /tmp/week10-rebuild/analysis
```

Then render the tables with
`ou_irregular.week10.plots.render("/tmp/week10-rebuild")`. Compare file digests
with the original tables and PNG/PDF inventory. PDF creation timestamps are
suppressed. Rendering depends on the recorded graphics environment.

## Running an explicitly requested new execution

Do this only when a new simulation is actually wanted; it is not a prerequisite
to reading or verifying the saved results. The fixed seed, full seven-cell
design and 2,000 replications per cell are defined by the protocol. New output
paths keep the completed package immutable.

```sh
python -m ou_irregular.week10.cli freeze --output /tmp/week10-new-freeze.json
python -m ou_irregular.week10.cli run --run-dir /tmp/week10-new-run \
  --freeze-file /tmp/week10-new-freeze.json --workers 4
python -m ou_irregular.week10.cli audit --run-dir /tmp/week10-new-run
```

The R=4 development run uses `--validation` with its separate seed. Resume is
explicit, verifies identical provenance and planned records, and is refused
after audit. Pipeline errors block successful audit. No rejection sampling or
replacement seed is used. A possible schedule with fewer than three observations
is rejected by the existing Sample contract and retained as a source failure;
the run would require diagnosis rather than being declared complete. The full
study's observation counts and absence of source failures are reported in its
source table and audit.

## Seal and documentation boundaries

The completion seal requires successful tests, the complete 14,000-source /
42,000-fit audit, fresh table and mechanism checks, original span readback,
current validation-file bindings and inspection records for every PNG/PDF pair.
The final validation record binds the reviewed scientific payload before seal.
Later manuscript revisions may cite these data; they must not rewrite the
sealed package or transfer an old figure review to changed pixels.

The original Draft 1 validation and subsequent wording-only revision record
remain historical. Draft 2 adds the supplement and has its own recorded hashes
and technical checks. Draft 3 corrects analytical interpretation outside that
seal; the [critical review](week10_analytical_review.md) governs its moment and
coverage qualifications. Academic review and the complete Week 11 manuscript
remain separate work.

## Reproducing the additive analytical review

With the same pinned environment prefix, run:

```sh
python -c 'from ou_irregular.week10.analytical_review import verify; print(verify("artifacts/week10_supplement", "artifacts/week10_analytical_review/calibration"))'
python -m pytest tests/week10 -q
```

The first command rederives four tables, verifies all 250 original sealed
files and reads all 42,000 fit records to compare every saved PFML center,
standard error and coverage indicator. It performs no source simulation or
model fitting. To rebuild the same tables in a new external directory, use
`build("artifacts/week10_supplement", "/tmp/week10-new-analytical-review")`
from the same module. An existing destination is refused.

The tests include deterministic boundary-pole examples, the observed-Hessian
mean identity, exact conditional mixture calculations, paired coverage
uncertainty and guard checks. They retain the original Week 10 test coverage.
The first-pass diagnostic output and its source are preserved separately in
the review package; a subsequent correction adds conservative paired bounds
for rare or zero discordances. Current `calibration/` is the reviewed output.
The original supplement and Draft 2 remain byte-identical.
