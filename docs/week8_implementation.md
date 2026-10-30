# Week 8 implementation checkpoints

The fixed exploratory protocol and gate were published at
`ff6872a0eed3220eb200ca39baf0b80e93efbd40` before implementation. Modules are
implemented, reviewed, tested and published in sequence. None of these module
checkpoints is a completed E4 scientific run.

## 1. Design and scientific identities

`ou_irregular/week8/design.py` defines the seven ordered cells, full R=2000
and separate fitted-validation R=4 plans, namespace 8, explicit independent
reference maps, and source/floor policy. `E4Cell` derives calibration and CV
from mechanism and ratio. Its complete flat scientific mapping is hashed;
root seed and replication count belong to the run identity. The strict root
parser accepts only the complete adopted plans. Internal CellStore partition
views remain possible but cannot pass root validation.

The focused suite passed 73 tests. An independent read-only mathematical
review passed 170 deterministic assertions using rational 2-by-2 resolvents,
independent hashes, and rehashed tampering of derived values. Root reviewed
the source and tests. These checks generated no schedules, OU observations
or fits and do not certify the later source generator or runner.

The first test invocation had one faulty tamper fixture: it replaced the
r=5 calm event probability with its existing value 0.5, so rejection was
incorrectly expected. The fixture was corrected to 0.4, and additional strict
type cases were added. Both test outputs are retained in
`artifacts/week8_development/design/`. This was a test-fixture correction,
not a change prompted by scientific outcomes.

Reviewed design source SHA-256:
`b3fbf27b4680a92a2c78e9647a8eb15eabdffaa175df5434a5e70424392f1ac7`.

At that checkpoint, remaining modules were source generation and analytical checks, guarded execution,
reporting, independent audit, plots, and runner/integration. The complete
source-only validation and fitted R=4 validation precede a reviewed full-run
freeze. No full E4 observations have been generated.

## 2. Source generator and fixed analytical-validation harness

The design checkpoint was published at
`6ce1a2ef2bac91ca12c234bfee892a131992364f` before this module began.
`E4Source` implements event-stationary MMPP gaps, independent full-marginal
controls and the regular reference. Its gap-only API draws no OU observations.
Its sample API uses the unchanged exact OU simulator with an independent
stationary initial value and innovation stream. Raw, transformed and stored-time
diagnostics are retained separately. See [the source contract](week8_source.md).

The first source test run passed 39 tests. Independent mathematical/source
review passed 238 deterministic assertions using rational moments, an 80-digit
Decimal floor reference, constructed race traces and pure namespace-8 addresses.
It allocated no real RNG generators and generated no study sources or fits.
Root reviewed the source, harness and tests.

Independent harness tests initially passed 22 cases and failed 15, exposing
two implementation gaps: incomplete schema/identity/count rejection and loss
of the terminal failed report if check construction itself raised. Both were
repaired before stochastic validation. Failed check construction now retains
all 1,536 planned sequence rows and 108 unavailable failed checks. The source
generator was unchanged. The repaired combined design/source/harness suite
passed **149 tests**; the independent harness suite separately passed **37**.
These counts overlap and must not be added as independent tests.

The initial harness source, tests and complete failing output are preserved in
`artifacts/week8_development/source/`, alongside final logs and review hashes.
Two initial test invocations inadvertently shared one log pathname; the retained
complete log is explicitly one output, not evidence for two separate invocations.
Both used the same source version and literal/mocked data. No research outputs
or stochastic-validation results were involved.

Reviewed source SHA-256:
`73adae7cf326984e43916a80993ac7cce541a0d406264eb0482a8c2849338ce5`.
Reviewed repaired harness SHA-256:
`9548f7c410a91d88374df86ace56b5e98759dbf28bcea50d3921a4b780bc46c5`.

The following fixed source-only validation is the next gate, after this
implementation checkpoint is published. It requires frozen dependencies and a
fresh output directory and exposes no seed, R, tolerance or overwrite override:

```bash
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLCONFIGDIR=/tmp/ou_recovery_mpl \
  .venv/bin/python -m scripts.week8_validate_source \
  --output artifacts/week8_source_validation
```

The later fitted R=4 validation and full E4 study remain pending. No full-seed
E4 observations have been generated.

### Observed source-only validation

The reviewed source/harness checkpoint was published as
`165c887bbbe5fb6a9ae2593a2e2f6766b7aea6bc` before execution. The exact command
above completed with exit code 0 in 8.53 seconds, with empty stderr.
`artifacts/week8_source_validation/validation.json` reports **passed**:
all 1,536 sequences and 108 fixed checks are retained, with zero source failures,
zero failed checks, zero OU observations and zero fits. Root checked complete
sequence keys, all check outcomes, the recorded freeze, and every declared
artifact digest. Native command/output records are preserved in
`artifacts/week8_development/source/validation_invocation/`.

The largest discrepancy was 2.280454 standard errors (MMPP r=1 raw adjacent
product), within the prespecified six-SE rule. At r=20, one MMPP gap and two
IID gaps were clipped among 255,744 gaps per cell. The other four cells have
no floor. These rare-event counts are descriptive and were not subjected to
an empirical normal acceptance interval. All source-only check results are
development evidence, not E4 estimator results or a confidence claim.

Source generation is now complete for its stated validation scope. Guarded
execution is the next module. Fitted validation, full data, summary/audit and
visual completion remain pending.

## 3. Guarded fitted execution

The source-validation evidence was published at
`a5dd51edb86999f7bff77b98115d24c0d8b7afa4` before this module began.
`ou_irregular/week8/execution.py` supplies pure preflight, guarded raw execution
and non-mutating saved-run readback. It accepts the exact full/R=4 designs,
uses the declared E4 source and unchanged inference services, and pins the
complete local import closure. See [execution semantics](week8_execution.md).

Independent review identified two imported package initializers absent from
the inherited dependency set. The initial focused suite passed 72 cases and
failed its dependency-closure case. E4's FILES set now explicitly includes
`ou_irregular/week5/__init__.py` and `ou_irregular/week6/__init__.py`; both are
currently docstring-only. This two-path change was the complete production
repair. Old source files and test expectations were unchanged.

The final focused suite passed **73 tests in 18.51 seconds**, using literal
records, mocked simulation/inference boundaries, real durable stores and
mocked process futures. Tests cover exact designs/provenance, freeze refusal,
closed-state handling, no-op cached resume, pipeline and statistical failures,
root drift and preserved incidents, corruption, orphan boundaries and exact
byte adoption, locks, worker identity and execution controls. Independent
read-only review verified all 38 imported repository modules are covered by
the 52-file dependency union. A separate pure preflight probe verified the
honest R=4 development allowance and refusal of an uncommitted full run,
without allocating sources or creating output.

Root reviewed production code and all test cases. `tests/week8/__init__.py`
gives these tests a package namespace alongside historical duplicate filenames;
a targeted collection check collected 227 tests without an import-name clash.
That collection check did not execute those tests. The focused execution tests
generated no schedules, OU paths or fits. Initial source and failure output,
final output, review hashes and collection output are retained under
`artifacts/week8_development/execution/`.

Reviewed execution SHA-256:
`5e7e0b80925b91f72d4940579c7e720c0dac6271ab4683a993593a6da88cc96a`.

The [prospective population-target note](week8_population_targets.md) was also
derived before fitted E4 observations. It binds to the actual existing PFML
and actual-gap Euler definitions, derives their expected-objective minimizers,
and explains matched-marginal interpretation. It changes no registered summary
target or estimator. Reporting is the next implementation module; fitted
validation and full E4 execution remain pending.
