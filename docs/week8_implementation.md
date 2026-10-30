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
