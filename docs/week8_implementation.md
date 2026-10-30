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

Remaining modules: source generation and analytical checks, guarded execution,
reporting, independent audit, plots, and runner/integration. The complete
source-only validation and fitted R=4 validation precede a reviewed full-run
freeze. No full E4 observations have been generated.
