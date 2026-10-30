# Auditor development validation plan — 31 October 2026

Before any fitted auditor fixture: verify literal summaries and the independent
schedule reconstructor against all 1,536 existing source-only sequences. The
full R=2000 seed is not used by these tests.

Then create exactly the fixed seven-cell R=4 development design (seed21261011,
28 paths/84 fits) once per test invocation, using the original source/fit services.
Retain the complete archive and invocation source snapshots outside temporary
pytest cleanup. The auditor reconstructs every source and replays replications0
and3 in all seven cells (42 additional fits). Corruption tests operate on copies,
not the retained archive. This is a native module-integration fixture, not the
later acceptance of the complete plotting/runner pipeline. Any subsequent
technical reproduction uses the same fixed development addresses, preserves
prior output and documents source changes. Neither development results nor
unit-test counts can provide full E4 research evidence.

This explicitly extends the earlier unit-only checkpoint evidence to test the
new auditor's real source/geometry/arithmetic composition. The full implementation
freeze and complete-runner R=4 validation still precede full E4 data. No seed,
replication count, statistical criterion, or inference gate is changed.
