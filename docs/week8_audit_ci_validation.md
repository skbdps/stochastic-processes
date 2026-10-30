# Week 8 auditor: clean-checkout validation

Recorded 31 October 2026. GitHub Actions run `41145059340`, validation revision `205cb5d6f363bf5be6a8eafba545563b042dc33c`, completed successfully. All **366 Week 8 tests passed in 171.74 seconds**. This combines the 296 existing module tests and 70 auditor tests; it is not 366 independent experimental replications.

## Exact publication

The validation branch started at scientific checkpoint `f82bbade0e94df2df3a6a4caa5937f403865cdd7`. Before running tests, CI restored 32 candidate files from a payload whose compressed SHA-256 is `39c18f3e7c49ef6f75712f68605373fde0de907e65fdc95ce3fc8440cd5dc10f`, and verified every restored file against its individual SHA-256. The payload contains source, tests, documentation and retained local development evidence, not a change to the scientific design.

After the tests, CI rechecked every candidate-file hash and uploaded those exact bytes together with its own complete native fixture and test evidence as Git objects. Its candidate publication tree is `8e1c9cae017e30b3118eeed5d170196cfade4d5e`, based on the scientific parent tree. No CI workflow or encoded transport payload is included in the scientific publication. CI did not update a branch; publication remains a separate guarded fast-forward.

The downloaded CI artifact was checked against the publication manifest: 73 file contents matched their expected SHA-256 values. The remaining entry is the canonical zero-byte `.execution.lock` Git blob `e69de29bb2d1d6434b8b29ae775ad8c2e48c5391`; GitHub's artifact ZIP omits hidden files, but the Git publication retains this empty lock marker. The local historical native archives remain in the separately supplied recovery bundle identified in `week8_audit_progress.md` and the committed archive inventory. They are not confused with this new CI fixture.

## Numerical scope

The CI native fixture uses the fixed validation design, R=4, seed 21261011: seven cells, 28 independent paths and 84 original fits. All 84 point estimates and intervals were available. The independent auditor reconstructed every source and replayed 42 complete fits at the fixed endpoint replications. All nine saved tables were checked. This is development evidence only; no full-seed E4 results are inferred.

The CI fixture and its full audit report are retained under `artifacts/week8_development/audit/ci/native_fixture/`. The exact command output, installed dependency list, Python version, CI head, candidate manifest and payload verification are retained alongside them. The runtime is CPython 3.12.14 on Ubuntu 24.04, with the repository-pinned direct requirements and one BLAS/OpenMP thread. Transitive dependencies are recorded rather than claimed identical to the historical lock file.

The auditor, its tests and its documentation now form a verified implementation checkpoint. Saved-table figures and guarded runner/integration remain the next modules. Complete R=4 pipeline validation and a published complete implementation freeze still precede the full R=2000 experiment.
