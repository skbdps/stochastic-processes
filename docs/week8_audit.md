# Week 8 independent E4 numerical audit

This module is a read-only verifier, not an experiment runner or a replacement
for visual inspection. It accepts only the complete fixed seven-cell full or
R=4 validation design. It uses no producer raw-record loader, table calculator,
reference-map implementation or source adapter to construct expected results.

## Checks and independence boundaries

Every planned source is reconstructed using independently derived namespace-8
addresses and a second arrival/switch race implementation. The original numeric
cell calibration is retained so random-call ordering and array-byte comparisons
are meaningful. MMPP retains the arrival phase; the IID control resets its
Palm phase at each gap. Flooring is applied after the raw race, with no hidden
chain advance through padding and no renormalization. Both source-array hashes,
both gap-vector hashes, the phase/switch trace, every source diagnostic and
available observed-data normalization triples are checked. Invalid source draws
cannot be replaced or resampled.

The calibrated cell model, analytical gap-moment functions and exact OU
primitive are shared dependencies, not a second independent derivation of those
components. Their earlier rational/spectral/source-only checks remain necessary.
The second race and raw-summary reductions provide an additional implementation
check without claiming complete independence of the entire numerical stack.

Strict parsing rejects duplicate keys, nonfinite constants/overflow, unexpected
fields, missing Cartesian keys, changed order/types and contradictory validity.
The reviewed independent core auditor checks saved inference geometry,
normalization, inverse-Hessian/delta covariance, numerical gates and unmodified
natural-scale endpoints. Coarse-Hessian reconstruction and complete fit replay
are performed only at the fixed endpoint replications, not claimed for every fit.

Core parameter, paired-estimator, spacing, outcome and numerical formulas are
reused from the independent auditor, not from reporting. E4-specific reference
selection, own-eligible-subset contrast arithmetic and source-path reductions
are separately implemented. Shared names and schemas define interoperability;
expected numeric values do not call the producer's summary functions. Every
field and planned key of all nine tables is verified, including an exact full
failure header when no failures occur. Failure rows retain their complete
candidates, reasons and nested diagnostics in canonical sorted export order.

## Replay, provenance and immutability

The full audit replays original fitting/inference services at replications 0
and R-1 in every cell: 14 paths and 42 fits for either fixed design. All raw and
nested fields, missingness and statuses are compared. These endpoint replays
share the original inference services and serializer; they are reproduction,
not a second estimator implementation.

Arithmetic tolerances remain rtol=2e-10 and atol=2e-12. Source hash comparisons
are exact. Saved covariance consistency uses the already reviewed core gate
rtol=2e-9 and atol=2e-12 times max(1, expected matrix Frobenius norm). None changes
an estimator, scientific Hessian acceptance gate or simulation design.

`replay=False` skips only endpoint fitting. Every source is still reconstructed,
and the report explicitly cannot certify full numerical completion. A closed
archive requires an existing bound full audit; a lighter readback rechecks its
invariant fields and endpoint descriptors against newly verified source arrays.

The auditor independently binds all raw partitions, cell manifests, resolved
design, source policy, table manifest and table binding. It verifies the exact
recorded dependency set/runtime and, for the full experiment, the committed
freeze. Full before/after byte and directory snapshots guard read-only behavior.
Closed-state audit and visual-review seals are verified without writing locks,
incident files, regenerated tables or metadata. Actual image legibility remains
a separate human-visible inspection, never inferred from hashes or filenames.

This contract describes checks; completion and test evidence are recorded
separately after actual execution.
