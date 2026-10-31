# Week 8 integrated runner validation and release gate

The complete runner is published at `edf0867494f1a040b6cd71d03d5feb6236343b40`. CI run 41198942670 accounted for all 2,061 collected test cases: 2,055 passed, six existing design-inapplicable skips, zero failures/errors. The six skips are floor-only checks instantiated for the non-paired follow-up fixture; none waives a scientific failure or strict historical replay. Exact names, reasons, commands and outputs are retained in `artifacts/week8_development/runner/ci/41198942670/`.

The fixed integrated R=4 CLI generated all seven cells, 28 independent paths and 84 fitted records using seed 21261011. Its independent audit reconstructed every source and replayed 42 endpoint fits. The pipeline had zero exceptions. Original CLI evidence remains unchanged. A separate copy is sealed in `artifacts/week8_validation/`, with its figure-renderer state alongside it.

Every one of its 18 PNG files is byte-for-byte identical to the images actually inspected earlier in the plotting-module review. The integrated visual record explicitly reuses that prior judgment with original reviewer/time and per-image hashes, while binding the new figure manifest. It does not pretend a new visual inspection occurred. No synthetic unit-test attestation is used. External replotting reproduced all 36 PNG/PDF files byte-for-byte without changing the sealed archive. See `artifacts/week8_evidence/validation_gate.json`.

The source-only validation remains 1,536 sequences and 108 passed checks, with every declared artifact rehashed at this gate. It did not generate fitted OU data.

## CI recovery and historical reproducibility

Earlier workflows and unsuccessful outputs remain on validation branches and in retained diagnostics. Workflow ordering, spawn-safe pytest entry and historical/current dependency environments were corrected; scientific code, archived outputs, seeds, thresholds and assertions were not changed to obtain a green run. The strict Week3 replay runs in its own process with the same exact input-hash and fitted-value assertions. All four historical isolation tests passed without skips. All other cases ran separately; unique test counts were checked against complete collection.

A separate diagnostic reconstructed all 180 archived Week3 sources exactly and reproduced the selected fitted values bit-for-bit under its default CPU dispatch. Disabling AVX-512 on that same host changed some observation hashes. This demonstrates hardware-dependent floating-point reproducibility, but does not prove the exact cause of the earlier combined-process failure, whose CPU state was not captured. The successful run records CPU flags, dispatch and NumPy/BLAS runtime alongside pinned versions. Exact replay elsewhere must verify hashes, not assume matching package versions suffice.

Primary references: https://numpy.org/doc/2.4/reference/simd/build-options.html and https://numpy.org/doc/stable/reference/random/compatibility.html .

This gate certifies the stated implementation and development validation only. Full E4 data had not been generated at the gate. Full R=2000 execution, its complete audit, actual new-figure inspection, interpretation and final publication remain required before calling the extension complete.
