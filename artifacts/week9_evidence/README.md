# Week 9 validation, review and failure evidence

The [completed publication package](../week9_publication/) is a derivative of
the existing studies. These records distinguish its actual checks from earlier
failures and unsupported-host probes. The renderer freeze is
`dad63eced437736cc32f7a481b687e226cb4bdaf`.

| Record | Actual outcome and scope |
| --- | --- |
| [week9_test_summary.json](week9_test_summary.json) and [final_week9_tests/](final_week9_tests/) | All 92 Week 9 tests passed under the recorded Python 3.12.10 renderer. Exact source/test snapshot, stdout/stderr, command and JUnit retained. |
| [final_build/](final_build/) and [clean_build/](clean_build/) | Both real C1-pinned builds succeeded from saved evidence, with zero new sources/fits. The second was a detached clean checkout. |
| [external_reproduction.json](external_reproduction.json) | All 141 immutable files matched byte-for-byte, including all numerical views and PNG/PDF outputs. |
| [actual_visual_review.json](actual_visual_review.json), [visual_inspection/](visual_inspection/) and [pdf_review/](pdf_review/) | Actual AI reviews, their original hashes/times, explicit identical-byte transfer, all 11 rendered PDF pages and font reports. No human-review claim. |
| [finalization_result.json](finalization_result.json) and [final_readback.json](final_readback.json) | Guarded completion against committed renderer bytes and subsequent read-only/idempotent verification. |
| [source_preservation.json](source_preservation.json) | Baseline tracked-file preservation; project README and artifact index are the two intentional navigation updates. |
| [development_build_01/](development_build_01/) | Failed real build: semantically equivalent reserialization changed the pinned catalog bytes. The verifier rejected it. Original logs/manifests and exact renderer sources retained; the full failed output also remains in local work evidence. Fixed by copying the original catalog bytes. |
| [cli_test_attempt01/](cli_test_attempt01/) | Initial harness produced 17 setup errors from importing a historical YAML dependency into the lean renderer. The original observed summary is retained. Full logs from a reproduction of the unchanged failed harness are explicitly labelled replay; they are not represented as the original invocation. |
| [cli_test_attempt02/](cli_test_attempt02/) | Actual run: one failure, 16 passes. It exposed dictionary-order-dependent claims CSV columns after JSON round-trip. The actual JUnit/observed failure and exact source snapshots are retained. Fixed by a stable column order. |
| [cli_test_attempt03/](cli_test_attempt03/) | Actual intermediate run: 21 passed. The final run adds complete figure-role/accounting negative checks and supersedes this smaller passing checkpoint. |
| [legacy_saved_tests/](legacy_saved_tests/) | Failed optional Mac probe: 78 passes, seven failures, 13 errors, three original inapplicable skips. Eighteen archived-environment mismatches and two extended-precision assumptions are documented without rewriting the outcomes or old tests. |

The final saved-data verification and passing Week 9 tests do not certify old
numerical replay on this Mac. Earlier Linux scientific/replay evidence remains
unchanged. The [completion report](../../docs/week9_completion.md) states that
boundary and the remaining mathematical/publication limitations.
