# Week 8 publication recovery

## Original execution and retained failure

Full E4 was executed once at scientific freeze `0ff57919f5d8345f52d8689732ff94448ca3ca1d` by [Actions run 41200456847](https://github.com/skbdps/stochastic-processes/actions/runs/41200456847), workflow execution commit `4deaa94d5735bcfff6c6cb47f8cff88df2fb2eb3`, job `124535057563`. The scientific CLI exited 0 with empty stderr, completed all seven R=2,000 cells (14,000 source paths and 42,000 fits), wrote the passed full audit, and generated 18 PNG/PDF pairs. Its original lifecycle remained `audited_pending_visual_review`. These execution facts do not themselves certify visual review or sealing.

The subsequent saved-evidence verification and summary creation succeeded. At 2026-10-31 12:49:35 UTC, the publication preparation failed at `POST repos/skbdps/stochastic-processes/git/trees`:

```text
gh: tree.sha 40af4433b49d2cef073384b77a2d2213963c5cb6 is not a valid blob (HTTP 422)
```

The retained workflow loops over all recovered artifact files, creates base64 Git blobs, then supplies the returned IDs as blob entries to one new tree based on the freeze. The traceback identifies this tree request as the failing operation. The retained log does not include the individual blob-creation responses or full submitted tree request. A later artifact upload succeeded; the scientific branch was not advanced by this failed publication operation. The failed workflow and original log remain evidence, rather than being replaced with an apparent success.

## Narrow object diagnosis

Computing Git's object digest `SHA1("blob " + decimal_byte_length + NUL + file_bytes)` over every original artifact member uniquely maps the rejected ID to:

| Item | Value |
|---|---|
| Original member | `week8_e4/figures/W8_source_span_sd.pdf` |
| Intended repository path | `artifacts/week8_e4/figures/W8_source_span_sd.pdf` |
| Byte length | 16,169 |
| Git blob SHA-1 | `40af4433b49d2cef073384b77a2d2213963c5cb6` |
| File SHA-256 | `09138fefc42b0d2b41b4f3aafe5480d13b7e92e69b622323d8fdea7dd7e6f7a9` |

A native read-only check reproduces the ID:

```sh
git hash-object --no-filters <original>/week8_e4/figures/W8_source_span_sd.pdf
```

It is not the empty-file object: the empty Git blob is `e69de29bb2d1d6434b8b29ae775ad8c2e48c5391`. Both the zero-byte lock file and zero-byte stderr have that different ID. The preserved PDF bytes are therefore the concrete object associated with the error. This identifies the affected object; it does not establish why GitHub rejected the object at that historical instant. No evidence here proves a scientific failure, permission failure, empty-lock-file defect, malformed PDF, or a particular GitHub backend fault. Connector attempts to read this binary object encountered UTF-8 decoding restrictions and do not establish its historical server state.

## Recovery boundary

The original artifact ZIP is preserved independently of the publishable copy. Its size is 35,558,419 bytes and SHA-256 is `e9a780206925c4085ff4bb47cf988c7ccec66ce40db038f3fdef82703fe52513`. Safe extraction recovered 124 regular-file members under `week8_e4/`, `week8_e4-figure-renderer/`, and `evidence/`; each extracted member was hashed. The untouched extraction is retained alongside the member-hash record. The original execution logs, metadata, numerical audit, tables and figure bytes are retained, including the renderer state beside the archive.

Publication repair does not rerun the experiment, change seeds or scientific source, replace saved numerical evidence, or relabel regenerated output as the original run. Only the existing guarded completion transition may add the actual visual review and completion seal to a verified working copy. The pre-seal root metadata remains in the untouched original.

## Selected publication alternative and guards

Recovery used a normal Git checkout and commit of the recovered, reviewed evidence instead of the failed handcrafted Git-tree preparation. [Actions 41218787133](https://github.com/skbdps/stochastic-processes/actions/runs/41218787133) passed the saved-evidence and external-renderer readbacks, verified every staged and committed blob, and published checkpoint `f1ad8cc5a7cd2831c6888f2840961178eeeceba7` on the isolated recovery branch. All 184 inventory file hashes were read back after publication. Inspect and fetch current remote heads first. Base work on the scientific branch after comparing any newer commits; retain the completion/orchestration and diagnostic branches without wholesale merges. Do not reset newer work, force-push, bypass protection, or edit frozen scientific modules to solve transport.

Before committing, compare every recovered source/raw/table/figure/renderer/evidence path to its original SHA-256. Allow only the documented guarded lifecycle changes and new review/recovery documents. Check all pre-existing tracked files against their base, with intentional documentation updates identified separately. Verify the completed archive and external figure reproduction before release. Use a fast-forward update only while the expected remote head still matches; if it has moved, inspect and integrate it before proceeding. A successful upload alone is not a completion certificate.

After publication, fetch/read back the remote main commit and key artifacts, checking the original run identity, raw and table hashes, all figure hashes, renderer state, actual review binding and completed seal. Record the resulting scientific/release commit IDs, final test accounting and readback hashes in the completion report and release evidence. The recovery remains unfinished until those checks actually pass.

## Final release record

The publication commit, remote main SHA, clean-checkout validation results, exact preservation comparison, completed-run verification and external replot readback are recorded by the release step in the Week 8 completion report and its accompanying evidence. This document explains the observed publication failure and repair boundary; it does not substitute a proposed procedure for those results.

The first recovery wrapper (run 41214670935) completed its full numerical re-audit, guarded seal, repeated read-only seal and external rendering, then looked for files one directory above the actual `external_replot/figures/` output. Its retained outputs supplied the corrected saved-byte continuation; the numerical audit and renderer were not rerun to repair that lookup. Forced staging of the exact artifact paths also preserves ignored Week 8 paths and empty lock bytes.

The original failed workflow remains at its pinned source branch. Full raw CI job logs remain in the original local handoff/evidence; public focused diagnostic records retain their SHA-256 and source links. Duplicate raw-log upload was declined by automatic approval review. The exact public scientific payload was subsequently explicitly approved by the owner. Original scientific CLI logs and runtime evidence remain part of the recovered artifact.
