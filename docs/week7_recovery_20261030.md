# Week 7 recovery after the interrupted publication

Recorded 30 October 2026 (India time and UTC).

The live `weeks5-8-sequential` branch was checked after reconnection and still
pointed to `7ed521b4567da5c973c409e47e1b1d9906fa3b47`, the completed Week 6
publication. Main remained `4e35ba23d79ec6dc16f10035df27729817e23678`.
The published Week 5 and Week 6 source, experiments and audits are intact.

Before the interruption, the session had built and reviewed nine Week 7 modules,
observed 971 passing tests, and sealed small validation runs. A local freeze was
created as `ccd91494218748f4789eb78f01dc35f8a6be285d`, with tree
`e3d422f50b12c166d30f56ae7c4753275881f983` (339 changed/added files,
7,399,273 bytes). Publication stalled before the branch update. The connector
later reported that the prepared tool call had been rejected because its catalog
changed. The remote tree was not found on the subsequent check.

After the runtime restart, the original checkout, its Git objects, the isolated
temporary outputs and their unpublished evidence files were unavailable in the
restored workspace. No recoverable copy was located in the designated workspace,
runtime caches or persistent file search. This does not establish whether any
provider-side backup exists. The recorded original source hashes and test results
remain historical conversation evidence; missing file bytes are not recreated
and represented as recovered originals.

**No full Week 7 selected-source replay, independent 10,000-replication study,
or full paired-floor study had started before the loss. Week 8 had not begun.**
The headline E1/E2 evidence is unaffected. Small development observations are
not headline or exploratory results.

## Recovery procedure

The published branch was cloned again with normal Git byte preservation. The
scientific implementation now proceeds in an isolated `/tmp` checkout, outside
the observed workspace synchronization path. The recorded Python/scientific
versions and pytest 9.1.1 are restored; notebook dependencies remain absent,
matching the headline provenance. The documented Week 5 identical-tree commit
object was restored without changing refs or artifacts.

Unpublished modules are recovered from exact available transcript text where
possible and otherwise reconstructed. New tests and independent reviews apply
to the actual reconstructed bytes. The old 971-test result is not a certificate
for this reconstruction. New development runs are explicitly disclosed recovery
validations using the already specified separate validation namespaces. Original
lost development outputs cannot be preserved as local files; that limitation
remains explicit rather than being hidden behind regenerated manifests.

The exploratory protocol's scientific choices are preserved: selection rules,
all fixed seeds/counts, diagnostic ranges/stencils, and paired epsilon coupling.
No new full scientific observations are generated before its reviewed code is
published. Intermediate module checkpoints are clearly labelled unfinished
Week 7 work and published with expected-head guards. Each completed week is
still audited, documented and published before implementation of the next week.

The prior development root-state incidents remain unattributed. The completed
Week 6 archive preserves its own incident evidence. The additional unpublished
Week 7 incident and isolated preview were recorded in the conversation but their
local forensic files were lost; this recovery does not fabricate replacements
or claim a causal repair. The Week 8 gate must retain this limitation.
