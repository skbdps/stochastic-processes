# Review record for Week 10 Draft 1

The reviewed document is [the first published methods manuscript](https://github.com/skbdps/stochastic-processes/blob/3b546268da1a290f07cf741185123192fad23e47/manuscript/week10_draft1.md).
Its reviewed SHA-256 is
`80a22e2eb97f48c0b720691da3d429778691e6a230e5f29466984b5bed256d35`.

Two content/source checks covered the technical methods and introduction;
their findings were incorporated into the draft and supporting documents.
Their scope excludes academic assessment, rendered-page review, publication
acceptance and certification of novelty. Local review times fall
on 1 November 2026 in Asia/Kolkata; UTC timestamps below retain the actual
recorded instants.

This summary describes the original reviewed revision. The
[editorial revision record](editorial_revision.json) identifies later wording
changes; the original [validation record](validation.json) retains its dated
hashes. The complete original review is preserved in the
[published release](https://github.com/skbdps/stochastic-processes/blob/3b546268da1a290f07cf741185123192fad23e47/artifacts/week10_draft1/review.md).

## Technical review

Review component: technical methods.

The complete initial technical review was recorded at
`2026-10-31T19:28:26.121400+00:00`, against manuscript digest
`120638b754ed2bfed4b1988ee97fecdba6b9670bd9296ac61d34784f99121948`.
The reviewer compared the transition law, all three estimators, covariance
normalization, numerical gates, population targets, designs, seeds, counts,
metrics, pairing and audit scope with the protocols and implemented files in
the [source map](../../manuscript/week10_source_map.md).

The initial disposition required three corrections:

1. Distinguish dimensional nominal/realized spans from their dimensionless
   counterparts.
2. Describe retained gap diagnostics, reconstruction identities and hashes
   accurately, without implying that full gap and observation arrays are
   stored in the fit partitions.
3. Name point availability as well as interval availability.

All were corrected. The draft also now defines the half-step gradient,
identifies the shared regular floor control, and states explicitly why the
overall practical classification is inapplicable to zero-truth mu.

At `2026-10-31T19:30:17.493100+00:00`, the reviewer checked the corrections and
substantive additions on the final manuscript digest above. The added PFML
population-mean versus realized-mean distinction was independently checked
in Aït-Sahalia and Mykland (2003), Section 4.3, equation 4.22, printed page 505.
**Final disposition: passed within the methods-review scope; no remaining
technical blocker identified.** This is not a new numerical audit or a rerun
of any research study.

## Introduction and literature review

Review component: introduction and literature.

The introduction, references and contribution scope were initially reviewed at
`2026-10-31T19:27:17+00:00`, against manuscript digest
`f85ba1510f448b7e95da6793e03d711e13a0996eaf8db9b128ed308c9aac70d6`.
Primary-source checks refreshed the three formal references, the PFML
attribution, the existence of finite-sample simulations in the 2003 paper,
and the IID assumption in the 2004 theory.

The review required the explicit distinction between the original theoretical
PFML mean-gap convention and this project's realized mean-gap procedure. It
also recommended narrower IID wording, neutral introductory phrasing and DOI
links. These changes are present in the final draft.

At `2026-10-31T19:29:17+00:00`, the reviewer read back the complete introduction,
formal references, declared scope and supporting pointers on the final digest
above. **Final disposition: passed within that review scope; no material
introduction, attribution or citation issue identified.** The wider related-work
section and pre-submission search remain later tasks.

| Reference | Checked material and limit |
| --- | --- |
| Vasicek (1977) | [Publisher record](https://www.sciencedirect.com/science/article/pii/0304405X77900162) confirms bibliographic metadata and short-rate diffusion motivation. Direct full-text retrieval was blocked; no bond-pricing formula or coefficient is imported into this draft. |
| Aït-Sahalia and Mykland (2003) | [Author-hosted published article](https://www.princeton.edu/~yacine/random.pdf), including the PFML convention and simulation section/Table V; [institutional metadata](https://collaborate.princeton.edu/en/publications/the-effects-of-random-and-discrete-sampling-when-estimating-conti/) confirms the citation. |
| Aït-Sahalia and Mykland (2004) | [Author-hosted published article](https://www.princeton.edu/~yacine/gmmrandom.pdf), Assumption 2 and OU nuisance conventions; [institutional metadata](https://collaborate.princeton.edu/en/publications/estimators-of-diffusions-with-randomly-spaced-discrete-observatio/) confirms the citation. |

The existing [Week 9 literature comparison](../../docs/week9_literature_comparison.md)
retains the broader comparison, recent nearby coverage work and outstanding
access limitations. This first draft does not claim exhaustive absence of
overlapping work or a matched dependent-gap covariance theorem.
