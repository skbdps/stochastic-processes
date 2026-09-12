# Expanded stochastic-process teaching notes

The current teaching edition is [`../../stochastic_processes_expanded_notes.html`](../../stochastic_processes_expanded_notes.html).
It preserves all 13 topics from the revised companion and expands the explanations
into 149 core subsections. The HTML uses native MathML and has no external runtime
assets. Download it and open it in a browser; no build is required to read it.

This is the exact HTML delivered as `stochastic_processes_expanded_notes.html` in
the conversation, not the earlier short revision or the alternate later attachment.

## Reproduce the delivered edition

Use a separate documentation environment with Python 3.13.5, Pandoc 3.1.11.1,
and the packages pinned in this directory's `requirements.txt`:

```bash
python -m pip install -r docs/expanded_notes/requirements.txt
python docs/expanded_notes/build_notes.py
python docs/expanded_notes/check_mathematics.py
python -m playwright install chromium
python docs/expanded_notes/check_browser.py
```

The builder joins the readable chapter sources without modifying their wording,
creates native MathML, applies the retained HTML templates, and extracts the
executable Python appendix. It checks the original Markdown, example-code and HTML
hashes before writing generated outputs. Intentional edits require reviewing and
updating the pinned edition hashes in `build_notes.py`; a changed output should not
be mistaken for the original delivered file.

Expected HTML SHA-256:
`022b791e5bfa4f8b91c4e261b21d77434d46ef2485d3c818f40b665b950b8044`

Expected HTML Git blob:
`aa516c81cf2eb8db3bf126dddd84fd40851a35f7`

The 47 test methods check worked-example arithmetic and selected identities,
including likelihood profiling. They are not formal theorem verification,
peer review, or certification of statistical coverage. Browser checks are scoped
to Chromium, internal navigation, native math structure, responsive widths,
print CSS, no-JavaScript reading and absence of external page resources.

The original learning notes and the earlier short revision remain unchanged.
Research implementation, results, milestones and the Week 3 branch are outside
this documentation update. Publishing to this private repository does not turn
on GitHub Pages or make the repository public.
