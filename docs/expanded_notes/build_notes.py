"""Rebuild the delivered teaching edition from readable Markdown and templates.

Pandoc 3.1.11.1 and beautifulsoup4 4.14.3 reproduce the original HTML exactly.
No network calls are made by this script. The expected hash pins the approved
edition; an intentional content edit requires reviewing and updating that pin.
"""
from __future__ import annotations

import hashlib
import html
import json
from pathlib import Path
import re
import subprocess

from bs4 import BeautifulSoup

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
EXPECTED_HTML = "022b791e5bfa4f8b91c4e261b21d77434d46ef2485d3c818f40b665b950b8044"
EXPECTED_MD = "1ecdad44a215d0270183ecd587472d81b2e8e90f1ff4ca21dd4af682d82822ba"
EXPECTED_CODE = "c0e4b92a8b6349e4ee9e8c8839b9449f8f6c538d4f12f0cad8340f3bfc4eee3d"


def digest(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def build() -> None:
    version = subprocess.run(["pandoc", "--version"], check=True, capture_output=True,
                             text=True).stdout.splitlines()[0]
    if version != "pandoc 3.1.11.1":
        raise RuntimeError(f"Use pandoc 3.1.11.1 for this exact edition, got {version}")
    sources = HERE / "chapter_sources"
    chapters = sorted(sources.glob("[0-9]*.md"))
    if len(chapters) != 13:
        raise ValueError("Expected all thirteen chapter sources")
    paths = [sources / "intro.md", *chapters, sources / "appendices.md"]
    markdown = "\n\n".join(path.read_text(encoding="utf-8") for path in paths)
    if digest(markdown) != EXPECTED_MD:
        raise ValueError("The assembled source differs from the delivered edition")
    code = (sources / "appendices.md").read_text(encoding="utf-8").split(
        "```python\n", 1)[1].rsplit("```", 1)[0]
    if digest(code) != EXPECTED_CODE:
        raise ValueError("The executable appendix differs from the delivered module")

    body = subprocess.run(
        ["pandoc", "--from", "markdown+tex_math_dollars+fenced_divs", "--to", "html5",
         "--mathml", "--wrap=none", "--no-highlight"],
        input=markdown, text=True, encoding="utf-8", check=True, capture_output=True
    ).stdout
    body = re.sub(r'(<math display="(inline|block)".*?</math>)',
                  lambda m: '<span class="math ' +
                  ('inline' if m[2] == 'inline' else 'display') + '">' + m[1] + '</span>',
                  body, flags=re.S)
    body = re.sub(r'(<table(?:\s[^>]*)?>.*?</table>)',
                  r'<div class="table-wrap">\1</div>', body, flags=re.S)
    soup = BeautifulSoup(body, "html.parser")
    for span in soup.select("span.math.display"):
        span["aria-label"] = "Displayed equation; scroll horizontally if needed"
        span["tabindex"] = "0"
        span.math["displaystyle"] = "true"
        for table in span.select("mtable"):
            table["displaystyle"] = "true"
    for wrap in soup.select("div.table-wrap"):
        wrap["aria-label"] = "Scrollable table"
        wrap["role"] = "region"
        wrap["tabindex"] = "0"

    sections = []
    for chunk in re.split(r'(?=<h2 id=")', str(soup)):
        if not chunk.strip():
            continue
        heading = re.match(r'<h2 id="([^"]+)">(.*?)</h2>', chunk, re.S)
        if heading is None:
            raise ValueError("Unexpected content outside the chapter structure")
        ident = heading[1]
        is_chapter = bool(re.match(r'\d+\.', BeautifulSoup(
            heading[2], "html.parser").get_text()))
        if is_chapter:
            items = []
            for h in re.finditer(r'<h3 id="([^"]+)">(.*?)</h3>', chunk, re.S):
                title = BeautifulSoup(h[2], "html.parser").get_text()
                items.append(f'<li><a href="#{h[1]}">{html.escape(title)}</a></li>')
            toc = ('<div class="chapter-map"><p class="map-label">IN THIS CHAPTER</p><ol>'
                   + ''.join(items) + '</ol></div>')
            chunk = chunk[:heading.end()] + toc + chunk[heading.end():]
            chunk += '<p class="return-link"><a href="#book-contents">Return to contents ↑</a></p>'
        if ident == "code-example":
            button = '<button class="js-only" id="copy-code" type="button">Copy Python example</button>'
            chunk = chunk[:heading.end()] + button + chunk[heading.end():]
        section_class = "chapter" if is_chapter else "appendix"
        sections.append(f'<section aria-labelledby="{ident}" class="{section_class}">' +
                        chunk + '</section>')
    prefix = (HERE / "templates/page_prefix.html").read_text(encoding="utf-8")
    suffix = (HERE / "templates/page_suffix.html").read_text(encoding="utf-8")
    rendered = prefix + ''.join(sections) + suffix
    actual = digest(rendered)
    if actual != EXPECTED_HTML:
        raise ValueError(f"HTML did not reproduce the approved edition: {actual}")

    output = ROOT / "stochastic_processes_expanded_notes.html"
    output.write_bytes(rendered.encode("utf-8"))
    (HERE / "stochastic_processes_expanded_notes.md").write_bytes(markdown.encode("utf-8"))
    (HERE / "ou_likelihood_examples.py").write_bytes(code.encode("utf-8"))
    report = {"html_sha256": actual, "html_bytes": output.stat().st_size,
              "markdown_sha256": digest(markdown), "example_sha256": digest(code),
              "chapters": len(chapters), "pandoc": version,
              "scope": "Byte-for-byte reproduction of the delivered edition; not a theorem proof."}
    (HERE / "validation").mkdir(exist_ok=True)
    (HERE / "validation/build_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    build()
