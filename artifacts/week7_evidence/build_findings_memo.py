"""Render the repository's one-page Week 7 memo; no scientific calculations."""
from pathlib import Path
import re
from urllib.parse import urljoin
from xml.sax.saxutils import escape

from reportlab import rl_config
from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import Paragraph, SimpleDocTemplate, Spacer
from pypdf import PdfReader

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "docs/week7_findings_memo.md"
OUTPUT = ROOT / "docs/week7_findings_memo.pdf"
BASE = "https://github.com/skbdps/stochastic-processes/blob/weeks5-8-sequential/docs/"
FONTDIR = Path("/usr/share/fonts/truetype/dejavu")
pdfmetrics.registerFont(TTFont("Memo", str(FONTDIR / "DejaVuSans.ttf")))
pdfmetrics.registerFont(TTFont("MemoBold", str(FONTDIR / "DejaVuSans-Bold.ttf")))
pdfmetrics.registerFontFamily("Memo", normal="Memo", bold="MemoBold")
rl_config.invariant = 1

body = ParagraphStyle("body", fontName="Memo", fontSize=9.4, leading=12.7,
                      textColor=colors.HexColor("#233247"), spaceAfter=7,
                      alignment=TA_LEFT)
title = ParagraphStyle("title", parent=body, fontName="MemoBold", fontSize=23,
                       leading=28, textColor=colors.HexColor("#12324D"), spaceAfter=7)
heading = ParagraphStyle("heading", parent=body, fontName="MemoBold", fontSize=11,
                         leading=14, textColor=colors.HexColor("#176B78"),
                         spaceBefore=6, spaceAfter=5)
evidence = ParagraphStyle("evidence", parent=body, fontSize=8.1, leading=10.5,
                          spaceAfter=3)


def inline(text):
    text = text.replace("\n", " ").replace("\u2013", "-").replace("\u2014", "-")
    text = escape(text)
    text = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", text)
    text = re.sub(r"`([^`]+)`", r"\1", text)
    return re.sub(r"\[([^\]]+)\]\(([^)]+)\)",
                  lambda m: '<link color="#176B78" href="' +
                  escape(urljoin(BASE, m.group(2)), {'"': '&quot;'}) +
                  '">' + m.group(1) + '</link>', text)


def footer(canvas, doc):
    canvas.saveState()
    canvas.setStrokeColor(colors.HexColor("#D7E1E8"))
    canvas.line(40, 32, A4[0] - 40, 32)
    canvas.setFont("Memo", 7.5)
    canvas.setFillColor(colors.HexColor("#536579"))
    canvas.drawString(40, 20, "OU estimation under irregular sampling | Week 7 exploratory evidence")
    canvas.drawRightString(A4[0] - 40, 20, str(doc.page))
    canvas.restoreState()


story = []
in_evidence = False
for block in SOURCE.read_text().strip().split("\n\n"):
    if block.startswith("# "):
        story.append(Paragraph(inline(block[2:]), title))
    elif block.startswith("## "):
        story.append(Paragraph(inline(block[3:]), heading))
    elif block.startswith("### "):
        in_evidence = True
        story.append(Paragraph(inline(block[4:]), heading))
    elif in_evidence:
        for line in block.splitlines():
            story.append(Paragraph(inline(line), evidence))
    else:
        story.append(Paragraph(inline(block), body))

doc = SimpleDocTemplate(str(OUTPUT), pagesize=A4, leftMargin=40, rightMargin=40,
                        topMargin=30, bottomMargin=43,
                        title="Week 7 findings memo", author="OU research project",
                        subject="Completed exploratory diagnostics and floor sensitivity")
doc.build(story, onFirstPage=footer, onLaterPages=footer)
reader = PdfReader(OUTPUT)
if len(reader.pages) != 1:
    raise RuntimeError(f"Memo must be one page; rendered {len(reader.pages)}")
text = reader.pages[0].extract_text()
for required in ("60,000", "72,000", "231", "0.095654", "1.00023261", "4,396"):
    if required not in text:
        raise RuntimeError("Missing memo content: " + required)
print(f"Created {OUTPUT}: one page, {len(text.split())} extracted words")
