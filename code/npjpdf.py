"""Shared reportlab helpers for the npj submission PDFs.

Derived from the journal-submission-prep skill kernel, with one correction:
Arial has no glyph for U+207B/U+2074/U+2075 (superscript minus and digits), so
unicode exponents in the markdown were silently dropped by the previous build.
`fix_exponents` rewrites those runs as reportlab <super> markup before escaping.
"""
import re
import os
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import mm as MM
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.enums import TA_JUSTIFY
from reportlab.lib import colors
from reportlab.platypus import Paragraph, Spacer, Image as RLImage, Table, TableStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfbase.pdfmetrics import registerFontFamily
from PIL import Image

SUPMAP = {"\u2070": "0", "\u00b9": "1", "\u00b2": "2", "\u00b3": "3", "\u2074": "4",
          "\u2075": "5", "\u2076": "6", "\u2077": "7", "\u2078": "8", "\u2079": "9",
          "\u207b": "-", "\u207a": "+"}
_SUPRE = re.compile("[" + "".join(SUPMAP) + "]+")


def fix_exponents(t):
    """Turn runs of unicode superscript characters into <super>...</super>.

    Applied only where the run follows a digit or a closing paren (i.e. a real
    exponent), so standalone U+00B2 in 'I2' / 'chi2' style tokens is preserved
    as a glyph Arial can render.
    """
    def rep(m):
        return "<super>" + "".join(SUPMAP[c] for c in m.group(0)) + "</super>"
    return _SUPRE.sub(rep, t)


def si_styles(font_dir="/System/Library/Fonts/Supplemental/", family="Arial"):
    faces = {"": f"{family}.ttf", "-B": f"{family} Bold.ttf",
             "-I": f"{family} Italic.ttf", "-BI": f"{family} Bold Italic.ttf"}
    for suf, fn in faces.items():
        path = os.path.join(font_dir, fn)
        if os.path.exists(path):
            pdfmetrics.registerFont(TTFont(family + suf, path))
    registerFontFamily(family, normal=family, bold=family + "-B",
                       italic=family + "-I", boldItalic=family + "-BI")
    INK = colors.HexColor("#1a1a1a")
    GREY = colors.HexColor("#666666")
    S = {
        "title": ParagraphStyle("title", fontName=family + "-B", fontSize=15, leading=19,
                                spaceAfter=8, textColor=INK),
        "h1": ParagraphStyle("h1", fontName=family + "-B", fontSize=13, leading=16,
                             spaceBefore=14, spaceAfter=7, textColor=INK),
        "h2": ParagraphStyle("h2", fontName=family + "-B", fontSize=10.5, leading=13,
                             spaceBefore=11, spaceAfter=5, textColor=INK),
        "body": ParagraphStyle("body", fontName=family, fontSize=8.6, leading=12.2,
                               spaceAfter=5, alignment=TA_JUSTIFY, textColor=INK),
        "cap": ParagraphStyle("cap", fontName=family, fontSize=8.0, leading=11.0,
                              spaceBefore=4, spaceAfter=9, alignment=TA_JUSTIFY, textColor=INK),
        "note": ParagraphStyle("note", fontName=family, fontSize=7.4, leading=10,
                               textColor=GREY, spaceAfter=6),
    }
    return S, A4[0] - 2 * 20 * MM


def md_inline_to_rl(t):
    import html as _h
    t = fix_exponents(t)
    t = _h.escape(t)
    t = t.replace("&lt;super&gt;", "<super>").replace("&lt;/super&gt;", "</super>")
    t = re.sub(r"\*\*\*(.+?)\*\*\*", r"<b><i>\1</i></b>", t)
    t = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", t)
    t = re.sub(r"(?<!\*)\*(?!\s)(.+?)(?<!\s)\*(?!\*)", r"<i>\1</i>", t)
    t = re.sub(r"`(.+?)`", r"<font face='Courier'>\1</font>", t)
    t = re.sub(r"\^((?:\\?[\d,\u2013\u2009*\u2020\u2021\u00a7\u2709])+)\^",
               lambda m: "<super>" + m.group(1).replace("\\", "") + "</super>", t)
    t = re.sub(r"\\([*_`^~\[\]()#+\-.!])", r"\1", t)
    return t


def rl_table(rows, widths=None, font_size=7.0, family="Arial"):
    INK = colors.HexColor("#1a1a1a")
    cs = ParagraphStyle("c", fontName=family, fontSize=font_size,
                        leading=font_size + 2.0, textColor=INK)
    hs = ParagraphStyle("h", fontName=family + "-B", fontSize=font_size,
                        leading=font_size + 2.0, textColor=colors.white)
    data = [[Paragraph(md_inline_to_rl(str(c)), hs if i == 0 else cs) for c in r]
            for i, r in enumerate(rows)]
    t = Table(data, colWidths=widths, repeatRows=1, hAlign="LEFT")
    t.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#3a4a5a")),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("TOPPADDING", (0, 0), (-1, -1), 2.5), ("BOTTOMPADDING", (0, 0), (-1, -1), 2.5),
        ("LEFTPADDING", (0, 0), (-1, -1), 4), ("RIGHTPADDING", (0, 0), (-1, -1), 4),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#f4f6f8")]),
        ("LINEBELOW", (0, -1), (-1, -1), 0.5, colors.HexColor("#cccccc")),
    ]))
    return t


def md_flowables(md, styles, frame_width, family="Arial"):
    out = []
    lines = md.split("\n")
    i = 0
    while i < len(lines):
        l = lines[i]
        if re.match(r"^\|.*\|\s*$", l):
            blk = []
            while i < len(lines) and re.match(r"^\|.*\|\s*$", lines[i]):
                blk.append(lines[i])
                i += 1
            rows = [[c.strip() for c in b.strip().strip("|").split("|")] for b in blk]
            rows = [r for r in rows if not all(set(c) <= set("-: ") for c in r)]
            if rows:
                nc = len(rows[0])
                out.append(rl_table(rows, [frame_width / nc] * nc, family=family))
                out.append(Spacer(1, 7))
            continue
        mi = re.match(r"^!\[.*?\]\((.+?)\)\s*$", l.strip())
        if mi and os.path.exists(mi.group(1)):
            im = Image.open(mi.group(1))
            w = frame_width
            h = w * im.height / im.width
            out.append(RLImage(mi.group(1), width=w, height=h))
            out.append(Spacer(1, 5))
            i += 1
            continue
        s = l.strip()
        if l.startswith("### "):
            out.append(Paragraph(md_inline_to_rl(l[4:]), styles["h2"]))
        elif l.startswith("## "):
            out.append(Paragraph(md_inline_to_rl(l[3:]), styles["h1"]))
        elif l.startswith("# "):
            out.append(Paragraph(md_inline_to_rl(l[2:]), styles["title"]))
        elif s == "---":
            out.append(Spacer(1, 5))
        elif s.startswith(("**Figure", "**Supplementary")):
            out.append(Paragraph(md_inline_to_rl(s), styles["cap"]))
        elif s.startswith(("- ", "* ")):
            out.append(Paragraph("\u2022 " + md_inline_to_rl(s[2:]), styles["body"]))
        elif s:
            out.append(Paragraph(md_inline_to_rl(s), styles["body"]))
        i += 1
    return out


def si_doc(path, title="Supplementary Information", header="Supplementary Information"):
    from reportlab.platypus import BaseDocTemplate, PageTemplate, Frame

    class _Doc(BaseDocTemplate):
        def afterPage(self):
            c = self.canv
            c.saveState()
            c.setFont("Arial", 7.2)
            c.setFillColor(colors.HexColor("#666666"))
            c.drawString(20 * MM, 12 * MM, header)
            c.drawRightString(A4[0] - 20 * MM, 12 * MM, str(c.getPageNumber()))
            c.setStrokeColor(colors.HexColor("#dddddd"))
            c.setLineWidth(0.4)
            c.line(20 * MM, 15 * MM, A4[0] - 20 * MM, 15 * MM)
            c.restoreState()

    d = _Doc(path, pagesize=A4, leftMargin=20 * MM, rightMargin=20 * MM,
             topMargin=18 * MM, bottomMargin=20 * MM, title=title)
    d.addPageTemplates([PageTemplate(id="all", frames=[
        Frame(d.leftMargin, d.bottomMargin, d.width, d.height, id="n")])])
    return d
