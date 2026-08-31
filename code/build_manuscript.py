"""Assemble Manuscript_with_figures.pdf for npj Precision Oncology submission."""
import re
import os
import pandas as pd
from PIL import Image
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import mm as MM
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import Paragraph, Spacer, Image as RLImage, PageBreak
from reportlab.lib import colors
import pypdf

from npjpdf import si_styles, md_inline_to_rl, md_flowables, rl_table, si_doc

mmN = open("manuscript_npj.md").read()
legN = open("figure_legends_npj.md").read()
refsN = open("references_npj.md").read()

TITLE = mmN.split("\n")[0].lstrip("# ").strip()
AUTHBLK = mmN.split("**Abstract**")[0].split("\n", 1)[1].strip()
ab = mmN.split("**Abstract**")[1].split("\n---")[0].strip()
after = mmN.split("**Abstract**")[1].split("\n---\n", 1)[1]
intro = after.split("## Results")[0].replace(ab, "").strip()
res = mmN.split("## Results")[1].split("## Discussion")[0]
dis = mmN.split("## Discussion")[1].split("## Methods")[0]
met = re.sub(r"\n---\s*$", "", mmN.split("## Methods")[1].split("## Data availability")[0].rstrip())
back = mmN[mmN.find("## Data availability"):mmN.find("## References")]
stat = mmN[mmN.find("## Acknowledgements"):]

figcaps, tabcaps = {}, {}
for blk in re.split(r"\n(?=\*\*Figure \d\.)", legN):
    m_ = re.match(r"\*\*Figure (\d)\.(.*?)\*\*(.*?)(?=\n\*\*(?:Figure|Table)|\Z)", blk, re.S)
    if m_:
        figcaps[int(m_.group(1))] = re.sub(
            r"\s*-{3,}\s*$", "",
            ("**Figure " + m_.group(1) + "." + m_.group(2) + "**" + m_.group(3)).strip())
for blk in re.split(r"\n(?=\*\*Table \d\.)", legN):
    m_ = re.match(r"\*\*Table (\d)\.(.*?)\*\*(.*?)(?=\n\*\*(?:Figure|Table)|\Z)", blk, re.S)
    if m_:
        tabcaps[int(m_.group(1))] = re.sub(
            r"\s*-{3,}\s*$", "",
            ("**Table " + m_.group(1) + "." + m_.group(2) + "**" + m_.group(3)).strip())

t1c = pd.read_csv("Table1_baseline_characteristics.csv").fillna("")
rows1 = [list(t1c.columns)] + t1c.astype(str).values.tolist()
t2d = pd.read_csv("Table2_interaction_summary.csv")
rows2 = [list(t2d.columns)] + t2d.astype(str).values.tolist()

S, FW = si_styles()
S["body"].fontSize = 9.6
S["body"].leading = 14.0
S["cap"].fontSize = 8.6
S["cap"].leading = 12.0
S["h1"].fontSize = 13.5
S["h2"].fontSize = 11.0
INK = colors.HexColor("#1a1a1a")
S["ptitle"] = ParagraphStyle("ptitle", fontName="Arial-B", fontSize=15.0, leading=20.0,
                             spaceAfter=8, textColor=INK)
S["auth"] = ParagraphStyle("auth", fontName="Arial", fontSize=9.2, leading=13.0,
                           spaceAfter=8, textColor=INK)
S["absh"] = ParagraphStyle("absh", fontName="Arial-B", fontSize=10.5, leading=14,
                           spaceAfter=4, textColor=INK)
S["ref"] = ParagraphStyle("ref", fontName="Arial", fontSize=8.2, leading=11.4,
                          spaceAfter=2.4, textColor=INK)


def build_ms(path="Manuscript_with_figures.pdf"):
    st = [Paragraph(md_inline_to_rl(TITLE), S["ptitle"])]
    for para in AUTHBLK.split("\n\n"):
        st.append(Paragraph(md_inline_to_rl(para.replace("\n", " ")), S["auth"]))
    st += [Paragraph("npj Precision Oncology \u2014 manuscript for submission", S["note"]),
           Spacer(1, 8),
           Paragraph("Abstract", S["absh"]),
           Paragraph(md_inline_to_rl(ab), S["body"]), Spacer(1, 10)]
    st += md_flowables(intro, S, FW)
    st.append(PageBreak())
    st.append(Paragraph("Results", S["h1"]))
    st += md_flowables(res, S, FW)
    st.append(Paragraph("Discussion", S["h1"]))
    st += md_flowables(dis, S, FW)
    for n in [1, 2, 3, 4]:
        st.append(PageBreak())
        p = f"fig/Figure{n}.png"
        im = Image.open(p)
        w = FW
        h = w * im.height / im.width
        maxh = A4[1] - 2 * 20 * MM - (len(figcaps[n].split()) / 13.0 * 12 + 26)
        if h > maxh:
            w = w * maxh / h
            h = maxh
        st += [RLImage(p, width=w, height=h), Spacer(1, 5),
               Paragraph(md_inline_to_rl(figcaps[n].replace("\n", " ")), S["cap"])]
    st.append(PageBreak())
    st += [Paragraph(md_inline_to_rl(tabcaps[1].replace("\n", " ")), S["cap"]), Spacer(1, 4),
           rl_table(rows1, [FW * x / 100 for x in [46, 19, 19, 16]], font_size=7.0),
           Spacer(1, 14),
           Paragraph(md_inline_to_rl(tabcaps[2].replace("\n", " ")), S["cap"]), Spacer(1, 4),
           rl_table(rows2, [FW * x / 100 for x in [16, 20, 19, 15, 8, 14, 8]], font_size=6.6)]
    st.append(PageBreak())
    st.append(Paragraph("Methods", S["h1"]))
    st += md_flowables(met, S, FW)
    st += md_flowables(back, S, FW)
    st += md_flowables(stat, S, FW)
    st.append(PageBreak())
    st.append(Paragraph("References", S["h1"]))
    for l in refsN.split("\n"):
        if re.match(r"^\d+\.\s", l.strip()):
            st.append(Paragraph(md_inline_to_rl(l.strip()), S["ref"]))
    si_doc(path, title=TITLE, header="Manuscript \u2014 npj Precision Oncology").build(st)
    return len(pypdf.PdfReader(path).pages)


npg = build_ms()
print("pages:", npg, "| size MB:", round(os.path.getsize("Manuscript_with_figures.pdf") / 1e6, 2))
