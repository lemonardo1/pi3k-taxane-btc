"""Render Cover_letter.pdf for the npj Precision Oncology submission."""
import os
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import mm as MM
from reportlab.platypus import Spacer
import pypdf

from npjpdf import si_styles, md_flowables, si_doc

cl = open("cover_letter_npj.md").read()


def build_cover(path="Cover_letter.pdf"):
    S, FW = si_styles()
    st = [Spacer(1, 4 * MM)]
    st += md_flowables(cl, S, FW)
    si_doc(path, title="Cover letter — npj Precision Oncology",
           header="Cover letter \u2014 npj Precision Oncology").build(st)
    return len(pypdf.PdfReader(path).pages)


npg = build_cover()
print("pages:", npg, "| size KB:", round(os.path.getsize("Cover_letter.pdf") / 1024, 1))
