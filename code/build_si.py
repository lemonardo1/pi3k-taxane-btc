"""Assemble Supplementary_Information.pdf for npj Precision Oncology submission."""
import re
import os
import pandas as pd
from PIL import Image
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import mm as MM
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import Paragraph, Spacer, Image as RLImage, PageBreak, KeepTogether
from reportlab.lib import colors
import pypdf

from npjpdf import si_styles, md_inline_to_rl, md_flowables, rl_table, si_doc

sfN = open("supplementary_figures_tables_npj.md").read()
si = open("supplementary_information_npj.md").read()
remark = open("REMARK_checklist.md").read()
smfull = open("supplementary_methods.md").read()

sm_body = re.sub(r"^## S1\.\d+ ", "### ", smfull.split("---", 1)[1].strip(),
                 flags=re.M).split("### Analyses attempted and not completed")[0].strip()
note1 = si.split("## Supplementary Note S1: analyses attempted and not completed")[1].split("---")[0].strip()
datinv = re.split(r"\n---\s*(?:\n|$)",
                  si.split("## Supplementary Data inventory", 1)[1], maxsplit=1)[0].strip()

figleg, tabcap_s = {}, {}
for blk in re.split(r"\n(?=\*\*Supplementary Figure S\d\.)", sfN):
    m_ = re.match(r"\*\*Supplementary Figure (S\d)\.(.*?)\*\*(.*?)(?=\n\*\*Supplementary|\Z)", blk, re.S)
    if m_:
        figleg[m_.group(1)] = re.sub(
            r"\s*-{3,}\s*$", "",
            ("**Supplementary Figure " + m_.group(1) + "." + m_.group(2) + "**" + m_.group(3)).strip())
for blk in re.split(r"\n(?=\*\*Supplementary Table S\d\.)", sfN):
    m_ = re.match(r"\*\*Supplementary Table (S\d)\.(.*?)\*\*(.*?)(?=\n\*\*Supplementary|\Z)", blk, re.S)
    if m_:
        tabcap_s[m_.group(1)] = re.sub(
            r"\s*-{3,}\s*$", "",
            ("**Supplementary Table " + m_.group(1) + "." + m_.group(2) + "**" + m_.group(3)).strip())

LBL = {"KRAS.1": "KRAS", "Cell_cycle": "Cell cycle", "RTK_RAS": "RTK\u2013RAS", "TGFB": "TGF-\u03b2",
       "ERBB2_AMP": "ERBB2 (amplification)", "ERBB2_SNV": "ERBB2 (mutation)",
       "PI3K_cont": "PI3K continuous", "PI3K_cont_full": "PI3K continuous, unshrunk",
       "PI3K_cont_k2": "PI3K continuous, k = 2", "PI3K_cont_k4": "PI3K continuous, k = 4",
       "PI3K(binary)": "PI3K", "A_full": "Fully adjusted", "B_reduced": "Reduced",
       "CHA_binary_only": "Unadjusted"}


def tidy_perm(d):
    key = "pathway" if "pathway" in d.columns else "gene"
    o = d[[key, "n_alt", "endpoint", "abs_z", "p_asymptotic", "p_perm_marginal", "p_bonferroni",
           "q_bh", "perm_fwer_p_maxT", "perm_fwer_p_minP", "sens_markerperm_maxT",
           "sens_ipw_maxT"]].copy()
    o.columns = ["Marker", "Altered n", "Endpoint", "|z|", "P asympt.", "P perm.", "P Bonf.",
                 "q BH", "FWER max-T", "FWER min-P", "Sens. marker-perm", "Sens. IPTW"]
    o["Marker"] = o["Marker"].replace(LBL)
    for c in o.columns[3:]:
        o[c] = o[c].map(lambda v: f"{v:.3g}" if pd.notna(v) else "")
    return o


def tidy_sens(d):
    o = d[["analysis", "spec", "score", "weighting", "n", "HR_interaction", "CI_low", "CI_high",
           "p"]].copy()
    o["analysis"] = o["analysis"].str.replace(r"^\((?:a|b|c|d|e)\)\s*", "", regex=True)
    o.columns = ["Analysis", "Specification", "Score", "Weighting", "n", "Interaction HR",
                 "CI low", "CI high", "P"]
    for c in ["Specification", "Score", "Weighting", "Analysis"]:
        o[c] = o[c].replace(LBL)
    o["95% CI"] = o.apply(lambda r: f"{r['CI low']:.2f}\u2013{r['CI high']:.2f}", axis=1)
    o["Interaction HR"] = o["Interaction HR"].map(lambda v: f"{v:.2f}")
    o["P"] = o["P"].map(lambda v: f"{v:.4g}")
    return o[["Analysis", "Specification", "Score", "Weighting", "n", "Interaction HR",
              "95% CI", "P"]]


tabs = {}
tabs["S1"] = pd.read_csv("TableS1_public_cohort_inventory.csv")
tabs["S2"] = pd.read_csv("TableS2_pathway_and_gene_screens.csv")
tabs["S3"] = pd.read_csv("TableS3_pik3ca_structural_metrics.csv")
tabs["S4"] = tidy_perm(pd.read_csv("permutation_pathway_screen.csv"))
tabs["S5"] = tidy_perm(pd.read_csv("permutation_gene_screen.csv"))
tabs["S6"] = tidy_sens(pd.read_csv("sensitivity_analyses.csv"))
tabs["S7"] = pd.read_csv("TableS7_baseline_by_treatment_arm.csv").fillna("")

d2 = tabs["S2"].copy()
d2["Feature"] = d2["Feature"].replace(LBL)
d2["95% CI"] = d2.apply(
    lambda r: f"{r['CI low']:.2f}\u2013{r['CI high']:.2f}" if pd.notna(r["CI low"]) else "", axis=1)
d2["Interaction HR"] = d2["Interaction HR"].map(lambda v: f"{v:.2f}")
d2["P"] = d2["P"].map(lambda v: f"{v:.4g}")
d2["q (BH)"] = d2["q (Benjamini\u2013Hochberg)"].map(lambda v: f"{v:.3g}" if pd.notna(v) else "")
tabs["S2"] = d2[["Level", "Feature", "Endpoint", "Altered patients", "Interaction HR", "95% CI",
                 "P", "q (BH)"]]

S, FW = si_styles()
S["body"].fontSize = 9.6
S["body"].leading = 14.0
S["cap"].fontSize = 8.6
S["cap"].leading = 12.0
S["h1"].fontSize = 13.5
S["h2"].fontSize = 11.0
INK = colors.HexColor("#1a1a1a")

TITLE = ("PI3K pathway alterations and attenuated taxane benefit in advanced biliary tract "
         "cancer")

toc = [["Section", "Contents"],
       ["Supplementary Figures",
        "S1 propensity overlap \u00b7 S2 cutpoint and discrimination \u00b7 S3 alteration landscape "
        "and gene weights \u00b7 S4 saturation background and interface geometry \u00b7 S5 "
        "permutation null distributions"],
       ["Supplementary Tables",
        "S1 public cohort inventory \u00b7 S2 pathway and gene screens \u00b7 S3 PIK3CA structural "
        "metrics \u00b7 S4\u2013S5 permutation-corrected screens \u00b7 S6 robustness and "
        "sensitivity analyses \u00b7 S7 baseline by treatment arm"],
       ["Supplementary Note S1", "Analyses attempted and not completed"],
       ["Supplementary Methods",
        "Cohorts, treatment annotation, variant scoring, structural and genomic axes, statistics"],
       ["REMARK checklist", "Twenty-item reporting checklist"],
       ["Supplementary Data inventory", "Derived files deposited at the repository"]]

WIDTHS = {
    "S1": [28, 20, 13, 17, 26, 13, 23],
    "S2": [13, 20, 14, 18, 18, 20, 15, 13],
    "S3": [15, 11, 17, 19, 21, 19, 21, 19, 21, 19, 19, 17],
    "S4": [18, 12, 13, 10, 15, 15, 14, 12, 15, 15, 17, 14],
    "S5": [18, 12, 13, 10, 15, 15, 14, 12, 15, 15, 17, 14],
    "S6": [40, 17, 21, 15, 9, 15, 17, 13],
    "S7": [22, 13, 38, 15, 15, 12],
}


def build_si(path="Supplementary_Information.pdf"):
    st = [
        Paragraph("Supplementary Information", S["title"]),
        Paragraph(md_inline_to_rl(TITLE), ParagraphStyle(
            "st", fontName="Arial", fontSize=10.5, leading=14, textColor=INK, spaceAfter=14)),
        rl_table(toc, [34 * MM, FW - 34 * MM], font_size=7.6),
        PageBreak(),
        Paragraph("Supplementary Figures", S["h1"]),
    ]
    for k in ["S1", "S2", "S3", "S4", "S5"]:
        p = f"fig/Figure{k}.png"
        im = Image.open(p)
        w = FW
        h = w * im.height / im.width
        if h > 105 * MM:
            w = w * 105 * MM / h
            h = 105 * MM
        st.append(KeepTogether([
            RLImage(p, width=w, height=h),
            Spacer(1, 4),
            Paragraph(md_inline_to_rl(figleg[k].replace("\n", " ")), S["cap"]),
        ]))
        st.append(Spacer(1, 10))

    st += [PageBreak(), Paragraph("Supplementary Tables", S["h1"])]
    for k in ["S1", "S2", "S3", "S4", "S5", "S6", "S7"]:
        d = tabs[k]
        st.append(Paragraph(md_inline_to_rl(tabcap_s[k].replace("\n", " ")), S["cap"]))
        st.append(Spacer(1, 3))
        rows = [list(d.columns)] + d.fillna("").astype(str).values.tolist()
        raw = WIDTHS[k]
        tot = sum(raw)
        fs = 5.5 if k in ("S3", "S4", "S5") else (6.2 if k == "S2" else 7.0)
        st.append(rl_table(rows, [FW * x / tot for x in raw], font_size=fs))
        st.append(Spacer(1, 11))

    st += [PageBreak(),
           Paragraph("Supplementary Note S1: analyses attempted and not completed", S["h1"])]
    st += md_flowables(note1, S, FW)
    st += [PageBreak(), Paragraph("Supplementary Methods", S["h1"])]
    st += md_flowables(sm_body, S, FW)
    st += [PageBreak()]
    st += md_flowables(remark, S, FW)
    st += [PageBreak(), Paragraph("Supplementary Data inventory", S["h1"])]
    st += md_flowables(datinv, S, FW)

    si_doc(path, title="Supplementary Information",
           header="Supplementary Information").build(st)
    return len(pypdf.PdfReader(path).pages)


nsi = build_si()
print("pages:", nsi, "| size MB:", round(os.path.getsize("Supplementary_Information.pdf") / 1e6, 2))
