"""Build the co-author review PDF v3: manuscript with ALL figures/tables embedded at
first-citation points, plus ALL supplementary tables S1-S8 rendered in full.
Pipeline: parse Genomics submission md -> move tables/figure legends inline ->
append supplementary CSVs as (split where wide) HTML tables -> markdown -> HTML -> WeasyPrint PDF.
Output: /workspace/P2-02_manuscript_review_v3.pdf
"""
import re
import markdown
import pandas as pd

SRC = "/workspace/P2-02_manuscript_v4.md"
OUT = "/workspace/P2-02_manuscript_review_v4.pdf"
IMG = "/mnt/results"  # base_url for figure PNGs
DATE = "2026-09-16"

FIG_FILES = {
    "Figure 1": "P2-02_Figure1_study_design_v2.png",
    "Figure 2": "P2-02_Figure2_MR_coloc_main_v6.png",
    "Figure 3": "P2-02_Figure3_locus_discrimination_v5.png",
    "Figure 4": "P2-02_Figure4_tissue_localization_v3.png",
    "Figure 5": "P2-02_Figure5_mechanism_hypothesis_v6.png",
    "Figure S1": "P2-02_FigureS1_singlecell_skin.png",
    "Figure S2": "P2-02_FigureS2_conditional_locus.png",
}

RES = "/mnt/results"
S_CSV = {k: f"{RES}/P2-02_SupplementaryTable{k}.csv" for k in [
    "S1_dilution_calculations", "S2_PheWAS_rs11635906", "S3_scan498",
    "S4_panel12_coloc", "S5_SMR_HEIDI", "S6_cisMVMR",
    "S7_immune_cell_eQTL_survey", "S8_PheWAS_trait_coloc",
    "S9_plasma_pQTL_coloc", "S10_locus_1Mb_gene_panel", "S11_pQTL_replication",
    "S12_FDR_hit_coloc", "S13_FIN_LD_sensitivity", "S14_extended_skin_endpoints",
    "S15_regulatory_annotation", "S16_ITPKA_vs_ITPKB", "S17_IL6ST_ANKRD55_discrimination"]}

lines = open(SRC, encoding="utf-8").read().split("\n")

# ---------- locate zones ----------
def find_line(prefix, start=0):
    for i in range(start, len(lines)):
        if lines[i].strip() == prefix:
            return i
    raise ValueError(prefix)

i_abs = find_line("## Abstract")
i_tables = find_line("## Tables")
i_legends = find_line("## Figure legends")
i_supp = find_line("## Supplementary tables")

title_lines = lines[:i_abs]
body_lines = lines[i_abs:i_tables]
tables_zone = lines[i_tables + 1:i_legends]
legends_zone = lines[i_legends + 1:i_supp]
supp_zone = lines[i_supp + 1:]

# ---------- parse the 4 main tables ----------
def parse_tables(zone):
    blocks, cur = [], None
    for ln in zone:
        s = ln.strip()
        if s.startswith("**Table "):
            if cur:
                blocks.append(cur)
            cur = {"title": s.strip("*"), "pipe": [], "note": []}
        elif cur is None:
            continue
        elif s.startswith("|"):
            cur["pipe"].append(ln.rstrip())
        elif s and s != "---":
            cur["note"].append(s)
    if cur:
        blocks.append(cur)
    return blocks

tables = {b["title"].split(".")[0].replace("Table ", ""): b for b in parse_tables(tables_zone)}
assert set(tables) == {"1", "2", "3", "4", "5", "6"}, tables.keys()

def table_md(num):
    b = tables[num]
    parts = [f'<p class="tbl-title"><strong>{b["title"]}</strong></p>', ""]
    parts += b["pipe"]
    note = " ".join(b["note"])
    if note:
        parts += ["", f'<p class="tbl-note">{note}</p>']
    return "\n".join(parts)

# ---------- parse figure legends ----------
def parse_legends(zone):
    out = {}
    for ln in zone:
        s = ln.strip()
        m = re.match(r"^\*\*(Figure (?:S?\d+)\..*?)\*\*\s*(.*)$", s)
        if m:
            out[f"Figure {m.group(1).split('.')[0].split(' ')[1]}"] = (m.group(1), m.group(2))
    return out

legends = parse_legends(legends_zone)
assert set(legends) == set(FIG_FILES), legends.keys()

def fig_html(key):
    head, rest = legends[key]
    return (f'<div class="figure"><img src="{FIG_FILES[key]}"/>'
            f'<p class="caption"><strong>{head}</strong> {rest}</p></div>')

# ---------- value formatting ----------
def fmt(v, kind=None):
    if v is None or (isinstance(v, float) and pd.isna(v)) or str(v) == "nan":
        return ""
    if kind == "p":
        try:
            x = float(v)
        except (TypeError, ValueError):
            return str(v)
        if x == 0:
            return "0"
        if abs(x) < 1e-3:
            mant, ex = f"{x:.1e}".split("e")
            return f"{mant}e{int(ex)}"
        return f"{x:.3f}".rstrip("0").rstrip(".")
    if kind == "f":
        try:
            return f"{float(v):.3g}"
        except (TypeError, ValueError):
            return str(v)
    if kind == "int":
        try:
            return str(int(float(v)))
        except (TypeError, ValueError):
            return str(v)
    return str(v)

def fmt_df(df, spec):
    out = df.copy()
    for c in out.columns:
        out[c] = out[c].map(lambda v, k=spec.get(c): fmt(v, k))
    return out

def tbl_html(df, cls="supp", weights=None):
    h = df.to_html(index=False, border=0, classes=cls, escape=True)
    if weights:  # fixed layout with proportional column widths
        tot = sum(weights)
        cols = "".join(f'<col style="width:{100*w/tot:.2f}%"/>' for w in weights)
        h = h.replace(">\n<thead", f'><colgroup>{cols}</colgroup>\n<thead', 1)
        h = h.replace('class="dataframe ' + cls + '"', 'class="dataframe ' + cls + ' fixed"')
    return h

# ---------- S1 (9x3, portrait) ----------
s1 = pd.read_csv(S_CSV["S1_dilution_calculations"])
s1_html = tbl_html(fmt_df(s1, {}))

# ---------- S2 (103x10, landscape) ----------
s2 = pd.read_csv(S_CSV["S2_PheWAS_rs11635906"])
s2_html = tbl_html(fmt_df(s2, {"eaf": "f", "beta": "f", "se": "f", "p": "p",
                               "n_cases": "int", "n_controls": "int"}), "supp s2")

# ---------- S3 (498x30 -> two landscape panels) ----------
s3 = pd.read_csv(S_CSV["S3_scan498"])
chr_s = pd.to_numeric(s3["chr_hg38"], errors="coerce").astype("Int64").astype(str)
pos_s = pd.to_numeric(s3["pos_hg38"], errors="coerce").astype("Int64").astype(str)
s3["locus"] = (chr_s + ":" + pos_s).str.replace("<NA>:<NA>", "", regex=False)
s3["OR (95% CI)"] = s3.apply(
    lambda r: (f"{fmt(r['or_mr'],'f')} ({fmt(r['or_lo'],'f')}–{fmt(r['or_hi'],'f')})"
               if pd.notna(r["or_mr"]) else ""), axis=1)
s3a_cols = ["pair_id", "gene", "exposure_id", "outcome_id", "outcome", "lead_rsid",
            "locus", "ea", "nea", "eaf_eqtl", "beta_eqtl", "se_eqtl", "p_eqtl"]
s3b_cols = ["pair_id", "beta_gwas", "se_gwas", "p_gwas", "af_outcome", "b_mr",
            "se_mr_fo", "p_mr", "OR (95% CI)", "se_mr_delta", "p_mr_delta",
            "q_mr_BH", "pass_bonferroni", "status", "coloc_triage"]
s3a = fmt_df(s3[s3a_cols], {"eaf_eqtl": "f", "beta_eqtl": "f", "se_eqtl": "f", "p_eqtl": "p"})
s3b = fmt_df(s3[s3b_cols], {"beta_gwas": "f", "se_gwas": "f", "p_gwas": "p",
                            "af_outcome": "f", "b_mr": "f", "se_mr_fo": "f", "p_mr": "p",
                            "se_mr_delta": "f", "p_mr_delta": "p", "q_mr_BH": "p"})
s3_html = ('<p class="tbl-note"><em>Part 1 of 2 — pair definition and eQTL instrument. '
           'locus = chr_hg38:pos_hg38.</em></p>'
           + tbl_html(s3a, "supp wide")
           + '<p class="tbl-note" style="page-break-before:always"><em>Part 2 of 2 — GWAS and '
             'MR results (same row order as part 1; pair_id is the join key). OR (95% CI) merges '
             'or_mr/or_lo/or_hi.</em></p>'
           + tbl_html(s3b, "supp wide"))

# ---------- S4 (152x24 -> landscape, p12/window ranges compressed) ----------
s4 = pd.read_csv(S_CSV["S4_panel12_coloc"])
P12 = ["PP.H4_p12_1e6", "PP.H4_p12_5e6", "PP.H4_p12_1e5", "PP.H4_p12_5e5", "PP.H4_p12_1e4"]
s4["PP.H4 p12=1e-6…1e-4"] = s4[P12].apply(lambda r: ", ".join(fmt(v, "f") for v in r), axis=1)
s4["PP.H4 ±250/±125kb"] = s4[["PP.H4_w250", "PP.H4_w125"]].apply(
    lambda r: ", ".join(fmt(v, "f") for v in r), axis=1)
s4_cols = ["panel_origin", "dataset", "gene", "endpoint", "nsnps", "N1",
           "PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4",
           "PP.H4 p12=1e-6…1e-4", "PP.H4 ±250/±125kb",
           "min_p_eqtl", "min_p_gwas", "lead_eqtl", "lead_gwas", "median_tpm", "note"]
s4f = fmt_df(s4[s4_cols], {"nsnps": "int", "N1": "int", "PP.H0": "f", "PP.H1": "f",
                           "PP.H2": "f", "PP.H3": "f", "PP.H4": "f",
                           "min_p_eqtl": "p", "min_p_gwas": "p", "median_tpm": "f"})
s4_html = ('<p class="tbl-note"><em>The five p12-prior and two window-sensitivity PP.H4 columns '
           'of the CSV are shown as ordered lists (p12 = 1e-6, 5e-6, 1e-5, 5e-5, 1e-4; '
           'windows ±250, ±125 kb).</em></p>'
           + tbl_html(s4f, "supp wide",
                      weights=[6, 8, 6, 5, 4, 4, 4.2, 4.2, 4.2, 4.2, 4.2, 9, 6, 5, 5, 7, 7, 4, 7]))

# ---------- S5 (24x20, landscape) ----------
s5 = pd.read_csv(S_CSV["S5_SMR_HEIDI"])
s5f = fmt_df(s5, {"Freq": "f", "b_GWAS": "f", "se_GWAS": "f", "p_GWAS": "p",
                  "b_eQTL": "f", "se_eQTL": "f", "p_eQTL": "p", "b_SMR": "f",
                  "se_SMR": "f", "p_SMR": "p", "p_HEIDI": "p", "nsnp_HEIDI": "int"})
s5_html = tbl_html(s5f, "supp wide")

# ---------- S6 (2x35 -> transposed portrait) ----------
s6 = pd.read_csv(S_CSV["S6_cisMVMR"])
s6t = s6.set_index("endpoint").T.reset_index().rename(columns={"index": "metric"})
s6t.columns.name = None
for c in s6t.columns[1:]:
    s6t[c] = s6t[c].map(lambda v: fmt(v, "p") if isinstance(v, float) and abs(v) < 1e-3 and v != 0
                        else (fmt(v, "f") if isinstance(v, float) else str(v)))
s6_html = tbl_html(s6t, "supp")

# ---------- S7 (143x31 -> two landscape panels) ----------
s7 = pd.read_csv(S_CSV["S7_immune_cell_eQTL_survey"])
s7["PP.H4 p12=1e-6…1e-4"] = s7[P12].apply(lambda r: ", ".join(fmt(v, "f") for v in r), axis=1)
s7["PP.H4 ±250/±125kb"] = s7[["PP.H4_w250", "PP.H4_w125"]].apply(
    lambda r: ", ".join(fmt(v, "f") for v in r), axis=1)
s7a_cols = ["gene", "study", "cell", "platform", "endpoint", "nsnps", "N1", "median_tpm",
            "probe", "note", "rs11635906_in", "rs11635906_eqtl_p", "rs11635906_eqtl_beta"]
s7b_cols = ["gene", "study", "cell", "endpoint", "PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4",
            "PP.H4 p12=1e-6…1e-4", "PP.H4 ±250/±125kb", "min_p_eqtl", "min_p_gwas",
            "lead_eqtl_rsid", "lead_eqtl_p", "lead_gwas_rsid", "lead_gwas_p"]
s7a = fmt_df(s7[s7a_cols], {"nsnps": "int", "N1": "int", "median_tpm": "f",
                            "rs11635906_eqtl_p": "p", "rs11635906_eqtl_beta": "f"})
s7b = fmt_df(s7[s7b_cols], {"PP.H0": "f", "PP.H1": "f", "PP.H2": "f", "PP.H3": "f", "PP.H4": "f",
                            "min_p_eqtl": "p", "min_p_gwas": "p",
                            "lead_eqtl_p": "p", "lead_gwas_p": "p"})
s7_html = ('<p class="tbl-note"><em>Part 1 of 2 — dataset measurability and rs11635906-specific '
           'association.</em></p>' + tbl_html(s7a, "supp wide")
           + '<p class="tbl-note" style="page-break-before:always"><em>Part 2 of 2 — coloc.abf '
             'results (same row order; gene/study/cell/endpoint are the join keys). p12-prior and '
             'window columns shown as ordered lists as in Table S4.</em></p>'
           + tbl_html(s7b, "supp wide"))

# ---------- S8 (33x15, landscape; rebuilt: 9 traits x 3 endpoints + cystatin C attribution) ----------
s8 = pd.read_csv(S_CSV["S8_PheWAS_trait_coloc"])
s8f = fmt_df(s8, {c: ("p" if c.startswith("min_p") else "f")
                  for c in s8.columns if c not in ("analysis", "trait", "endpoint", "note")})
s8_html = tbl_html(s8f, "supp wide")

# ---------- S9-S17 (new in v4) ----------
s9 = pd.read_csv(S_CSV["S9_plasma_pQTL_coloc"])
s9f = fmt_df(s9, {c: ("p" if c.startswith("min_p") else "f")
                  for c in s9.columns if c.startswith(("PP.", "H4_", "min_p"))})
s9_html = tbl_html(s9f, "supp wide")

s10 = pd.read_csv(S_CSV["S10_locus_1Mb_gene_panel"])
s10f = fmt_df(s10, {"start": "int", "end": "int", "rs11635906_eqtl_beta": "f",
                    "rs11635906_eqtl_se": "f", "rs11635906_eqtl_p": "p", "rs11635906_eaf": "f",
                    **{c: "f" for c in s10.columns if c.startswith("coloc_")}})
s10_html = tbl_html(s10f, "supp wide")

s11 = pd.read_csv(S_CSV["S11_pQTL_replication"])
s11_html = tbl_html(fmt_df(s11, {"n": "int", "rs11635906_beta": "f", "rs11635906_se": "f",
                                 "rs11635906_p": "p", "window_min_p": "p", "window_nsnps": "int"}))

s12 = pd.read_csv(S_CSV["S12_FDR_hit_coloc"])
s12f = fmt_df(s12, {c: ("p" if c.startswith("min_p") else "f")
                    for c in s12.columns if c.startswith(("PP.", "H4_", "min_p"))})
s12_html = tbl_html(s12f, "supp wide")

s13 = pd.read_csv(S_CSV["S13_FIN_LD_sensitivity"])
s13f = fmt_df(s13, {c: "f" for c in s13.columns if c.startswith("PP.H")})
s13_html = tbl_html(s13f, "supp wide")

s14 = pd.read_csv(S_CSV["S14_extended_skin_endpoints"])
s14_html = tbl_html(fmt_df(s14, {"n_cases": "int", "beta_G": "f", "se": "f", "p": "p"}))

s15 = pd.read_csv(S_CSV["S15_regulatory_annotation"])
s15_html = tbl_html(fmt_df(s15, {}), "supp")

s16 = pd.read_csv(S_CSV["S16_ITPKA_vs_ITPKB"])
s16_html = tbl_html(fmt_df(s16, {"n_cells": "int", "n_pos": "int", "pct_pos": "f",
                                 "mean_expr": "f", "mean_pos": "f"}))

s17 = pd.read_csv(S_CSV["S17_IL6ST_ANKRD55_discrimination"])
s17_html = tbl_html(fmt_df(s17, {"nsnp": "int", "H3": "f", "H4": "f",
                                 "rs7731626_eq_p": "p", "rs7731626_gw_p": "p"}))

# ---------- walk body, insert blocks after anchor paragraphs ----------
ANCHORS = [
    ("An overview is shown in Figure 1", fig_html("Figure 1")),
    ("(Table 1, Figure 2)", table_md("1") + "\n\n" + fig_html("Figure 2")),
    ("PP.H4=0.962 for dermatitis/eczema (Table 2)", table_md("2")),
    ("twelve named genes are cis-regulated", fig_html("Figure 3") + "\n\n" + table_md("3")),
    ("(Figure S2)", fig_html("Figure S2")),
    ("beyond FinnGen (Table 4)", table_md("4")),
    ("(Table 5, Supplementary Table S6)", table_md("5")),
    ("Table 6, Supplementary Table S7)", table_md("6")),
    ("(Figure 4A)", fig_html("Figure 4")),
    ("(Figure S1)", fig_html("Figure S1")),
]

out, used = [], set()
for ln in body_lines:
    if ln.strip().startswith("### 4.4"):
        out.append(fig_html("Figure 5"))
        out.append("")
        used.add("fig5")
    out.append(ln)
    for j, (anchor, block) in enumerate(ANCHORS):
        if j not in used and anchor in ln:
            out += ["", block, ""]
            used.add(j)
assert len(used) == 11, used

# ---------- supplementary zone: text descriptions only (tables in separate folder) ----------
supp_out = ["## Supplementary tables", ""]
for ln in supp_zone:
    s = ln.strip()
    if not s or s == "---":
        continue
    supp_out.append(s)
supp_out += ["", "*All supplementary tables (S1–S17) are provided as separate CSV files in the review package folder; they are rendered in full in the appendix below.*", ""]

body_md = "\n".join(out) + "\n\n---\n\n" + "\n".join(supp_out)

# ---------- appendix: full supplementary table rendering ----------
SUPP_APP = [
    ("Supplementary Table S1. Quantitative dilution calculations", s1_html, False),
    ("Supplementary Table S2. Phenome-wide associations of rs11635906", s2_html, True),
    ("Supplementary Table S3. Full scan results (498 pairs)", s3_html, True),
    ("Supplementary Table S4. Extended twelve-gene panel colocalization", s4_html, True),
    ("Supplementary Table S5. SMR-HEIDI results", s5_html, True),
    ("Supplementary Table S6. cis multivariable MR results", s6_html, False),
    ("Supplementary Table S7. Purified immune-cell eQTL survey", s7_html, True),
    ("Supplementary Table S8. Non-atopic and blood-cell trait colocalization", s8_html, True),
    ("Supplementary Table S9. Plasma pQTL colocalization", s9_html, True),
    ("Supplementary Table S10. Exhaustive +/-1 Mb locus gene panel", s10_html, True),
    ("Supplementary Table S11. rs11635906 pQTL replication", s11_html, False),
    ("Supplementary Table S12. FDR-level scan hits: colocalization", s12_html, True),
    ("Supplementary Table S13. Finnish-matched LD sensitivity", s13_html, True),
    ("Supplementary Table S14. Extended FinnGen skin endpoints", s14_html, False),
    ("Supplementary Table S15. Variant-level regulatory annotation", s15_html, False),
    ("Supplementary Table S16. ITPKA versus ITPKB expression", s16_html, False),
    ("Supplementary Table S17. Positive-control locus discrimination at 5q11.2", s17_html, False),
]
app_parts = []
for t, h, land in SUPP_APP:
    if land:
        app_parts.append(f'<section class="landscape"><h3>{t}</h3>{h}</section>')
    else:
        app_parts.append(f'<h3 style="page-break-before:always">{t}</h3>{h}')
appendix_html = "\n".join(app_parts)
body_md = re.sub(r"(\[[A-Z][^\]\n]{3,}?\])", r"<mark>\1</mark>", body_md)

# ---------- title page ----------
title = title_lines[0].lstrip("# ").strip()
meta = [re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", t.strip()) for t in title_lines[1:] if t.strip()]
meta = [re.sub(r"(\[[^\]]+\])", r"<mark>\1</mark>", m) for m in meta]
titlepage = (
    '<div class="titlepage"><p class="draft-tag">DRAFT FOR CO-AUTHOR REVIEW<br/>'
    f'<span class="draft-sub">{DATE} &middot; not for distribution &middot; '
    'please send comments by page &amp; section number</span></p>'
    f'<h1 class="tp-title">{title}</h1>'
    + "".join(f'<p class="tp-meta">{m}</p>' for m in meta)
    + "</div>"
)

# ---------- render ----------
body_html = markdown.markdown(body_md, extensions=["tables"])
html_prefix = f"""<!DOCTYPE html><html><head><meta charset="utf-8"/><style>
@page {{ size: A4; margin: 2.2cm 2.2cm 2.4cm 2.2cm;
  @top-center {{ content: "Draft for co-author review \\2014 {DATE} (not for distribution)";
    font-family: 'Liberation Sans'; font-size: 8pt; color: #999; }}
  @bottom-center {{ content: "Page " counter(page) " of " counter(pages);
    font-family: 'Liberation Sans'; font-size: 8pt; color: #999; }} }}
@page :first {{ @top-center {{ content: none; }} }}
@page landscape {{ size: A4 landscape; }}
section.landscape {{ page: landscape; break-after: page; }}
body {{ font-family: 'Liberation Serif', serif; font-size: 11pt; line-height: 1.4; color: #111; }}
h1 {{ font-family: 'Liberation Sans'; font-size: 14pt; margin: 14pt 0 8pt; }}
h2 {{ font-family: 'Liberation Sans'; font-size: 12.5pt; margin: 14pt 0 6pt; page-break-after: avoid; }}
h3 {{ font-family: 'Liberation Sans'; font-size: 11pt; margin: 10pt 0 4pt; page-break-after: avoid; }}
p {{ margin: 5pt 0; text-align: justify; }}
.figure {{ text-align: center; margin: 12pt 0; page-break-inside: avoid; }}
.figure img {{ max-width: 100%; max-height: 21cm; }}
.caption {{ font-size: 9pt; color: #444; text-align: left; margin-top: 4pt; }}
table {{ border-collapse: collapse; width: 100%; font-size: 9pt; margin: 6pt 0; font-family: 'Liberation Sans'; }}
th, td {{ border: 0.5pt solid #999; padding: 2.5pt 4pt; text-align: left; vertical-align: top; }}
thead {{ display: table-header-group; }}
tr {{ page-break-inside: avoid; }}
.tbl-title {{ font-size: 10pt; margin: 10pt 0 2pt; page-break-after: avoid; }}
.tbl-note {{ font-size: 8.5pt; color: #444; text-align: left; }}
table.supp {{ font-size: 8pt; }}
table.supp.s2 {{ font-size: 7.5pt; }}
table.supp.wide {{ font-size: 6pt; }}
table.supp.wide th, table.supp.wide td {{ padding: 1.5pt 2.5pt; }}
table.supp.fixed {{ table-layout: fixed; }}
table.supp.fixed th, table.supp.fixed td {{ word-break: break-all; overflow-wrap: anywhere; }}
mark {{ background: #fff3a3; padding: 0 1pt; }}
.titlepage {{ page-break-after: always; text-align: center; padding-top: 150pt; }}
.draft-tag {{ display: inline-block; border: 1.2pt solid #b00; color: #b00; font-family: 'Liberation Sans';
  font-weight: bold; font-size: 12pt; padding: 8pt 16pt; letter-spacing: 1pt; }}
.draft-sub {{ font-size: 8.5pt; font-weight: normal; letter-spacing: 0; color: #833; }}
.tp-title {{ font-size: 17pt; line-height: 1.35; margin: 36pt 24pt 18pt; }}
.tp-meta {{ font-size: 10.5pt; color: #333; text-align: center; }}
</style></head><body>""" 
html = html_prefix + titlepage + body_html + appendix_html + "</body></html>"

open("/workspace/review_v3.html", "w", encoding="utf-8").write(html)

from weasyprint import HTML
HTML(string=html, base_url=IMG + "/").write_pdf(OUT)
print("saved", OUT)

import pypdf
r = pypdf.PdfReader(OUT)
print("pages:", len(r.pages))
