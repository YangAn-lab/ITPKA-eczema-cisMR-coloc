#!/usr/bin/env python3
"""P2-02 Figure 4 (v3): tissue-localization evidence summary.
Panel A: ITPKA log2FC forest plot (unchanged from v2).
Panel B: evidence layers as a colour-coded heat strip (replaces embedded table).
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import numpy as np
from scipy.stats import norm

matplotlib.rcParams['font.family'] = ['Liberation Sans', 'Arimo', 'DejaVu Sans']
matplotlib.rcParams['svg.fonttype'] = 'none'

contrasts = [
    ("GSE121212  Lesional vs non-lesional (paired, n=48)",      0.013, 0.958, "GSE121212"),
    ("GSE121212  Lesional vs healthy control (49 vs 38)",       0.196, 0.449, "GSE121212"),
    ("GSE157194  Lesional vs non-lesional, month 0 (paired, n=54)", -0.068, 0.648, "GSE157194"),
    ("GSE157194  Dupilumab: lesional month 3 vs 0 (paired, n=21)",  -0.177, 0.483, "GSE157194"),
    ("GSE157194  Ciclosporin: lesional month 3 vs 0 (paired, n=8)",  0.130, 0.768, "GSE157194"),
]
rows = []
for label, lfc, p, ds in contrasts:
    z = norm.isf(p / 2.0)
    se = abs(lfc) / z if z > 0 else np.nan
    rows.append(dict(label=label, lfc=lfc, p=p, ds=ds, se=se,
                     lo=lfc - 1.96 * se, hi=lfc + 1.96 * se))
colors = {"GSE121212": "#0072B2", "GSE157194": "#E69F00"}

fig = plt.figure(figsize=(15.5, 5.4))

# ---------------- Panel A (unchanged) ----------------
axA = fig.add_axes([0.265, 0.13, 0.385, 0.74])
ys = np.arange(len(rows))[::-1] * 1.0
for y, r in zip(ys, rows):
    c = colors[r["ds"]]
    axA.plot([r["lo"], r["hi"]], [y, y], color=c, lw=1.8, zorder=2)
    for xcap in (r["lo"], r["hi"]):
        axA.plot([xcap, xcap], [y - 0.13, y + 0.13], color=c, lw=1.8, zorder=2)
    axA.scatter([r["lfc"]], [y], s=95, color=c, zorder=3)
x_pcol = 1.18
axA.text(x_pcol, len(rows) - 1 + 0.62, "p value", fontsize=10, fontweight="bold",
         ha="left", va="center", color="#333333")
for y, r in zip(ys, rows):
    axA.text(x_pcol, y, f"{r['p']:.3f}", va="center", ha="left", fontsize=10, color="#333333")
axA.axvline(0, color="black", ls="--", lw=1.1, zorder=1)
axA.set_yticks(ys)
axA.set_yticklabels([r["label"] for r in rows], fontsize=10.5)
axA.set_xlim(-0.95, 1.55)
axA.set_xticks([-0.75, -0.50, -0.25, 0, 0.25, 0.50, 0.75, 1.00])
axA.set_xticklabels(["-0.75", "-0.50", "-0.25", "0", "0.25", "0.50", "0.75", "1.00"])
axA.set_xlabel("ITPKA log2 fold-change (95% CI)", fontsize=11.5)
axA.spines[["top", "right"]].set_visible(False)
axA.set_ylim(-0.7, len(rows) - 1 + 1.15)
axA.text(-0.95, len(rows) - 1 + 0.62,
         "GSE121212 paired Wilcoxon cross-check: p = 0.471 (n.s.)",
         fontsize=9.5, style="italic", color="#555555", ha="left", va="center",
         bbox=dict(facecolor="white", edgecolor="none", pad=1.5), zorder=5)
axA.text(0.30, -0.58, "●", fontsize=10, color=colors["GSE121212"], ha="left", va="center")
axA.text(0.36, -0.58, "GSE121212 (Tsoi 2019)", fontsize=9, color="#555555", ha="left", va="center")
axA.text(0.95, -0.58, "●", fontsize=10, color=colors["GSE157194"], ha="left", va="center")
axA.text(1.01, -0.58, "GSE157194 (Möbus 2021)", fontsize=9, color="#555555", ha="left", va="center")
axA.set_title("A", loc="left", fontsize=17, fontweight="bold", x=-0.02)

# ---------------- Panel B: evidence heat strip ----------------
axB = fig.add_axes([0.685, 0.13, 0.305, 0.74])
axB.axis("off")

# verdict: +1 supportive of blood/immune mechanism (blue), 0 negative/absent in skin (warm grey)
layers = [
    ("Blood eQTL colocalization", "eQTLGen + GTEx v8",        "PP.H4 > 0.96 in both",        +1),
    ("Immune-cell enrichment",    "HPA",                      "T-cell lineage / MAIT",       +1),
    ("Single-cell (skin)",        "CELLxGENE (672k cells)",   "Memory/helper T; absent in KC", +1),
    ("Skin cis-eQTL",             "GTEx v8/v10 skin",         "Absent; PP.H4 ≤ 0.21",        0),
    ("Skin RNA",                  "HPA / GTEx",               "0.1 nTPM; 0.74–0.83 TPM",     0),
    ("Skin protein",              "HPA",                      "Not detected (4 cell types)", 0),
    ("Lesional expression",       "GSE121212 + GSE157194",    "Null in 5/5 contrasts",       0),
]
C_POS, C_NEG = "#0279EE", "#C9C4BC"
n = len(layers)
axB.set_xlim(0, 1); axB.set_ylim(-1.75, n + 0.6)
for i, (layer, res, result, v) in enumerate(layers):
    y = n - 1 - i
    c = C_POS if v == +1 else C_NEG
    axB.add_patch(FancyBboxPatch((0.0, y + 0.12), 0.045, 0.76, boxstyle="round,pad=0.008",
                                 facecolor=c, edgecolor="none"))
    axB.text(0.075, y + 0.62, layer, fontsize=9.6, fontweight="bold", va="center", color="#1a1a1a")
    axB.text(0.075, y + 0.18, f"{res} — {result}", fontsize=8.4, va="center", color="#444444")
# legend
axB.add_patch(FancyBboxPatch((0.0, -0.72), 0.045, 0.5, boxstyle="round,pad=0.008",
                             facecolor=C_POS, edgecolor="none"))
axB.text(0.075, -0.47, "supports circulating immune-cell mechanism", fontsize=8.4, va="center", color="#333333")
axB.add_patch(FancyBboxPatch((0.0, -1.42), 0.045, 0.5, boxstyle="round,pad=0.008",
                             facecolor=C_NEG, edgecolor="none"))
axB.text(0.075, -1.17, "no skin-intrinsic signal detected", fontsize=8.4, va="center", color="#333333")
axB.set_title("B", loc="left", fontsize=17, fontweight="bold", x=-0.02)

fig.savefig("/workspace/P2-02_Figure4_tissue_localization_v3.png", dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig("/workspace/P2-02_Figure4_tissue_localization_v3.svg", bbox_inches="tight", facecolor="white")
print("saved figure 4 v3")
