import matplotlib
matplotlib.rcParams['font.family'] = ['Liberation Sans', 'Arimo', 'DejaVu Sans']
matplotlib.rcParams['svg.fonttype'] = 'none'
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

BLUE, ORANGE, GREY = '#0279EE', '#FF9400', '#666666'

# Panel A data (from P2-02_GEO皮损表达验证_结果表.csv; Wilcoxon row has no log2FC)
rows = [
    ("GSE121212  Lesional vs non-lesional (paired, n=48)", 0.013, 0.958, BLUE),
    ("GSE121212  Lesional vs healthy control (49 vs 38)",   0.196, 0.449, BLUE),
    ("GSE157194  Lesional vs non-lesional, month 0 (paired, n=54)", -0.068, 0.648, ORANGE),
    ("GSE157194  Dupilumab: lesional month 3 vs 0 (paired, n=21)",  -0.177, 0.483, ORANGE),
    ("GSE157194  Ciclosporin: lesional month 3 vs 0 (paired, n=8)",   0.130, 0.768, ORANGE),
]

fig = plt.figure(figsize=(13.5, 5.4))
gs = GridSpec(1, 2, width_ratios=[1.05, 1.35], wspace=0.30)

# ---------- Panel A ----------
ax = fig.add_subplot(gs[0])
ys = list(range(len(rows), 0, -1))
for (lab, fc, p, c), y in zip(rows, ys):
    ax.scatter(fc, y, s=64, color=c, zorder=3)
    ax.annotate(f"p = {p:.3f}", (fc, y), xytext=(8, 0), textcoords='offset points',
                va='center', fontsize=8.5, color='#333333')
ax.axvline(0, color='black', lw=0.9, ls='--', zorder=1)
ax.set_yticks(ys)
ax.set_yticklabels([r[0] for r in rows], fontsize=8.8)
ax.set_xlim(-0.42, 0.42)
ax.set_xlabel("ITPKA log2 fold-change", fontsize=10)
ax.set_ylim(0.3, 6.4)
ax.text(0.0, 6.05, "GSE121212 paired Wilcoxon cross-check: p = 0.471 (n.s.)",
        ha='center', fontsize=8.2, color=GREY, style='italic')
ax.spines[['top', 'right']].set_visible(False)
ax.set_title("A", loc='left', fontsize=14, fontweight='bold', x=-0.28)

# ---------- Panel B: evidence summary table ----------
axb = fig.add_subplot(gs[1])
axb.axis('off')
tbl_rows = [
    ["Blood eQTL colocalization", "eQTLGen + GTEx v8 blood", "PP.H4 > 0.96 in both resources"],
    ["Immune-cell enrichment", "Human Protein Atlas", "T-cell lineage / MAIT enriched"],
    ["Skin cis-eQTL (rs11635906)", "GTEx v8/v10, sun-exp. + non-sun-exp.", "Absent; coloc PP.H4 \u2264 0.21"],
    ["Skin RNA", "HPA / GTEx", "Very low (0.1 nTPM / 0.74\u20130.83 TPM)"],
    ["Skin protein", "Human Protein Atlas", "Not detected in 4 skin cell types"],
    ["Lesional-skin expression", "GSE121212 + GSE157194", "Null in 6/6 contrasts"],
]
tbl = axb.table(cellText=[["Evidence layer", "Resource", "Result"]] + tbl_rows,
                colWidths=[0.34, 0.36, 0.44], loc='center', cellLoc='left')
tbl.auto_set_font_size(False)
tbl.set_fontsize(9.2)
tbl.scale(1, 1.75)
for j in range(3):
    c = tbl[0, j]
    c.set_facecolor('#ECE9E2'); c.set_text_props(fontweight='bold')
for i in range(1, 7):
    for j in range(3):
        tbl[i, j].set_facecolor('#FAF9F3' if i % 2 else 'white')
axb.set_title("B", loc='left', fontsize=14, fontweight='bold', x=-0.06)

fig.savefig('/mnt/results/P2-02_fig4_组织定位证据汇总.png', dpi=300, bbox_inches='tight')
fig.savefig('/mnt/results/P2-02_fig4_组织定位证据汇总.svg', bbox_inches='tight')
plt.show()
print("saved")