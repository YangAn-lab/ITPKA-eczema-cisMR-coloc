#!/usr/bin/env python3
"""P2-02 Figure 1 v2: study-design flowchart, pure matplotlib (Elsevier-compliant).
Reflects the corrected scan: 4 Bonferroni hits, IL6ST positive control, ITPKA novel nomination."""
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['font.family'] = ['Liberation Sans', 'Arimo', 'DejaVu Sans']
matplotlib.rcParams['svg.fonttype'] = 'none'
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

BLUE, GREEN, ORANGE = '#0279EE', '#75A025', '#FF9400'
NEUT, NEUT_BG = '#444444', '#FAF9F3'

fig, ax = plt.subplots(figsize=(10.2, 12.0))
ax.set_xlim(0, 10); ax.set_ylim(0, 12); ax.axis('off')

def box(x, y, w, h, title, body, fc=NEUT_BG, ec=NEUT, title_fs=10.5, body_fs=8.6,
        title_color=None, lw=1.4, bold_title=True):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.09,rounding_size=0.12',
                                fc=fc, ec=ec, lw=lw, zorder=2))
    tc = title_color or ec
    ax.text(x + w / 2, y + h - 0.24, title, ha='center', va='top', fontsize=title_fs,
            fontweight='bold' if bold_title else 'normal', color=tc, zorder=3)
    ax.text(x + w / 2, y + h - 0.56, body, ha='center', va='top', fontsize=body_fs,
            color='#333333', zorder=3, linespacing=1.45)

def arrow(x1, y1, x2, y2, color=NEUT, lw=1.6, style='-|>', connectionstyle='arc3,rad=0'):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle=style, mutation_scale=14,
                                 color=color, lw=lw, zorder=1,
                                 connectionstyle=connectionstyle))

# 1. Scan
box(1.3, 10.55, 7.4, 1.28, 'cis-MR scan: 498 SASP-gene \u00d7 skin-disease pairs',
    '83 SenMayo genes \u00d7 6 FinnGen R12 skin endpoints\n'
    'instrument: lead blood cis-eQTL per gene (eQTLGen, p < 5\u00d710\u207b\u2078; 480 pairs testable)')

# 2. Threshold
box(1.3, 8.95, 7.4, 1.18, 'Scan-wide Bonferroni threshold (0.05/498 = 1.0\u00d710\u207b\u2074)',
    '4 pairs pass:  ITPKA\u00d7eczema (p=2.9\u00d710\u207b\u2078), ITPKA\u00d7AD (p=1.4\u00d710\u207b\u2077),\n'
    'IL6ST\u00d7AD (p=2.1\u00d710\u207b\u2076), IL6ST\u00d7eczema (p=2.8\u00d710\u207b\u2076)')

# 3. Triage
box(1.3, 7.50, 7.4, 1.02, 'Colocalization triage (coloc.abf)',
    'all four pairs retained (ITPKA PP.H4 = 0.962 / 0.550; IL6ST PP.H4 = 0.994 / 0.997)')

arrow(5, 10.55, 5, 10.13); arrow(5, 8.95, 5, 8.52)

# 4. Split: IL6ST positive control (left, terminal) vs ITPKA novel (right)
arrow(5, 7.50, 2.55, 6.98, connectionstyle='arc3,rad=0.18')
arrow(5, 7.50, 7.45, 6.98, connectionstyle='arc3,rad=-0.18')
box(0.45, 5.72, 4.2, 1.26, 'IL6ST (5q11.2) \u2014 positive control',
    'established AD risk locus (rs7731626 genome-wide\n'
    'significant in Oliva et al.); recovery validates\n'
    'scan sensitivity \u2014 not pursued further',
    fc='#F2F7EA', ec=GREEN, title_color=GREEN)
box(5.35, 5.72, 4.2, 1.26, 'ITPKA (15q15.1) \u2014 novel nomination',
    'no prior link to skin or inflammatory skin disease;\n'
    'index variant rs11635906 (blood eQTL lead)\n'
    'taken forward for full characterization',
    fc='#EAF2FD', ec=BLUE, title_color=BLUE, lw=1.8)

# 5. Deep-dive panel (2 cols x 3 rows) under ITPKA: bracket connector with clear gap
arrow(7.45, 5.72, 7.45, 5.52, color=BLUE, lw=1.6)
ax.plot([2.76, 7.73], [5.50, 5.50], color=BLUE, lw=1.4, zorder=1)
for xc in (2.76, 7.73):
    ax.plot([xc, xc], [5.50, 5.34], color=BLUE, lw=1.4, zorder=1)
dd = [
    ('Colocalization robustness', 'window / p12-prior sensitivity;\ncoloc.susie relaxation\n(eczema PP.H4 = 0.962)'),
    ('Locus-wide discrimination', '12 genes \u00d7 6 eQTL resources;\nconditional analysis; SMR-HEIDI;\ncis multivariable MR'),
    ('Cross-cohort replication', 'Oliva & EAGLE lookups;\noverlap-aware meta-analysis\n(p = 2.2\u00d710\u207b\u2078, I\u00b2 = 0%)'),
    ('Phenome-wide scan', '64 independent phenotypes;\ntrait-level colocalization \u2192\natopic-spectrum restricted'),
    ('Tissue localization', 'skin eQTL, lesional transcriptome,\nsingle-cell, HPA protein \u2192\nskin-intrinsic evidence null'),
    ('Immune-cell attribution', '33 purified immune-cell eQTL\ndatasets; expression\nT-cell-lineage enriched'),
]
bw, bh, gx, gy = 4.42, 1.12, 0.55, 0.28
x0, y0 = 0.55, 1.38
for k, (t, b) in enumerate(dd):
    r, c = divmod(k, 2)
    bx = x0 + c * (bw + gx); by = y0 + (2 - r) * (bh + gy)
    box(bx, by, bw, bh, t, b, fc='white', ec='#9CC4E4', title_fs=9.3, body_fs=7.9,
        title_color=BLUE)

# 6. Conclusion
arrow(5.0, 1.36, 5.0, 1.12, color=ORANGE, lw=2.0)
box(1.3, 0.08, 7.4, 1.00, '',
    'Blood ITPKA up-regulation is a genetically supported, immune-mediated\n'
    'causal risk factor for dermatitis/eczema and atopic dermatitis',
    fc='#FDF3E7', ec=ORANGE, body_fs=9.6)

plt.savefig('/workspace/P2-02_Figure1_study_design_v2.png', dpi=300, bbox_inches='tight')
plt.savefig('/workspace/P2-02_Figure1_study_design_v2.svg', bbox_inches='tight')
print('Figure 1 v2 saved')
