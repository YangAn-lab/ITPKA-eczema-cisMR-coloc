import pandas as pd, numpy as np
import matplotlib
matplotlib.rcParams['font.family'] = ['Liberation Sans', 'Arimo', 'DejaVu Sans']
matplotlib.rcParams['svg.fonttype'] = 'none'
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

DIR = '/workspace/coloc_panel'
univ = pd.read_csv(f'{DIR}/snp_universe.csv').dropna(subset=['rsid'])
eq = pd.read_csv(f'{DIR}/eqtlgen_ITPKA.csv').sort_values('p').drop_duplicates('rsid')
eqm = eq.merge(univ[['rsid','position']].rename(columns={'position':'pos38'}), on='rsid')
eqm['mlogp_e'] = -np.log10(eqm.p)
gw = pd.read_csv(f'{DIR}/gwas_finngen_L12_DERMATITISECZEMA.tsv', sep='\t', header=None,
    names=['chromosome','position','ref','alt','rsids','nearest_genes','pval','mlogp','beta','sebeta','af_alt','af_alt_cases','af_alt_controls'])
VAR = 41487062

fig = plt.figure(figsize=(10.5, 8.5))
gs = GridSpec(3, 2, height_ratios=[2.2, 1, 1], hspace=0.42, wspace=0.28)

# ---- Panel A: 双轨 locus 图 ----
ax1 = fig.add_subplot(gs[0, :])
ax1.scatter(eqm.pos38/1e6, eqm.mlogp_e, s=6, c='#0279EE', alpha=0.6, linewidths=0)
ax1.axhline(-np.log10(5e-8), color='grey', ls=':', lw=0.8)
ax1.axvline(VAR/1e6, color='red', ls='--', lw=1)
ax1.text(VAR/1e6+0.012, 26, 'rs11635906', color='red', fontsize=8)
ax1.set_ylabel('ITPKA blood eQTL (eQTLGen)\n$-\log_{10}(p)$', fontsize=9)
ax1.set_xlim(40.99, 41.99)
ax1.spines[['top','right']].set_visible(False)
ax1.set_title('A', loc='left', fontsize=12, fontweight='bold')
plt.setp(ax1.get_xticklabels(), visible=False)

ax2 = fig.add_subplot(gs[1, :], sharex=ax1)
ax2.scatter(gw.position/1e6, gw.mlogp, s=6, c='#FF9400', alpha=0.6, linewidths=0)
ax2.axhline(-np.log10(5e-8), color='grey', ls=':', lw=0.8)
ax2.plot([VAR/1e6, VAR/1e6], [0, ax2.get_ylim()[1]], color='red', ls='--', lw=1)
ax2.set_ylabel('Dermatitis/eczema GWAS\n(FinnGen R12) $-\log_{10}(p)$', fontsize=9)
ax2.set_xlabel('chr15 position (Mb, hg38)', fontsize=9)
ax2.spines[['top','right']].set_visible(False)
genes = [('RTF1', 41.408408, 41.483563, -2.6, -4.3), ('ITPKA', 41.493360, 41.503554, -2.6, -6.1), ('RPAP1', 41.516557, 41.544281, -2.6, -4.3)]
for g, s, e, ybar, ytxt in genes:
    ax2.plot([s, e], [ybar, ybar], color='black', lw=3, solid_capstyle='butt', clip_on=False)
    ax2.text((s+e)/2, ytxt, g, ha='center', fontsize=8, style='italic')
ax2.set_ylim(bottom=-7.2)

# ---- Panel B: MR 森林图（OR 文本移到须线上方，避免与 Panel C 轴签重叠）----
ax3 = fig.add_subplot(gs[2, 0])
mr = pd.DataFrame({'ep': ['Dermatitis/eczema', 'Atopic dermatitis'], 'b': [0.237, 0.323], 'se': [0.0428, 0.0613]})
mr['or_'] = np.exp(mr.b); mr['lo'] = np.exp(mr.b-1.96*mr.se); mr['hi'] = np.exp(mr.b+1.96*mr.se)
y = [1, 0]
ax3.errorbar(mr['or_'], y, xerr=[mr['or_']-mr.lo, mr.hi-mr['or_']], fmt='o', color='#0279EE', capsize=4, ms=6)
for yi, (_, r) in zip(y, mr.iterrows()):
    ax3.text(r['or_'], yi+0.30, f"OR {r['or_']:.2f} ({r.lo:.2f}\u2013{r.hi:.2f})",
             va='bottom', ha='center', fontsize=8.5)
ax3.axvline(1, color='grey', ls='--', lw=0.8)
ax3.set_yticks(y); ax3.set_yticklabels(mr.ep, fontsize=9)
ax3.set_ylim(-0.55, 1.62)
ax3.set_xlim(0.9, 1.82); ax3.set_xlabel('Odds ratio per 1-unit higher blood ITPKA', fontsize=9)
ax3.set_title('B', loc='left', fontsize=12, fontweight='bold')
ax3.spines[['top','right']].set_visible(False)

# ---- Panel C: coloc PP.H4 + p12 范围 ----
ax4 = fig.add_subplot(gs[2, 1])
cc = pd.read_csv(f'{DIR}/coloc_eqtlgen_results.csv')
it = cc[cc.gene=='ITPKA']
x = [0, 1]
for xi, ep in zip(x, ['Eczema','AD']):
    r = it[it.endpoint==ep].iloc[0]
    ax4.bar(xi, r['PP.H4'], width=0.5, color=['#0279EE','#75A025'][xi], edgecolor='white')
    ax4.plot([xi, xi], [r['PP.H4_p12_1e6'], r['PP.H4_p12_1e4']], color='black', lw=1.5)
    ax4.text(xi, r['PP.H4']+0.03, f"{r['PP.H4']:.3f}", ha='center', fontsize=9)
ax4.axhline(0.75, color='red', ls='--', lw=1)
ax4.text(1.38, 0.70, 'PP.H4 = 0.75', color='red', fontsize=8)
ax4.set_xticks(x); ax4.set_xticklabels(['Dermatitis/eczema', 'Atopic dermatitis'], fontsize=9)
ax4.set_ylim(0, 1.12); ax4.set_ylabel('PP.H4 (coloc.abf)', fontsize=9)
ax4.set_title('C', loc='left', fontsize=12, fontweight='bold')
ax4.spines[['top','right']].set_visible(False)

plt.savefig('/mnt/results/P2-02_Figure2_MR_coloc_main.png', dpi=300, bbox_inches='tight')
plt.savefig('/mnt/results/P2-02_Figure2_MR_coloc_main.svg', bbox_inches='tight')
plt.show()
print('Figure2 saved with new name')