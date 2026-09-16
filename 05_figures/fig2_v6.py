#!/usr/bin/env python3
"""P2-02 G7 v4 图件：Figure 2 v4（LD 着色+Oliva 第三柱）、Figure 3 v4（热图）、Figure S2（条件分析 locus 图）
数据源：coloc_panel 修正版结果 + coloc_oliva 分析K结果 + ld_eur 面板（hg38, alt 定向）
"""
import pandas as pd, numpy as np, gzip
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['font.family'] = ['Liberation Sans', 'Arimo', 'DejaVu Sans']
matplotlib.rcParams['svg.fonttype'] = 'none'
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import norm

DIR = '/workspace/coloc_panel'; OL = '/workspace/coloc_oliva'; OUT = '/workspace'
IDX = 41487062  # rs11635906 hg38

# ---------- 共享数据 ----------
ld = pd.read_csv(f'{DIR}/ld_eur_snps.csv')
ld['snp'] = (ld.chromosome.astype(str) + ':' + ld.position.astype(str) + ':' + ld.ref + ':' + ld.alt)
D = np.loadtxt(gzip.open(f'{DIR}/ld_eur.csv.gz', 'rt'), delimiter=',')
snp_to_i = {s: i for i, s in enumerate(ld.snp)}
i_idx = snp_to_i[f'15:{IDX}:A:G']
r_with_idx = D[:, i_idx]                      # 有符号 r（alt 定向）
i_top = snp_to_i['15:41490486:A:C']           # Oliva block top rs12440045
r_with_top = D[:, i_top]

univ = pd.read_csv(f'{DIR}/snp_universe.csv').dropna(subset=['rsid'])

def load_eqtl(sym):
    eq = pd.read_csv(f'{DIR}/eqtlgen_{sym}.csv').sort_values('p').drop_duplicates('rsid')
    m = eq.merge(univ[['rsid', 'chromosome', 'position', 'ref', 'alt']], on='rsid',
                 suffixes=('_eq19', ''))   # position_eq19=hg19；position=univ hg38
    exact = (m.ea == m.alt) & (m.nea == m.ref); swap = (m.ea == m.ref) & (m.nea == m.alt)
    m = m[exact | swap].copy(); sw = swap[exact | swap].values
    m['z_alt'] = np.where(sw, -1.0, 1.0) * (m.beta / m.se)
    m['snp'] = (m.chromosome.astype(str) + ':' + m.position.astype(str) + ':' + m.ref + ':' + m.alt)
    return m

def load_finngen(ep):
    gw = pd.read_csv(f'{DIR}/gwas_finngen_{ep}.tsv', sep='\t', header=None,
        names=['chromosome','position','ref','alt','rsids','nearest_genes','pval','mlogp',
               'beta','sebeta','af_alt','af_alt_cases','af_alt_controls'])
    gw['z_alt'] = gw.beta / gw.sebeta
    gw['snp'] = (gw.chromosome.astype(str) + ':' + gw.position.astype(str) + ':' + gw.ref + ':' + gw.alt)
    return gw

def load_oliva():
    gw = pd.read_csv(f'{OL}/oliva_eur_locus.tsv', sep='\t', header=None,
        names=['chromosome','position','ea','oa','beta','sebeta','eaf','pval','vid','rsid','z','nstudy','n','eff','conv','code'])
    # 定向到 LD 面板 alt：面板 ref/alt vs Oliva oa/ea
    gw['snp'] = ('15:' + gw.position.astype(str) + ':' + gw.oa + ':' + gw.ea)
    gw['z_alt'] = gw.beta / gw.sebeta
    return gw

def add_ld(df):
    """给含 snp 列的 df 加 r2_idx（与 index 的 r²）与面板内标记"""
    ii = df.snp.map(snp_to_i)
    df = df.assign(in_ld=ii.notna())
    df['r2_idx'] = np.where(df.in_ld, r_with_idx[ii.fillna(0).astype(int)] ** 2, np.nan)
    return df

def cond_track(df, r_col_snp):
    """返回 df 中 LD 面板内 SNP 的 (pos, -log10 p_cond)，相对指定条件 SNP；剔除 |r|>0.99"""
    ii = df.snp.map(snp_to_i); have = ii.notna()
    sub = df[have].copy(); ii = ii[have].astype(int)
    r = D[ii, r_col_snp]
    zc = df.loc[df.snp == ld.snp[r_col_snp], 'z_alt']
    zc = zc.iloc[0] if len(zc) else None
    if zc is None:  # 条件 SNP 不在 df 中（如 Oliva 轨道缺 index 行不会发生，但保险）
        raise ValueError('conditioning SNP missing')
    keep = np.abs(r) < 0.99
    z_cond = (sub.z_alt.values[keep] - r[keep] * zc) / np.sqrt(1 - r[keep] ** 2)
    p_cond = 2 * norm.sf(np.abs(z_cond))
    return sub.position.values[keep] / 1e6, -np.log10(p_cond)

# ================= Figure 2 v4 =================
eqm = add_ld(load_eqtl('ITPKA'))
gw = add_ld(load_finngen('L12_DERMATISECZEMA') if False else load_finngen('L12_DERMATITISECZEMA'))

cmap_ld = plt.cm.viridis
def ld_colors(df):
    c = np.full((len(df), 4), 0.82)  # 灰
    c[:, 3] = 1.0
    m = df.in_ld.values
    c[m] = cmap_ld(df.r2_idx[m].values)
    return c

fig = plt.figure(figsize=(10.5, 9.6))
gs = GridSpec(4, 2, height_ratios=[2.3, 1.05, 0.30, 1.05], hspace=0.42, wspace=0.28)

ax1 = fig.add_subplot(gs[0, :])
ax1.scatter(eqm.position / 1e6, -np.log10(eqm.p), s=8, c=ld_colors(eqm), linewidths=0)
ax1.scatter([IDX / 1e6], [-np.log10(eqm.loc[eqm.rsid == 'rs11635906', 'p'].iloc[0])],
            marker='D', s=55, facecolor='none', edgecolor='red', linewidths=1.6, zorder=6)
ax1.axhline(-np.log10(5e-8), color='grey', ls=':', lw=0.8)
ax1.text(IDX / 1e6 + 0.012, 26, 'rs11635906', color='red', fontsize=8)
# signal-B block leads (within 5 kb of index; mark on top spine)
for xpos, lab, c in [(41.482225, 'rs1942', 'purple'), (41.490486, 'rs12440045', '#7E57C2')]:
    ax1.scatter([xpos], [1.005], marker='v', s=38, color=c, clip_on=False,
                transform=ax1.get_xaxis_transform(), zorder=7)
ax1.text(0.012, 0.965, '\u25bc rs1942 (RTF1 eQTL lead)', color='purple', fontsize=7.5,
         transform=ax1.transAxes, va='top')
ax1.text(0.012, 0.895, '\u25bc rs12440045 (block-B lead)', color='#7E57C2', fontsize=7.5,
         transform=ax1.transAxes, va='top')
ax1.set_ylabel('ITPKA blood eQTL (eQTLGen)\n$-\\log_{10}(p)$', fontsize=9)
ax1.set_xlim(40.99, 41.99)
ax1.spines[['top', 'right']].set_visible(False)
ax1.set_title('A', loc='left', fontsize=12, fontweight='bold')
plt.setp(ax1.get_xticklabels(), visible=False)
sm = plt.cm.ScalarMappable(cmap=cmap_ld, norm=plt.Normalize(0, 1))
cb = fig.colorbar(sm, ax=ax1, pad=0.01, fraction=0.025)
cb.set_label('LD $r^2$ with rs11635906', fontsize=8); cb.ax.tick_params(labelsize=7)

ax2 = fig.add_subplot(gs[1, :], sharex=ax1)
ax2.scatter(gw.position / 1e6, gw.mlogp, s=8, c=ld_colors(gw), linewidths=0)
gw_idx = gw[gw.position == IDX]
ax2.scatter([IDX / 1e6], [gw_idx.mlogp.iloc[0]], marker='D', s=55,
            facecolor='none', edgecolor='red', linewidths=1.6, zorder=6)
ax2.axhline(-np.log10(5e-8), color='grey', ls=':', lw=0.8)
ax2.set_ylabel('Dermatitis/eczema GWAS\n(FinnGen R12) $-\\log_{10}(p)$', fontsize=9)
ax2.spines[['top', 'right']].set_visible(False)
for xpos, c in [(41.482225, 'purple'), (41.490486, '#7E57C2')]:
    ax2.axvline(xpos, color=c, ls=':', lw=0.7, alpha=0.55)
ax2.set_ylim(bottom=0)
plt.setp(ax2.get_xticklabels(), visible=False)

# gene models on a dedicated strip axis directly below the GWAS track
pos2 = ax2.get_position()
axg = fig.add_axes([pos2.x0, pos2.y0 - 0.062, pos2.width, 0.062], sharex=ax1)
for g, s, e, yl in [('RTF1', 41.408408, 41.483563, 0), ('ITPKA', 41.493360, 41.503554, 0),
                    ('RPAP1', 41.516557, 41.544281, 0)]:
    axg.plot([s, e], [0, 0], color='black', lw=3, solid_capstyle='butt')
# stagger labels vertically to avoid collision (ITPKA/RPAP1 are close)
axg.text(41.446, -0.70, 'RTF1', ha='center', va='top', fontsize=8, style='italic')
axg.text(41.4985, -0.70, 'ITPKA', ha='center', va='top', fontsize=8, style='italic')
axg.text(41.5305, -1.45, 'RPAP1', ha='center', va='top', fontsize=8, style='italic')
axg.plot([41.5305, 41.5305], [-0.15, -1.35], color='#888888', lw=0.6)
axg.set_ylim(-2.6, 0.45); axg.set_yticks([])
axg.set_ylabel('genes', fontsize=8, rotation=0, ha='right', va='center')
axg.set_xlabel('chr15 position (Mb, hg38)', fontsize=9)
axg.spines[['top', 'right', 'left']].set_visible(False)

ax3 = fig.add_subplot(gs[3, 0])
mr = pd.DataFrame({'ep': ['Dermatitis/eczema', 'Atopic dermatitis'], 'b': [0.237, 0.323], 'se': [0.0428, 0.0613]})
mr['or_'] = np.exp(mr.b); mr['lo'] = np.exp(mr.b - 1.96 * mr.se); mr['hi'] = np.exp(mr.b + 1.96 * mr.se)
y = [1, 0]
ax3.errorbar(mr['or_'], y, xerr=[mr['or_'] - mr.lo, mr.hi - mr['or_']], fmt='o', color='#0279EE', capsize=4, ms=6)
for yi, (_, r) in zip(y, mr.iterrows()):
    ax3.text(r['or_'], yi + 0.30, f"OR {r['or_']:.2f} ({r.lo:.2f}\u2013{r.hi:.2f})", va='bottom', ha='center', fontsize=8.5)
ax3.axvline(1, color='grey', ls='--', lw=0.8)
ax3.set_yticks(y); ax3.set_yticklabels(mr.ep, fontsize=9)
ax3.set_ylim(-0.55, 1.62); ax3.set_xlim(0.9, 1.82)
ax3.set_xlabel('Odds ratio per 1-unit higher blood ITPKA', fontsize=9)
ax3.set_title('B', loc='left', fontsize=12, fontweight='bold')
ax3.spines[['top', 'right']].set_visible(False)

ax4 = fig.add_subplot(gs[3, 1])
cc = pd.read_csv(f'{DIR}/coloc_eqtlgen_results_fixed.csv')
co = pd.read_csv(f'{OL}/coloc_oliva_eur_results.csv')
rows_c = [
    ('Dermatitis/eczema\nFinnGen (67,474)', cc[(cc.gene == 'ITPKA') & (cc.endpoint == 'Eczema')].iloc[0], '#0279EE'),
    ('Atopic dermatitis\nFinnGen (31,245)', cc[(cc.gene == 'ITPKA') & (cc.endpoint == 'AD')].iloc[0], '#75A025'),
    ('Atopic dermatitis\nOliva EUR (42,963)', co[co.gene == 'ITPKA'].iloc[0], '#FF9400'),
]
for xi, (lab, r, c) in enumerate(rows_c):
    ax4.bar(xi, r['PP.H4'], width=0.55, color=c, edgecolor='white')
    ax4.plot([xi, xi], [r['PP.H4_p12_1e6'], r['PP.H4_p12_1e4']], color='black', lw=1.5)
    ax4.text(xi, max(r['PP.H4'], r['PP.H4_p12_1e4']) + 0.035, f"{r['PP.H4']:.3f}", ha='center', fontsize=9)
ax4.axhline(0.75, color='red', ls='--', lw=1)
ax4.text(0.5, 0.77, 'PP.H4 = 0.75', color='red', fontsize=8, ha='center', va='bottom')
ax4.set_xticks(range(3)); ax4.set_xticklabels([r[0] for r in rows_c], fontsize=7.8)
ax4.set_ylim(0, 1.12); ax4.set_ylabel('PP.H4 (coloc.abf)', fontsize=9)
ax4.set_title('C', loc='left', fontsize=12, fontweight='bold')
ax4.spines[['top', 'right']].set_visible(False)

plt.savefig(f'{OUT}/P2-02_Figure2_MR_coloc_main_v6.png', dpi=300, bbox_inches='tight')
plt.savefig(f'{OUT}/P2-02_Figure2_MR_coloc_main_v6.svg', bbox_inches='tight')
plt.close()
print('Figure 2 v5 saved')