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

plt.savefig(f'{OUT}/P2-02_Figure2_MR_coloc_main_v5.png', dpi=300, bbox_inches='tight')
plt.savefig(f'{OUT}/P2-02_Figure2_MR_coloc_main_v5.svg', bbox_inches='tight')
plt.close()
print('Figure 2 v5 saved')

# ================= Figure 3 v5: 12-gene heatmap from S4 =================
s4 = pd.read_csv('/mnt/results/P2-02_SupplementaryTableS4_panel12_coloc.csv')

# gene order: ITPKA first, then original panel, then extension genes (by eQTLGen x Eczema PP.H4)
ecz = s4[(s4.dataset == 'eQTLGen_blood') & (s4.endpoint == 'Eczema')].set_index('gene')['PP.H4']
orig5 = ['NDUFAF1', 'RTF1', 'CHP1', 'RPAP1', 'OIP5-AS1']
new6 = ecz.drop(index=['ITPKA'] + orig5).sort_values(ascending=False).index.tolist()
gene_order = ['ITPKA'] + orig5 + new6
gene_labels = list(gene_order)

row_spec = [
    ('eQTLGen_blood', 'Eczema', 'eQTLGen blood (n~30.7k) \u00d7 Eczema FinnGen'),
    ('eQTLGen_blood', 'AD', 'eQTLGen blood \u00d7 AD FinnGen'),
    ('eQTLGen_blood', 'AD_OlivaEUR', 'eQTLGen blood \u00d7 AD Oliva EUR'),
    ('GTExv8_WholeBlood', 'Eczema', 'GTEx v8 whole blood (n=670) \u00d7 Eczema'),
    ('GTExv8_WholeBlood', 'AD', 'GTEx v8 whole blood \u00d7 AD FinnGen'),
    ('GTExv8_SkinSE', 'Eczema', 'GTEx v8 skin SE (n=605) \u00d7 Eczema'),
    ('GTExv8_SkinSE', 'AD', 'GTEx v8 skin SE \u00d7 AD'),
    ('GTExv8_SkinNSE', 'Eczema', 'GTEx v8 skin NSE (n=517) \u00d7 Eczema'),
    ('GTExv8_SkinNSE', 'AD', 'GTEx v8 skin NSE \u00d7 AD'),
    ('BLUEPRINT_Neutro', 'Eczema', 'BLUEPRINT neutrophil (n=196) \u00d7 Eczema'),
    ('BLUEPRINT_Neutro', 'AD', 'BLUEPRINT neutrophil \u00d7 AD'),
    ('BLUEPRINT_Mono', 'Eczema', 'BLUEPRINT monocyte (n=191) \u00d7 Eczema'),
    ('BLUEPRINT_Mono', 'AD', 'BLUEPRINT monocyte \u00d7 AD'),
]
M = np.full((len(row_spec), len(gene_order)), np.nan)
for i, (ds, ep, _) in enumerate(row_spec):
    for j, g in enumerate(gene_order):
        hit = s4[(s4.dataset == ds) & (s4.endpoint == ep) & (s4.gene == g)]
        if len(hit) and pd.notna(hit.iloc[0]['PP.H4']):
            M[i, j] = hit.iloc[0]['PP.H4']

cmap_h = LinearSegmentedColormap.from_list('pp', ['#FAF9F3', '#9CC4E4', '#0279EE'])
cmap_h.set_bad('#DDDDDD')
fig, ax = plt.subplots(figsize=(9.8, 6.9))
im = ax.imshow(M, cmap=cmap_h, vmin=0, vmax=1, aspect='auto')
ax.set_xticks(range(len(gene_order)))
ax.set_xticklabels(gene_labels, fontsize=9, style='italic', rotation=38, ha='right',
                   rotation_mode='anchor')
# ITPKA column header emphasis
ax.get_xticklabels()[0].set_fontweight('bold')
ax.set_yticks(range(len(row_spec))); ax.set_yticklabels([r[2] for r in row_spec], fontsize=8.5)
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        v = M[i, j]
        if np.isnan(v):
            ax.text(j, i, 'n.t.', ha='center', va='center', fontsize=7.5, color='#444444')
        else:
            strong = v >= 0.75
            ax.text(j, i, f'{v:.2f}' + ('*' if strong else ''), ha='center', va='center',
                    fontsize=8, fontweight='bold' if strong else 'normal',
                    color='white' if v > 0.65 else '#222222')
# separator between original panel and extension genes
ax.axvline(5.5, color='#888888', lw=1.0, ls='-')
ax.text(5.5, -0.95, 'original panel | added by extension', ha='center', va='bottom', fontsize=8, color='#444444')
ax.set_xticks(np.arange(-0.5, len(gene_order)), minor=True)
ax.set_yticks(np.arange(-0.5, len(row_spec)), minor=True)
ax.grid(which='minor', color='white', lw=1.4)
ax.tick_params(which='both', length=0)
cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
cb.set_label('PP.H4 (coloc.abf)', fontsize=9); cb.ax.tick_params(labelsize=8)

plt.tight_layout()
plt.savefig(f'{OUT}/P2-02_Figure3_locus_discrimination_v5.png', dpi=300, bbox_inches='tight')
plt.savefig(f'{OUT}/P2-02_Figure3_locus_discrimination_v5.svg', bbox_inches='tight')
plt.close()
print('Figure 3 v5 saved')

# ================= Figure S2 条件分析 =================
eq_it = load_eqtl('ITPKA'); eq_rt = load_eqtl('RTF1')
gw_ecz = load_finngen('L12_DERMATITISECZEMA'); gw_ol = load_oliva()

fig, axes = plt.subplots(2, 2, figsize=(11.5, 7.6), sharex=True)
xlim = (41.20, 41.60)

def plot_cond(ax, df_all, label_before, label_after, r_col_snp, color_after, vlines=True):
    ax.scatter(df_all.position / 1e6, -np.log10(df_all.pval if 'pval' in df_all else df_all.p),
               s=5, c='#BBBBBB', linewidths=0, label=label_before)
    xc, yc = cond_track(df_all, r_col_snp)
    ax.scatter(xc, yc, s=7, c=color_after, linewidths=0, label=label_after)
    if vlines:
        ax.axvline(IDX / 1e6, color='red', ls='--', lw=0.9)
        ax.axvline(41482225 / 1e6, color='purple', ls=':', lw=0.9)
    ax.set_xlim(xlim); ax.spines[['top', 'right']].set_visible(False)
    ax.legend(fontsize=7.5, frameon=False, loc='upper right')

plot_cond(axes[0, 0], gw_ecz, 'unconditioned', 'conditioned on rs11635906', i_idx, '#FF9400')
axes[0, 0].set_ylabel('Eczema GWAS (FinnGen)\n$-\\log_{10}(p)$', fontsize=9)
axes[0, 0].set_title('A', loc='left', fontsize=12, fontweight='bold')

plot_cond(axes[0, 1], eq_it, 'unconditioned', 'conditioned on rs11635906', i_idx, '#0279EE')
axes[0, 1].set_ylabel('ITPKA eQTL (eQTLGen)\n$-\\log_{10}(p)$', fontsize=9)
axes[0, 1].set_title('B', loc='left', fontsize=12, fontweight='bold')

plot_cond(axes[1, 0], eq_rt, 'unconditioned', 'conditioned on rs11635906', i_idx, '#75A025')
axes[1, 0].set_ylabel('RTF1 eQTL (eQTLGen)\n$-\\log_{10}(p)$', fontsize=9)
axes[1, 0].set_xlabel('chr15 position (Mb, hg38)', fontsize=9)
axes[1, 0].set_title('C', loc='left', fontsize=12, fontweight='bold')

# Panel D：Oliva AD 双向条件
ax = axes[1, 1]
ax.scatter(gw_ol.position / 1e6, -np.log10(gw_ol.pval), s=5, c='#BBBBBB', linewidths=0, label='unconditioned')
xc1, yc1 = cond_track(gw_ol, i_idx)
ax.scatter(xc1, yc1, s=7, c='#FF9400', linewidths=0, label='conditioned on rs11635906')
xc2, yc2 = cond_track(gw_ol, i_top)
ax.scatter(xc2, yc2, s=7, c='#7E57C2', linewidths=0, label='conditioned on rs12440045')
ax.axvline(IDX / 1e6, color='red', ls='--', lw=0.9)
ax.axvline(41482225 / 1e6, color='purple', ls=':', lw=0.9)
ax.set_xlim(xlim); ax.spines[['top', 'right']].set_visible(False)
ax.legend(fontsize=7.5, frameon=False, loc='upper right')
ax.set_ylabel('AD GWAS (Oliva EUR)\n$-\\log_{10}(p)$', fontsize=9)
ax.set_xlabel('chr15 position (Mb, hg38)', fontsize=9)
ax.set_title('D', loc='left', fontsize=12, fontweight='bold')

fig.text(0.52, 0.005, 'red dashed = rs11635906 (index);  purple dotted = rs1942 (RTF1 eQTL lead)',
         ha='center', fontsize=8, color='#444444')
plt.tight_layout(rect=[0, 0.02, 1, 1])
plt.savefig(f'{OUT}/P2-02_FigureS2_conditional_locus.png', dpi=300, bbox_inches='tight')
plt.savefig(f'{OUT}/P2-02_FigureS2_conditional_locus.svg', bbox_inches='tight')
plt.close()
print('Figure S2 saved')
