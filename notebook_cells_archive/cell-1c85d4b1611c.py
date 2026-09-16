# FigureS1: 正常皮肤单细胞 ITPKA 表达谱（CELLxGENE Census, 672,457 细胞）
import pandas as pd, numpy as np
import matplotlib
matplotlib.rcParams['font.family'] = ['Liberation Sans', 'Arimo', 'DejaVu Sans']
matplotlib.rcParams['svg.fonttype'] = 'none'
import matplotlib.pyplot as plt

bc = pd.read_csv('/mnt/results/P2-02_单细胞_正常皮肤ITPKA_大类汇总.csv')
bc = bc[bc['broad_class'] != 'Other'].sort_values('pct_expressing')
# 加入 T 细胞亚型细分（来自细类表）
ct = pd.read_csv('/mnt/results/P2-02_单细胞_正常皮肤ITPKA_细胞类型谱.csv')
tsub = ct[ct['cell_type'].isin(['CD8-positive, alpha-beta memory T cell, CD45RO-positive',
                                'helper T cell','regulatory T cell',
                                'naive thymus-derived CD4-positive, alpha-beta T cell'])].copy()
tsub['broad_class'] = ['  CD8 memory T','  helper T','  regulatory T','  naive CD4 T']
tsub = tsub.rename(columns={})

fig, ax = plt.subplots(figsize=(8.5, 5.2))
# 主图: 大类
labels = list(bc['broad_class'])
vals = list(bc['pct_expressing'])
ns = list(bc['n_cells'])
ypos = np.arange(len(labels))
cols = ['#0279EE' if l in ('T cell','Dendritic cell','Macrophage') else '#999999' for l in labels]
ax.barh(ypos, vals, color=cols, height=0.62)
for y, v, n in zip(ypos, vals, ns):
    ax.text(v+0.12, y, f'{v:.1f}%  (n={n:,})', va='center', fontsize=9)
ax.set_yticks(ypos); ax.set_yticklabels(labels, fontsize=10)
ax.set_xlabel('Cells with detectable ITPKA (%)', fontsize=11)
ax.set_xlim(0, 17)
ax.spines[['top','right']].set_visible(False)
ax.set_title('ITPKA expression across normal-skin cell classes\n(CELLxGENE Census 2025-01-30; 672,457 cells, 233 donors)',
             fontsize=11, loc='left')
# 注记 T 细胞亚型
ax.text(16.8, len(labels)-0.5,
        'Within T cells: CD8 memory 11.2%, helper 7.4%, Treg 6.0%, naive CD4/CD8 ≤0.15%',
        fontsize=8.5, ha='right', style='italic', color='#0279EE')
plt.tight_layout()
plt.savefig('/workspace/P2-02_FigureS1_singlecell_skin.png', dpi=300)
plt.savefig('/workspace/P2-02_FigureS1_singlecell_skin.svg')
import shutil
shutil.copy('/workspace/P2-02_FigureS1_singlecell_skin.png', '/mnt/results/P2-02_FigureS1_singlecell_skin.png')
shutil.copy('/workspace/P2-02_FigureS1_singlecell_skin.svg', '/mnt/results/P2-02_FigureS1_singlecell_skin.svg')
print('FigureS1 saved')