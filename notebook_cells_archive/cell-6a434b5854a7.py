import pandas as pd
import matplotlib
matplotlib.rcParams['font.family'] = ['Liberation Sans', 'Arimo', 'DejaVu Sans']
matplotlib.rcParams['svg.fonttype'] = 'none'
import matplotlib.pyplot as plt
import numpy as np

panel = pd.read_csv('/workspace/coloc_panel/coloc_panel_results.csv')
eqg   = pd.read_csv('/workspace/coloc_panel/coloc_eqtlgen_results.csv')
df = pd.concat([panel[panel.note=='ok'], eqg[eqg.note=='ok']], ignore_index=True)

ds_order = ['eQTLGen_blood', 'GTExv8_WholeBlood', 'GTExv8_SkinSE', 'GTExv8_SkinNSE', 'BLUEPRINT_Neutro', 'BLUEPRINT_Mono']
ds_label = {'eQTLGen_blood': 'eQTLGen blood (n~30.4k)',
            'GTExv8_WholeBlood': 'GTEx v8 whole blood (n=670)',
            'GTExv8_SkinSE': 'GTEx v8 skin, sun-exposed (n=605)',
            'GTExv8_SkinNSE': 'GTEx v8 skin, not sun-exposed (n=517)',
            'BLUEPRINT_Neutro': 'BLUEPRINT neutrophil (n=196)',
            'BLUEPRINT_Mono': 'BLUEPRINT monocyte (n=191)'}
gene_order = ['ITPKA', 'NDUFAF1', 'RTF1', 'CHP1', 'RPAP1', 'OIP5-AS1']
colors = {'eQTLGen_blood': '#E9ED4C', 'GTExv8_WholeBlood': '#0279EE', 'GTExv8_SkinSE': '#FF9400',
          'GTExv8_SkinNSE': '#F5C98A', 'BLUEPRINT_Neutro': '#75A025', 'BLUEPRINT_Mono': '#B5D48A'}

fig, axes = plt.subplots(1, 2, figsize=(11.5, 6.8), sharex=True)
for ax, ep, title in zip(axes, ['Eczema', 'AD'],
                         ['FinnGen R12 dermatitis/eczema (67,474 cases)',
                          'FinnGen R12 atopic dermatitis (31,245 cases)']):
    sub = df[df.endpoint == ep]
    ypos, ylabels = [], []
    y = 0
    for ds in ds_order:
        for g in gene_order:
            row = sub[(sub.dataset == ds) & (sub.gene == g)]
            if len(row) == 0:
                continue
            r = row.iloc[0]
            ax.barh(y, r['PP.H4'], color=colors[ds], height=0.72, edgecolor='white', linewidth=0.4)
            ax.plot([r['PP.H4_p12_1e6'], r['PP.H4_p12_1e4']], [y, y], color='black', linewidth=1.0, zorder=5)
            ax.text(min(r['PP.H4'] + 0.02, 1.02), y, f"{r['PP.H4']:.2f}", va='center', fontsize=7.5)
            ypos.append(y); ylabels.append(g)
            y += 1
        y += 0.8
    ax.axvline(0.75, color='red', linestyle='--', linewidth=1)
    ax.text(0.76, -0.9, 'PP.H4 = 0.75', color='red', fontsize=8, va='bottom', ha='left')
    ax.set_yticks(ypos); ax.set_yticklabels(ylabels, fontsize=8.5)
    ax.invert_yaxis()
    ax.set_xlim(0, 1.12); ax.set_xlabel('PP.H4 (coloc.abf)')
    ax.set_title(title, fontsize=10)
    ax.spines[['top', 'right']].set_visible(False)
    ax.set_ylim(y - 0.2, -1.6)

handles = [plt.Rectangle((0, 0), 1, 1, color=colors[d]) for d in ds_order]
fig.legend(handles, [ds_label[d] for d in ds_order], loc='lower center', ncol=3,
           fontsize=8, frameon=False, bbox_to_anchor=(0.5, -0.02))
plt.tight_layout(rect=[0, 0.06, 1, 1])
plt.savefig('/mnt/results/P2-02_fig_位点内多基因coloc鉴别_v2.png', dpi=300, bbox_inches='tight')
plt.savefig('/mnt/results/P2-02_fig_位点内多基因coloc鉴别_v2.svg', bbox_inches='tight')
plt.show()
print('figure v2 saved')