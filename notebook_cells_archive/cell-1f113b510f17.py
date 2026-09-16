import pandas as pd
import numpy as np

# ---- 合并 GTEx v8 / BLUEPRINT 面板结果 + eQTLGen 结果 → 对比表 ----
panel = pd.read_csv('/workspace/coloc_panel/coloc_panel_results.csv')
eqg   = pd.read_csv('/workspace/coloc_panel/coloc_eqtlgen_results.csv')
ws    = pd.read_csv('/workspace/coloc_panel/coloc_eqtlgen_window_sensitivity.csv')

panel_ok = panel[panel.note == 'ok'].copy()
eqg_ok   = eqg[eqg.note == 'ok'].copy()

ds_order = ['eQTLGen_blood', 'GTExv8_WholeBlood', 'GTExv8_SkinSE', 'GTExv8_SkinNSE',
            'BLUEPRINT_Neutro', 'BLUEPRINT_Mono']
ds_label = {'eQTLGen_blood': 'eQTLGen blood (n≈30.4k)',
            'GTExv8_WholeBlood': 'GTEx v8 whole blood (n=670)',
            'GTExv8_SkinSE': 'GTEx v8 skin, sun-exposed (n=605)',
            'GTExv8_SkinNSE': 'GTEx v8 skin, not sun-exposed (n=517)',
            'BLUEPRINT_Neutro': 'BLUEPRINT neutrophil (n=196)',
            'BLUEPRINT_Mono': 'BLUEPRINT monocyte (n=191)'}
gene_order = ['ITPKA', 'NDUFAF1', 'RTF1', 'CHP1', 'RPAP1', 'OIP5-AS1']

allres = pd.concat([panel_ok, eqg_ok], ignore_index=True)

rows = []
for g in gene_order:
    for ep in ['Eczema', 'AD']:
        row = {'gene': g, 'endpoint': ep}
        for ds in ds_order:
            r = allres[(allres.dataset == ds) & (allres.gene == g) & (allres.endpoint == ep)]
            row[ds_label[ds]] = round(float(r['PP.H4'].iloc[0]), 3) if len(r) else np.nan
        # eQTLGen 窗口稳定性
        w = ws[(ws.gene == g) & (ws.endpoint == ep)]
        row['eQTLGen PP.H4 @±250kb'] = round(float(w[w.window_kb == 250]['PP.H4'].iloc[0]), 3) if len(w[w.window_kb == 250]) else np.nan
        row['eQTLGen PP.H4 @±125kb'] = round(float(w[w.window_kb == 125]['PP.H4'].iloc[0]), 3) if len(w[w.window_kb == 125]) else np.nan
        rows.append(row)
comp = pd.DataFrame(rows)
comp.to_csv('/mnt/results/P2-02_eQTLGen敏感性_对比表.csv', index=False)
print(comp.to_string(index=False))