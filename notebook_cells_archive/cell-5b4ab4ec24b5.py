# 分析D: 6 基因在 GTEx v8 全血中的 eQTL 可测性/信号强度
import pandas as pd, numpy as np
from scipy.stats import norm

cols_v8 = ["variant","r2","pvalue","mto_id","molecular_trait_id","maf","gene_id",
           "median_tpm","beta","se","an","ac","chromosome","position","ref","alt","type","rsid"]
wb = pd.read_csv('/workspace/coloc_panel/eqtl_gtexv8_whole_blood.tsv', sep='\t', header=None, names=cols_v8,
                 dtype={'pvalue':float,'beta':float,'se':float,'maf':float})
print('rows:', len(wb))

panel = {'ENSG00000137825':'ITPKA','ENSG00000137806':'NDUFAF1','ENSG00000187446':'CHP1',
         'ENSG00000103932':'RPAP1','ENSG00000137815':'RTF1','ENSG00000247556':'OIP5-AS1'}

rows = []
for gid, sym in panel.items():
    sub = wb[wb['gene_id'].str.startswith(gid)]
    if len(sub)==0:
        rows.append({'gene':sym,'n_variants_tested':0}); continue
    sub = sub.copy()
    sub['z'] = (sub['beta']/sub['se']).abs()
    lead = sub.loc[sub['pvalue'].idxmin()]
    at_iv = sub[sub['position']==41487062]
    rows.append({
        'gene': sym,
        'median_tpm': round(sub['median_tpm'].iloc[0],3),
        'n_variants_tested': len(sub),
        'lead_variant': lead['rsid'] if pd.notna(lead['rsid']) else lead['variant'],
        'lead_p': f"{lead['pvalue']:.2e}",
        'lead_abs_z': round(lead['z'],2),
        'p_at_rs11635906': (f"{at_iv['pvalue'].iloc[0]:.2e}" if len(at_iv) else 'not tested'),
        'beta_at_rs11635906': (round(at_iv['beta'].iloc[0],3) if len(at_iv) else '-'),
    })

meas = pd.DataFrame(rows)
print(meas.to_string(index=False))
meas.to_csv('/mnt/results/P2-02_GTExv8全血_6基因eQTL可测性.csv', index=False)
print('\nsaved: P2-02_GTExv8全血_6基因eQTL可测性.csv')