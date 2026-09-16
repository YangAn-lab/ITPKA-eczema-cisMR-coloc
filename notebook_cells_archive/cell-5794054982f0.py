# 分析H 第1步：He 2020 (GSE147424) 每样本 T 细胞占比（CD3D/E/G 任一>0）
# 目的：定量 AD 皮损中 T 细胞扩增倍数 → 稀释效应估算
import gzip, os, re
import numpy as np, pandas as pd

group = {'sample1':('AD','Lesional'),'sample2':('AD','Lesional'),'sample3':('AD','Non-lesional'),
         'sample4':('Healthy','Healthy'),'sample5':('AD','Lesional'),'sample6':('Healthy','Healthy'),
         'sample7':('AD','Lesional'),'sample8':('Healthy','Healthy'),'sample9':('Healthy','Healthy'),
         'sample10':('Healthy','Healthy'),'sample11':('AD','Non-lesional'),'sample12':('Healthy','Healthy'),
         'sample13':('Healthy','Healthy'),'sample14':('AD','Non-lesional'),'sample15':('AD','Non-lesional'),
         'sample16':('AD','Non-lesional'),'sample17':('Healthy','Healthy')}

MARKERS = {'"CD3D"':'CD3D','"CD3E"':'CD3E','"CD3G"':'CD3G','"CD8A"':'CD8A','"CD4"':'CD4','"ITPKA"':'ITPKA'}
rows = []
base = '/workspace/gse147424'
for f in sorted(os.listdir(base)):
    if not f.endswith('.gz'): continue
    sname = re.search(r'(sample\d+)', f).group(1)
    vals = {}
    with gzip.open(os.path.join(base, f), 'rt') as fh:
        header = fh.readline()
        n_cells = len(header.rstrip('\n').split(',')) - 1
        for line in fh:
            g = line.split(',', 1)[0]
            if g in MARKERS:
                vals[MARKERS[g]] = np.array([float(v) for v in line.rstrip('\n').split(',')[1:] if v != ''])
    cd3 = np.maximum.reduce([vals.get('CD3D', np.zeros(n_cells)),
                             vals.get('CD3E', np.zeros(n_cells)),
                             vals.get('CD3G', np.zeros(n_cells))])
    tcell = cd3 > 0
    itpka = vals.get('ITPKA', np.zeros(n_cells))
    dis, tt = group[sname]
    rows.append({'sample': sname, 'disease': dis, 'tissue_type': tt, 'n_cells': n_cells,
                 'n_T': int(tcell.sum()), 'pct_T': round(100*tcell.mean(), 2),
                 'n_CD8A': int((vals.get('CD8A', np.zeros(n_cells)) > 0).sum()),
                 'n_ITPKA': int((itpka > 0).sum()),
                 'n_ITPKA_in_T': int(((itpka > 0) & tcell).sum())})

tfrac = pd.DataFrame(rows)
print(tfrac.to_string(index=False))
print('\n=== 按组汇总 ===')
g = tfrac.groupby('tissue_type').agg(n_samples=('sample','size'), cells=('n_cells','sum'),
                                     T=('n_T','sum'), ITPKA=('n_ITPKA','sum'), ITPKA_in_T=('n_ITPKA_in_T','sum'))
g['pct_T'] = (100*g['T']/g['cells']).round(2)
g['pct_ITPKA_all'] = (100*g['ITPKA']/g['cells']).round(4)
g['pct_ITPKA_in_T'] = (100*g['ITPKA_in_T']/g['T']).round(3)
print(g.to_string())
tfrac.to_csv('/workspace/coloc_panel/he2020_tcell_fractions.csv', index=False)
