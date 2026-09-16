# 分析A 补充: He 2020 (GSE147424) AD 皮损 scRNA 中 ITPKA 检出率（修正文件名解析）
import gzip, os, re
import numpy as np, pandas as pd

group = {'sample1':('AD','Lesional'),'sample2':('AD','Lesional'),'sample3':('AD','Non-lesional'),
         'sample4':('Healthy','Healthy'),'sample5':('AD','Lesional'),'sample6':('Healthy','Healthy'),
         'sample7':('AD','Lesional'),'sample8':('Healthy','Healthy'),'sample9':('Healthy','Healthy'),
         'sample10':('Healthy','Healthy'),'sample11':('AD','Non-lesional'),'sample12':('Healthy','Healthy'),
         'sample13':('Healthy','Healthy'),'sample14':('AD','Non-lesional'),'sample15':('AD','Non-lesional'),
         'sample16':('AD','Non-lesional'),'sample17':('Healthy','Healthy')}

res = []
base = '/workspace/gse147424'
for f in sorted(os.listdir(base)):
    if not f.endswith('.gz'): continue
    gsm = f.split('_')[0]
    sname = re.search(r'(sample\d+)', f).group(1)
    path = os.path.join(base, f)
    with gzip.open(path, 'rt') as fh:
        header = fh.readline().rstrip('\n').split(',')
        n_cells = len(header) - 1
        itpka_vals = None
        for line in fh:
            if line.startswith('"ITPKA"'):
                parts = line.rstrip('\n').split(',')
                itpka_vals = np.array([float(v) for v in parts[1:] if v != ''])
                break
    dis, tt = group[sname]
    if itpka_vals is None:
        res.append({'sample': sname, 'gsm': gsm, 'disease': dis, 'tissue_type': tt,
                    'n_cells': n_cells, 'n_expr': 0, 'pct_expr': 0.0, 'mean_expr_in_expr': 0.0,
                    'note': 'ITPKA filtered out (below detection)'})
    else:
        nz = itpka_vals[itpka_vals > 0]
        res.append({'sample': sname, 'gsm': gsm, 'disease': dis, 'tissue_type': tt,
                    'n_cells': len(itpka_vals), 'n_expr': len(nz),
                    'pct_expr': round(100*len(nz)/len(itpka_vals), 3),
                    'mean_expr_in_expr': round(nz.mean(), 3) if len(nz) else 0.0,
                    'note': ''})

he = pd.DataFrame(res)
print(he.to_string(index=False))
print('\n=== 按组汇总 ===')
summ = he.groupby(['disease','tissue_type']).agg(
    n_samples=('sample','size'), total_cells=('n_cells','sum'),
    total_expr=('n_expr','sum')).reset_index()
summ['pct_expr_pooled'] = (100*summ['total_expr']/summ['total_cells']).round(3)
print(summ.to_string(index=False))
he.to_csv('/mnt/results/P2-02_单细胞_He2020_AD皮损ITPKA检出.csv', index=False)
print('\nsaved: P2-02_单细胞_He2020_AD皮损ITPKA检出.csv')