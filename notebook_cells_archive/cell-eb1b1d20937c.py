# 分析D 补充: 各基因 GTEx lead 变异与 rs11635906 的 LD（1000G EUR 本地矩阵）
import pandas as pd, numpy as np

snps = pd.read_csv('/workspace/coloc_panel/ld_eur_snps.csv')
ld = np.load('/workspace/coloc_panel/ld_eur.npy')
print('LD matrix:', ld.shape, '| snps:', len(snps), '| cols:', list(snps.columns))

# 找到 rs11635906 与各 lead 的索引
snps = snps.reset_index(drop=True)
def find_idx(rsid):
    m = snps[snps.apply(lambda r: rsid in str(r.values), axis=1)]
    return m.index.tolist()

targets = ['rs11635906','rs12440111','rs146883891','rs13329233','rs2297379','rs9302109','rs11633799']
idx_map = {}
for t in targets:
    ix = find_idx(t)
    idx_map[t] = ix[0] if ix else None
    print(t, '-> idx', idx_map[t])

if idx_map.get('rs11635906') is not None:
    i0 = idx_map['rs11635906']
    print('\n=== r 与 rs11635906 ===')
    for gene, lead in [('ITPKA','rs12440111'),('NDUFAF1','rs146883891'),('CHP1','rs13329233'),
                       ('RPAP1','rs2297379'),('RTF1','rs9302109'),('OIP5-AS1','rs11633799')]:
        j = idx_map.get(lead)
        if j is not None:
            r = ld[i0, j]
            print(f'{gene:10s} lead {lead:14s} r={r:+.3f}  r2={r**2:.3f}')
        else:
            print(f'{gene:10s} lead {lead:14s} 不在 LD 矩阵中')