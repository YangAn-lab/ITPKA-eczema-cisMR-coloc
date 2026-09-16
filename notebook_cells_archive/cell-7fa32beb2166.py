from scipy.stats import norm
import numpy as np

# 锁定输入值
# GWAS 侧（FinnGen R12 tabix 实测）
gwas = {
    'eczema': dict(beta=0.0356413, se=0.0064273),
    'AD':     dict(beta=0.0484656, se=0.00919896),
}
# 暴露侧（GTEx Portal API 实测 NES 与 p；SE 由 p 反推）
expo = {
    'v10': dict(nes=0.2464, p=3.314e-08),
    'v8':  dict(nes=0.2164, p=6.827e-06),
}

rows = []
for ev, e in expo.items():
    z_e = norm.isf(e['p']/2)          # 双侧 p -> z
    se_e = e['nes']/z_e
    for out, g in gwas.items():
        b = g['beta']/e['nes']
        se1 = g['se']/e['nes']                       # 一阶 SE
        p1 = 2*norm.sf(abs(b/se1))
        # delta 法（二阶）：Var(b) ≈ Var(G)/E² + G²·Var(E)/E⁴
        se2 = np.sqrt(g['se']**2/e['nes']**2 + g['beta']**2*se_e**2/e['nes']**4)
        p2 = 2*norm.sf(abs(b/se2))
        rows.append((ev, out, round(b,3), round(se1,4), f'{p1:.1e}', round(se2,4), f'{p2:.1e}'))

import pandas as pd
df = pd.DataFrame(rows, columns=['exposure','outcome','b','SE_1st','p_1st','SE_delta','p_delta'])
print(df.to_string(index=False))
print('\ndelta 法最大 p 值:', max(float(r[6]) for r in rows))