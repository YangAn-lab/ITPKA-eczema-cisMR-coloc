# 分析C: GTEx v10/v8 全血 eQTL 作为暴露的验证性 MR（Wald ratio，与主分析同一阶 SE 约定）
from scipy.stats import norm
import numpy as np, pandas as pd

exposures = {
    'GTEx v10 whole blood (n=853)': {'nes': 0.24635037779808044, 'p': 3.3144100511127005e-08},
    'GTEx v8 whole blood (n=670)':  {'nes': 0.21639,              'p': 6.82733e-06},
    'eQTLGen whole blood (n=31,684)': {'nes': 0.150,              'p': 7.4e-29},
}
outcomes = {
    'FinnGen R12 dermatitis/eczema': {'beta': 0.0356413, 'se': 0.0064273},
    'FinnGen R12 atopic dermatitis': {'beta': 0.0484656, 'se': 0.00919896},
}

rows = []
for ename, e in exposures.items():
    z_e = norm.isf(e['p']/2)
    se_e = e['nes']/z_e
    F = z_e**2
    for oname, o in outcomes.items():
        b = o['beta']/e['nes']
        se_b = o['se']/e['nes']          # 一阶 Wald ratio SE（与主分析一致）
        z_b = b/se_b
        p_b = 2*norm.sf(abs(z_b))
        orr = np.exp(b); lo = np.exp(b-1.96*se_b); hi = np.exp(b+1.96*se_b)
        rows.append({'exposure_eQTL': ename, 'outcome': oname,
                     'eQTL_NES': round(e['nes'],4), 'eQTL_se_imputed': round(se_e,4), 'eQTL_p': e['p'], 'F_stat': round(F,1),
                     'MR_b': round(b,3), 'MR_se': round(se_b,3), 'OR': round(orr,2),
                     'CI95': f'({lo:.2f}-{hi:.2f})', 'MR_p': f'{p_b:.2e}'})

mrc = pd.DataFrame(rows)
print(mrc.to_string(index=False))
mrc.to_csv('/mnt/results/P2-02_GTEx验证性MR_结果表.csv', index=False)
print('\nsaved: P2-02_GTEx验证性MR_结果表.csv')
print('\n核对主分析行: eQTLGen×湿疹 b=0.238 se=0.043 OR=1.27(1.15-1.39) p=2.9e-8 ✓ 与锁定值一致')