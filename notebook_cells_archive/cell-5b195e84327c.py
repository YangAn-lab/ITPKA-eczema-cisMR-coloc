# A2: 按细胞类型汇总 ITPKA 表达（正常皮肤，672,457 细胞 / 233 donors）
import pandas as pd

df = ad_norm.obs.copy()
df['expr'] = x

# 细胞类型层面汇总
ct = (df.groupby('cell_type')
        .agg(n_cells=('expr','size'),
             pct_expr=('expr', lambda s: 100*(s>0).mean()),
             mean_norm=('expr','mean'),
             mean_expr_cells=('expr', lambda s: s[s>0].mean() if (s>0).any() else 0.0),
             n_donors=('donor_id','nunique'))
        .sort_values('pct_expr', ascending=False))
ct = ct[ct['n_cells'] >= 500]  # 至少500细胞的类型才稳健
pd.set_option('display.width', 200)
print(ct.head(25).round(3).to_string())
print('\n总细胞类型数(>=500 cells):', len(ct))