import cellxgene_census, pandas as pd, numpy as np

with cellxgene_census.open_soma(census_version='2025-01-30') as census:
    ad = cellxgene_census.get_anndata(
        census, organism='homo_sapiens',
        obs_value_filter="tissue_general == 'skin of body' and cell_type == 'chondrocyte' and disease == 'normal' and is_primary_data == True",
        var_value_filter="feature_name == 'ITPKA'",
        obs_column_names=['donor_id','tissue'])
X = ad.X.toarray().flatten() if hasattr(ad.X,'toarray') else np.asarray(ad.X).flatten()
print('ITPKA+ chondrocytes: %d/%d (%.2f%%), mean UMI=%.2f' % ((X>0).sum(), len(X), 100*(X>0).mean(), X.mean()))
print('占总皮肤细胞比例: 2139/672457 = %.3f%%' % (100*2139/672457))
# 按解剖部位拆分 ITPKA 表达
df = pd.DataFrame({'tissue': ad.obs['tissue'].values, 'expr': X})
print(df.groupby('tissue')['expr'].agg(n='count', n_pos=lambda s:(s>0).sum(), pct=lambda s:100*(s>0).mean()).round(2).sort_values('n', ascending=False).head(6))