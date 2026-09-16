# A3: 检查表达值量级 + 用 raw counts 重算（normalized 层值过小，疑似部分数据集无 normalized 层回退）
print('nonzero values stats:')
nz = x[x>0]
print('n nonzero:', len(nz), '| min/median/max:', nz.min(), np.median(nz), nz.max())

# 用 raw 层重取，自行 log1p 标准化以便跨数据集可比
ad_raw = cellxgene_census.get_anndata(
    census, organism="Homo sapiens",
    obs_value_filter="tissue_general == 'skin of body' and disease == 'normal' and is_primary_data == True",
    var_value_filter=f"feature_id == '{ITPKA_ENS}'",
    obs_column_names=["cell_type", "donor_id", "dataset_id"],
    X_name="raw",
)
xr = np.asarray(ad_raw.X.todense()).ravel()
print('\nraw layer: nonzero n=%d, max=%.1f' % ((xr>0).sum(), xr.max()))
print('raw pct>0 = %.2f%%' % (100*(xr>0).mean()))