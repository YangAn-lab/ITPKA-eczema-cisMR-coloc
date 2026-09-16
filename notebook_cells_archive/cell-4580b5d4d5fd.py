# A1: Census 正常皮肤 ITPKA 各细胞类型表达（修正：X_name 而非 layer）
import numpy as np

ITPKA_ENS = "ENSG00000137825"
# 先确认 var 中的 feature_id
var = census["census_data"]["homo_sapiens"]["ms"]["RNA"].var.read().concat().to_pandas()
hit = var[var['feature_name'] == 'ITPKA']
print(hit[['feature_id','feature_name']])

norm_skin = skin[skin['disease'] == 'normal']
print('normal skin cells:', len(norm_skin), '| donors:', norm_skin['donor_id'].nunique())

ad_norm = cellxgene_census.get_anndata(
    census, organism="Homo sapiens",
    obs_value_filter="tissue_general == 'skin of body' and disease == 'normal' and is_primary_data == True",
    var_value_filter=f"feature_id == '{ITPKA_ENS}'",
    obs_column_names=["cell_type", "donor_id", "dataset_id"],
    X_name="normalized",
)
print(ad_norm.shape)
x = np.asarray(ad_norm.X.todense()).ravel()
print('X stats: mean=%.4f, pct>0=%.2f%%' % (x.mean(), 100*(x>0).mean()))