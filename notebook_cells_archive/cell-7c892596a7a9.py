import cellxgene_census, pandas as pd

with cellxgene_census.open_soma(census_version='2025-01-30') as census:
    obs = cellxgene_census.get_obs(
        census, 'homo_sapiens',
        value_filter="tissue_general == 'skin of body' and cell_type == 'chondrocyte' and disease == 'normal' and is_primary_data == True",
        column_names=['soma_joinid','dataset_id','donor_id','cell_type','tissue','tissue_general',
                      'assay','suspension_type','disease','development_stage','sex'])
print('chondrocytes (skin of body, normal, primary):', len(obs))
print('\n=== tissue 解剖部位 ==='); print(obs['tissue'].value_counts())
print('\n=== dataset 分布 ==='); print(obs['dataset_id'].value_counts())
print('\n=== assay ==='); print(obs['assay'].value_counts())
print('\n=== development_stage ==='); print(obs['development_stage'].value_counts())
print('\nn_donors:', obs['donor_id'].nunique())
obs.to_csv('/workspace/coloc_oliva/census_chondrocyte_obs.csv', index=False)