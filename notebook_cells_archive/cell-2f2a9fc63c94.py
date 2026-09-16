import cellxgene_census, pandas as pd

with cellxgene_census.open_soma(census_version='2025-01-30') as census:
    # 软骨细胞：皮肤组织内的全部 chondrocyte，取 obs + ITPKA 表达
    obs = cellxgene_census.get_obs(
        census, 'homo_sapiens',
        value_filter="tissue_general == 'skin' and cell_type == 'chondrocyte'",
        column_names=['soma_joinid','dataset_id','donor_id','cell_type','tissue','tissue_general',
                      'assay','suspension_type','disease','development_stage','sex','self_reported_ethnicity'])
    print('chondrocytes in skin:', len(obs))
    print('\n=== disease ==='); print(obs['disease'].value_counts())
    print('\n=== tissue (具体解剖部位) ==='); print(obs['tissue'].value_counts())
    print('\n=== dataset 分布 ==='); print(obs['dataset_id'].value_counts())
    print('\n=== assay ==='); print(obs['assay'].value_counts())
    print('\n=== development_stage ==='); print(obs['development_stage'].value_counts().head(8))
    print('\nn_donors:', obs['donor_id'].nunique())