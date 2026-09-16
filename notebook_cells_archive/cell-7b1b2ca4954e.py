skin = obs.read(value_filter="tissue_general == 'skin of body' and is_primary_data == True",
                column_names=["disease", "cell_type", "donor_id", "dataset_id"]).concat().to_pandas()
print('skin cells:', len(skin))
print(skin['disease'].value_counts().head(12))
print('---AD 细胞类型 top20---')
ad = skin[skin['disease']=='atopic dermatitis']
print('AD cells:', len(ad), '| donors:', ad['donor_id'].nunique(), '| datasets:', ad['dataset_id'].nunique())
print(ad['cell_type'].value_counts().head(20))