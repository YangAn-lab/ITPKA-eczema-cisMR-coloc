# tissue_general 里找皮肤相关标签
all_tissues = obs.read(column_names=["tissue_general"]).concat().to_pandas()
t = all_tissues['tissue_general'].value_counts()
print([x for x in t.index if 'skin' in x.lower()])