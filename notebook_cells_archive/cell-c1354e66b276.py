# 全局搜 dermatitis/eczema 标签（不限组织）
alld = obs.read(column_names=["disease"]).concat().to_pandas()
labels = sorted(alld['disease'].unique())
print([x for x in labels if any(k in x.lower() for k in ['dermatitis','eczema','atopic'])])