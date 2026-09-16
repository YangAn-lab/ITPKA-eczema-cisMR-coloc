import cellxgene_census
import pandas as pd

CENSUS_VERSION = "2025-01-30"
census = cellxgene_census.open_soma(census_version=CENSUS_VERSION)
obs = census["census_data"]["homo_sapiens"].obs

# 皮肤细胞的疾病标签与细胞量
skin = obs.read(value_filter="tissue_general == 'skin' and is_primary_data == True",
                column_names=["disease", "cell_type", "donor_id", "dataset_id"]).concat().to_pandas()
print('skin cells:', len(skin))
print(skin['disease'].value_counts().head(15))