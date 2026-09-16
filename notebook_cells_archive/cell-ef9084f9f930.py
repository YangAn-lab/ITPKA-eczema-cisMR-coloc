import requests, json, pandas as pd

# ±1 Mb window around rs11635906 (hg38 chr15:41,487,062)
lo, hi = 40487062, 42487062
r = requests.get(f'https://rest.ensembl.org/overlap/region/human/15:{lo}-{hi}?feature=gene',
                 headers={'Content-Type': 'application/json'}, timeout=60)
r.raise_for_status()
genes = r.json()
gdf = pd.DataFrame(genes)[['id', 'external_name', 'start', 'end', 'strand', 'biotype']]
gdf = gdf.sort_values('start').reset_index(drop=True)
print(f'{len(gdf)} genes in chr15:{lo}-{hi} (hg38, ±1 Mb)')
print(gdf.to_string())
gdf.to_csv('/workspace/stageA_window_genes.csv', index=False)