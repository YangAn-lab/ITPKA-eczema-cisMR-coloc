import pandas as pd, requests, time

gdf = pd.read_csv('/workspace/stageA_window_genes.csv')
cand = gdf[gdf['biotype'].isin(['protein_coding', 'lncRNA']) & gdf['external_name'].notna()].copy()
dsids = [f'eqtl-a-{g}' for g in cand['id']]

r = requests.post(f'{API}/gwasinfo', headers=H, json={'id': dsids}, timeout=60)
existing = {d['id'] for d in r.json()}
print(f'{len(existing)}/{len(dsids)} eQTLGen datasets exist')

r = requests.post(f'{API}/associations', headers=H,
                  json={'variant': ['rs11635906'], 'id': sorted(existing)}, timeout=120)
recs = r.json()
scr = pd.DataFrame(recs)[['id', 'beta', 'se', 'p', 'eaf', 'n']]
scr['ensg'] = scr['id'].str.replace('eqtl-a-', '')
scr = scr.merge(cand[['id', 'external_name', 'start', 'end']].rename(columns={'id': 'ensg'}), on='ensg')
scr = scr.sort_values('p').reset_index(drop=True)
scr.to_csv('/workspace/stageA_eqtlgen_screen.csv', index=False)
pd.set_option('display.float_format', lambda x: f'{x:.3g}')
print(scr[['external_name', 'beta', 'p', 'eaf', 'n']].to_string(index=False))