import requests, time, pandas as pd

JWT = open('/workspace/.opengwas_jwt').read().strip()
H = {'Authorization': f'Bearer {JWT}'}
API = 'https://api.opengwas.io/api'

gdf = pd.read_csv('/workspace/stageA_window_genes.csv')
cand = gdf[gdf['biotype'].isin(['protein_coding', 'lncRNA']) & gdf['external_name'].notna()].copy()
print(f'{len(cand)} named protein-coding/lncRNA genes to screen')

# Step 1: does an eQTLGen dataset exist in OpenGWAS? Step 2: rs11635906 association p-value
rows = []
for _, g in cand.iterrows():
    dsid = f"eqtl-a-{g['id']}"
    try:
        gi = requests.get(f'{API}/gwasinfo/{dsid}', headers=H, timeout=30)
        exists = gi.status_code == 200 and isinstance(gi.json(), list) and len(gi.json()) > 0
    except Exception:
        exists = False
    pval, beta = None, None
    if exists:
        try:
            a = requests.get(f'{API}/associations/rs11635906/{dsid}', headers=H, timeout=30)
            if a.status_code == 200:
                recs = a.json()
                if recs:
                    pval, beta = recs[0].get('p'), recs[0].get('beta')
        except Exception:
            pass
    rows.append({'gene_id': g['id'], 'symbol': g['external_name'], 'biotype': g['biotype'],
                 'opengwas_exists': exists, 'rs11635906_p': pval, 'rs11635906_beta': beta})
    time.sleep(0.15)

scr = pd.DataFrame(rows)
scr.to_csv('/workspace/stageA_eqtlgen_screen.csv', index=False)
have = scr[scr['opengwas_exists']]
print(f'\n{len(have)} genes have eQTLGen datasets in OpenGWAS')
print(have.sort_values('rs11635906_p')[['symbol','rs11635906_p','rs11635906_beta']].to_string(index=False))