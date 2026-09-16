#!/usr/bin/env python3
# P2-02 Stage B2: pull window summary stats for top non-atopic PheWAS traits from OpenGWAS.
# POST /associations, constraint N(id) x N(variant) <= 64 -> 1 id x 64 variants per call.
import json, time, urllib.request
import pandas as pd

DIR = '/workspace/coloc_panel'
JWT = open('/workspace/.opengwas_jwt').read().strip()
rsids = pd.read_csv(f'{DIR}/snp_universe.csv')['rsid'].dropna().unique().tolist()
print(f'{len(rsids)} rsids in universe')

TRAITS = {
    'height_GCST90029008':      'ebi-a-GCST90029008',   # p=3.2e-20, n=673878
    'leanmass_GCST90000025':    'ebi-a-GCST90000025',   # p=9.4e-17, n=450243
    'menopause_ukb-b-17422':    'ukb-b-17422',          # p=1.6e-15, n=143819
    'cystatinC_GCST90014003':   'ebi-a-GCST90014003',   # p=8.1e-15, n=389834
    'egfr_GCST90026654':        'ebi-a-GCST90026654',   # p=4.3e-12, n=1159871
    'neutrophil_GCST90002351':  'ebi-a-GCST90002351',   # p=1.1e-13, n=519288 (blood control)
}

def post_batch(variants, dsid, retries=3):
    body = json.dumps({'variant': variants, 'id': [dsid]}).encode()
    req = urllib.request.Request('https://api.opengwas.io/api/associations', data=body,
        headers={'Authorization': f'Bearer {JWT}', 'Content-Type': 'application/json'})
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                return json.loads(r.read())
        except Exception as e:
            if attempt == retries - 1:
                print(f'  FAILED batch ({dsid}, {len(variants)} variants): {e}')
                return []
            time.sleep(2 * (attempt + 1))

for name, dsid in TRAITS.items():
    out = []
    for i in range(0, len(rsids), 64):
        batch = rsids[i:i+64]
        recs = post_batch(batch, dsid)
        for r in recs:
            out.append({k: r.get(k) for k in ['id','trait','chr','position','rsid','ea','nea','eaf','beta','se','p','n']})
        time.sleep(0.15)
    df = pd.DataFrame(out)
    df.to_csv(f'{DIR}/phewas_trait_{name}.csv', index=False)
    hit = df[df.rsid == 'rs11635906']
    print(f'{name}: {len(df)}/{len(rsids)} variants | rs11635906: ' +
          (f"beta={hit.iloc[0].beta:+.5f} p={hit.iloc[0].p:.1e}" if len(hit) else 'NOT RETURNED'))
