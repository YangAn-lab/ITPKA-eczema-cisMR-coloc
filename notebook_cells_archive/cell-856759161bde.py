import requests, time, json
import pandas as pd

JWT = open('/workspace/.opengwas_jwt').read().strip()
HEAD = {'Authorization': f'Bearer {JWT}', 'Content-Type': 'application/json'}
snps = pd.read_csv('/workspace/coloc_panel/snp_universe.csv')
rsids = snps.rsid.tolist()

genes = {'ITPKA': 'ENSG00000137825', 'NDUFAF1': 'ENSG00000137806', 'CHP1': 'ENSG00000187446',
         'RPAP1': 'ENSG00000103932', 'RTF1': 'ENSG00000137815', 'OIP5-AS1': 'ENSG00000247556'}

def fetch_gene(sym, gid, chunk=64):
    ds = f'eqtl-a-{gid}'
    out, fails = [], 0
    for i in range(0, len(rsids), chunk):
        vs = rsids[i:i+chunk]
        ok = False
        for attempt in range(4):
            try:
                r = requests.post('https://api.opengwas.io/api/associations', headers=HEAD,
                                  json={'variant': vs, 'id': [ds]}, timeout=120)
                if r.status_code == 200:
                    out.extend(r.json()); ok = True; break
                elif r.status_code == 429:
                    time.sleep(8 * (attempt + 1))
                else:
                    print(sym, 'chunk', i, 'HTTP', r.status_code, r.text[:120]); break
            except Exception as e:
                time.sleep(4)
        if not ok: fails += 1
        time.sleep(0.25)
    df = pd.DataFrame(out)
    df.to_csv(f'/workspace/coloc_panel/eqtlgen_{sym}.csv', index=False)
    return df, fails

for sym, gid in genes.items():
    df, fails = fetch_gene(sym, gid)
    if len(df):
        hit = df[df.rsid == 'rs11635906']
        b = f"beta={hit.beta.values[0]:.4f}, p={hit.p.values[0]:.2e}" if len(hit) else 'rs11635906 缺失'
        print(f'{sym}: {len(df)} 行 (失败块 {fails}), {b}')
    else:
        print(f'{sym}: 0 行 (失败块 {fails})')