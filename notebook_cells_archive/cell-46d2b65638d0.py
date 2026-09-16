import pandas as pd, requests, time

univ = pd.read_csv('/workspace/coloc_panel/snp_universe.csv')
rsids = univ['rsid'].dropna().unique().tolist()
print(f'{len(rsids)} rsids in SNP universe')

NEW_GENES = {  # symbol -> ensembl
    'NUSAP1': 'ENSG00000137804', 'OIP5': 'ENSG00000104147', 'PLA2G4B': 'ENSG00000243708',
    'EHD4': 'ENSG00000103966', 'LINC_260814': 'ENSG00000260814', 'JMJD7': 'ENSG00000243789'}

def fetch_dataset(ds_id, rsids, batch=400):
    out = []
    for i in range(0, len(rsids), batch):
        chunk = rsids[i:i+batch]
        for attempt in range(3):
            try:
                r = requests.post(f'{API}/associations', headers=H,
                                  json={'variant': chunk, 'id': [ds_id]}, timeout=180)
                if r.status_code == 200:
                    out.extend(r.json()); break
            except Exception:
                pass
            time.sleep(3 * (attempt + 1))
        else:
            print(f'  WARN batch {i} failed for {ds_id}')
        time.sleep(0.3)
    return pd.DataFrame(out)

for sym, ensg in NEW_GENES.items():
    df = fetch_dataset(f'eqtl-a-{ensg}', rsids)
    df.to_csv(f'/workspace/coloc_panel/eqtlgen_{sym}.csv', index=False)
    print(sym, len(df), 'rows, min p =', df['p'].min() if len(df) else 'NA')