def fetch_dataset3(ds_id, rsids, batch=64):
    out = []
    for i in range(0, len(rsids), batch):
        chunk = rsids[i:i+batch]
        for attempt in range(4):
            try:
                r = requests.post(f'{API}/associations', headers=H,
                                  json={'variant': chunk, 'id': [ds_id]}, timeout=60)
                if r.status_code == 200:
                    out.extend(r.json()); break
            except Exception:
                pass
            time.sleep(1.5 * (attempt + 1))
        else:
            print(f'  FAIL batch {i} for {ds_id}')
        time.sleep(0.15)
    return pd.DataFrame(out)

for sym, ensg in NEW_GENES.items():
    df = fetch_dataset3(f'eqtl-a-{ensg}', rsids)
    hit = df[df['rsid'] == 'rs11635906']
    df.to_csv(f'/workspace/coloc_panel/eqtlgen_{sym}.csv', index=False)
    print(sym, len(df), 'rows | rs11635906 p =', hit['p'].iloc[0] if len(hit) else 'MISSING')