def fetch_dataset2(ds_id, rsids, batch=100):
    out = []
    for i in range(0, len(rsids), batch):
        chunk = rsids[i:i+batch]
        ok = False
        for attempt in range(4):
            try:
                r = requests.post(f'{API}/associations', headers=H,
                                  json={'variant': chunk, 'id': [ds_id]}, timeout=120)
                if r.status_code == 200:
                    out.extend(r.json()); ok = True; break
            except Exception:
                pass
            time.sleep(2 * (attempt + 1))
        if not ok:
            print(f'  FAIL batch {i} for {ds_id}')
        time.sleep(0.25)
    return pd.DataFrame(out)

for sym, ensg in NEW_GENES.items():
    df = fetch_dataset2(f'eqtl-a-{ensg}', rsids)
    # sanity: rs11635906 must be present with the screened p-value
    hit = df[df['rsid'] == 'rs11635906']
    df.to_csv(f'/workspace/coloc_panel/eqtlgen_{sym}.csv', index=False)
    print(sym, len(df), 'rows | rs11635906 p =', hit['p'].iloc[0] if len(hit) else 'MISSING')