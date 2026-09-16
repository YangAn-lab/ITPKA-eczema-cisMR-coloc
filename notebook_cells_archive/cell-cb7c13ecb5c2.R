# debug: small batch + inspect what the 15-row fallback was
test_rsids = rsids[:10] + ['rs11635906']
r = requests.post(f'{API}/associations', headers=H,
                  json={'variant': test_rsids, 'id': ['eqtl-a-ENSG00000137804']}, timeout=60)
print('batch of 11 ->', r.status_code, 'rows:', len(r.json()) if r.status_code == 200 else r.text[:200])
if r.status_code == 200:
    d = pd.DataFrame(r.json())
    print(d[['rsid', 'p']].to_string(index=False))

# what was in the failed NUSAP1 pull?
prev = pd.read_csv('/workspace/coloc_panel/eqtlgen_NUSAP1.csv')
print('\nprevious 15-row file:'); print(prev[['rsid', 'position', 'p']].head(15).to_string(index=False))