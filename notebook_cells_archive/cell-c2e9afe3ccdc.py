r2 = requests.post(f'{API}/associations', headers=H,
                   json={'variant': ['rs11635906'], 'id': ['eqtl-a-ENSG00000137825']}, timeout=60)
print('POST associations ->', r2.status_code)
d = r2.json()
print(json.dumps(d[0] if isinstance(d, list) else d, indent=1)[:600])