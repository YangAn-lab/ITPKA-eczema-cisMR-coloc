chunk = rsids[0:100]
r = requests.post(f'{API}/associations', headers=H, json={'variant': chunk, 'id': ['eqtl-a-ENSG00000137804']}, timeout=120)
print('status:', r.status_code)
print('error body:', r.text[:500])
# bisect: try first 50
r50 = requests.post(f'{API}/associations', headers=H, json={'variant': chunk[:50], 'id': ['eqtl-a-ENSG00000137804']}, timeout=120)
print('first 50:', r50.status_code, len(r50.json()) if r50.status_code == 200 else r50.text[:150])
r50b = requests.post(f'{API}/associations', headers=H, json={'variant': chunk[50:], 'id': ['eqtl-a-ENSG00000137804']}, timeout=120)
print('last 50:', r50b.status_code, len(r50b.json()) if r50b.status_code == 200 else r50b.text[:150])