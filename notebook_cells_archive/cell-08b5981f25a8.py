import requests, json
JWT = open('/workspace/.opengwas_jwt').read().strip()
H = {'Authorization': f'Bearer {JWT}', 'Content-Type': 'application/json'}
API = 'https://api.opengwas.io/api'

# correct endpoints: POST /gwasinfo and POST /associations
r = requests.post(f'{API}/gwasinfo', headers=H, json={'id': ['eqtl-a-ENSG00000137825']}, timeout=30)
print('POST gwasinfo ->', r.status_code, '|', r.text[:250])
r2 = requests.post(f'{API}/associations', headers=H,
                   json={'variant': ['rs11635906'], 'dataset': ['eqtl-a-ENSG00000137825']}, timeout=60)
print('POST associations ->', r2.status_code, '|', r2.text[:350])