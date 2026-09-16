import requests
JWT = open('/workspace/.opengwas_jwt').read().strip()
H = {'Authorization': f'Bearer {JWT}'}
for url in ['https://api.opengwas.io/api/gwasinfo/eqtl-a-ENSG00000137825',
            'https://api.opengwas.io/api/associations/rs11635906/eqtl-a-ENSG00000137825']:
    r = requests.get(url, headers=H, timeout=30)
    print(url.split('/api/')[1], '->', r.status_code, '|', r.text[:200].replace('\n', ' '))
# also check whoami/token status
r = requests.get('https://api.opengwas.io/api/profile', headers=H, timeout=30)
print('profile ->', r.status_code, '|', r.text[:200])