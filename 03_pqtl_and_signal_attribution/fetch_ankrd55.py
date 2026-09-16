#!/usr/bin/env python3
# 从 OpenGWAS 拉取 ANKRD55 (ENSG00000164512) 在 IL6ST 窗口的 eQTLGen 关联
import requests, json, time, pandas as pd

JWT = open('/workspace/.opengwas_jwt').read().strip()
H = {'Authorization': f'Bearer {JWT}', 'Content-Type': 'application/json'}
URL = 'https://api.opengwas.io/api/associations'

il6st = pd.read_csv('/workspace/scan498/eqtlgen_IL6ST.csv')
variants = il6st['rsid'].dropna().unique().tolist()
print(f"variants to fetch: {len(variants)}")

rows = []
B = 60  # N(id)xN(variant) <= 64 上限留余量
for i in range(0, len(variants), B):
    batch = variants[i:i+B]
    for attempt in range(4):
        r = requests.post(URL, headers=H, json={'id': ['eqtl-a-ENSG00000164512'], 'variant': batch}, timeout=60)
        if r.status_code == 200:
            break
        time.sleep(3 * (attempt + 1))
    if r.status_code != 200:
        print(f"batch {i}: HTTP {r.status_code} {r.text[:100]}")
        continue
    rows.extend(r.json())
    if (i // B) % 10 == 0:
        print(f"{i}/{len(variants)} rows={len(rows)}")
    time.sleep(0.4)

d = pd.DataFrame(rows)
d.to_csv('/workspace/fable_fix/pqtl/eqtlgen_ANKRD55_il6st_window.csv', index=False)
print("saved", d.shape)
if len(d):
    print(d[['rsid','beta','se','p','n']].assign(absz=(d['beta']/d['se']).abs()).sort_values('absz', ascending=False).head(3).to_string(index=False))
