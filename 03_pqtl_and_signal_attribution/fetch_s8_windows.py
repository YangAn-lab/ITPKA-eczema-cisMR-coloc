#!/usr/bin/env python3
"""C1-2: 从 OpenGWAS 批量拉取 S8/血细胞性状在位点窗口内的关联（universe rsid 列表）"""
import requests, json, time, sys
import pandas as pd

JWT = open('/workspace/.opengwas_jwt').read().strip()
HDR = {"Authorization": f"Bearer {JWT}", "Content-Type": "application/json"}
URL = "https://api.opengwas.io/api/associations"

univ = pd.read_csv('/workspace/coloc_panel/snp_universe.csv')
rsids = univ['rsid'].dropna().unique().tolist()
print(f"universe rsids: {len(rsids)}")

TRAITS = {
    "ebi-a-GCST90014003": "cystatinC",
    "ebi-a-GCST90029008": "height",
    "ebi-a-GCST90000025": "leanmass",
    "ukb-b-17422":        "menopause",
    "ebi-a-GCST90026654": "egfr",
    "ebi-a-GCST90002351": "neutrophil",
    "ebi-a-GCST90002389": "lymphocyte_pct",
    "ebi-a-GCST90028992": "eosinophil_count",
    "ebi-a-GCST90002382": "eosinophil_pct",
}

BATCH = 64
for tid, name in TRAITS.items():
    rows = []
    for i in range(0, len(rsids), BATCH):
        chunk = rsids[i:i+BATCH]
        for attempt in range(4):
            try:
                r = requests.post(URL, headers=HDR, json={"variant": chunk, "id": [tid]}, timeout=180)
                if r.status_code == 200:
                    rows.extend(r.json()); break
                elif r.status_code == 429:
                    wait = 20 * (attempt + 1); print(f"  429 限流, {wait}s 后重试"); time.sleep(wait)
                else:
                    print(f"  HTTP {r.status_code}: {r.text[:150]}"); time.sleep(5)
            except Exception as e:
                print(f"  异常 {e}; 重试"); time.sleep(10)
        time.sleep(0.4)
    df = pd.DataFrame(rows)
    out = f"/workspace/fable_fix/pqtl/s8win_{name}.csv"
    df.to_csv(out, index=False)
    print(f"{name} ({tid}): {len(df)} 行 -> {out}")
print("ALL DONE")
