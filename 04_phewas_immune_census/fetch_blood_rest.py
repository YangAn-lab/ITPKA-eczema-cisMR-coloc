import requests, time
import pandas as pd
JWT = open('/workspace/.opengwas_jwt').read().strip()
HDR = {"Authorization": f"Bearer {JWT}", "Content-Type": "application/json"}
URL = "https://api.opengwas.io/api/associations"
univ = pd.read_csv('/workspace/coloc_panel/snp_universe.csv')
rsids = univ['rsid'].dropna().unique().tolist()
TRAITS = {
    "ebi-a-GCST90002389": "lymphocyte_pct",
    "ebi-a-GCST90028992": "eosinophil_count",
    "ebi-a-GCST90002382": "eosinophil_pct",
}
for tid, name in TRAITS.items():
    rows = []
    for i in range(0, len(rsids), 64):
        chunk = rsids[i:i+64]
        for attempt in range(6):
            try:
                r = requests.post(URL, headers=HDR, json={"variant": chunk, "id": [tid]}, timeout=180)
                if r.status_code == 200:
                    rows.extend(r.json()); break
                elif r.status_code == 429:
                    time.sleep(30 * (attempt + 1))
                else:
                    print(f"  HTTP {r.status_code}: {r.text[:120]}", flush=True); time.sleep(8)
            except Exception as e:
                print(f"  异常 {e}", flush=True); time.sleep(15)
        time.sleep(0.8)
    pd.DataFrame(rows).to_csv(f"/workspace/fable_fix/pqtl/s8win_{name}.csv", index=False)
    print(f"{name}: {len(rows)} 行", flush=True)
print("ALL DONE", flush=True)
