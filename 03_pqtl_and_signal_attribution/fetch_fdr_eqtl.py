import requests, time
import pandas as pd
JWT = open('/workspace/.opengwas_jwt').read().strip()
HDR = {"Authorization": f"Bearer {JWT}", "Content-Type": "application/json"}
URL = "https://api.opengwas.io/api/associations"
LOCI = [
    ("FGF2",   "eqtl-a-ENSG00000138685", "fg_FGF2_atopic.tsv"),
    ("CSF2RB", "eqtl-a-ENSG00000100368", "fg_CSF2RB_eczema.tsv"),
    ("EDN1",   "eqtl-a-ENSG00000078401", "fg_EDN1_eczema.tsv"),
]
for gene, eid, gwas_file in LOCI:
    gw = pd.read_csv(gwas_file, sep='\t', header=None,
                     names=['chr','pos','ref','alt','rsid','gene','p','mlogp','beta','se','af','af_c','af_ct'])
    rsids = gw['rsid'].dropna().unique().tolist()
    rows = []
    for i in range(0, len(rsids), 64):
        chunk = rsids[i:i+64]
        for attempt in range(6):
            try:
                r = requests.post(URL, headers=HDR, json={"variant": chunk, "id": [eid]}, timeout=180)
                if r.status_code == 200:
                    rows.extend(r.json()); break
                elif r.status_code == 429:
                    time.sleep(30*(attempt+1))
                else:
                    print(f"  HTTP {r.status_code}: {r.text[:100]}", flush=True); time.sleep(8)
            except Exception as e:
                print(f"  异常 {e}", flush=True); time.sleep(15)
        time.sleep(0.5)
        if (i//64) % 20 == 0:
            print(f"  {gene}: {i}/{len(rsids)}", flush=True)
    pd.DataFrame(rows).to_csv(f"eqtl_{gene}_window.csv", index=False)
    print(f"{gene}: {len(rows)} 行", flush=True)
print("ALL DONE", flush=True)
