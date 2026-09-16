#!/usr/bin/env python3
"""Step 1: fetch eQTLGen clumped tophits (p<5e-8) for the 83 SenMayo genes from OpenGWAS,
take lead variant per gene, liftOver hg19->hg38. Faithful to manuscript Methods 2.2."""
import csv, json, time, urllib.request, os

JWT = open('/workspace/.opengwas_jwt').read().strip()
CFG = '/workspace/pipeline_pack/mr_pipeline/config/scan_sasp_skin.csv'
OUT = '/workspace/scan498'

# parse config: unique genes
genes = {}
with open(CFG) as f:
    for row in csv.DictReader(f):
        sym = row['pair_id'].split('_')[1]
        genes[row['exposure_id']] = sym
print(f"{len(genes)} unique genes")

def post(url, payload, tries=4):
    body = json.dumps(payload).encode()
    for t in range(tries):
        try:
            req = urllib.request.Request(url, data=body, headers={
                'Authorization': f'Bearer {JWT}', 'Content-Type': 'application/json'})
            with urllib.request.urlopen(req, timeout=120) as r:
                return json.loads(r.read())
        except Exception as e:
            print(f"  retry {t+1}: {e}")
            time.sleep(5 * (t + 1))
    return None

# batch tophits, 10 ids per call
ids = list(genes)
hits = []
for i in range(0, len(ids), 10):
    batch = ids[i:i+10]
    d = post('https://api.opengwas.io/api/tophits',
             {'id': batch, 'pval': 5e-8, 'clump': 1, 'r2': 0.001, 'kb': 10000})
    if isinstance(d, list):
        hits.extend(d)
    else:
        print('ERROR batch', batch, d)
    time.sleep(1)
print(f"{len(hits)} clumped tophits rows")

# lead variant per gene (min p)
by_gene = {}
for h in hits:
    gid = h.get('id') if 'id' in h else None
    # OpenGWAS tophits rows carry 'id' = dataset id
    if gid is None:
        continue
    if gid not in by_gene or h['p'] < by_gene[gid]['p']:
        by_gene[gid] = h
print(f"{len(by_gene)} genes with >=1 tophit at 5e-8")

rows = []
for gid, sym in genes.items():
    h = by_gene.get(gid)
    if h is None:
        rows.append({'gene': sym, 'exposure_id': gid, 'rsid': '', 'chr_hg19': '', 'pos_hg19': '',
                     'ea': '', 'nea': '', 'eaf': '', 'beta_eqtl': '', 'se_eqtl': '', 'p_eqtl': '',
                     'n_eqtl': '', 'status': 'no_instrument_p5e-8'})
    else:
        rows.append({'gene': sym, 'exposure_id': gid, 'rsid': h['rsid'], 'chr_hg19': h['chr'],
                     'pos_hg19': h['position'], 'ea': h['ea'].upper(), 'nea': h['nea'].upper(),
                     'eaf': h.get('eaf', ''), 'beta_eqtl': h['beta'], 'se_eqtl': h['se'],
                     'p_eqtl': h['p'], 'n_eqtl': h.get('n', ''), 'status': 'ok'})
with open(f'{OUT}/lead_variants_hg19.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)
print(f"wrote {OUT}/lead_variants_hg19.csv; genes without instrument: "
      f"{sum(1 for r in rows if r['status'] != 'ok')}")
