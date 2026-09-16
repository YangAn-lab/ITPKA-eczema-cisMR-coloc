#!/usr/bin/env python3
"""Pull eQTLGen + FinnGen windows for IL1B / IL6ST loci for coloc triage of the
corrected-scan Bonferroni hits. Universe = FinnGen variants within +/-500 kb of the lead
(hg38); eQTLGen queried by chr:pos (hg19) after liftOver."""
import csv, json, time, subprocess, urllib.request, os
from pyliftover import LiftOver

JWT = open('/workspace/.opengwas_jwt').read().strip()
OUT = '/workspace/scan498'
BASE = 'https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release'
WIN = 500000

LOCI = {
    'IL1B': {'ensg': 'ENSG00000125538', 'chr': '17', 'pos38': 40010626},
    'IL6ST': {'ensg': 'ENSG00000134352', 'chr': '5', 'pos38': 56148856},
}
ENDPOINTS = ['L12_DERMATITISECZEMA', 'L12_ATOPIC', 'L12_ACTINKERA']

# reverse chain for hg38 -> hg19
if not os.path.exists(f'{OUT}/hg38ToHg19.over.chain.gz'):
    urllib.request.urlretrieve(
        'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz',
        f'{OUT}/hg38ToHg19.over.chain.gz')
lo = LiftOver(f'{OUT}/hg38ToHg19.over.chain.gz')

def tabix(ep, chrom, lo_, hi):
    url = f'{BASE}/finngen_R12_{ep}.gz'
    out = subprocess.run(['tabix', url, f'{chrom}:{lo_}-{hi}'],
                         capture_output=True, text=True, timeout=300).stdout
    rows = []
    for ln in out.splitlines():
        if ln.startswith('#') or not ln.strip():
            continue
        f = ln.split('\t')
        if len(f) >= 13:
            rows.append(f)
    return rows

def post_assoc(ds, variants, tries=4):
    body = json.dumps({'variant': variants, 'id': [ds]}).encode()
    for t in range(tries):
        try:
            req = urllib.request.Request('https://api.opengwas.io/api/associations',
                                         data=body, headers={
                'Authorization': f'Bearer {JWT}', 'Content-Type': 'application/json'})
            with urllib.request.urlopen(req, timeout=120) as r:
                return json.loads(r.read())
        except Exception as e:
            time.sleep(3 * (t + 1))
    print(f'  FAIL batch ({len(variants)} variants) for {ds}')
    return []

for sym, L in LOCI.items():
    chrom, c = L['chr'], L['pos38']
    print(f'=== {sym} chr{chrom}:{c} ===')
    univ_rows = tabix('L12_DERMATITISECZEMA', chrom, c - WIN, c + WIN)
    print(f'  universe (FinnGen eczema window): {len(univ_rows)} variants')
    # universe with hg19 mapping
    univ = []
    for f in univ_rows:
        pos = int(f[1])
        c19 = lo.convert_coordinate('chr' + chrom, pos)
        univ.append({'chr38': chrom, 'pos38': pos, 'ref': f[2], 'alt': f[3],
                     'pos19': c19[0][1] if c19 else None})
    mapped = [u for u in univ if u['pos19']]
    print(f'  lifted to hg19: {len(mapped)}/{len(univ)}')
    # eQTLGen pull by chr:pos19 in batches
    ds = f"eqtl-a-{L['ensg']}"
    q = [f"{chrom}:{u['pos19']}" for u in mapped]
    out = []
    for i in range(0, len(q), 64):
        out.extend(post_assoc(ds, q[i:i+64]))
        time.sleep(0.15)
        if (i // 64) % 20 == 0:
            print(f'    eQTLGen batch {i//64}/{(len(q)+63)//64}, rows so far {len(out)}')
    with open(f'{OUT}/eqtlgen_{sym}.csv', 'w', newline='') as f:
        if out:
            w = csv.DictWriter(f, fieldnames=list(out[0].keys()))
            w.writeheader(); w.writerows(out)
    print(f'  eQTLGen rows: {len(out)}')
    # universe file
    with open(f'{OUT}/universe_{sym}.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['chr38', 'pos38', 'ref', 'alt', 'pos19'])
        w.writeheader(); w.writerows(univ)
    # endpoint windows
    for ep in ENDPOINTS:
        rows = tabix(ep, chrom, c - WIN, c + WIN)
        with open(f'{OUT}/gwas_{ep}_{sym}.tsv', 'w') as f:
            for r in rows:
                f.write('\t'.join(r) + '\n')
        print(f'  {ep}: {len(rows)} rows')
print('done')
