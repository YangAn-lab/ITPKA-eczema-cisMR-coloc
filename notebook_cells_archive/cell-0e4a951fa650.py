import requests, time
import pandas as pd

JWT = open('/workspace/.opengwas_jwt').read().strip()
HEAD = {'Authorization': f'Bearer {JWT}', 'Content-Type': 'application/json'}
snps = pd.read_csv('/workspace/coloc_panel/snp_universe.csv')
rsids = set(snps.rsid)
genes = {'ITPKA': 'ENSG00000137825', 'NDUFAF1': 'ENSG00000137806', 'CHP1': 'ENSG00000187446',
         'RPAP1': 'ENSG00000103932', 'RTF1': 'ENSG00000137815', 'OIP5-AS1': 'ENSG00000247556'}

for sym, gid in genes.items():
    f = f'/workspace/coloc_panel/eqtlgen_{sym}.csv'
    have = pd.read_csv(f)
    missing = [r for r in rsids - set(have.rsid)]
    if not missing:
        print(sym, '完整'); continue
    # 重试缺失 (部分缺失是 eQTLGen 本身不含该 SNP, 重试一次以区分)
    out = []
    for i in range(0, len(missing), 64):
        vs = missing[i:i+64]
        for attempt in range(3):
            r = requests.post('https://api.opengwas.io/api/associations', headers=HEAD,
                              json={'variant': vs, 'id': [f'eqtl-a-{gid}']}, timeout=120)
            if r.status_code == 200:
                out.extend(r.json()); break
            time.sleep(6)
        time.sleep(0.25)
    add = pd.DataFrame(out)
    if len(add):
        have = pd.concat([have, add], ignore_index=True).drop_duplicates(subset=['rsid'])
        have.to_csv(f, index=False)
    print(f'{sym}: 补拉 {len(missing)} 个缺失 rsid, 返回 {len(add)} 行, 合计 {len(have)} 行')