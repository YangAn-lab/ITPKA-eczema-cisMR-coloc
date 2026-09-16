import json, csv

d = json.load(open('/workspace/coloc_oliva/phewas_opengwas.json'))
print('total records p<1e-5:', len(d))
rows = sorted(d, key=lambda r: r['p'])
for r in rows[:45]:
    n = r.get('n', '?')
    print("{:.2e}  beta={:+.4f}  n={:>8}  {:<30} {}".format(r['p'], r['beta'], n, r['id'], r['trait'][:70]))
keys = ['id', 'trait', 'chr', 'position', 'rsid', 'ea', 'nea', 'eaf', 'beta', 'se', 'p', 'n']
with open('/workspace/coloc_oliva/phewas_rs11635906_full.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=keys, extrasaction='ignore')
    w.writeheader()
    for r in rows:
        w.writerow(r)
print('saved phewas_rs11635906_full.csv')
