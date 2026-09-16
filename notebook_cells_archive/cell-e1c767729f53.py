cd /workspace/coloc_oliva && python3 << 'PYEOF'
import csv
rows = list(csv.DictReader(open('phewas_rs11635906_full.csv')))
print('--- records 40-103 ---')
for r in rows[40:]:
    print('{:>10}  beta={:+.4f}  {:<28} {}'.format(r['p'], float(r['beta']), r['id'], r['trait'][:65]))
print()
print('--- 疾病/免疫相关关键词检索 ---')
kws = ['colitis','crohn','asthma','arthrit','osteoa','psoria','eczema','dermat','allerg','rhinit','autoimm','lupus','diabetes type 1','celiac','sclerosis','hypothyroid','thyroid','vitiligo','alopecia','IBD','inflamm']
for r in rows:
    t = r['trait'].lower()
    if any(k.lower() in t for k in kws):
        print('{:>10}  beta={:+.4f}  {:<28} {}'.format(r['p'], float(r['beta']), r['id'], r['trait']))
PYEOF