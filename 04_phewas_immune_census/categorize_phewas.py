import csv

rows = list(csv.DictReader(open('/workspace/coloc_oliva/phewas_rs11635906_full.csv')))

def cat(r):
    t = r['trait'].lower(); i = r['id']
    if i.startswith('eqtl-'): return 'cis-eQTL (eQTLGen blood)'
    if any(k in t for k in ['eczema','dermat','atopic','hayfever','rhinitis','allerg']): return 'Atopic spectrum (eczema/rhinitis/allergy)'
    if 'asthma' in t: return 'Asthma'
    if 'colitis' in t or 'crohn' in t: return 'IBD'
    if any(k in t for k in ['eosinophil','neutrophil','lymphocyte','basophil','monocyte','white blood','platelet','red blood','reticulocyte']): return 'Blood cell traits'
    if any(k in t for k in ['height','lean mass','fat percentage','body fat','weight','bmi','hip','waist']): return 'Anthropometric/body composition'
    if any(k in t for k in ['creatinine','cystatin','egfr','glomerular']): return 'Kidney function'
    if 'menopause' in t: return 'Reproductive aging'
    return 'Other'

for r in rows:
    r['category'] = cat(r)

# 每类取最显著代表
from collections import defaultdict
best = {}
for r in rows:
    c = r['category']
    if c not in best or float(r['p']) < float(best[c]['p']):
        best[c] = r

print('=== 分类汇总（每类最显著）===')
for c, r in sorted(best.items(), key=lambda kv: float(kv[1]['p'])):
    n_in = sum(1 for x in rows if x['category'] == c)
    print("{:<38} n={:<3} best: {:<45} p={:.2e} beta={:+.4f}".format(c, n_in, r['trait'][:45], float(r['p']), float(r['beta'])))

# 基因符号标注
sym = {'ENSG00000103932':'RPAP1','ENSG00000137804':'NUSAP1','ENSG00000137825':'ITPKA',
       'ENSG00000137815':'RTF1','ENSG00000104147':'OIP5'}
for r in rows:
    if r['trait'] in sym: r['trait'] = sym[r['trait']] + ' expression (eQTLGen)'

keys = ['category','id','trait','ea','nea','eaf','beta','se','p','n']
with open('/workspace/coloc_oliva/phewas_rs11635906_categorized.csv','w',newline='') as f:
    w = csv.DictWriter(f, fieldnames=keys, extrasaction='ignore'); w.writeheader()
    for r in rows: w.writerow(r)
print('\nsaved phewas_rs11635906_categorized.csv, total', len(rows))

# 关键阴性检查：自身免疫/关节
neg = ['crohn','rheumatoid','osteoarthritis','psoriasis','lupus','celiac','type 1 diabetes','vitiligo','alopecia','multiple sclerosis','ankylosing']
hits = [r['trait'] for r in rows if any(k in r['trait'].lower() for k in neg)]
print('自身免疫/关节疾病命中(p<1e-5):', hits if hits else '无')
