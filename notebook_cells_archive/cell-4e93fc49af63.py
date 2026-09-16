from scipy.stats import norm
import math

# FinnGen R12 exact values at rs11635906 (G allele)
gwas = {
    'eczema': dict(b=0.0356413, se=0.0064273, p=2.93461e-08),
    'AD':     dict(b=0.0484656, se=0.00919896, p=1.37467e-07),
}
# Real GTEx whole-blood eQTL values (v10: GTEx_Analysis_v10_eQTL.tar signif_pairs;
# v8: local eqtl_gtexv8_whole_blood.tsv)
eqtl = {
    'v10': dict(b=0.24635, se=0.044116, n=853),
    'v8':  dict(b=0.21639, se=0.0476722, n=670),
}

print(f"{'exposure':<22}{'outcome':<9}{'b_MR':>8}{'se1':>8}{'se_delta':>9}{'z':>7}{'p_delta':>10}{'F':>7}")
results = {}
for ename, e in eqtl.items():
    F = (e['b']/e['se'])**2
    for oname, g in gwas.items():
        b = g['b']/e['b']
        se1 = g['se']/abs(e['b'])                       # first-order
        se_d = math.sqrt(g['se']**2/e['b']**2 + g['b']**2*e['se']**2/e['b']**4)  # delta method
        z = b/se_d
        p = 2*norm.sf(abs(z))
        results[(ename,oname)] = dict(b=b, se1=se1, se_d=se_d, z=z, p=p, F=F)
        print(f"GTEx {ename} WB (n={e['n']})  {oname:<9}{b:8.4f}{se1:8.4f}{se_d:9.4f}{z:7.3f}{p:10.2e}{F:7.1f}")

pmax = max(r['p'] for r in results.values())
print(f"\nmax delta-method p across 4 GTEx rows = {pmax:.2e}")
print(f"F: v10={results[('v10','eczema')]['F']:.1f}, v8={results[('v8','eczema')]['F']:.1f}")
# cross-check first-order rows reproduce Table 1 printed values
for k,r in results.items():
    print(k, 'OR=', round(math.exp(r['b']),3), 'se1=', round(r['se1'],3))