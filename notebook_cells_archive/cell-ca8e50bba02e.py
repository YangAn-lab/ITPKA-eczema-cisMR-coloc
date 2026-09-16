import pandas as pd

cols_v8 = ["variant","r2","pvalue","mto_id","molecular_trait_id","maf","gene_id",
           "median_tpm","beta","se","an","ac","chromosome","position","ref","alt","type","rsid"]
gt = pd.read_csv('/workspace/coloc_panel/eqtl_gtexv8_whole_blood.tsv', sep='\t', header=None, names=cols_v8)
snps = gt[['chromosome','position','ref','alt','rsid','maf']].drop_duplicates(subset=['position','ref','alt'])
snps['rsid'] = snps.rsid.fillna('')

fg = pd.read_csv('/workspace/coloc_panel/gwas_finngen_L12_DERMATITISECZEMA.tsv', sep='\t', header=None,
                 names=['chromosome','position','ref','alt','rsids','nearest_genes','pval','mlogp','beta','sebeta','af_alt','af_alt_cases','af_alt_controls'])
fg_pos = set(fg.position)
snps['in_finngen'] = snps.position.isin(fg_pos)
shared = snps[snps.in_finngen].copy()
has_rsid = shared.rsid.str.startswith('rs')
print('GTEx v8 窗口唯一变异:', len(snps), '| 与 FinnGen 共有:', len(shared))
print('rsid 覆盖: 有效', int(has_rsid.sum()), '| 缺失:', int((~has_rsid).sum()))
shared[has_rsid].to_csv('/workspace/coloc_panel/snp_universe.csv', index=False)
print('已存 snp_universe.csv,', int(has_rsid.sum()), '行')