# A4: 用 raw counts 出正式细胞类型表达谱表
dfr = ad_raw.obs.copy()
dfr['expr'] = xr

ct = (dfr.groupby('cell_type')
        .agg(n_cells=('expr','size'),
             pct_expr=('expr', lambda s: 100*(s>0).mean()),
             mean_umi_expr=('expr', lambda s: s[s>0].mean() if (s>0).any() else 0.0),
             n_donors=('donor_id','nunique'))
        .sort_values('pct_expr', ascending=False))
ct = ct[ct['n_cells'] >= 500]
print(ct.head(20).round(3).to_string())

# 归并到大类（与论文叙事对应）
def broad(c):
    c = c.lower()
    if 'keratinocyte' in c: return 'Keratinocyte'
    if 't cell' in c or c=='t cell' or 'thymocyte' in c: return 'T cell'
    if 'fibroblast' in c: return 'Fibroblast'
    if 'melanocyte' in c: return 'Melanocyte'
    if 'langerhans' in c: return 'Langerhans cell'
    if 'dendritic' in c: return 'Dendritic cell'
    if 'macrophage' in c: return 'Macrophage'
    if 'mast' in c: return 'Mast cell'
    if 'endothelial' in c: return 'Endothelial cell'
    if 'chondrocyte' in c: return 'Chondrocyte'
    if 'secretory' in c or 'sweat' in c or 'sebaceous' in c: return 'Glandular epithelial'
    return 'Other'

dfr['broad'] = dfr['cell_type'].map(broad)
bc = (dfr.groupby('broad')
        .agg(n_cells=('expr','size'),
             pct_expr=('expr', lambda s: 100*(s>0).mean()),
             mean_umi_expr=('expr', lambda s: s[s>0].mean() if (s>0).any() else 0.0),
             n_donors=('donor_id','nunique'))
        .sort_values('pct_expr', ascending=False))
print('\n=== 大类汇总 ===')
print(bc.round(3).to_string())