# A5: T 细胞亚型细看 + 导出 CSV
pd.set_option('display.max_rows', 100)
tsub = ct[ct.index.str.lower().str.contains('t cell|thymocyte|nk|lymphocyte')]
print('=== T/NK 亚型 ===')
print(tsub.round(3).to_string())

# 合并导出：细类 + 大类
ct_out = ct.reset_index()
ct_out.columns = ['cell_type','n_cells','pct_expressing','mean_umi_in_expressing','n_donors']
bc_out = bc.reset_index()
bc_out.columns = ['broad_class','n_cells','pct_expressing','mean_umi_in_expressing','n_donors']

ct_out.to_csv('/mnt/results/P2-02_单细胞_正常皮肤ITPKA_细胞类型谱.csv', index=False)
bc_out.to_csv('/mnt/results/P2-02_单细胞_正常皮肤ITPKA_大类汇总.csv', index=False)
print('\nsaved CSVs')
print('\n关键数字：T cell 5.12%% (86,752 cells, 45 donors) | DC 4.27%% | Macrophage 2.65%% | Keratinocyte 0.66%% | Fibroblast 1.32%% | Langerhans ~0%%')