import pandas as pd
summary = pd.DataFrame({
    'item': ['chondrocyte_cells_total','n_donors','pct_of_all_skin_cells','ITPKA_pos_cells','ITPKA_pos_pct',
             'mean_UMI_all','main_dataset_id','main_dataset_title','main_dataset_cells','main_dataset_pct',
             'top_anatomical_site','top_site_cells','top_site_pct','donor_age_note','second_dataset'],
    'value': [2139, 21, 0.318, 312, 14.59, 0.17,
              'e7c5ba7e-4be0-4696-af50-a07d0dc7aac5',
              'BCC and Normal - scRNA-seq (Multi-scale spatial mapping of cell populations across anatomical sites in healthy human skin and BCC)',
              2121, 99.2, 'skin of external ear', 1395, 65.2,
              'donors predominantly 80-90 years old (86yo:1348, 85yo:437, 90yo:50, 80yo:62, 59yo:56)',
              'fd072bc3 (developmental immune atlas, 18 cells)']
})
summary.to_csv('/workspace/coloc_oliva/chondrocyte_check_summary.csv', index=False)
import shutil
shutil.copy('/workspace/coloc_oliva/chondrocyte_check_summary.csv', '/mnt/results/P2-02_软骨细胞核查_Census.csv')
shutil.copy('/workspace/coloc_oliva/census_chondrocyte_obs.csv', '/mnt/results/P2-02_软骨细胞核查_细胞明细.csv')
print(summary.to_string(index=False))