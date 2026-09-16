# 修复 S1 CSV：规范引号重写（A v4 指出的列偏移问题）
import csv

rows = [
    ["quantity", "value", "note"],
    ["ITPKA+ cells as fraction of normal skin cells (Census)", "0.66%", "86,752 T cells / 672,457 cells x 5.12% ITPKA+ among T cells"],
    ["ITPKA+ cells as fraction of lesional AD skin cells (He 2020)", "0.067%", "7/10,433 cells"],
    ["ITPKA+ cells as fraction of healthy skin cells (He 2020)", "0.029%", "8/27,864 cells"],
    ["T-cell fraction: lesional vs healthy skin (He 2020)", "2.91x", "16.2% vs 5.6%"],
    ["T-cell fraction: lesional vs non-lesional skin (He 2020)", "3.00x", "16.2% vs 5.4%"],
    ["Lesional/healthy ratio of ITPKA+ cell fractions", "2.34 (95% CI 0.85-6.44)", "7 vs 8 detected cells; CI spans 1 - subset expansion statistically unresolvable"],
    ["Genetic effect diluted to bulk skin (Census basis)", "beta_bulk ~ 0.0010", "0.150 x 0.66% -> log2FC ~ 0.0014 per allele"],
    ["Genetic effect diluted to bulk skin (He 2020 basis)", "beta_bulk ~ 0.00010", "0.150 x 0.067% - undetectable"],
    ["Bulk skin absolute abundance (HPA)", "0.1 nTPM", "even a 3x compositional change -> ~0.3 nTPM; below reliable bulk quantification limit"],
]
path = "/mnt/results/P2-02_SupplementaryTableS1_dilution_calculations.csv"
with open(path, "w", newline="", encoding="utf-8") as f:
    csv.writer(f, quoting=csv.QUOTE_MINIMAL).writerows(rows)

# 验证：pandas 必须无错解析且形状正确
import pandas as pd
df = pd.read_csv(path)
print("S1 parsed shape:", df.shape)
assert df.shape == (9, 3)
print(df.to_string(max_colwidth=70))