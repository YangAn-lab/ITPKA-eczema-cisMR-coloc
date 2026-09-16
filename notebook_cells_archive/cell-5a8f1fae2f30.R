# 分析H 第2步：稀释效应定量估算
# 输入：Census 正常皮肤（T 细胞 86,752/672,457 = 12.9%；ITPKA+ 占 T 5.12%）
#       He 2020 实测 T 细胞占比（LS 16.24% / NL 5.42% / H 5.58%）
#       HPA 皮肤 bulk ITPKA = 0.1 nTPM；eQTLGen β=0.150
import numpy as np, pandas as pd
from scipy import stats

# 1) ITPKA+ 细胞占全部皮肤细胞比例
census_t_frac = 86752/672457
census_itpka_in_t = 0.051226
f_census = census_t_frac * census_itpka_in_t          # 正常皮肤
f_he_LS, f_he_H = 0.000671, 0.000287                   # He 2020 实测

# 2) T 细胞扩增倍数（He 2020 实测）
t_LS, t_NL, t_H = 0.1624, 0.0542, 0.0558

# 3) ITPKA+ 细胞分数比 LS vs H 的 Poisson 精确 CI（7/10433 vs 8/27864）
# 用条件二项精确检验构造率比 CI（近似：对数法）
import math
x1, n1, x2, n2 = 7, 10433, 8, 27864
rr = (x1/n1)/(x2/n2)
se_logrr = math.sqrt(1/x1 - 1/n1 + 1/x2 - 1/n2)
rr_lo, rr_hi = math.exp(math.log(rr)-1.96*se_logrr), math.exp(math.log(rr)+1.96*se_logrr)

# 4) 遗传效应稀释：bulk 每等位效应 ≈ β_eqtl × 表达细胞占比
beta_eqtl = 0.150
dil_census = beta_eqtl * f_census
dil_he = beta_eqtl * f_he_LS

# 5) 绝对丰度论证：0.1 nTPM 基线 × 组成变化上限 ~3× → 0.3 nTPM（< 可靠定量限 ~1 TPM）
rows = [
 ('ITPKA+ 细胞占正常皮肤细胞（Census）', f'{100*f_census:.2f}%', '86,752 T/672,457 细胞 × 5.12%'),
 ('ITPKA+ 细胞占 AD 皮损细胞（He 2020）', f'{100*f_he_LS:.3f}%', '7/10,433 细胞'),
 ('ITPKA+ 细胞占健康皮肤细胞（He 2020）', f'{100*f_he_H:.3f}%', '8/27,864 细胞'),
 ('T 细胞占比：皮损 vs 健康（He 2020）', f'{t_LS/t_H:.2f}×', f'{100*t_LS:.1f}% vs {100*t_H:.1f}%'),
 ('T 细胞占比：皮损 vs 非皮损（He 2020）', f'{t_LS/t_NL:.2f}×', f'{100*t_LS:.1f}% vs {100*t_NL:.1f}%'),
 ('ITPKA+ 细胞分数比 LS/H', f'{rr:.2f} (95%CI {rr_lo:.2f}–{rr_hi:.2f})', '7 vs 8 个细胞，CI 跨 1，无法确认亚群扩增'),
 ('遗传效应 bulk 稀释（Census 口径）', f'β_bulk ≈ {dil_census:.4f}', '0.150 × 0.66% → log2FC≈0.0014/等位'),
 ('遗传效应 bulk 稀释（He 2020 口径）', f'β_bulk ≈ {dil_he:.5f}', '0.150 × 0.067% → 不可检出'),
 ('bulk 绝对丰度（HPA 皮肤）', '0.1 nTPM', '≤3× 组成变化 → 0.3 nTPM，低于 bulk 可靠定量限'),
]
df = pd.DataFrame(rows, columns=['quantity','value','note'])
print(df.to_string(index=False))
df.to_csv('/workspace/coloc_panel/dilution_estimates.csv', index=False)
print('\nlog2FC 上限（组成变化）:', round(math.log2(t_LS/t_H),2), '| 遗传稀释 log2FC:', round(dil_census/math.log(2),4))
