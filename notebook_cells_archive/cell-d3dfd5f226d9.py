# He 2020 组间 Fisher 精确检验（检出细胞数极少，如实报告功效不足）
from scipy.stats import fisher_exact
# LS vs H
or1, p1 = fisher_exact([[7, 10433-7],[8, 27864-8]])
# LS vs NL
or2, p2 = fisher_exact([[7, 10433-7],[0, 7035-0]])
print(f'LS vs Healthy: OR={or1:.2f}, p={p1:.3f}')
print(f'LS vs Non-lesional: OR={or2:.2f}, p={p2:.3f}')
print('结论: ITPKA 在 AD scRNA 中检出过于稀疏（全数据集 15/45,332 细胞），不足以做差异分析；与 bulk RNA-seq 零结果及皮肤低表达一致')