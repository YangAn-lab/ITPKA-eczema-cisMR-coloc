# 分析B: 跨队列验证表（rs11635906 G 等位对疾病风险的 GWAS 层面效应）
import pandas as pd

rows = [
    # cohort, ancestry, cases, N_total, beta, se, p, note
    ('FinnGen R12 L12_DERMATITISECZEMA', 'European (Finnish)', '67,474', '500,348', 0.0356413, 0.0064273, 2.93e-08, 'discovery endpoint (dermatitis/eczema)'),
    ('FinnGen R12 L12_ATOPIC', 'European (Finnish)', '31,245', '464,119', 0.0484656, 0.00919896, 1.37e-07, 'discovery endpoint (atopic dermatitis)'),
    ('Oliva 2025 multi-ancestry', 'Multi-ancestry', '56,146', '652,156', 0.0393, 0.00843, 3.15e-06, 'independent replication'),
    ('Oliva 2025 EUR', 'European', '42,963', '445,165', 0.0385, 0.00897, 1.79e-05, 'independent replication'),
    ('Oliva 2025 ASN', 'East Asian', '5,014', '176,149', 0.0597, 0.0647, 0.356, 'underpowered; same direction'),
    ('Oliva 2025 AFR', 'African', '7,063', '22,942', 0.0134, 0.0299, 0.653, 'underpowered; same direction'),
    ('EAGLE 2015 (European analysis)', 'European', 'subset of 21,399*', '34,565', 0.0398, 0.0211, 0.0596, 'predominantly childhood-onset eczema; same direction & magnitude'),
]
xval = pd.DataFrame(rows, columns=['cohort','ancestry','n_cases','n_total_this_snp','beta_G','se','p_value','note'])
xval['direction'] = '+'
print(xval.to_string(index=False))
xval.to_csv('/mnt/results/P2-02_跨队列验证_rs11635906.csv', index=False)
print('\nsaved: P2-02_跨队列验证_rs11635906.csv')
print('*EAGLE 总病例 21,399（多祖先）；European_N=34,565 为该 SNP 欧洲分析实际样本量')
print('7/7 队列方向一致; 4/7 名义显著; EUR 效应量聚于 +0.036~+0.048')