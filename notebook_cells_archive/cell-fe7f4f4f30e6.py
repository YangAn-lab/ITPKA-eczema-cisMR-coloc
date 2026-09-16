import pandas as pd, csv, io

# A v4 指出 S1 CSV 可能有解析问题（未加引号的逗号致列偏移）——逐文件验证
for path in ['/mnt/results/P2-02_SupplementaryTableS1_dilution_calculations.csv',
             '/mnt/results/P2-02_SupplementaryTableS2_PheWAS_rs11635906.csv']:
    print('='*20, path.split('/')[-1])
    raw = open(path, encoding='utf-8').read()
    # 每行字段数（朴素按逗号切分）是否一致
    counts = [len(l.split(',')) for l in raw.strip().split('\n')]
    print('naive comma-split field counts per line:', counts)
    try:
        df = pd.read_csv(path)
        print('pandas shape:', df.shape)
        print(df.head(15).to_string(max_colwidth=60))
    except Exception as e:
        print('pandas parse error:', e)
    print()