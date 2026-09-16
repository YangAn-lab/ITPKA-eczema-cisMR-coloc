#!/usr/bin/env python3
# P2-02 审计: 从 1000G hg38 EUR 无关个体计算窗口 LD 矩阵
# 剂量定向到 GTEx v8 alt 等位基因尺度 (与 coloc 数据 beta 方向一致)
import pysam, numpy as np, pandas as pd

DIR = "/workspace/coloc_panel"
univ = pd.read_csv(f"{DIR}/snp_universe.csv")  # chromosome position ref alt rsid maf
univ = univ.dropna(subset=["rsid"]).reset_index(drop=True)
key2idx = {}
for i, r in univ.iterrows():
    key2idx[(int(r.position), r.ref, r.alt)] = i  # GTEx 尺度

eur = set(open(f"{DIR}/eur_unrel.txt").read().split())
vcf = pysam.VariantFile(f"{DIR}/kgp_chr15_window.vcf")
samples = list(vcf.header.samples)
keep = [s for s in samples if s in eur]
print(f"VCF samples={len(samples)}, EUR unrel kept={len(keep)}")

n = len(univ)
G = np.full((n, len(keep)), np.nan, dtype=np.float32)
found = 0
for rec in vcf:  # 文件已是窗口切片, 顺序遍历
    if len(rec.alts or []) != 1:
        continue
    ref, alt = rec.ref, rec.alts[0]
    if len(ref) != 1 or len(alt) != 1:
        continue
    k1 = (rec.pos, ref, alt)
    k2 = (rec.pos, alt, ref)
    if k1 in key2idx:
        i, flip = key2idx[k1], False
    elif k2 in key2idx:
        i, flip = key2idx[k2], True   # 1000G alt == GTEx ref → 翻转剂量
    else:
        continue
    gt = np.array([rec.samples[s]["GT"] for s in keep], dtype=np.float32)  # (m,2)
    if (gt < 0).any():
        continue
    dos = gt.sum(axis=1)
    if flip:
        dos = 2.0 - dos
    G[i] = dos
    found += 1
print(f"universe SNPs={n}, matched in 1000G={found}")

ok = ~np.isnan(G).all(axis=1)
print(f"rows with genotypes: {ok.sum()}")
# 方差>0 的才能算相关
Gok = G[ok]
var = Gok.var(axis=1)
mono = var == 0
print(f"monomorphic in EUR: {mono.sum()}")
Gok = Gok[~mono]
idx = np.where(ok)[0][~mono]
D = np.corrcoef(Gok)
np.save(f"{DIR}/ld_eur.npy", D)
univ.iloc[idx][["chromosome", "position", "ref", "alt", "rsid"]].to_csv(
    f"{DIR}/ld_eur_snps.csv", index=False)
print("LD matrix:", D.shape, "saved")
