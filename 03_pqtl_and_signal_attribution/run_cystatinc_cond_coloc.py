#!/usr/bin/env python3
# 胱抑素 C 第三信号对主共定位的影响：条件化/区域剔除敏感性
# 1) 从 1000G VCF 计算 rs13329240/rs7165675 与全部 LD 面板变异的 r（EUR 无关）
# 2) 近似条件化 ITPKA eQTL 与湿疹 GWAS → 重跑 coloc.abf（R 脚本完成）
# 3) 区域剔除敏感性：剔除 chr15:41.15-41.35Mb 后重跑
import pysam, numpy as np, pandas as pd

DIR = "/workspace/coloc_panel"
OUT = "/workspace/fable_fix/pqtl"
univ = pd.read_csv(f"{DIR}/snp_universe.csv").dropna(subset=["rsid"]).reset_index(drop=True)
ld_snps = pd.read_csv(f"{DIR}/ld_eur_snps.csv")["rsid"].tolist()

# 目标变异在 universe 中的键
targets = {}
for v in ["rs13329240", "rs7165675"]:
    r = univ[univ.rsid == v].iloc[0]
    targets[v] = (int(r.position), r.ref, r.alt)

eur = set(open(f"{DIR}/eur_unrel.txt").read().split())
vcf = pysam.VariantFile(f"{DIR}/kgp_chr15_window.vcf")
keep = [s for s in vcf.header.samples if s in eur]

# 提取目标变异剂量（定向到 universe ref→alt 尺度）
dos_t = {}
for rec in vcf:
    if len(rec.alts or []) != 1:
        continue
    ref, alt = rec.ref, rec.alts[0]
    if len(ref) != 1 or len(alt) != 1:
        continue
    for v, (pos, uref, ualt) in targets.items():
        if rec.pos != pos:
            continue
        if (ref, alt) == (uref, ualt):
            flip = False
        elif (ref, alt) == (ualt, uref):
            flip = True
        else:
            continue
        gt = np.array([rec.samples[s]["GT"] for s in keep], dtype=np.float32)
        if (gt < 0).any():
            continue
        d = gt.sum(axis=1)
        dos_t[v] = (2 - d) if flip else d
        print(f"{v}: pos={pos} flip={flip} EUR AF={d.mean()/2:.4f}")

# 面板变异剂量（复用 make_ld 逻辑，仅取 LD 面板内变异）
ld_set = set(ld_snps)
key2idx = {}
for i, r in univ.iterrows():
    if r.rsid in ld_set:
        key2idx[(int(r.position), r.ref, r.alt)] = r.rsid
vcf.close()
vcf = pysam.VariantFile(f"{DIR}/kgp_chr15_window.vcf")  # 重开（非索引流 reset 不可靠）
G = {}
for rec in vcf:
    if len(rec.alts or []) != 1:
        continue
    ref, alt = rec.ref, rec.alts[0]
    if len(ref) != 1 or len(alt) != 1:
        continue
    k1 = (rec.pos, ref, alt); k2 = (rec.pos, alt, ref)
    if k1 in key2idx:
        rsid, flip = key2idx[k1], False
    elif k2 in key2idx:
        rsid, flip = key2idx[k2], True
    else:
        continue
    gt = np.array([rec.samples[s]["GT"] for s in keep], dtype=np.float32)
    if (gt < 0).any():
        continue
    d = gt.sum(axis=1)
    G[rsid] = (2 - d) if flip else d

# 计算 r
rows = []
for v, dt in dos_t.items():
    dt_c = dt - dt.mean()
    for rsid in ld_snps:
        if rsid not in G:
            continue
        g = G[rsid]; g_c = g - g.mean()
        denom = np.sqrt((g_c**2).sum() * (dt_c**2).sum())
        if denom == 0:
            continue
        r = float((g_c * dt_c).sum() / denom)
        rows.append({"cond_variant": v, "rsid": rsid, "r": r})
out = pd.DataFrame(rows)
out.to_csv(f"{OUT}/cystatinc_cluster_ld.csv", index=False)
print(out.groupby("cond_variant")["r"].apply(lambda s: s.abs().max()))
print("saved cystatinc_cluster_ld.csv", out.shape)
