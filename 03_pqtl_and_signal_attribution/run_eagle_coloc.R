#!/usr/bin/env Rscript
# P2-02 C1-3：EAGLE AD（hg19, n≈34.5-40.5k）× ITPKA eQTL coloc + 双向条件分析
suppressMessages({library(data.table); library(coloc)})
DIR <- "/workspace/fable_fix/pqtl"; CP <- "/workspace/coloc_panel"

D <- as.matrix(fread(file.path(CP, "ld_eur.csv.gz"), header = FALSE))
ldsnps <- fread(file.path(CP, "ld_eur_snps.csv"))
ldsnps[, snp := paste(chromosome, position, ref, alt, sep = ":")]
rownames(D) <- colnames(D) <- ldsnps$snp
univ <- fread(file.path(CP, "snp_universe.csv"))  # rsid -> hg38 position/ref/alt

# EAGLE 窗口（hg19）
ea <- fread(file.path(DIR, "EAGLE_window.tsv"))
ea <- ea[!is.na(beta) & !is.na(se) & se > 0]
setorder(ea, p.value); ea <- ea[!duplicated(rsID)]
ea <- ea[, .(rsid = rsID, pos19 = position, ref_e = reference_allele, alt_e = other_allele,
             eaf, beta_e = beta, se_e = se, p_e = p.value, N_e = European_N)]

# eQTLGen ITPKA（hg19）
eq <- fread(file.path(CP, "eqtlgen_ITPKA.csv"))
setorder(eq, p); eq <- eq[!duplicated(rsid)]
eq <- eq[, .(rsid, pos19_eq = position, ea, nea, eaf_eq = eaf, beta_eq = beta, se_eq = se, p_eq = p, n_eq = n)]

# 按 rsid 合并（hg19 坐标应一致；校验）
m <- merge(ea, eq, by = "rsid")
m[, pos_diff := abs(pos19 - pos19_eq)]
cat("合并", nrow(m), "变异；坐标不一致(>100bp):", sum(m$pos_diff > 100), "\n")
# 等位对齐：EAGLE beta 相对 reference_allele（C0 审计已验证：rs11635906 ref=G eaf=0.2657 β=+0.0398 与 meta 一致）
# 统一到 eQTLGen ea 方向：ref_e==ea 时方向已一致；alt_e==ea 时需翻转
exact <- m$ref_e == m$nea & m$alt_e == m$ea      # EAGLE alt == eQTL ea -> beta(ref) 需翻转
swap <- m$ref_e == m$ea & m$alt_e == m$nea       # EAGLE ref == eQTL ea -> 方向已一致
m <- m[exact | swap]; sw <- swap[exact | swap]
m[, beta_ea := ifelse(sw, beta_e, -beta_e)]      # sw: ref==ea 不翻；exact: alt==ea 翻转
m[, eaf_ea := ifelse(sw, eaf, 1 - eaf)]          # eaf 相对 ref_e
cat("等位匹配后:", nrow(m), "\n")

# 映射到 hg38（universe）以接 LD 面板
m2 <- merge(m, univ[, .(rsid, chromosome, position, ref, alt, maf)], by = "rsid")
m2[, snp := paste(chromosome, position, ref, alt, sep = ":")]
cat("映射到 LD 面板:", sum(m2$snp %in% ldsnps$snp), "/", nrow(m2), "\n")

# ---- coloc.abf ----
run_abf <- function(mg, p12 = 1e-5) {
  maf <- pmin(mg$eaf_ea, 1 - mg$eaf_ea)
  keep <- !is.na(maf) & maf > 1e-6 & maf < 0.999999
  mg <- mg[keep]; maf <- maf[keep]
  N1 <- round(median(mg$n_eq, na.rm = TRUE)); N2 <- round(median(mg$N_e, na.rm = TRUE))
  d1 <- list(beta = mg$beta_eq, varbeta = mg$se_eq^2, N = N1, type = "quant", MAF = maf, snp = mg$rsid)
  d2 <- list(beta = mg$beta_ea, varbeta = mg$se_e^2, N = N2, type = "cc", s = 0.5, snp = mg$rsid)
  # EAGLE: European_N 为总样本（case+control 混合）；s 未知，取 0.5 保守（Paternoster 2015 EUR: 5,606 cases / 20,565 controls 在 discovery；此处为全 EUR meta）
  coloc.abf(d1, d2, p12 = p12)$summary
}
pp <- run_abf(m2)
sens <- sapply(c(1e-6, 5e-6, 1e-5, 5e-5, 1e-4), function(p12) run_abf(m2, p12)["PP.H4.abf"])
cat(sprintf("\n=== ITPKA eQTL x EAGLE AD coloc.abf === nsnp=%d\n", nrow(m2)))
cat(sprintf("H0=%.4f H1=%.4f H2=%.4f H3=%.4f H4=%.4f\n", pp["PP.H0.abf"], pp["PP.H1.abf"], pp["PP.H2.abf"], pp["PP.H3.abf"], pp["PP.H4.abf"]))
cat("p12 敏感性 H4:", paste(round(sens, 5), collapse = " "), "\n")

# ---- coloc.susie ----
mg <- m2[snp %in% ldsnps$snp]
ord <- match(mg$snp, ldsnps$snp); Dsub <- D[ord, ord, drop = FALSE]
maf <- pmin(mg$eaf_ea, 1 - mg$eaf_ea); keep <- !is.na(maf) & maf > 1e-6 & maf < 0.999999
mg <- mg[keep]; maf <- maf[keep]; Dsub <- Dsub[keep, keep]
d1 <- list(beta = mg$beta_eq, varbeta = mg$se_eq^2, N = round(median(mg$n_eq)), type = "quant",
           MAF = maf, snp = mg$snp, LD = Dsub, position = mg$position)
d2 <- list(beta = mg$beta_ea, varbeta = mg$se_e^2, N = round(median(mg$N_e)), type = "cc",
           s = 0.5, snp = mg$snp, LD = Dsub, position = mg$position)
res <- tryCatch(coloc.susie(d1, d2), error = function(e) e)
if (inherits(res, "error")) cat("susie ERROR:", conditionMessage(res), "\n") else {
  s <- as.data.table(res$summary)
  if (nrow(s) > 0 && "idx1" %in% names(s)) print(s[, .(idx1, hit1, idx2, hit2, PP.H4 = round(PP.H4.abf,4), PP.H3 = round(PP.H3.abf,4))])
  else cat("susie: 无 CS 对（EAGLE 信号弱，符合预期）\n")
}

# ---- 双向条件分析（LD 近似）----
cat("\n=== EAGLE 条件分析 ===\n")
z <- mg$beta_ea / mg$se_e
for (cond_rs in c("rs11635906", "rs12440045")) {
  i0 <- which(mg$rsid == cond_rs)
  if (!length(i0)) { cat(cond_rs, "不在面板\n"); next }
  r <- Dsub[, i0]
  zc <- (z - r * z[i0]) / sqrt(pmax(1 - r^2, 1e-6)); zc[i0] <- NA
  cat(sprintf("条件化 %s (z=%.2f): 残留 |z|>3.9 (p<1e-4) = %d, 最大残留 |z| = %.2f\n",
              cond_rs, z[i0], sum(abs(zc) > 3.9, na.rm = TRUE), max(abs(zc), na.rm = TRUE)))
}
# index 变异在 EAGLE 的 z
i0 <- which(mg$rsid == "rs11635906")
cat(sprintf("\nrs11635906 in EAGLE: z=%.2f (p=%.3f), 在合并面板中排名 %d/%d\n",
            z[i0], 2*pnorm(-abs(z[i0])), rank(-abs(z))[i0], length(z)))
fwrite(m2, file.path(DIR, "eagle_merged_panel.csv"))
cat("saved eagle_merged_panel.csv\n")
