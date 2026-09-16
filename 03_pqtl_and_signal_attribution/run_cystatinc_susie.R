#!/usr/bin/env Rscript
# C1-2 核心：胱抑素 C 信号归属 —— coloc.susie（CS 级）+ 条件分析
suppressMessages({library(data.table); library(coloc)})
DIR <- "/workspace/fable_fix/pqtl"; CP <- "/workspace/coloc_panel"

D <- as.matrix(fread(file.path(CP, "ld_eur.csv.gz"), header = FALSE))
ldsnps <- fread(file.path(CP, "ld_eur_snps.csv"))
ldsnps[, snp := paste(chromosome, position, ref, alt, sep = ":")]
rownames(D) <- colnames(D) <- ldsnps$snp

univ <- fread(file.path(CP, "snp_universe.csv")); univ[, key := paste(chromosome, position, sep = ":")]
read_gwas <- function(ep) {
  dt <- fread(file.path(CP, paste0("gwas_finngen_", ep, ".tsv")), header = FALSE)
  setnames(dt, c("chromosome","position","ref","alt","rsids","nearest_genes","pval","mlogp",
                 "beta","sebeta","af_alt","af_alt_cases","af_alt_controls"))
  dt[, key := paste(chromosome, position, sep = ":")]; dt
}
gwas_list <- list(
  Eczema_FinnGen = list(dt = read_gwas("L12_DERMATITISECZEMA"), cases = 67474, total = 500348),
  AD_FinnGen     = list(dt = read_gwas("L12_ATOPIC"),           cases = 31245, total = 500348))

# ITPKA eQTL（hg38 key）
eq <- fread(file.path(CP, "eqtlgen_ITPKA.csv")); setorder(eq, p); eq <- eq[!duplicated(rsid)]
eqm <- merge(eq, univ, by = "rsid")
ex0 <- eqm$ea == eqm$alt & eqm$nea == eqm$ref; sw0 <- eqm$ea == eqm$ref & eqm$nea == eqm$alt
eqm <- eqm[ex0 | sw0]; swq <- sw0[ex0 | sw0]
eqm[, beta_q := ifelse(swq, -beta, beta)][, eaf_q := ifelse(swq, 1 - eaf, eaf)]
pos38 <- if ("position.y" %in% names(eqm)) eqm$position.y else eqm$position
eqm[, key := paste(eqm$chromosome, pos38, sep = ":")]
eqtl_dt <- eqm[, .(key, rsid, beta_q, se_q = se, eaf_q, p_q = p)]

# 胱抑素 C
s8 <- fread(file.path(DIR, "s8win_cystatinC.csv"))
s8 <- s8[!is.na(rsid) & !is.na(beta) & !is.na(se) & se > 0]
setorder(s8, p); s8 <- s8[!duplicated(rsid)]
s8 <- s8[, .(rsid, ea = toupper(ea), nea = toupper(nea), beta, se, p)]

merge_univ <- function(s8) {
  m <- merge(s8, univ, by = "rsid")
  exact <- m$ea == m$alt & m$nea == m$ref; swap <- m$ea == m$ref & m$nea == m$alt
  m <- m[exact | swap]; sw <- swap[exact | swap]
  m[, beta_t := ifelse(sw, -beta, beta)]
  m
}
m0 <- merge_univ(s8)

run_susie_pair <- function(mg, N1, d2_meta, label) {
  maf <- mg$maf
  keep <- !is.na(maf) & maf > 1e-6 & maf < 0.999999
  mg <- mg[keep]; maf <- maf[keep]
  ord <- match(mg$snp, ldsnps$snp); ok <- !is.na(ord)
  mg <- mg[ok]; maf <- maf[ok]; ord <- ord[ok]
  Dsub <- D[ord, ord, drop = FALSE]
  d1 <- list(beta = mg$beta_t, varbeta = mg$se^2, N = N1, type = "quant",
             MAF = maf, snp = mg$snp, LD = Dsub, position = mg$position)
  if (d2_meta$type == "cc") {
    d2 <- list(beta = mg$beta2, varbeta = mg$varbeta2, N = d2_meta$total, type = "cc",
               s = d2_meta$s, snp = mg$snp, LD = Dsub, position = mg$position)
  } else {
    d2 <- list(beta = mg$beta2, varbeta = mg$varbeta2, N = d2_meta$total, type = "quant",
               MAF = maf, snp = mg$snp, LD = Dsub, position = mg$position)
  }
  res <- tryCatch(coloc.susie(d1, d2), error = function(e) e)
  cat(sprintf("\n===== %s ===== nsnp=%d\n", label, nrow(mg)))
  if (inherits(res, "error")) { cat("ERROR:", conditionMessage(res), "\n"); return(NULL) }
  s <- as.data.table(res$summary)
  if (nrow(s) == 0 || !"idx1" %in% names(s)) { cat("无 CS 对（单侧无信号）\n"); return(res) }
  print(s[, .(idx1, hit1, idx2, hit2, PP.H4 = round(PP.H4.abf, 4), PP.H3 = round(PP.H3.abf, 4))])
  res
}

## ① 胱抑素 C × FinnGen 湿疹/AD（GWAS 侧）
for (ep in names(gwas_list)) {
  gw <- gwas_list[[ep]]
  mg <- merge(m0, gw$dt, by = "key", suffixes = c("", "_g"))
  ex <- mg$ref == mg$ref_g & mg$alt == mg$alt_g; swp <- mg$ref == mg$alt_g & mg$alt == mg$ref_g
  mg <- mg[ex | swp]; s2 <- swp[ex | swp]
  mg[, beta2 := ifelse(s2, -beta_g, beta_g)][, varbeta2 := sebeta^2]
  mg[, snp := paste(chromosome, position, ref, alt, sep = ":")]
  setorder(mg, pval); mg <- mg[!duplicated(snp)]
  run_susie_pair(mg, 389834, list(type = "cc", total = gw$total, s = gw$cases/gw$total),
                 paste0("cystatinC x ", ep))
}

## ② 胱抑素 C × ITPKA eQTL
mg <- merge(m0, eqtl_dt, by = "key")
mg[, beta2 := beta_q][, varbeta2 := se_q^2]
mg[, snp := paste(chromosome, position, ref, alt, sep = ":")]
setorder(mg, p_q); mg <- mg[!duplicated(snp)]
run_susie_pair(mg, 389834, list(type = "quant", total = 30744), "cystatinC x ITPKA_eQTL")

## ③ 条件分析：胱抑素 C 信号在条件化 index(rs11635906) 后是否残留（GCTA 风格近似：用 LD 计算条件 z）
cat("\n===== 条件分析（LD 近似）：胱抑素 C 条件化 rs11635906 前后 =====\n")
mg <- merge(m0, univ[, .(rsid)], by = "rsid")  # 已在 universe
mg[, snp := paste(chromosome, position, ref, alt, sep = ":")]
mg <- mg[snp %in% ldsnps$snp]
ord <- match(mg$snp, ldsnps$snp); Dsub <- D[ord, ord]
z <- mg$beta_t / mg$se
idx_pos <- 41487062
i0 <- which(mg$position == idx_pos)
if (length(i0) == 1) {
  r <- Dsub[, i0]
  z_cond <- (z - r * z[i0]) / sqrt(pmax(1 - r^2, 1e-6))
  res <- data.table(rsid = mg$rsid, position = mg$position, z_raw = round(z, 2),
                    z_cond_index = round(z_cond, 2))
  setorder(res, -abs(z_cond_index))
  cat("条件化 rs11635906 后 |z| 最大的 10 个变异：\n")
  print(head(res, 10))
  cat("\nrs11635906 自身: z_raw =", round(z[i0], 2), "\n")
  cat("条件化后 p<5e-8 残留变异数:", sum(abs(z_cond) > 5.45 & (1:length(z)) != i0), "\n")
}
