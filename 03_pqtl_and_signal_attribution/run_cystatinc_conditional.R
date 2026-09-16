#!/usr/bin/env Rscript
# C1-2：胱抑素 C 信号归属 —— 双向 LD 条件分析
# 问题：胱抑素 C 关联跟随信号 A（rs11635906/ITPKA eQTL）还是独立信号 B？
suppressMessages({library(data.table)})
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

# 胱抑素 C（OpenGWAS hg19 -> universe hg38）
s8 <- fread(file.path(DIR, "s8win_cystatinC.csv"))
s8 <- s8[!is.na(rsid) & !is.na(beta) & !is.na(se) & se > 0]
setorder(s8, p); s8 <- s8[!duplicated(rsid)]
s8 <- s8[, .(rsid, ea = toupper(ea), nea = toupper(nea), beta_cc = beta, se_cc = se, p_cc = p)]

# ITPKA eQTL
eq <- fread(file.path(CP, "eqtlgen_ITPKA.csv")); setorder(eq, p); eq <- eq[!duplicated(rsid)]
eq <- eq[, .(rsid, beta_eq = beta, se_eq = se, p_eq = p, ea_eq = ea, nea_eq = nea)]

# 合并：胱抑素 C + eQTL + 湿疹 GWAS + AD GWAS，以 universe 对齐等位（全部翻到 alt 方向）
m <- merge(s8, univ, by = "rsid")
ex <- m$ea == m$alt & m$nea == m$ref; sw <- m$ea == m$ref & m$nea == m$alt
m <- m[ex | sw]; s1 <- sw[ex | sw]
m[, beta_cc := ifelse(s1, -beta_cc, beta_cc)]
m <- m[, .(rsid, chromosome, position, ref, alt, maf, beta_cc, se_cc, p_cc)]

for (ep in c("L12_DERMATITISECZEMA", "L12_ATOPIC")) {
  gw <- read_gwas(ep)
  tag <- ifelse(ep == "L12_DERMATITISECZEMA", "eczema", "ad")
  g <- gw[, .(rsid = rsids, position, ref, alt, beta_g = beta, se_g = sebeta, p_g = pval)]
  mm <- merge(m, g, by = c("rsid", "position"), suffixes = c("", "_g"))
  ex <- mm$ref == mm$ref_g & mm$alt == mm$alt_g; swp <- mm$ref == mm$alt_g & mm$alt == mm$ref_g
  mm <- mm[ex | swp]; s2 <- swp[ex | swp]
  mm[, beta_g := ifelse(s2, -beta_g, beta_g)]
  m[, paste0("beta_", tag) := mm$beta_g[match(rsid, mm$rsid)]]
  m[, paste0("se_", tag) := mm$se_g[match(rsid, mm$rsid)]]
}
# eQTL 并入
mme <- merge(m, eq, by = "rsid")
ex <- mme$ea_eq == mme$alt & mme$nea_eq == mme$ref; swp <- mme$ea_eq == mme$ref & mme$nea_eq == mme$alt
mme <- mme[ex | swp]; s3 <- swp[ex | swp]
mme[, beta_eq := ifelse(s3, -beta_eq, beta_eq)]
m[, beta_eq := mme$beta_eq[match(rsid, mme$rsid)]]
m[, se_eq := mme$se_eq[match(rsid, mme$rsid)]]

m[, snp := paste(chromosome, position, ref, alt, sep = ":")]
m <- m[snp %in% ldsnps$snp]
ord <- match(m$snp, ldsnps$snp); Dsub <- D[ord, ord]
cat("合并后面板:", nrow(m), "变异\n")

cond_z <- function(z, D, i0) {
  r <- D[, i0]
  (z - r * z[i0]) / sqrt(pmax(1 - r^2, 1e-6))
}
report <- function(m, zcol, label, cond_rsid) {
  z <- m[[zcol[1]]] / m[[zcol[2]]]
  i0 <- which(m$rsid == cond_rsid)
  if (length(i0) != 1) { cat(label, ": 条件变异", cond_rsid, "不在面板\n"); return(invisible(NULL)) }
  zc <- cond_z(z, Dsub, i0)
  zc[i0] <- NA
  dt <- data.table(rsid = m$rsid, position = m$position, z_raw = round(z, 2), z_cond = round(zc, 2))
  dt[, absz := abs(z_cond)]
  setorder(dt, -absz)
  cat(sprintf("\n--- %s | 条件化 %s (z=%.2f) ---\n", label, cond_rsid, z[i0]))
  print(head(dt[, .(rsid, position, z_raw, z_cond)], 8))
  cat("条件化后 |z|>5.45 (p<5e-8) 残留:", sum(abs(zc) > 5.45, na.rm = TRUE),
      "| |z|>4.7 (p<2.6e-6):", sum(abs(zc) > 4.7, na.rm = TRUE),
      "| |z|>3.9 (p<1e-4):", sum(abs(zc) > 3.9, na.rm = TRUE), "\n")
  invisible(dt)
}

# 关键变异 z 值一览
key_snps <- c("rs11635906","rs13329240","rs7165675","rs55782852","rs316618","rs1942","rs12440045")
cat("\n=== 关键变异 z 值（胱抑素 C / 湿疹 / AD / eQTL）===\n")
for (rs in key_snps) {
  i <- which(m$rsid == rs)
  if (!length(i)) { cat(rs, "缺失\n"); next }
  cat(sprintf("%-12s pos=%d  cc z=%6.2f  ecz z=%6.2f  ad z=%6.2f  eqtl z=%7.2f\n",
              rs, m$position[i],
              m$beta_cc[i]/m$se_cc[i], m$beta_eczema[i]/m$se_eczema[i],
              m$beta_ad[i]/m$se_ad[i], m$beta_eq[i]/m$se_eq[i]))
}

# 归属检验
cat("\n########## A. 胱抑素 C 条件化 rs11635906（信号 A index）##########")
report(m, c("beta_cc","se_cc"), "胱抑素 C", "rs11635906")
cat("\n########## B. 胱抑素 C 条件化自身 lead rs7165675 ##########")
report(m, c("beta_cc","se_cc"), "胱抑素 C", "rs7165675")
cat("\n########## C. 湿疹 GWAS 条件化胱抑素 C lead rs7165675 ##########")
report(m, c("beta_eczema","se_eczema"), "湿疹 GWAS", "rs7165675")
cat("\n########## D. 湿疹 GWAS 条件化 rs11635906 ##########")
report(m, c("beta_eczema","se_eczema"), "湿疹 GWAS", "rs11635906")
cat("\n########## E. ITPKA eQTL 条件化胱抑素 C lead rs7165675 ##########")
report(m, c("beta_eq","se_eq"), "ITPKA eQTL", "rs7165675")

fwrite(m, file.path(DIR, "cystatinc_merged_panel.csv"))
cat("\nsaved cystatinc_merged_panel.csv\n")
