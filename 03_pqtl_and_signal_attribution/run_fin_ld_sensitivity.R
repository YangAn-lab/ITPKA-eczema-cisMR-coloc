#!/usr/bin/env Rscript
# P2-02 C1-7：1000G FIN 子集 LD 敏感性 —— ITPKA eQTL × FinnGen 湿疹 coloc.susie 重跑
suppressMessages({library(data.table); library(coloc)})
CP <- "/workspace/coloc_panel"; OUT <- "/workspace/fable_fix/pqtl"

# EUR LD（既有）与 FIN LD（新算）
De <- as.matrix(fread(file.path(CP, "ld_eur.csv.gz"), header = FALSE))
se <- fread(file.path(CP, "ld_eur_snps.csv")); se[, snp := paste(chromosome, position, ref, alt, sep = ":")]
rownames(De) <- colnames(De) <- se$snp
Df <- as.matrix(fread(file.path(CP, "ld_fin.csv.gz"), header = FALSE))
sf <- fread(file.path(CP, "ld_fin_snps.csv")); sf[, snp := paste(chromosome, position, ref, alt, sep = ":")]
rownames(Df) <- colnames(Df) <- sf$snp
cat("EUR LD:", nrow(se), "SNPs | FIN LD:", nrow(sf), "SNPs | 共有:", length(intersect(se$snp, sf$snp)), "\n")

# LD 矩阵直接比较（共有 SNP 子集）
common <- intersect(se$snp, sf$snp)
ie <- match(common, se$snp); if_ <- match(common, sf$snp)
De_c <- De[ie, ie]; Df_c <- Df[if_, if_]
offdiag <- upper.tri(De_c)
cat("EUR vs FIN LD 相关（上三角 r 值相关）:", cor(De_c[offdiag], Df_c[offdiag]), "\n")
# index 变异的 LD 谱比较
idx_snp <- se[rsid == "rs11635906", snp]
if (idx_snp %in% common) {
  i0 <- which(common == idx_snp)
  r_e <- De_c[i0, ]; r_f <- Df_c[i0, ]
  cat("rs11635906 LD 谱: cor(r_EUR, r_FIN) =", cor(r_e, r_f), "\n")
  cat("  |r|>0.5 的 SNP 数: EUR =", sum(abs(r_e) > 0.5), ", FIN =", sum(abs(r_f) > 0.5), "\n")
}

# coloc.susie 主分析重跑（FIN LD）
univ <- fread(file.path(CP, "snp_universe.csv")); univ[, key := paste(chromosome, position, sep = ":")]
gw <- fread(file.path(CP, "gwas_finngen_L12_DERMATITISECZEMA.tsv"), header = FALSE)
setnames(gw, c("chromosome","position","ref","alt","rsids","nearest_genes","pval","mlogp",
               "beta","sebeta","af_alt","af_alt_cases","af_alt_controls"))
gw[, key := paste(chromosome, position, sep = ":")]
eq <- fread(file.path(CP, "eqtlgen_ITPKA.csv")); setorder(eq, p); eq <- eq[!duplicated(rsid)]
m <- merge(eq, univ, by = "rsid")
exact <- m$ea == m$alt & m$nea == m$ref; swap <- m$ea == m$ref & m$nea == m$alt
m <- m[exact | swap]; sw <- swap[exact | swap]
m[, beta_e := ifelse(sw, -beta, beta)][, eaf_alt := ifelse(sw, 1 - eaf, eaf)]
mg <- merge(m, gw, by = "key", suffixes = c("", "_g"))
ex <- mg$ref == mg$ref_g & mg$alt == mg$alt_g; swp <- mg$ref == mg$alt_g & mg$alt == mg$ref_g
mg <- mg[ex | swp]; s2 <- swp[ex | swp]
mg[, beta_gwas := ifelse(s2, -beta_g, beta_g)]
mg[, snp := paste(chromosome, position, ref, alt, sep = ":")]
setorder(mg, pval); mg <- mg[!duplicated(snp)]

run_susie_with_ld <- function(mg, ldsnps, D, label) {
  sub <- mg[snp %in% ldsnps$snp]
  ord <- match(sub$snp, ldsnps$snp)
  Dsub <- D[ord, ord, drop = FALSE]
  maf <- pmin(sub$eaf_alt, 1 - sub$eaf_alt)
  keep <- !is.na(maf) & maf > 1e-6 & maf < 0.999999
  sub <- sub[keep]; maf <- maf[keep]; Dsub <- Dsub[keep, keep]
  d1 <- list(beta = sub$beta_e, varbeta = sub$se^2, N = round(median(sub$n)), type = "quant",
             MAF = maf, snp = sub$snp, LD = Dsub, position = sub$position)
  d2 <- list(beta = sub$beta_gwas, varbeta = sub$sebeta^2, N = 500348, type = "cc",
             s = 67474/500348, snp = sub$snp, LD = Dsub, position = sub$position)
  res <- tryCatch(coloc.susie(d1, d2), error = function(e) e)
  cat(sprintf("\n=== ITPKA x Eczema coloc.susie [%s LD] === nsnp=%d\n", label, nrow(sub)))
  if (inherits(res, "error")) { cat("ERROR:", conditionMessage(res), "\n"); return(NULL) }
  s <- as.data.table(res$summary)
  if (nrow(s) > 0 && "idx1" %in% names(s)) {
    print(s[, .(idx1, hit1, idx2, hit2, PP.H4 = round(PP.H4.abf, 4), PP.H3 = round(PP.H3.abf, 4))])
  } else cat("无 CS 对\n")
  res
}

r_eur <- run_susie_with_ld(mg, se, De, "EUR")
r_fin <- run_susie_with_ld(mg, sf, Df, "FIN")
saveRDS(list(EUR = r_eur, FIN = r_fin), file.path(OUT, "fin_ld_sensitivity.rds"))
cat("\nsaved fin_ld_sensitivity.rds\n")
