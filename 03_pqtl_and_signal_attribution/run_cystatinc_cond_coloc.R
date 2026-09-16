#!/usr/bin/env Rscript
# 胱抑素 C 第三信号条件化后的主共定位稳健性
suppressMessages({library(coloc); library(data.table)})
DIR <- "/workspace/coloc_panel"; OUT <- "/workspace/fable_fix/pqtl"

eq <- fread(file.path(DIR, "eqtlgen_ITPKA.csv"))
gw <- fread(file.path(DIR, "gwas_finngen_L12_DERMATITISECZEMA.tsv"), header=FALSE)
setnames(gw, c("chromosome","position","ref","alt","rsids","nearest_genes","pval","mlogp","beta","sebeta","af_alt","af_alt_cases","af_alt_controls"))
univ <- fread(file.path(DIR, "snp_universe.csv"))
ld  <- fread(file.path(OUT, "cystatinc_cluster_ld.csv"))

# GWAS 展开 rsids（多 rsid 分号分隔取第一个匹配 universe 的）
gw[, rsid := sapply(strsplit(rsids, ";"), function(x) x[1])]
m <- merge(eq, gw, by.x="rsid", by.y="rsid", suffixes=c("_eq","_gw"))
m <- merge(m, univ[, .(rsid, maf)], by="rsid")
m <- m[maf > 1e-6 & maf < 0.999999]
cat("merged nsnps:", nrow(m), "\n")

run_coloc <- function(b1, se1, b2, se2, maf, N1, N2, label) {
  d1 <- list(beta=b1, varbeta=se1^2, N=N1, type="quant", MAF=maf, snp=seq_along(b1))
  d2 <- list(beta=b2, varbeta=se2^2, N=N2, type="cc", MAF=maf, snp=seq_along(b2))
  r <- coloc.abf(d1, d2)
  cat(sprintf("%-42s H3=%.4f H4=%.4f\n", label, r$summary["PP.H3.abf"], r$summary["PP.H4.abf"]))
  invisible(r$summary)
}

res <- list()
res[["unconditioned"]] <- run_coloc(m$beta_eq, m$se, m$beta_gw, m$sebeta, m$maf, 30744, 500348, "unconditioned eQTL x eczema")

# 近似条件化：z_{j|c} = (z_j - r z_c)/sqrt(1-r^2)；转回 beta 尺度用 se_c = se/sqrt(1-r^2)
for (cv in c("rs13329240","rs7165675")) {
  rr <- ld[cond_variant == cv, .(rsid, r)]
  mm <- merge(m, rr, by="rsid")
  mm <- mm[abs(r) < 0.99]
  # eQTL 侧
  zc_eq <- mm[mm$rsid == cv, beta_eq/se]; if (length(zc_eq)==0) { # 条件变异不在 merged 集则用其最近代理
    prox <- mm[which.max(abs(r))][1]
    zc_eq <- prox$beta_eq/prox$se; cat(cv, "eQTL 侧用代理", prox$rsid, "r=", prox$r, "\n") }
  # GWAS 侧
  zc_gw <- mm[mm$rsid == cv, beta_gw/sebeta]; if (length(zc_gw)==0) {
    prox <- mm[which.max(abs(r))][1]
    zc_gw <- prox$beta_gw/prox$sebeta; cat(cv, "GWAS 侧用代理", prox$rsid, "r=", prox$r, "\n") }
  f <- sqrt(1 - mm$r^2)
  b1c <- (mm$beta_eq - mm$r * zc_eq * mm$se) / f^2   # beta 尺度条件化
  se1c <- mm$se / f
  b2c <- (mm$beta_gw - mm$r * zc_gw * mm$sebeta) / f^2
  se2c <- mm$sebeta / f
  res[[paste0("cond_", cv)]] <- run_coloc(b1c, se1c, b2c, se2c, mm$maf, 30744, 500348, paste0("conditioned on ", cv, " (n=", nrow(mm), ")"))
}

# 区域剔除：剔除胱抑素 C 峰区 41.15-41.35 Mb（hg38）
pos <- univ$position[match(m$rsid, univ$rsid)]
keep <- !(pos >= 41150000 & pos <= 41350000)
res[["drop_region"]] <- run_coloc(m$beta_eq[keep], m$se[keep], m$beta_gw[keep], m$sebeta[keep], m$maf[keep], 30744, 500348, paste0("region 41.15-41.35Mb dropped (n=", sum(keep), ")"))

out <- data.frame(analysis=names(res),
                  H3=sapply(res, `[`, "PP.H3.abf"), H4=sapply(res, `[`, "PP.H4.abf"))
fwrite(out, file.path(OUT, "cystatinc_cond_coloc_results.csv"))
print(out)
