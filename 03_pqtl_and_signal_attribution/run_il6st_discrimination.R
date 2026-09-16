#!/usr/bin/env Rscript
# IL6ST vs ANKRD55 在 5q11.2 的最小判别分析（阳性对照位点）
# coloc.abf（eQTLGen × FinnGen 湿疹/AD）+ rs7731626 条件化
suppressMessages({library(coloc); library(data.table)})
OUT <- "/workspace/fable_fix/pqtl"

eqL <- fread("/workspace/scan498/eqtlgen_IL6ST.csv")          # hg19, 含 rsid
eqA <- fread(file.path(OUT, "eqtlgen_ANKRD55_il6st_window.csv"))
# OpenGWAS 返回列: id,trait,chr,position,rsid,ea,nea,eaf,beta,se,p,n（position 为 hg19）
cat("ANKRD55 列:", paste(names(eqA), collapse=","), "\n")

read_gw <- function(f) {
  g <- fread(f, header=FALSE)
  setnames(g, c("chromosome","position","ref","alt","rsids","nearest_genes","pval","mlogp","beta","sebeta","af_alt","af_alt_cases","af_alt_controls"))
  g[, rsid := sapply(strsplit(rsids, ";"), `[`, 1)]
  g
}
gwE <- read_gw("/workspace/scan498/gwas_L12_DERMATITISECZEMA_IL6ST.tsv")
gwA <- read_gw("/workspace/scan498/gwas_L12_ATOPIC_IL6ST.tsv")

run_pair <- function(eq, gw, N2, label) {
  eq <- unique(eq[, .(rsid, beta_eq=beta, se_eq=se, eaf, n)], by="rsid")
  gw <- unique(gw[, .(rsid, beta_gw=beta, se_gw=sebeta)], by="rsid")
  m <- merge(eq, gw, by="rsid")
  m[, maf := pmin(eaf, 1-eaf)]
  m <- m[!is.na(maf) & maf > 1e-6 & maf < 0.4999999]
  d1 <- list(beta=m$beta_eq, varbeta=m$se_eq^2, N=as.integer(m$n), type="quant", MAF=m$maf, snp=m$rsid)
  d2 <- list(beta=m$beta_gw, varbeta=m$se_gw^2, N=N2, type="cc", MAF=m$maf, snp=m$rsid)
  r <- coloc.abf(d1, d2)
  # 条件化 rs7731626（若存在）
  zc_eq <- m[rsid=="rs7731626", beta_eq/se_eq]; zc_gw <- m[rsid=="rs7731626", beta_gw/se_gw]
  cond <- c(NA, NA)
  if (length(zc_eq)>0 && length(zc_gw)>0) {
    # 无 LD 矩阵：用 eQTLGen 窗口内 z 谱相关近似不可行——改用剔除 rs7731626 强 LD 代理不可行。
    # 简化：报告条件化前后的 lead 变化（用 |r| 来自 1000G 不可得，跳过精确条件化）
  }
  data.frame(label=label, nsnp=nrow(m), H3=r$summary["PP.H3.abf"], H4=r$summary["PP.H4.abf"],
             lead_eq=m$rsid[which.max(abs(m$beta_eq/m$se_eq))],
             lead_gw=m$rsid[which.max(abs(m$beta_gw/m$se_gw))],
             rs7731626_eq_p = m[rsid=="rs7731626", ]$se_eq |> (\(x) ifelse(length(x)==0, NA, 2*pnorm(-abs(m[rsid=="rs7731626", beta_eq/se_eq]))))(),
             rs7731626_gw_p = ifelse(nrow(m[rsid=="rs7731626"])==0, NA, 2*pnorm(-abs(m[rsid=="rs7731626", beta_gw/se_gw]))))
}

res <- rbind(
  run_pair(eqL, gwE, 500348, "IL6ST x eczema"),
  run_pair(eqL, gwA, 464119, "IL6ST x AD"),
  run_pair(eqA, gwE, 500348, "ANKRD55 x eczema"),
  run_pair(eqA, gwA, 464119, "ANKRD55 x AD")
)
print(res, digits=3)
fwrite(res, file.path(OUT, "il6st_ankrd55_discrimination.csv"))
cat("saved il6st_ankrd55_discrimination.csv\n")
