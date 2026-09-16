#!/usr/bin/env Rscript
# P2-02 G6.6 修正版：eQTLGen×FinnGen coloc 全面重跑（修复 merge 列名覆盖 bug）
# 修正要点：merge 后显式使用 GWAS 侧 beta 列（beta_g）做等位翻转，杜绝 eQTL beta 覆盖
# 输出：① 6基因×2终点 coloc.abf + p12 敏感性  ② ITPKA 窗口敏感性  ③ ITPKA coloc.susie
suppressMessages({library(data.table); library(coloc)})
DIR <- "/workspace/coloc_panel"

panel <- data.table(
  gene_id = c("ENSG00000137825","ENSG00000137806","ENSG00000187446","ENSG00000103932","ENSG00000137815","ENSG00000247556"),
  symbol  = c("ITPKA","NDUFAF1","CHP1","RPAP1","RTF1","OIP5-AS1"))

univ <- fread(file.path(DIR, "snp_universe.csv")); univ[, key := paste(chromosome, position, sep = ":")]
read_gwas <- function(ep) {
  dt <- fread(file.path(DIR, paste0("gwas_finngen_", ep, ".tsv")), header = FALSE)
  setnames(dt, c("chromosome","position","ref","alt","rsids","nearest_genes","pval","mlogp",
                 "beta","sebeta","af_alt","af_alt_cases","af_alt_controls"))
  dt[, key := paste(chromosome, position, sep = ":")]; dt
}
gwas_list <- list(
  Eczema = list(dt = read_gwas("L12_DERMATITISECZEMA"), cases = 67474, total = 500348),
  AD     = list(dt = read_gwas("L12_ATOPIC"),           cases = 31245, total = 500348))

# 修正版合并：显式列名，杜绝歧义
build_merged <- function(sym, gw) {
  eq <- fread(file.path(DIR, paste0("eqtlgen_", sym, ".csv")))
  setorder(eq, p); eq <- eq[!duplicated(rsid)]
  m <- merge(eq, univ, by = "rsid")
  exact <- m$ea == m$alt & m$nea == m$ref; swap <- m$ea == m$ref & m$nea == m$alt
  m <- m[exact | swap]; sw <- swap[exact | swap]
  m[, beta_e := ifelse(sw, -beta, beta)][, eaf_alt := ifelse(sw, 1 - eaf, eaf)]
  mg <- merge(m, gw, by = "key", suffixes = c("", "_g"))
  ex <- mg$ref == mg$ref_g & mg$alt == mg$alt_g; swp <- mg$ref == mg$alt_g & mg$alt == mg$ref_g
  mg <- mg[ex | swp]; s2 <- swp[ex | swp]
  # 修正：显式使用 GWAS 侧 beta（merge 后名为 beta_g）
  stopifnot("beta_g" %in% names(mg))
  mg[, beta_gwas := ifelse(s2, -beta_g, beta_g)]
  mg[, snp := paste(chromosome, position, ref, alt, sep = ":")]
  mg
}

run_abf <- function(mg, gw_meta, p12 = 1e-5) {
  n1 <- round(median(mg$n, na.rm = TRUE)); maf <- pmin(mg$eaf_alt, 1 - mg$eaf_alt)
  d1 <- list(beta = mg$beta_e, varbeta = mg$se^2, N = n1, type = "quant", MAF = maf, snp = mg$snp)
  d2 <- list(beta = mg$beta_gwas, varbeta = mg$sebeta^2, N = gw_meta$total, type = "cc",
             s = gw_meta$cases / gw_meta$total, snp = mg$snp)
  coloc.abf(d1, d2, p12 = p12)$summary
}

## ① 6 基因 × 2 终点 + p12 敏感性
results <- list(); k <- 1
for (sym in panel$symbol) {
  for (ep_name in names(gwas_list)) {
    mg <- build_merged(sym, gwas_list[[ep_name]]$dt)
    pp <- run_abf(mg, gwas_list[[ep_name]])
    sens <- sapply(c(1e-6, 5e-6, 1e-5, 5e-5, 1e-4), function(p12) run_abf(mg, gwas_list[[ep_name]], p12)["PP.H4.abf"])
    results[[k]] <- data.table(
      dataset = "eQTLGen_blood", gene = sym, endpoint = ep_name, nsnps = nrow(mg),
      N1 = round(median(mg$n, na.rm = TRUE)),
      PP.H0 = pp["PP.H0.abf"], PP.H1 = pp["PP.H1.abf"], PP.H2 = pp["PP.H2.abf"],
      PP.H3 = pp["PP.H3.abf"], PP.H4 = pp["PP.H4.abf"],
      PP.H4_p12_1e6 = sens[1], PP.H4_p12_5e6 = sens[2], PP.H4_p12_1e5 = sens[3],
      PP.H4_p12_5e5 = sens[4], PP.H4_p12_1e4 = sens[5],
      min_p_eqtl = min(mg$p), min_p_gwas = min(mg$pval),
      lead_eqtl_rsid = mg[which.min(p), rsid], lead_gwas_pos = mg[which.min(pval), position],
      note = "fixed")
    k <- k + 1
  }
  cat(sym, "done\n")
}
res <- rbindlist(results)
fwrite(res, file.path(DIR, "coloc_eqtlgen_results_fixed.csv"))
print(res[, .(gene, endpoint, nsnps, PP.H3 = round(PP.H3, 4), PP.H4 = round(PP.H4, 4),
              p12_lo = round(PP.H4_p12_1e6, 4), p12_hi = round(PP.H4_p12_1e4, 4))])

## ② ITPKA 窗口敏感性（±250/125 kb）
idx_pos <- 41487062
ws <- list(); k <- 1
for (ep_name in names(gwas_list)) {
  mg_full <- build_merged("ITPKA", gwas_list[[ep_name]]$dt)
  for (w in c(500000, 250000, 125000)) {
    mg <- mg_full[abs(position - idx_pos) <= w]
    pp <- run_abf(mg, gwas_list[[ep_name]])
    ws[[k]] <- data.table(endpoint = ep_name, window_kb = w/1000, nsnps = nrow(mg),
                          PP.H3 = pp["PP.H3.abf"], PP.H4 = pp["PP.H4.abf"])
    k <- k + 1
  }
}
ws <- rbindlist(ws)
fwrite(ws, file.path(DIR, "window_sensitivity_fixed.csv"))
print(ws[, .(endpoint, window_kb, nsnps, PP.H3 = round(PP.H3, 4), PP.H4 = round(PP.H4, 4))])
cat("\nsaved coloc_eqtlgen_results_fixed.csv + window_sensitivity_fixed.csv\n")
