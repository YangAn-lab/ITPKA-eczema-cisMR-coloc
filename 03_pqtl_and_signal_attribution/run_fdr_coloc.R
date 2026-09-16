#!/usr/bin/env Rscript
# P2-02 C1-6：FGF2×AD、CSF2RB×湿疹、EDN1×湿疹 coloc（abf + p12 敏感性）
# eQTL: OpenGWAS eqtl-a（hg19，按 rsid 合并）；GWAS: FinnGen R12 窗口（hg38）
suppressMessages({library(data.table); library(coloc)})
DIR <- "/workspace/fable_fix/pqtl"

LOCI <- list(
  list(gene = "FGF2",   eqtl_file = "eqtl_FGF2_window.csv",   gwas_file = "fg_FGF2_atopic.tsv",
       endpoint = "AD_FinnGen", cases = 31245, total = 500348),
  list(gene = "CSF2RB", eqtl_file = "eqtl_CSF2RB_window.csv", gwas_file = "fg_CSF2RB_eczema.tsv",
       endpoint = "Eczema_FinnGen", cases = 67474, total = 500348),
  list(gene = "EDN1",   eqtl_file = "eqtl_EDN1_window.csv",   gwas_file = "fg_EDN1_eczema.tsv",
       endpoint = "Eczema_FinnGen", cases = 67474, total = 500348))

run_one <- function(L) {
  eq <- fread(file.path(DIR, L$eqtl_file))
  eq <- eq[!is.na(rsid) & !is.na(beta) & !is.na(se) & se > 0]
  setorder(eq, p); eq <- eq[!duplicated(rsid)]
  eq <- eq[, .(rsid, ea = toupper(ea), nea = toupper(nea), beta_e = beta, se_e = se, p_e = p, n_e = n)]
  gw <- fread(file.path(DIR, L$gwas_file), header = FALSE,
              col.names = c("chr","pos","ref","alt","rsid","gene","pval","mlogp","beta","sebeta","af","af_c","af_ct"))
  gw <- gw[!is.na(rsid)]
  mg <- merge(eq, gw, by = "rsid")
  ex <- mg$ea == mg$alt & mg$nea == mg$ref; sw <- mg$ea == mg$ref & mg$nea == mg$alt
  mg <- mg[ex | sw]; s1 <- sw[ex | sw]
  mg[, beta_eq := ifelse(s1, -beta_e, beta_e)]   # 翻到 GWAS alt 方向
  mg[, af_alt := af]
  setorder(mg, pval); mg <- mg[!duplicated(rsid)]
  maf <- pmin(mg$af_alt, 1 - mg$af_alt)
  keep <- !is.na(maf) & maf > 1e-6 & maf < 0.999999
  mg <- mg[keep]; maf <- maf[keep]
  N1 <- round(median(mg$n_e, na.rm = TRUE))
  run_abf <- function(p12) {
    d1 <- list(beta = mg$beta_eq, varbeta = mg$se_e^2, N = N1, type = "quant", MAF = maf, snp = mg$rsid)
    d2 <- list(beta = mg$beta, varbeta = mg$sebeta^2, N = L$total, type = "cc",
               s = L$cases / L$total, snp = mg$rsid)
    coloc.abf(d1, d2, p12 = p12)$summary
  }
  pp <- run_abf(1e-5)
  sens <- sapply(c(1e-6, 5e-6, 1e-5, 5e-5, 1e-4), function(p) run_abf(p)["PP.H4.abf"])
  lead_eq <- mg[which.min(p_e), rsid]; lead_gw <- mg[which.min(pval), rsid]
  cat(sprintf("\n=== %s eQTL x %s === nsnp=%d N1=%d\n", L$gene, L$endpoint, nrow(mg), N1))
  cat(sprintf("H0=%.4f H1=%.4f H2=%.4f H3=%.4f H4=%.4f\n",
              pp["PP.H0.abf"], pp["PP.H1.abf"], pp["PP.H2.abf"], pp["PP.H3.abf"], pp["PP.H4.abf"]))
  cat("p12 敏感性 H4:", paste(round(sens, 4), collapse = " "), "\n")
  cat(sprintf("eQTL lead: %s (p=%.2e) | GWAS lead: %s (p=%.2e) | 同一变异: %s\n",
              lead_eq, min(mg$p_e), lead_gw, min(mg$pval), lead_eq == lead_gw))
  data.table(gene = L$gene, endpoint = L$endpoint, nsnps = nrow(mg), N1 = N1,
             PP.H0 = pp["PP.H0.abf"], PP.H1 = pp["PP.H1.abf"], PP.H2 = pp["PP.H2.abf"],
             PP.H3 = pp["PP.H3.abf"], PP.H4 = pp["PP.H4.abf"],
             H4_p12_1e6 = sens[1], H4_p12_5e6 = sens[2], H4_p12_5e5 = sens[4], H4_p12_1e4 = sens[5],
             lead_eqtl = lead_eq, lead_gwas = lead_gw, same_lead = lead_eq == lead_gw,
             min_p_eqtl = min(mg$p_e), min_p_gwas = min(mg$pval))
}

res <- rbindlist(lapply(LOCI, run_one))
fwrite(res, file.path(DIR, "fdr_hit_coloc_results.csv"))
cat("\nsaved fdr_hit_coloc_results.csv\n")
