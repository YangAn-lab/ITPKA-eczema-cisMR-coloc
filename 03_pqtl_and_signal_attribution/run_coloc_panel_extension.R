#!/usr/bin/env Rscript
# P2-02 C1-4：面板扩展 coloc —— 6 个未纳入基因 × 3 终点（复用 G6.6 修正版管线）
suppressMessages({library(data.table); library(coloc)})
DIR <- "/workspace/coloc_panel"; OL <- "/workspace/coloc_oliva"
OUT <- "/workspace/fable_fix/pqtl"

new_genes <- c("EHD4","JMJD7","LINC_260814","NUSAP1","OIP5","PLA2G4B")

univ <- fread(file.path(DIR, "snp_universe.csv")); univ[, key := paste(chromosome, position, sep = ":")]
read_gwas <- function(ep) {
  dt <- fread(file.path(DIR, paste0("gwas_finngen_", ep, ".tsv")), header = FALSE)
  setnames(dt, c("chromosome","position","ref","alt","rsids","nearest_genes","pval","mlogp",
                 "beta","sebeta","af_alt","af_alt_cases","af_alt_controls"))
  dt[, key := paste(chromosome, position, sep = ":")]; dt
}
gw_ol <- fread(file.path(OL, "oliva_eur_locus.tsv"), header = FALSE)
setnames(gw_ol, c("chromosome","base_pair_location","effect_allele","other_allele","beta",
                  "standard_error","effect_allele_frequency","p_value","variant_id","rsid",
                  "z","nstudy","n","effects","hm_coordinate_conversion","hm_code"))
gw_ol <- gw_ol[, .(chromosome, position = base_pair_location, ref = other_allele, alt = effect_allele,
                   beta, sebeta = standard_error, pval = p_value, af_alt = effect_allele_frequency, rsid)]
gw_ol[, key := paste(chromosome, position, sep = ":")]

gwas_list <- list(
  Eczema_FinnGen = list(dt = read_gwas("L12_DERMATITISECZEMA"), cases = 67474, total = 500348),
  AD_FinnGen     = list(dt = read_gwas("L12_ATOPIC"),           cases = 31245, total = 500348),
  AD_OlivaEUR    = list(dt = gw_ol,                             cases = 42963, total = 451435))

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
  stopifnot("beta_g" %in% names(mg))
  mg[, beta_gwas := ifelse(s2, -beta_g, beta_g)]
  mg[, snp := paste(chromosome, position, ref, alt, sep = ":")]
  setorder(mg, pval); mg <- mg[!duplicated(snp)]
  mg
}

run_abf <- function(mg, gw_meta, p12 = 1e-5) {
  n1 <- round(median(mg$n, na.rm = TRUE)); maf <- pmin(mg$eaf_alt, 1 - mg$eaf_alt)
  keep <- !is.na(maf) & maf > 1e-6 & maf < 0.999999 & mg$se > 0 & mg$sebeta > 0
  mg <- mg[keep]; maf <- maf[keep]
  d1 <- list(beta = mg$beta_e, varbeta = mg$se^2, N = n1, type = "quant", MAF = maf, snp = mg$snp)
  d2 <- list(beta = mg$beta_gwas, varbeta = mg$sebeta^2, N = gw_meta$total, type = "cc",
             s = gw_meta$cases / gw_meta$total, snp = mg$snp)
  coloc.abf(d1, d2, p12 = p12)$summary
}

results <- list(); k <- 1
for (sym in new_genes) {
  f <- file.path(DIR, paste0("eqtlgen_", sym, ".csv"))
  if (!file.exists(f)) { cat(sym, "暴露文件缺失\n"); next }
  for (ep in names(gwas_list)) {
    mg <- build_merged(sym, gwas_list[[ep]]$dt)
    if (nrow(mg) < 50) { cat(sym, ep, "SNP 太少:", nrow(mg), "\n"); next }
    pp <- run_abf(mg, gwas_list[[ep]])
    sens <- sapply(c(1e-6, 1e-4), function(p12) run_abf(mg, gwas_list[[ep]], p12)["PP.H4.abf"])
    results[[k]] <- data.table(
      gene = sym, endpoint = ep, nsnps = nrow(mg), N1 = round(median(mg$n, na.rm = TRUE)),
      PP.H0 = pp["PP.H0.abf"], PP.H1 = pp["PP.H1.abf"], PP.H2 = pp["PP.H2.abf"],
      PP.H3 = pp["PP.H3.abf"], PP.H4 = pp["PP.H4.abf"],
      H4_p12_1e6 = sens[1], H4_p12_1e4 = sens[2],
      min_p_eqtl = min(mg$p), min_p_gwas = min(mg$pval))
    k <- k + 1
    cat(sym, "x", ep, ": H3=", round(pp["PP.H3.abf"],4), "H4=", round(pp["PP.H4.abf"],4), "\n")
  }
}
res <- rbindlist(results)
fwrite(res, file.path(OUT, "coloc_panel_extension_results.csv"))
cat("\nsaved coloc_panel_extension_results.csv\n")
