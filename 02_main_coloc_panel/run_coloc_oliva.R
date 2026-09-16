#!/usr/bin/env Rscript
# P2-02 G7 分析K：eQTLGen × Oliva 2025 EUR AD（42,963 cases / 408,472 controls, GCST90503109）coloc
# 复用 G6.6 修正版管线（显式 beta_g 防覆盖）；Oliva 列映射：ref=other_allele, alt=effect_allele
suppressMessages({library(data.table); library(coloc)})
DIR <- "/workspace/coloc_panel"
OL  <- "/workspace/coloc_oliva"

panel <- data.table(
  gene_id = c("ENSG00000137825","ENSG00000137806","ENSG00000187446","ENSG00000103932","ENSG00000137815","ENSG00000247556"),
  symbol  = c("ITPKA","NDUFAF1","CHP1","RPAP1","RTF1","OIP5-AS1"))

univ <- fread(file.path(DIR, "snp_universe.csv")); univ[, key := paste(chromosome, position, sep = ":")]

# Oliva EUR locus（tabix 切片 15:41.2-41.6Mb, hg38 harmonised）
gw <- fread(file.path(OL, "oliva_eur_locus.tsv"), header = FALSE)
setnames(gw, c("chromosome","base_pair_location","effect_allele","other_allele","beta",
               "standard_error","effect_allele_frequency","p_value","variant_id","rsid",
               "z","nstudy","n","effects","hm_coordinate_conversion","hm_code"))
gw <- gw[, .(chromosome, position = base_pair_location, ref = other_allele, alt = effect_allele,
             beta, sebeta = standard_error, pval = p_value, af_alt = effect_allele_frequency, rsid, n)]
gw[, key := paste(chromosome, position, sep = ":")]
gw_meta <- list(dt = gw, cases = 42963, total = 451435)   # s = 0.09517

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
  setorder(mg, pval); mg <- mg[!duplicated(snp)]   # Oliva 多等位重复行去重
  mg
}

run_abf <- function(mg, gw_meta, p12 = 1e-5) {
  n1 <- round(median(mg$n, na.rm = TRUE)); maf <- pmin(mg$eaf_alt, 1 - mg$eaf_alt)
  d1 <- list(beta = mg$beta_e, varbeta = mg$se^2, N = n1, type = "quant", MAF = maf, snp = mg$snp)
  d2 <- list(beta = mg$beta_gwas, varbeta = mg$sebeta^2, N = gw_meta$total, type = "cc",
             s = gw_meta$cases / gw_meta$total, snp = mg$snp)
  coloc.abf(d1, d2, p12 = p12)$summary
}

## ① 6 基因 × Oliva AD + p12 敏感性
results <- list(); k <- 1
for (sym in panel$symbol) {
  mg <- build_merged(sym, gw_meta$dt)
  pp <- run_abf(mg, gw_meta)
  sens <- sapply(c(1e-6, 5e-6, 1e-5, 5e-5, 1e-4), function(p12) run_abf(mg, gw_meta, p12)["PP.H4.abf"])
  results[[k]] <- data.table(
    dataset = "eQTLGen_blood", gene = sym, endpoint = "AD_OlivaEUR", nsnps = nrow(mg),
    N1 = round(median(mg$n, na.rm = TRUE)),
    PP.H0 = pp["PP.H0.abf"], PP.H1 = pp["PP.H1.abf"], PP.H2 = pp["PP.H2.abf"],
    PP.H3 = pp["PP.H3.abf"], PP.H4 = pp["PP.H4.abf"],
    PP.H4_p12_1e6 = sens[1], PP.H4_p12_5e6 = sens[2], PP.H4_p12_1e5 = sens[3],
    PP.H4_p12_5e5 = sens[4], PP.H4_p12_1e4 = sens[5],
    min_p_eqtl = min(mg$p), min_p_gwas = min(mg$pval),
    lead_eqtl_rsid = mg[which.min(p), rsid], lead_gwas_pos = mg[which.min(pval), position],
    note = "oliva_eur_fixed")
  k <- k + 1
  cat(sym, "done\n")
}
res <- rbindlist(results)
fwrite(res, file.path(OL, "coloc_oliva_eur_results.csv"))
print(res[, .(gene, nsnps, PP.H3 = round(PP.H3, 4), PP.H4 = round(PP.H4, 4),
              p12_lo = round(PP.H4_p12_1e6, 4), p12_hi = round(PP.H4_p12_1e4, 4))])

## ② ITPKA 窗口敏感性
idx_pos <- 41487062
mg_full <- build_merged("ITPKA", gw_meta$dt)
ws <- list(); k <- 1
for (w in c(500000, 250000, 125000)) {
  mg <- mg_full[abs(position - idx_pos) <= w]
  pp <- run_abf(mg, gw_meta)
  ws[[k]] <- data.table(endpoint = "AD_OlivaEUR", window_kb = w/1000, nsnps = nrow(mg),
                        PP.H3 = pp["PP.H3.abf"], PP.H4 = pp["PP.H4.abf"])
  k <- k + 1
}
ws <- rbindlist(ws)
fwrite(ws, file.path(OL, "window_sensitivity_oliva.csv"))
print(ws[, .(window_kb, nsnps, PP.H3 = round(PP.H3, 4), PP.H4 = round(PP.H4, 4))])

## ③ ITPKA coloc.susie（LD 子集）
D <- as.matrix(fread(file.path(DIR, "ld_eur.csv.gz"), header = FALSE))
ldsnps <- fread(file.path(DIR, "ld_eur_snps.csv"))
ldsnps[, snp := paste(chromosome, position, ref, alt, sep = ":")]
rownames(D) <- colnames(D) <- ldsnps$snp
mg <- mg_full[snp %in% ldsnps$snp]
ord <- match(mg$snp, ldsnps$snp)
Dsub <- D[ord, ord, drop = FALSE]
n1 <- round(median(mg$n, na.rm = TRUE)); maf <- pmin(mg$eaf_alt, 1 - mg$eaf_alt)
d1 <- list(beta = mg$beta_e, varbeta = mg$se^2, N = n1, type = "quant",
           MAF = maf, snp = mg$snp, LD = Dsub, position = mg$position)
d2 <- list(beta = mg$beta_gwas, varbeta = mg$sebeta^2, N = gw_meta$total, type = "cc",
           s = gw_meta$cases / gw_meta$total, snp = mg$snp, LD = Dsub, position = mg$position)
res_s <- tryCatch(coloc.susie(d1, d2), error = function(e) e)
if (inherits(res_s, "error")) {
  cat("susie ERROR:", conditionMessage(res_s), "\n")
} else {
  s <- as.data.table(res_s$summary)
  cat(sprintf("\n=== ITPKA x Oliva EUR susie === nsnp=%d, eQTL CS=%d, GWAS CS=%d\n",
              nrow(mg), uniqueN(s$idx1), uniqueN(s$idx2)))
  print(s[, .(idx1, hit1, idx2, hit2, PP.H4 = round(PP.H4.abf, 4), PP.H3 = round(PP.H3.abf, 4))])
  fwrite(s, file.path(OL, "coloc_susie_oliva.csv"))
}
cat("\nsaved coloc_oliva_eur_results.csv + window_sensitivity_oliva.csv\n")
