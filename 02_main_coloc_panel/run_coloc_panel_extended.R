#!/usr/bin/env Rscript
# P2-02 Stage A1: extended 12-gene panel coloc (eQTLGen blood x FinnGen eczema/AD)
# Panel = original 6 (GTEx v8 criterion) + 6 new (eQTLGen rs11635906 p<1e-4 / FDR-sig):
# NUSAP1, OIP5, PLA2G4B, EHD4, LINC_260814 (ENSG00000260814), JMJD7
# Full PP.H0-H4 + p12 sensitivity + window sensitivity for new genes.
suppressMessages({library(data.table); library(coloc)})
DIR <- "/workspace/coloc_panel"

panel <- data.table(
  symbol = c("ITPKA","NDUFAF1","CHP1","RPAP1","RTF1","OIP5-AS1",
             "NUSAP1","OIP5","PLA2G4B","EHD4","LINC_260814","JMJD7"),
  origin = c(rep("original", 6), rep("new", 6)))

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

results <- list(); k <- 1
for (i in seq_len(nrow(panel))) {
  sym <- panel$symbol[i]
  for (ep_name in names(gwas_list)) {
    mg <- build_merged(sym, gwas_list[[ep_name]]$dt)
    if (nrow(mg) < 50) {
      results[[k]] <- data.table(gene = sym, origin = panel$origin[i], endpoint = ep_name,
                                 nsnps = nrow(mg), note = "<50 shared SNPs"); k <- k + 1; next
    }
    pp <- run_abf(mg, gwas_list[[ep_name]])
    sens <- sapply(c(1e-6, 5e-6, 1e-5, 5e-5, 1e-4), function(p12) run_abf(mg, gwas_list[[ep_name]], p12)["PP.H4.abf"])
    # window sensitivity for new genes
    ws <- sapply(c(250000, 125000), function(w) {
      msub <- mg[abs(position - 41487062) <= w]
      if (nrow(msub) < 50) return(NA_real_)
      run_abf(msub, gwas_list[[ep_name]])["PP.H4.abf"]
    })
    results[[k]] <- data.table(
      gene = sym, origin = panel$origin[i], endpoint = ep_name, nsnps = nrow(mg),
      N1 = round(median(mg$n, na.rm = TRUE)),
      PP.H0 = pp["PP.H0.abf"], PP.H1 = pp["PP.H1.abf"], PP.H2 = pp["PP.H2.abf"],
      PP.H3 = pp["PP.H3.abf"], PP.H4 = pp["PP.H4.abf"],
      PP.H4_p12_1e6 = sens[1], PP.H4_p12_5e6 = sens[2], PP.H4_p12_1e5 = sens[3],
      PP.H4_p12_5e5 = sens[4], PP.H4_p12_1e4 = sens[5],
      PP.H4_w250 = ws[1], PP.H4_w125 = ws[2],
      min_p_eqtl = min(mg$p), min_p_gwas = min(mg$pval),
      lead_eqtl_rsid = mg[which.min(p), rsid],
      beta_lead_eqtl = mg[which.min(p), beta_e],
      note = "ok")
    k <- k + 1
  }
  cat(sym, "done\n")
}
res <- rbindlist(results, fill = TRUE)
fwrite(res, file.path(DIR, "coloc_panel_extended_results.csv"))
print(res[, .(gene, origin, endpoint, nsnps, PP.H3 = round(PP.H3, 4), PP.H4 = round(PP.H4, 4),
              p12_lo = round(PP.H4_p12_1e6, 4), p12_hi = round(PP.H4_p12_1e4, 4),
              w250 = round(PP.H4_w250, 4), w125 = round(PP.H4_w125, 4))])
