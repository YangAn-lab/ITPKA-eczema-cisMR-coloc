#!/usr/bin/env Rscript
# Coloc triage for corrected-scan Bonferroni hits: IL6ST x AD, IL6ST x eczema
# Same fixed-pipeline conventions as the panel coloc (position+allele matching, coloc.abf).
suppressMessages({library(data.table); library(coloc)})
DIR <- "/workspace/scan498"

univ <- fread(file.path(DIR, "universe_IL6ST.csv"))  # chr38,pos38,ref,alt,pos19
univ <- univ[!is.na(pos19)]
univ[, key38 := paste(chr38, pos38, sep = ":")]
univ[, key19 := paste(chr38, pos19, sep = ":")]

eq <- fread(file.path(DIR, "eqtlgen_IL6ST.csv"))
eq[, key19 := paste(chr, position, sep = ":")]
setorder(eq, p); eq <- eq[!duplicated(key19)]
m <- merge(eq, univ, by = "key19", suffixes = c("", "_u"))
exact <- m$ea == m$alt & m$nea == m$ref
swap  <- m$ea == m$ref & m$nea == m$alt
m <- m[exact | swap]; sw <- swap[exact | swap]
m[, beta_e := ifelse(sw, -beta, beta)]
m[, eaf_alt := ifelse(sw, 1 - eaf, eaf)]
cat("eQTLGen ∩ universe (allele-checked):", nrow(m), "\n")

read_gwas <- function(ep) {
  dt <- fread(file.path(DIR, paste0("gwas_", ep, "_IL6ST.tsv")), header = FALSE,
              col.names = c("chromosome","position","ref","alt","rsids","nearest_genes","pval",
                            "mlogp","beta","sebeta","af_alt","af_alt_cases","af_alt_controls"))
  dt[, key38 := paste(chromosome, position, sep = ":")]
  dt
}
gwas_meta <- list(
  L12_ATOPIC = list(cases = 31245, total = 500348, label = "AD"),
  L12_DERMATITISECZEMA = list(cases = 67474, total = 500348, label = "Eczema"))

out_list <- list(); k <- 1
for (ep in names(gwas_meta)) {
  gw <- read_gwas(ep)
  mg <- merge(m, gw, by = "key38", suffixes = c("", "_g"))
  ex <- mg$ref == mg$ref_g & mg$alt == mg$alt_g
  swp <- mg$ref == mg$alt_g & mg$alt == mg$ref_g
  mg <- mg[ex | swp]; s2 <- swp[ex | swp]
  mg[, beta_gwas := ifelse(s2, -beta_g, beta_g)]
  mg[, snp := paste(chromosome, position, ref, alt, sep = ":")]
  cat(ep, "merged:", nrow(mg), "\n")
  n1 <- round(median(mg$n, na.rm = TRUE))
  maf <- pmin(mg$eaf_alt, 1 - mg$eaf_alt)
  meta <- gwas_meta[[ep]]
  d1 <- list(beta = mg$beta_e, varbeta = mg$se^2, N = n1, type = "quant", MAF = maf, snp = mg$snp)
  d2 <- list(beta = mg$beta_gwas, varbeta = mg$sebeta^2, N = meta$total, type = "cc",
             s = meta$cases / meta$total, snp = mg$snp)
  pp <- coloc.abf(d1, d2)$summary
  sens <- sapply(c(1e-6, 5e-6, 1e-5, 5e-5, 1e-4), function(p12)
    coloc.abf(d1, d2, p12 = p12)$summary["PP.H4.abf"])
  out_list[[k]] <- data.table(
    gene = "IL6ST", endpoint = meta$label, nsnps = nrow(mg), N1 = n1,
    PP.H0 = pp["PP.H0.abf"], PP.H1 = pp["PP.H1.abf"], PP.H2 = pp["PP.H2.abf"],
    PP.H3 = pp["PP.H3.abf"], PP.H4 = pp["PP.H4.abf"],
    PP.H4_p12_1e6 = sens[1], PP.H4_p12_5e6 = sens[2], PP.H4_p12_1e5 = sens[3],
    PP.H4_p12_5e5 = sens[4], PP.H4_p12_1e4 = sens[5],
    min_p_eqtl = min(mg$p), min_p_gwas = min(mg$pval),
    lead_eqtl_rsid = mg[which.min(p), rsid],
    lead_gwas_rsid = mg[which.min(pval), rsids])
  k <- k + 1
}
res <- rbindlist(out_list)
fwrite(res, file.path(DIR, "coloc_triage_IL6ST.csv"))
print(res[, .(gene, endpoint, nsnps, PP.H3 = round(PP.H3, 4), PP.H4 = round(PP.H4, 4),
              p12_lo = round(PP.H4_p12_1e6, 4), p12_hi = round(PP.H4_p12_1e4, 4),
              min_p_eqtl, min_p_gwas, lead_eqtl_rsid, lead_gwas_rsid)])
