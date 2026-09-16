#!/usr/bin/env Rscript
# P2-02 G6.6 修正版：coloc.susie（ITPKA × 湿疹/AD，正确 GWAS 数据）
suppressMessages({library(data.table); library(coloc)})
DIR <- "/workspace/coloc_panel"

D <- as.matrix(fread(file.path(DIR, "ld_eur.csv.gz"), header = FALSE))
ldsnps <- fread(file.path(DIR, "ld_eur_snps.csv"))
ldsnps[, snp := paste(chromosome, position, ref, alt, sep = ":")]
rownames(D) <- colnames(D) <- ldsnps$snp

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
  stopifnot("beta_g" %in% names(mg))
  mg[, beta_gwas := ifelse(s2, -beta_g, beta_g)]   # 修正：显式 GWAS 列
  mg[, snp := paste(chromosome, position, ref, alt, sep = ":")]
  mg[snp %in% ldsnps$snp]
}

out <- list()
for (ep_name in names(gwas_list)) {
  gw <- gwas_list[[ep_name]]
  mg <- build_merged("ITPKA", gw$dt)
  ord <- match(mg$snp, ldsnps$snp)
  Dsub <- D[ord, ord, drop = FALSE]
  n1 <- round(median(mg$n, na.rm = TRUE)); maf <- pmin(mg$eaf_alt, 1 - mg$eaf_alt)
  d1 <- list(beta = mg$beta_e, varbeta = mg$se^2, N = n1, type = "quant",
             MAF = maf, snp = mg$snp, LD = Dsub, position = mg$position)
  d2 <- list(beta = mg$beta_gwas, varbeta = mg$sebeta^2, N = gw$total, type = "cc",
             s = gw$cases / gw$total, snp = mg$snp, LD = Dsub, position = mg$position)
  res <- tryCatch(coloc.susie(d1, d2), error = function(e) e)
  if (inherits(res, "error")) { cat(ep_name, "ERROR:", conditionMessage(res), "\n"); next }
  s <- as.data.table(res$summary)
  cat(sprintf("\n=== ITPKA x %s（修正版）=== nsnp=%d, eQTL 可信集=%d, GWAS 可信集=%d\n",
              ep_name, nrow(mg), uniqueN(s$idx1), uniqueN(s$idx2)))
  print(s[, .(idx1, hit1, idx2, hit2, PP.H4 = round(PP.H4.abf, 4), PP.H3 = round(PP.H3.abf, 4))])
  out[[ep_name]] <- list(summary = s, nsnp = nrow(mg))
}
saveRDS(out, file.path(DIR, "coloc_susie_fixed_results.rds"))
cat("\nsaved coloc_susie_fixed_results.rds\n")
