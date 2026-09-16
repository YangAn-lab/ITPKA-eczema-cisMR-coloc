#!/usr/bin/env Rscript
# P2-02 审计补充: RPAP1/RTF1 coloc.susie @ ±125kb 窗口 (快速版)
suppressMessages({library(data.table); library(coloc)})
DIR <- "/workspace/coloc_panel"
CENTER <- 41487062L; W <- 125000

D <- as.matrix(fread(file.path(DIR, "ld_eur.csv.gz"), header = FALSE))
ldsnps <- fread(file.path(DIR, "ld_eur_snps.csv"))
ldsnps[, snp := paste(chromosome, position, ref, alt, sep = ":")]

univ <- fread(file.path(DIR, "snp_universe.csv"))
univ[, key := paste(chromosome, position, sep = ":")]

read_gwas <- function(ep) {
  dt <- fread(file.path(DIR, paste0("gwas_finngen_", ep, ".tsv")), header = FALSE)
  setnames(dt, c("chromosome","position","ref","alt","rsids","nearest_genes","pval","mlogp",
                 "beta","sebeta","af_alt","af_alt_cases","af_alt_controls"))
  dt[, key := paste(chromosome, position, sep = ":")]
  dt
}
gwas_list <- list(
  Eczema = list(dt = read_gwas("L12_DERMATITISECZEMA"), cases = 67474, total = 500348),
  AD     = list(dt = read_gwas("L12_ATOPIC"),           cases = 31245, total = 500348))

run_susie <- function(sym, ep_name) {
  gw <- gwas_list[[ep_name]]
  eq <- fread(file.path(DIR, paste0("eqtlgen_", sym, ".csv")))
  setorder(eq, p); eq <- eq[!duplicated(rsid)]
  m <- merge(eq, univ, by = "rsid")
  exact <- m$ea == m$alt & m$nea == m$ref
  swap  <- m$ea == m$ref & m$nea == m$alt
  m <- m[exact | swap]
  sw <- swap[exact | swap]
  m[, beta_e := ifelse(sw, -beta, beta)]
  m[, eaf_alt := ifelse(sw, 1 - eaf, eaf)]
  mg <- merge(m, gw$dt, by = "key", suffixes = c("", "_g"))
  ex <- mg$ref == mg$ref_g & mg$alt == mg$alt_g
  swp <- mg$ref == mg$alt_g & mg$alt == mg$ref_g
  mg <- mg[ex | swp]
  s2 <- swp[ex | swp]
  mg[, beta_g := ifelse(s2, -beta, beta)]
  mg[, snp := paste(chromosome, position, ref, alt, sep = ":")]
  mg <- mg[abs(position - CENTER) <= W]
  mg <- mg[snp %in% ldsnps$snp]
  ord <- match(mg$snp, ldsnps$snp)
  Dsub <- D[ord, ord, drop = FALSE]
  rownames(Dsub) <- colnames(Dsub) <- mg$snp
  n1 <- round(median(mg$n, na.rm = TRUE))
  maf <- pmin(mg$eaf_alt, 1 - mg$eaf_alt)
  d1 <- list(beta = mg$beta_e, varbeta = mg$se^2, N = n1, type = "quant",
             MAF = maf, snp = mg$snp, LD = Dsub, position = mg$position)
  d2 <- list(beta = mg$beta_g, varbeta = mg$sebeta^2, N = gw$total, type = "cc",
             s = gw$cases / gw$total, snp = mg$snp, LD = Dsub, position = mg$position)
  res <- tryCatch(coloc.susie(d1, d2), error = function(e) e)
  if (inherits(res, "error")) { cat(sym, ep_name, "ERROR:", conditionMessage(res), "\n"); return(NULL) }
  s <- as.data.table(res$summary)
  best <- s[which.max(PP.H4.abf)]
  cat(sym, ep_name, "nsnp=", nrow(mg), " CS_eqtl=", uniqueN(s$idx1), " CS_gwas=", uniqueN(s$idx2),
      " bestPP.H4=", round(best$PP.H4.abf, 4),
      " (eqtlCS", best$idx1, "hit", best$hit1, "x gwasCS", best$idx2, "hit", best$hit2, ")",
      " bestPP.H3=", round(best$PP.H3.abf, 4), "\n")
  list(summary = s, best_PP.H4 = best$PP.H4.abf)
}

out <- list()
for (sym in c("RPAP1", "RTF1"))
  for (ep in c("Eczema", "AD"))
    out[[paste(sym, ep)]] <- run_susie(sym, ep)
saveRDS(out, file.path(DIR, "coloc_susie_narrow_results.rds"))
cat("done\n")
