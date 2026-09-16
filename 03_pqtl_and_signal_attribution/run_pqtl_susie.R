#!/usr/bin/env Rscript
# C1-1 补充：coloc.susie 确认 TYRO3 pQTL × 湿疹/AD 为独立信号（CS 级）
suppressMessages({library(data.table); library(coloc)})
DIR <- "/workspace/fable_fix/pqtl"; CP <- "/workspace/coloc_panel"

D <- as.matrix(fread(file.path(CP, "ld_eur.csv.gz"), header = FALSE))
ldsnps <- fread(file.path(CP, "ld_eur_snps.csv"))
ldsnps[, snp := paste(chromosome, position, ref, alt, sep = ":")]
rownames(D) <- colnames(D) <- ldsnps$snp

univ <- fread(file.path(CP, "snp_universe.csv")); univ[, key := paste(chromosome, position, sep = ":")]
read_gwas <- function(ep) {
  dt <- fread(file.path(CP, paste0("gwas_finngen_", ep, ".tsv")), header = FALSE)
  setnames(dt, c("chromosome","position","ref","alt","rsids","nearest_genes","pval","mlogp",
                 "beta","sebeta","af_alt","af_alt_cases","af_alt_controls"))
  dt[, key := paste(chromosome, position, sep = ":")]; dt
}
gwas_list <- list(
  Eczema_FinnGen = list(dt = read_gwas("L12_DERMATITISECZEMA"), cases = 67474, total = 500348),
  AD_FinnGen     = list(dt = read_gwas("L12_ATOPIC"),           cases = 31245, total = 500348))

read_pqtl <- function(tag) {
  d <- fread(file.path(DIR, paste0(tag, "_window.tsv")))
  eaf <- d$hm_effect_allele_frequency
  if ("effect_allele_frequency" %in% names(d)) eaf <- ifelse(is.na(eaf), d$effect_allele_frequency, eaf)
  d[, .(rsid = hm_rsid, ea = hm_effect_allele, oa = hm_other_allele,
        beta = hm_beta, se = standard_error, p = p_value, eaf)]
}
pqtl_N <- list(ICE_TYRO3 = 5365, HELIC_TYRO3 = 1313)

build_merged <- function(tag, gw) {
  pq <- read_pqtl(tag); pq <- pq[!is.na(rsid)]
  setorder(pq, p); pq <- pq[!duplicated(rsid)]
  m <- merge(pq, univ, by = "rsid")
  ea_u <- toupper(m$ea); oa_u <- toupper(m$oa)
  exact <- ea_u == m$alt & oa_u == m$ref; swap <- ea_u == m$ref & oa_u == m$alt
  m <- m[exact | swap]; sw <- swap[exact | swap]
  m[, beta_p := ifelse(sw, -beta, beta)][, eaf_alt := ifelse(sw, 1 - eaf, eaf)]
  m[, eaf_alt := ifelse(is.na(eaf_alt), maf, eaf_alt)]
  mg <- merge(m, gw, by = "key", suffixes = c("", "_g"))
  ex <- mg$ref == mg$ref_g & mg$alt == mg$alt_g; swp <- mg$ref == mg$alt_g & mg$alt == mg$ref_g
  mg <- mg[ex | swp]; s2 <- swp[ex | swp]
  mg[, beta_gwas := ifelse(s2, -beta_g, beta_g)]
  mg[, snp := paste(chromosome, position, ref, alt, sep = ":")]
  setorder(mg, pval); mg <- mg[!duplicated(snp)]
  mg[snp %in% ldsnps$snp]
}

out <- list()
for (tag in names(pqtl_N)) {
  for (ep in names(gwas_list)) {
    gw <- gwas_list[[ep]]
    mg <- build_merged(tag, gw$dt)
    maf <- pmin(mg$eaf_alt, 1 - mg$eaf_alt)
    keep <- !is.na(maf) & maf > 1e-6 & maf < 0.999999 & mg$se > 0 & mg$sebeta > 0
    mg <- mg[keep]; maf <- maf[keep]
    ord <- match(mg$snp, ldsnps$snp)
    Dsub <- D[ord, ord, drop = FALSE]
    d1 <- list(beta = mg$beta_p, varbeta = mg$se^2, N = pqtl_N[[tag]], type = "quant",
               MAF = maf, snp = mg$snp, LD = Dsub, position = mg$position)
    d2 <- list(beta = mg$beta_gwas, varbeta = mg$sebeta^2, N = gw$total, type = "cc",
               s = gw$cases / gw$total, snp = mg$snp, LD = Dsub, position = mg$position)
    res <- tryCatch(coloc.susie(d1, d2), error = function(e) e)
    if (inherits(res, "error")) { cat(tag, ep, "ERROR:", conditionMessage(res), "\n"); next }
    s <- as.data.table(res$summary)
    if (nrow(s) == 0 || !"idx1" %in% names(s)) {
      cat(sprintf("\n=== %s x %s === nsnp=%d, 无可信集对（susie 未找到 CS），abf 结果为准\n", tag, ep, nrow(mg)))
      out[[paste(tag, ep, sep = " x ")]] <- list(summary = s, nsnp = nrow(mg)); next
    }
    cat(sprintf("\n=== %s x %s === nsnp=%d, pQTL CS=%d, GWAS CS=%d\n",
                tag, ep, nrow(mg), uniqueN(s$idx1), uniqueN(s$idx2)))
    print(s[, .(idx1, hit1, idx2, hit2, PP.H4 = round(PP.H4.abf, 4), PP.H3 = round(PP.H3.abf, 4))])
    out[[paste(tag, ep, sep = " x ")]] <- list(summary = s, nsnp = nrow(mg))
  }
}
saveRDS(out, file.path(DIR, "pqtl_coloc_susie_results.rds"))
cat("\nsaved pqtl_coloc_susie_results.rds\n")
