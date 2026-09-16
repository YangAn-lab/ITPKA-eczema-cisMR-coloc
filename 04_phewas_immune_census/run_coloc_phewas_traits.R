#!/usr/bin/env Rscript
# P2-02 Stage B2 (Claude M5): do the strongest non-atopic PheWAS associations at
# rs11635906 share the ITPKA cis-eQTL causal variant (horizontal pleiotropy via
# ITPKA expression) or are they independent (LD hitchhiking)?
# coloc.abf: trait GWAS x ITPKA eQTLGen blood eQTL; secondary: trait x FinnGen eczema.
suppressMessages({library(data.table); library(coloc)})
DIR <- "/workspace/coloc_panel"

traits <- data.table(
  file = c("phewas_trait_height_GCST90029008.csv", "phewas_trait_leanmass_GCST90000025.csv",
           "phewas_trait_menopause_ukb-b-17422.csv", "phewas_trait_cystatinC_GCST90014003.csv",
           "phewas_trait_egfr_GCST90026654.csv", "phewas_trait_neutrophil_GCST90002351.csv"),
  label = c("Height", "Appendicular lean mass", "Age at menopause", "Cystatin C",
            "eGFR (creatinine)", "Neutrophil count"),
  N = c(673878, 450243, 143819, 389834, 1159871, 519288))

univ <- fread(file.path(DIR, "snp_universe.csv"))
eq <- fread(file.path(DIR, "eqtlgen_ITPKA.csv")); setorder(eq, p); eq <- eq[!duplicated(rsid)]
gw <- fread(file.path(DIR, "gwas_finngen_L12_DERMATITISECZEMA.tsv"), header = FALSE)
setnames(gw, c("chromosome","position","ref","alt","rsids","nearest_genes","pval","mlogp",
               "beta","sebeta","af_alt","af_alt_cases","af_alt_controls"))

run_abf <- function(b1, v1, n1, maf, b2, v2, d2meta, p12 = 1e-5) {
  d1 <- list(beta = b1, varbeta = v1, N = n1, type = "quant", MAF = maf, snp = seq_along(b1))
  if (d2meta$type == "quant") {
    d2 <- list(beta = b2, varbeta = v2, N = d2meta$N, type = "quant", MAF = maf, snp = seq_along(b1))
  } else {
    d2 <- list(beta = b2, varbeta = v2, N = d2meta$N, type = "cc", s = d2meta$s, snp = seq_along(b1))
  }
  coloc.abf(d1, d2, p12 = p12)$summary
}

results <- list(); k <- 1
for (ti in seq_len(nrow(traits))) {
  tr <- fread(file.path(DIR, traits$file[ti]))
  tr <- tr[!is.na(beta) & !is.na(se) & se > 0]
  setorder(tr, p); tr <- tr[!duplicated(rsid)]
  # merge trait <-> eQTL by rsid, harmonize trait beta to eQTL ea
  m <- merge(eq, tr[, .(rsid, ea_t = ea, nea_t = nea, beta_t = beta, se_t = se, p_t = p)],
             by = "rsid")
  m <- merge(m, univ[, .(rsid, maf_u = maf)], by = "rsid")
  exact <- m$ea == m$ea_t & m$nea == m$nea_t
  swap  <- m$ea == m$nea_t & m$nea == m$ea_t
  m <- m[exact | swap]
  sw <- m$ea == m$nea_t & m$nea == m$ea_t
  m[, beta_t := ifelse(sw, -beta_t, beta_t)]
  m[, maf := pmin(maf_u, 1 - maf_u)]
  # (1) trait x ITPKA eQTL
  pp <- run_abf(m$beta, m$se^2, round(median(m$n, na.rm = TRUE)), m$maf,
                m$beta_t, m$se_t^2, list(type = "quant", N = traits$N[ti]))
  sens <- sapply(c(1e-6, 1e-5, 1e-4), function(p12)
    run_abf(m$beta, m$se^2, round(median(m$n, na.rm = TRUE)), m$maf,
            m$beta_t, m$se_t^2, list(type = "quant", N = traits$N[ti]), p12)["PP.H4.abf"])
  # (2) trait x FinnGen eczema GWAS (harmonize trait to GWAS alt via rsid merge)
  g2 <- merge(gw, tr[, .(rsid, ea_t = ea, nea_t = nea, beta_t = beta, se_t = se, p_t = p)],
              by.x = "rsids", by.y = "rsid")
  g2 <- merge(g2, univ[, .(rsid, maf_u = maf)], by.x = "rsids", by.y = "rsid")
  ex2 <- g2$alt == g2$ea_t & g2$ref == g2$nea_t
  sw2 <- g2$alt == g2$nea_t & g2$ref == g2$ea_t
  g2 <- g2[ex2 | sw2]
  s2m <- g2$alt == g2$nea_t & g2$ref == g2$ea_t
  g2[, beta_t := ifelse(s2m, -beta_t, beta_t)]
  g2[, maf := pmin(maf_u, 1 - maf_u)]
  pp2 <- run_abf(g2$beta, g2$sebeta^2, 0, g2$maf, g2$beta_t, g2$se_t^2,
                 list(type = "quant", N = traits$N[ti]))
  # note: run_abf arg order is (trait1=d1 ...); here d1=eczema GWAS cc
  pp2 <- coloc.abf(
    list(beta = g2$beta, varbeta = g2$sebeta^2, N = 500348, type = "cc", s = 67474/500348, snp = seq_len(nrow(g2))),
    list(beta = g2$beta_t, varbeta = g2$se_t^2, N = traits$N[ti], type = "quant", MAF = g2$maf, snp = seq_len(nrow(g2))))$summary
  r116 <- m[rsid == "rs11635906"]
  results[[k]] <- data.table(
    trait = traits$label[ti], nsnps_eqtl_merge = nrow(m),
    lead_trait_p = min(m$p_t), rs11635906_trait_p = if (nrow(r116)) r116$p_t else NA_real_,
    rs11635906_trait_beta = if (nrow(r116)) r116$beta_t else NA_real_,
    eqtl_PP.H3 = pp["PP.H3.abf"], eqtl_PP.H4 = pp["PP.H4.abf"],
    eqtl_H4_p12_1e6 = sens[1], eqtl_H4_p12_1e4 = sens[3],
    gwas_nsnps = nrow(g2), gwas_PP.H3 = pp2["PP.H3.abf"], gwas_PP.H4 = pp2["PP.H4.abf"])
  k <- k + 1
  cat(traits$label[ti], "done\n")
}
res <- rbindlist(results)
fwrite(res, file.path(DIR, "coloc_phewas_traits_results.csv"))
print(res[, .(trait, nsnps_eqtl_merge, lead_trait_p = signif(lead_trait_p, 2),
              eqtl_PP.H3 = round(eqtl_PP.H3, 3), eqtl_PP.H4 = round(eqtl_PP.H4, 3),
              H4_p12_lo = round(eqtl_H4_p12_1e6, 3), H4_p12_hi = round(eqtl_H4_p12_1e4, 3),
              gwas_PP.H3 = round(gwas_PP.H3, 3), gwas_PP.H4 = round(gwas_PP.H4, 3))])
