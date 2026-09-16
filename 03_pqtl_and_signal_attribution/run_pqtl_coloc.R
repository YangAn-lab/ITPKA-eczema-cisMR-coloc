#!/usr/bin/env Rscript
# P2-02 C1-1：pQTL × 湿疹/AD GWAS coloc
# 5 个 pQTL 数据集 × 3 个 GWAS 终点；abf + p12 敏感性；susie 按需
# 合并策略：pQTL rsid -> universe(hg38 pos/ref/alt) -> GWAS key；等位翻转同 G6.6 修正版
suppressMessages({library(data.table); library(coloc)})
DIR <- "/workspace/fable_fix/pqtl"
CP  <- "/workspace/coloc_panel"
OL  <- "/workspace/coloc_oliva"

univ <- fread(file.path(CP, "snp_universe.csv")); univ[, key := paste(chromosome, position, sep = ":")]

read_gwas <- function(ep) {
  dt <- fread(file.path(CP, paste0("gwas_finngen_", ep, ".tsv")), header = FALSE)
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

# ---- pQTL 窗口读取（统一为 rsid/ea/oa/beta/se/p/eaf）----
read_pqtl <- function(tag) {
  if (tag == "INTERVAL_TYRO3") {
    d <- fread(file.path(DIR, "INTERVAL_TYRO3_window_hg38.csv"))
    d[, .(rsid, ea, oa, beta, se, p, eaf = NA_real_)]
  } else {
    d <- fread(file.path(DIR, paste0(tag, "_window.tsv")))
    eaf <- if ("hm_effect_allele_frequency" %in% names(d)) d$hm_effect_allele_frequency else rep(NA_real_, nrow(d))
    if ("effect_allele_frequency" %in% names(d)) eaf <- ifelse(is.na(eaf), d$effect_allele_frequency, eaf)
    d[, .(rsid = hm_rsid, ea = hm_effect_allele, oa = hm_other_allele,
          beta = hm_beta, se = standard_error, p = p_value, eaf)]
  }
}
pqtl_meta <- list(
  INTERVAL_ITPKA = list(protein = "ITPKA", N = 3301, cohort = "INTERVAL"),
  ICE_ITPKA      = list(protein = "ITPKA", N = 5359, cohort = "Icelandic"),
  INTERVAL_TYRO3 = list(protein = "TYRO3", N = 3301, cohort = "INTERVAL"),
  ICE_TYRO3      = list(protein = "TYRO3", N = 5365, cohort = "Icelandic"),
  HELIC_TYRO3    = list(protein = "TYRO3", N = 1313, cohort = "HELIC"))

build_merged <- function(tag, gw) {
  pq <- read_pqtl(tag)
  pq <- pq[!is.na(rsid)]
  setorder(pq, p); pq <- pq[!duplicated(rsid)]
  m <- merge(pq, univ, by = "rsid")                      # universe: hg38 position/ref/alt/maf
  ea_u <- toupper(m$ea); oa_u <- toupper(m$oa)
  exact <- ea_u == m$alt & oa_u == m$ref; swap <- ea_u == m$ref & oa_u == m$alt
  m <- m[exact | swap]; sw <- swap[exact | swap]
  m[, beta_p := ifelse(sw, -beta, beta)]
  m[, eaf_alt := ifelse(sw, 1 - eaf, eaf)]
  m[, eaf_alt := ifelse(is.na(eaf_alt), maf, eaf_alt)]   # 缺失 EAF 用 1000G EUR MAF 近似
  mg <- merge(m, gw, by = "key", suffixes = c("", "_g"))
  ex <- mg$ref == mg$ref_g & mg$alt == mg$alt_g; swp <- mg$ref == mg$alt_g & mg$alt == mg$ref_g
  mg <- mg[ex | swp]; s2 <- swp[ex | swp]
  stopifnot("beta_g" %in% names(mg))
  mg[, beta_gwas := ifelse(s2, -beta_g, beta_g)]
  mg[, snp := paste(chromosome, position, ref, alt, sep = ":")]
  setorder(mg, pval); mg <- mg[!duplicated(snp)]
  mg
}

run_abf <- function(mg, N1, gw_meta, p12 = 1e-5) {
  maf <- pmin(mg$eaf_alt, 1 - mg$eaf_alt)
  keep <- !is.na(maf) & maf > 1e-6 & maf < 0.999999 &
          !is.na(mg$beta_p) & !is.na(mg$beta_gwas) & mg$se > 0 & mg$sebeta > 0
  mg <- mg[keep]; maf <- maf[keep]
  d1 <- list(beta = mg$beta_p, varbeta = mg$se^2, N = N1, type = "quant", MAF = maf, snp = mg$snp)
  d2 <- list(beta = mg$beta_gwas, varbeta = mg$sebeta^2, N = gw_meta$total, type = "cc",
             s = gw_meta$cases / gw_meta$total, snp = mg$snp)
  coloc.abf(d1, d2, p12 = p12)$summary
}

results <- list(); k <- 1
for (tag in names(pqtl_meta)) {
  for (ep in names(gwas_list)) {
    mg <- build_merged(tag, gwas_list[[ep]]$dt)
    if (nrow(mg) < 50) { cat(tag, ep, "SNP 太少:", nrow(mg), "\n"); next }
    pp <- run_abf(mg, pqtl_meta[[tag]]$N, gwas_list[[ep]])
    sens <- sapply(c(1e-6, 5e-6, 1e-5, 5e-5, 1e-4),
                   function(p12) run_abf(mg, pqtl_meta[[tag]]$N, gwas_list[[ep]], p12)["PP.H4.abf"])
    results[[k]] <- data.table(
      pqtl = tag, protein = pqtl_meta[[tag]]$protein, cohort = pqtl_meta[[tag]]$cohort,
      N_pqtl = pqtl_meta[[tag]]$N, endpoint = ep, nsnps = nrow(mg),
      PP.H0 = pp["PP.H0.abf"], PP.H1 = pp["PP.H1.abf"], PP.H2 = pp["PP.H2.abf"],
      PP.H3 = pp["PP.H3.abf"], PP.H4 = pp["PP.H4.abf"],
      H4_p12_1e6 = sens[1], H4_p12_5e6 = sens[2], H4_p12_5e5 = sens[4], H4_p12_1e4 = sens[5],
      min_p_pqtl = min(mg$p), min_p_gwas = min(mg$pval),
      lead_pqtl = mg[which.min(p), rsid], lead_gwas = mg[which.min(pval), rsid])
    k <- k + 1
    cat(tag, "x", ep, "done (", nrow(mg), "snps )\n")
  }
}
res <- rbindlist(results)
fwrite(res, file.path(DIR, "pqtl_coloc_abf_results.csv"))
print(res[, .(pqtl, endpoint, nsnps, H3 = round(PP.H3, 4), H4 = round(PP.H4, 4),
              minp_pqtl = signif(min_p_pqtl, 2), minp_gwas = signif(min_p_gwas, 2))])
cat("\nsaved pqtl_coloc_abf_results.csv\n")
