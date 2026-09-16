#!/usr/bin/env Rscript
# P2-02 C1-2：S8 性状（胱抑素 C 等）× 湿疹/AD GWAS 与 × ITPKA eQTL 的 coloc
# abf + p12 敏感性 + susie（CS 级归属）；OpenGWAS 窗口数据（hg19，按 rsid 合并）
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

# ITPKA eQTL 作为第三"终点"（quant）
eq <- fread(file.path(CP, "eqtlgen_ITPKA.csv"))
setorder(eq, p); eq <- eq[!duplicated(rsid)]
eqm <- merge(eq, univ, by = "rsid")
ex0 <- eqm$ea == eqm$alt & eqm$nea == eqm$ref; sw0 <- eqm$ea == eqm$ref & eqm$nea == eqm$alt
eqm <- eqm[ex0 | sw0]; swq <- sw0[ex0 | sw0]
eqm[, beta_q := ifelse(swq, -beta, beta)][, eaf_q := ifelse(swq, 1 - eaf, eaf)]
# eq 自带 hg19 position -> merge 后 universe 的 hg38 坐标在 position.y
pos38 <- if ("position.y" %in% names(eqm)) eqm$position.y else eqm$position
eqm[, key := paste(eqm$chromosome, pos38, sep = ":")]
eqtl_itpka <- list(dt = eqm[, .(key, rsid, beta_q, se_q = se, eaf_q, p_q = p, n_q = n)],
                   N = 30744)

# S8 性状窗口
s8_files <- list(
  cystatinC      = list(file = "s8win_cystatinC.csv",      N = 389834),
  egfr           = list(file = "s8win_egfr.csv",           N = 1159871),
  height         = list(file = "s8win_height.csv",         N = 673878),
  leanmass       = list(file = "s8win_leanmass.csv",       N = 450243),
  menopause      = list(file = "s8win_menopause.csv",      N = 143819),
  neutrophil     = list(file = "s8win_neutrophil.csv",     N = 519288),
  lymphocyte_pct = list(file = "s8win_lymphocyte_pct.csv", N = 408112),
  eosinophil_count = list(file = "s8win_eosinophil_count.csv", N = NA),
  eosinophil_pct   = list(file = "s8win_eosinophil_pct.csv",   N = NA))

read_s8 <- function(f) {
  d <- fread(file.path(DIR, f))
  d <- d[!is.na(rsid) & !is.na(beta) & !is.na(se) & se > 0]
  setorder(d, p); d <- d[!duplicated(rsid)]
  d[, .(rsid, ea = toupper(ea), nea = toupper(nea), beta, se, p, n)]
}

# 合并 S8 性状与 universe（取 hg38 坐标/等位），再与 GWAS 或 eQTL 合并
build_merged_gwas <- function(s8, gw) {
  m <- merge(s8, univ, by = "rsid")
  exact <- m$ea == m$alt & m$nea == m$ref; swap <- m$ea == m$ref & m$nea == m$alt
  m <- m[exact | swap]; sw <- swap[exact | swap]
  m[, beta_t := ifelse(sw, -beta, beta)]
  mg <- merge(m, gw, by = "key", suffixes = c("", "_g"))
  ex <- mg$ref == mg$ref_g & mg$alt == mg$alt_g; swp <- mg$ref == mg$alt_g & mg$alt == mg$ref_g
  mg <- mg[ex | swp]; s2 <- swp[ex | swp]
  mg[, beta_gwas := ifelse(s2, -beta_g, beta_g)]
  mg[, snp := paste(chromosome, position, ref, alt, sep = ":")]
  setorder(mg, pval); mg <- mg[!duplicated(snp)]
  mg
}
build_merged_eqtl <- function(s8, eq) {
  m <- merge(s8, univ, by = "rsid")
  exact <- m$ea == m$alt & m$nea == m$ref; swap <- m$ea == m$ref & m$nea == m$alt
  m <- m[exact | swap]; sw <- swap[exact | swap]
  m[, beta_t := ifelse(sw, -beta, beta)]
  mg <- merge(m, eq$dt, by = "key")
  mg[, snp := paste(chromosome, position, ref, alt, sep = ":")]
  setorder(mg, p_q); mg <- mg[!duplicated(snp)]
  mg
}

run_abf <- function(b1, v1, N1, maf, snp, gw_meta, p12 = 1e-5) {
  d1 <- list(beta = b1, varbeta = v1, N = N1, type = "quant", MAF = maf, snp = snp)
  d2 <- list(beta = gw_meta$beta2, varbeta = gw_meta$varbeta2, N = gw_meta$total,
             type = gw_meta$type, snp = snp)
  if (gw_meta$type == "cc") d2$s <- gw_meta$s else d2$MAF <- maf
  coloc.abf(d1, d2, p12 = p12)$summary
}

results <- list(); k <- 1
for (tr in names(s8_files)) {
  f <- file.path(DIR, s8_files[[tr]]$file)
  if (!file.exists(f)) { cat(tr, "文件缺失，跳过\n"); next }
  s8 <- read_s8(s8_files[[tr]]$file)
  N1 <- s8_files[[tr]]$N
  if (is.na(N1)) N1 <- round(median(s8$n, na.rm = TRUE))
  for (ep in c(names(gwas_list), "ITPKA_eQTL")) {
    if (ep == "ITPKA_eQTL") {
      mg <- build_merged_eqtl(s8, eqtl_itpka)
      maf <- pmin(pmin(mg$eaf_q, 1 - mg$eaf_q), 0.5)
      gwm <- list(beta2 = mg$beta_q, varbeta2 = mg$se_q^2, total = eqtl_itpka$N, type = "quant", s = NA)
      keep <- !is.na(maf) & maf > 1e-6 & mg$se > 0 & mg$se_q > 0
      mg <- mg[keep]; maf <- maf[keep]
      pp <- run_abf(mg$beta_t, mg$se^2, N1, maf, mg$snp, gwm)
      sens <- sapply(c(1e-6, 1e-4), function(p12) run_abf(mg$beta_t, mg$se^2, N1, maf, mg$snp, gwm, p12)["PP.H4.abf"])
      minp2 <- min(mg$p_q)
    } else {
      gw <- gwas_list[[ep]]
      mg <- build_merged_gwas(s8, gw$dt)
      maf <- mg$maf
      gwm <- list(beta2 = mg$beta_gwas, varbeta2 = mg$sebeta^2, total = gw$total, type = "cc",
                  s = gw$cases / gw$total)
      keep <- !is.na(maf) & maf > 1e-6 & mg$se > 0 & mg$sebeta > 0
      mg <- mg[keep]; maf <- maf[keep]
      pp <- run_abf(mg$beta_t, mg$se^2, N1, maf, mg$snp, gwm)
      sens <- sapply(c(1e-6, 1e-4), function(p12) run_abf(mg$beta_t, mg$se^2, N1, maf, mg$snp, gwm, p12)["PP.H4.abf"])
      minp2 <- min(mg$pval)
    }
    results[[k]] <- data.table(
      trait = tr, endpoint = ep, nsnps = nrow(mg), N1 = N1,
      PP.H0 = pp["PP.H0.abf"], PP.H1 = pp["PP.H1.abf"], PP.H2 = pp["PP.H2.abf"],
      PP.H3 = pp["PP.H3.abf"], PP.H4 = pp["PP.H4.abf"],
      H4_p12_1e6 = sens[1], H4_p12_1e4 = sens[2],
      min_p_trait = min(mg$p), min_p_endpoint = minp2)
    k <- k + 1
    cat(tr, "x", ep, "done (", nrow(mg), "snps ) H3=", round(pp["PP.H3.abf"],3),
        "H4=", round(pp["PP.H4.abf"],3), "\n")
  }
}
res <- rbindlist(results)
fwrite(res, file.path(DIR, "s8_coloc_abf_results.csv"))
cat("\nsaved s8_coloc_abf_results.csv\n")
