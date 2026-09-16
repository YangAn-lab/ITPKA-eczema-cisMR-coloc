#!/usr/bin/env Rscript
# P2-02 Stage A1b: new-panel genes in GTEx v8 whole blood x FinnGen eczema/AD
suppressMessages({library(data.table); library(coloc)})
DIR <- "/workspace/coloc_panel"

new_panel <- data.table(
  gene_id = c("ENSG00000137804","ENSG00000104147","ENSG00000243708","ENSG00000103966","ENSG00000260814","ENSG00000243789"),
  symbol  = c("NUSAP1","OIP5","PLA2G4B","EHD4","LINC_260814","JMJD7"))

cols_v8 <- c("variant","r2","pvalue","mto_id","molecular_trait_id","maf","gene_id",
             "median_tpm","beta","se","an","ac","chromosome","position","ref","alt","type","rsid")
eq_all <- fread(file.path(DIR, "eqtl_gtexv8_whole_blood.tsv"), header = FALSE)
setnames(eq_all, cols_v8)
eq_all[, `:=`(position = as.integer(position), beta = as.numeric(beta), se = as.numeric(se),
              maf = as.numeric(maf), pvalue = as.numeric(pvalue), an = as.numeric(an))]
n1 <- round(unique(eq_all$an)[1] / 2)

read_gwas <- function(ep) {
  dt <- fread(file.path(DIR, paste0("gwas_finngen_", ep, ".tsv")), header = FALSE)
  setnames(dt, c("chromosome","position","ref","alt","rsids","nearest_genes","pval","mlogp",
                 "beta","sebeta","af_alt","af_alt_cases","af_alt_controls"))
  dt[, key := paste(chromosome, position, sep = ":")]; dt
}
gwas_list <- list(
  Eczema = list(dt = read_gwas("L12_DERMATITISECZEMA"), cases = 67474, total = 500348),
  AD     = list(dt = read_gwas("L12_ATOPIC"),           cases = 31245, total = 500348))

harmonise <- function(eq, gw) {
  eq[, key := paste(chromosome, position, sep = ":")]
  m <- merge(eq, gw, by = "key", suffixes = c("_e", "_g"))
  exact <- m$ref_e == m$ref_g & m$alt_e == m$alt_g
  swap  <- m$ref_e == m$alt_g & m$alt_e == m$ref_g
  m <- m[exact | swap]; sw <- swap[exact | swap]
  m[, beta_g := ifelse(sw, -beta_g, beta_g)]
  m[, af_g := ifelse(sw, 1 - af_alt, af_alt)]
  setorder(m, pvalue); m <- m[!duplicated(paste(key, molecular_trait_id))]
  m
}

results <- list(); k <- 1
for (i in seq_len(nrow(new_panel))) {
  gid <- new_panel$gene_id[i]; sym <- new_panel$symbol[i]
  eq <- eq_all[gene_id == gid]
  if (nrow(eq) < 50) {
    results[[k]] <- data.table(gene = sym, endpoint = NA, nsnps = nrow(eq),
                               note = "gene absent or <50 SNPs in GTEx v8 whole blood"); k <- k + 1; next
  }
  for (ep_name in names(gwas_list)) {
    m <- harmonise(eq, gwas_list[[ep_name]]$dt)
    if (nrow(m) < 50) { results[[k]] <- data.table(gene = sym, endpoint = ep_name, nsnps = nrow(m), note = "<50 shared"); k <- k + 1; next }
    m[, snp := paste(chromosome_e, position_e, ref_e, alt_e, sep = ":")]
    d1 <- list(beta = m$beta_e, varbeta = m$se^2, N = n1, type = "quant",
               MAF = pmin(m$maf, 1 - m$maf), snp = m$snp)
    d2 <- list(beta = m$beta_g, varbeta = m$sebeta^2, N = gwas_list[[ep_name]]$total,
               type = "cc", s = gwas_list[[ep_name]]$cases / gwas_list[[ep_name]]$total, snp = m$snp)
    pp <- coloc.abf(d1, d2)$summary
    sens <- sapply(c(1e-6, 1e-4), function(p12) coloc.abf(d1, d2, p12 = p12)$summary["PP.H4.abf"])
    results[[k]] <- data.table(
      dataset = "GTExv8_WholeBlood", gene = sym, endpoint = ep_name, nsnps = nrow(m), N1 = n1,
      PP.H0 = pp["PP.H0.abf"], PP.H1 = pp["PP.H1.abf"], PP.H2 = pp["PP.H2.abf"],
      PP.H3 = pp["PP.H3.abf"], PP.H4 = pp["PP.H4.abf"],
      PP.H4_p12_1e6 = sens[1], PP.H4_p12_1e4 = sens[2],
      min_p_eqtl = min(m$pvalue), lead_eqtl = m[which.min(pvalue), snp],
      median_tpm = unique(m$median_tpm)[1], note = "ok")
    k <- k + 1
  }
  cat(sym, "done\n")
}
res <- rbindlist(results, fill = TRUE)
fwrite(res, file.path(DIR, "coloc_gtexv8wb_extended_results.csv"))
print(res[, .(gene, endpoint, nsnps, PP.H3 = round(PP.H3,3), PP.H4 = round(PP.H4,4),
              p12_lo = round(PP.H4_p12_1e6,4), p12_hi = round(PP.H4_p12_1e4,4),
              min_p_eqtl, median_tpm, note)])
