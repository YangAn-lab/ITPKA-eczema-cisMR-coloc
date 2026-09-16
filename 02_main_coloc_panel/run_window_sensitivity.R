#!/usr/bin/env Rscript
# P2-02 审计: eQTLGen coloc 窗口宽度敏感性 (±500/250/125 kb, 以 rs11635906 b38 41487062 为中心)
# 对 eQTLGen 中 PP.H4>0.75 的三个基因 (ITPKA/RPAP1/RTF1) × 两端点重跑
suppressMessages({library(data.table); library(coloc)})
DIR <- "/workspace/coloc_panel"
CENTER <- 41487062L

panel <- data.table(
  gene_id = c("ENSG00000137825","ENSG00000103932","ENSG00000137815"),
  symbol  = c("ITPKA","RPAP1","RTF1"))

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

results <- list(); k <- 1
for (i in seq_len(nrow(panel))) {
  sym <- panel$symbol[i]
  eq <- fread(file.path(DIR, paste0("eqtlgen_", sym, ".csv")))
  setorder(eq, p); eq <- eq[!duplicated(rsid)]
  m <- merge(eq, univ, by = "rsid")
  exact <- m$ea == m$alt & m$nea == m$ref
  swap  <- m$ea == m$ref & m$nea == m$alt
  m <- m[exact | swap]
  sw <- swap[exact | swap]
  m[, beta_e := ifelse(sw, -beta, beta)]
  m[, eaf_alt := ifelse(sw, 1 - eaf, eaf)]
  for (ep_name in names(gwas_list)) {
    gw <- gwas_list[[ep_name]]$dt
    mg <- merge(m, gw, by = "key", suffixes = c("", "_g"))
    ex <- mg$ref == mg$ref_g & mg$alt == mg$alt_g
    swp <- mg$ref == mg$alt_g & mg$alt == mg$ref_g
    mg <- mg[ex | swp]
    s2 <- swp[ex | swp]
    mg[, beta_g := ifelse(s2, -beta, beta)]
    for (w in c(500000, 250000, 125000)) {
      sub <- mg[abs(position - CENTER) <= w]
      if (nrow(sub) < 50) { results[[k]] <- data.table(gene=sym, endpoint=ep_name, window_kb=w/1000, nsnps=nrow(sub), note="<50"); k<-k+1; next }
      sub[, snp := paste(chromosome, position, ref, alt, sep = ":")]
      n1 <- round(median(sub$n, na.rm = TRUE))
      maf <- pmin(sub$eaf_alt, 1 - sub$eaf_alt)
      d1 <- list(beta=sub$beta_e, varbeta=sub$se^2, N=n1, type="quant", MAF=maf, snp=sub$snp)
      d2 <- list(beta=sub$beta_g, varbeta=sub$sebeta^2, N=gwas_list[[ep_name]]$total,
                 type="cc", s=gwas_list[[ep_name]]$cases/gwas_list[[ep_name]]$total, snp=sub$snp)
      pp <- coloc.abf(d1, d2)$summary
      results[[k]] <- data.table(gene=sym, endpoint=ep_name, window_kb=w/1000, nsnps=nrow(sub),
        PP.H3=pp["PP.H3.abf"], PP.H4=pp["PP.H4.abf"], note="ok")
      k <- k + 1
    }
  }
  cat(sym, "done\n")
}
res <- rbindlist(results, fill = TRUE)
fwrite(res, file.path(DIR, "coloc_eqtlgen_window_sensitivity.csv"))
print(res[, .(gene, endpoint, window_kb, nsnps, PP.H3=round(PP.H3,4), PP.H4=round(PP.H4,4))])
