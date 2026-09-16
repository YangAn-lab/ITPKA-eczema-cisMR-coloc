#!/usr/bin/env Rscript
# Extend fill: 6 extended-panel genes x {GTEx skin SE/NSE, BLUEPRINT neutro/mono} x FinnGen AD
suppressMessages({library(data.table); library(coloc)})
DIR <- "/workspace/coloc_panel"
newg <- data.table(symbol = c("NUSAP1","OIP5","PLA2G4B","EHD4","LINC_260814","JMJD7"),
                   ensg = c("ENSG00000137804","ENSG00000104147","ENSG00000243708",
                            "ENSG00000103966","ENSG00000260814","ENSG00000243789"))
gw <- fread(file.path(DIR, "gwas_finngen_L12_ATOPIC.tsv"), header = FALSE)
setnames(gw, c("chromosome","position","ref","alt","rsids","nearest_genes","pval","mlogp",
               "beta","sebeta","af_alt","af_alt_cases","af_alt_controls"))
gw[, key := paste(chromosome, position, sep = ":")]
AD <- list(dt = gw, cases = 31245, total = 500348)

run_abf <- function(b_e, se_e, n1, maf, b_g, se_g, gm, p12 = 1e-5) {
  d1 <- list(beta = b_e, varbeta = se_e^2, N = n1, type = "quant", MAF = maf, snp = seq_along(b_e))
  d2 <- list(beta = b_g, varbeta = se_g^2, N = gm$total, type = "cc", s = gm$cases/gm$total, snp = seq_along(b_e))
  coloc.abf(d1, d2, p12 = p12)$summary
}
merge_ecat <- function(file, ensg, gw) {
  dt <- fread(file.path(DIR, file), header = FALSE)
  if (ncol(dt) == 19) {
    setnames(dt, c("molecular_trait_id","chromosome","position","ref","alt","variant",
                   "ma_samples","maf","pvalue","beta","se","type","ac","an","r2",
                   "mt_object_id","gene_id","median_tpm","rsid"))
  } else if (ncol(dt) == 18) {
    setnames(dt, c("variant","r2","pvalue","mto_id","molecular_trait_id","maf","gene_id",
                   "median_tpm","beta","se","an","ac","chromosome","position","ref","alt",
                   "type","rsid"))
  } else stop("unexpected column count in ", file)
  dt[, `:=`(position = as.integer(position), beta = as.numeric(beta), se = as.numeric(se),
            maf = as.numeric(maf), pvalue = as.numeric(pvalue), an = as.numeric(an))]
  dt <- dt[gene_id == ensg]
  if (nrow(dt) == 0) return(NULL)
  setorder(dt, pvalue); dt <- dt[!duplicated(variant)]
  dt[, key := paste(chromosome, position, sep = ":")]
  mg <- merge(dt, gw, by = "key", suffixes = c("", "_g"))
  ex <- mg$ref == mg$ref_g & mg$alt == mg$alt_g; swp <- mg$ref == mg$alt_g & mg$alt == mg$ref_g
  mg <- mg[ex | swp]; s2 <- swp[ex | swp]
  mg[, beta_gwas := ifelse(s2, -beta_g, beta_g)]
  setorder(mg, pval); mg <- mg[!duplicated(paste(chromosome, position, ref, alt))]
  mg
}
results <- list(); k <- 1
for (i in seq_len(nrow(newg))) {
  sym <- newg$symbol[i]; ensg <- newg$ensg[i]
  for (ds in c("eqtl_blueprint_neutrophil.tsv", "eqtl_blueprint_monocyte.tsv",
               "eqtl_gtexv8_skin_se.tsv", "eqtl_gtexv8_skin_nse.tsv")) {
    mg2 <- merge_ecat(ds, ensg, AD$dt)
    dsn <- sub("eqtl_|\\.tsv", "", ds)
    if (is.null(mg2) || nrow(mg2) < 50) {
      results[[k]] <- data.table(gene = sym, dataset = dsn, endpoint = "AD",
        nsnps = ifelse(is.null(mg2), 0L, nrow(mg2)), PP.H3 = NA, PP.H4 = NA,
        note = ifelse(is.null(mg2), "gene not tested", "<50 shared SNPs")); k <- k + 1; next
    }
    pp2 <- run_abf(mg2$beta, mg2$se, round(median(mg2$an/2, na.rm=TRUE)), pmin(mg2$maf,1-mg2$maf),
                   mg2$beta_gwas, mg2$sebeta, AD)
    results[[k]] <- data.table(gene = sym, dataset = dsn, endpoint = "AD",
      nsnps = nrow(mg2), PP.H3 = pp2["PP.H3.abf"], PP.H4 = pp2["PP.H4.abf"], note = "ok"); k <- k + 1
  }
  cat(sym, "done\n")
}
res <- rbindlist(results, fill = TRUE)
fwrite(res, file.path(DIR, "coloc_panel_fill_ad_results.csv"))
print(res[, .(gene, dataset, nsnps, PP.H3 = round(PP.H3, 4), PP.H4 = round(PP.H4, 4), note)])
