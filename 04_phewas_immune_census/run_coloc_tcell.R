#!/usr/bin/env Rscript
# P2-02 Stage A2: T-cell (and comparator immune cell) eQTL coloc
# Question (Claude M2): does the ITPKA whole-blood eQTL signal replicate/colocalize
# in purified immune cell datasets? Measurability documented per dataset.
# Window: chr15 40,987,062-41,987,062 (hg38, +/-500kb around rs11635906).
# eQTL Catalogue files are GRCh38 -> direct merge with FinnGen R12 on chrom:pos.
suppressMessages({library(data.table); library(coloc)})
DIR <- "/workspace/coloc_panel"

meta <- data.table(
  file = c("eqtl_tcell_BLUEPRINT_CD4T.tsv",
           "eqtl_tcell_OneK1K_CD4_CTL.tsv","eqtl_tcell_OneK1K_CD4_Naive.tsv",
           "eqtl_tcell_OneK1K_CD4_TCM.tsv","eqtl_tcell_OneK1K_CD4_TEM.tsv",
           "eqtl_tcell_OneK1K_CD8_Naive.tsv","eqtl_tcell_OneK1K_CD8_TCM.tsv",
           "eqtl_tcell_OneK1K_CD8_TEM.tsv","eqtl_tcell_OneK1K_MAIT.tsv",
           "eqtl_tcell_OneK1K_NK.tsv","eqtl_tcell_OneK1K_Treg.tsv","eqtl_tcell_OneK1K_gdT.tsv",
           "eqtl_tcell_DICE_Tfh_memory.tsv","eqtl_tcell_DICE_Th17_memory.tsv",
           "eqtl_tcell_DICE_Th1_memory.tsv","eqtl_tcell_DICE_Th2_memory.tsv",
           "eqtl_tcell_DICE_Th1-17_memory.tsv","eqtl_tcell_DICE_Treg_memory.tsv",
           "eqtl_tcell_DICE_Treg_naive.tsv","eqtl_tcell_DICE_B-cell_naive.tsv",
           "eqtl_tcell_DICE_CD4_T-cell_naive.tsv","eqtl_tcell_DICE_CD4_T-cell_anti-CD3-CD28.tsv",
           "eqtl_tcell_DICE_CD8_T-cell_naive.tsv","eqtl_tcell_DICE_CD8_T-cell_anti-CD3-CD28.tsv",
           "eqtl_tcell_DICE_monocyte_CD16_naive.tsv","eqtl_tcell_DICE_monocyte_naive.tsv",
           "eqtl_tcell_DICE_NK-cell_naive.tsv",
           "eqtl_tcell_GENCORD_T.tsv","eqtl_tcell_Bossini_Treg.tsv",
           "eqtl_tcell_CEDAR_CD4.tsv","eqtl_tcell_CEDAR_CD8.tsv",
           "eqtl_tcell_Kasela_CD4.tsv","eqtl_tcell_Kasela_CD8.tsv"),
  study = c("BLUEPRINT", rep("OneK1K", 11), rep("DICE", 15), "GENCORD", "BossiniCastillo",
            "CEDAR", "CEDAR", "Kasela2017", "Kasela2017"),
  cell = c("CD4+ T",
           "CD4 CTL","CD4 naive","CD4 TCM","CD4 TEM","CD8 naive","CD8 TCM","CD8 TEM",
           "MAIT","NK","Treg","gdT",
           "Tfh memory","Th17 memory","Th1 memory","Th2 memory","Th1-17 memory",
           "Treg memory","Treg naive","B naive","CD4 naive","CD4 anti-CD3-CD28",
           "CD8 naive","CD8 anti-CD3-CD28","CD16 monocyte","monocyte","NK",
           "T (cord blood)","Treg","CD4+ T","CD8+ T","CD4+ T","CD8+ T"),
  N = c(167, 371,949,954,782,729,457,946,186,951,745,453,
        89,89,82,89,88,89,89,91,88,89,89,88,90,91,90,
        184,119, 290,277,280,269),
  platform = c(rep("RNAseq", 29), rep("microarray", 4)))

genes <- data.table(symbol = c("ITPKA","NUSAP1","RTF1"),
                    ensg   = c("ENSG00000137825","ENSG00000137804","ENSG00000137815"))

read_gwas <- function(ep) {
  dt <- fread(file.path(DIR, paste0("gwas_finngen_", ep, ".tsv")), header = FALSE)
  setnames(dt, c("chromosome","position","ref","alt","rsids","nearest_genes","pval","mlogp",
                 "beta","sebeta","af_alt","af_alt_cases","af_alt_controls"))
  dt[, key := paste(chromosome, position, sep = ":")]; dt
}
gwas_list <- list(
  Eczema = list(dt = read_gwas("L12_DERMATITISECZEMA"), cases = 67474, total = 500348),
  AD     = list(dt = read_gwas("L12_ATOPIC"),           cases = 31245, total = 500348))

read_eqtl_gene <- function(file, ensg) {
  dt <- fread(file.path(DIR, file), header = FALSE)
  setnames(dt, c("molecular_trait_id","chromosome","position","ref","alt","variant",
                 "ma_samples","maf","pvalue","beta","se","type","ac","an","r2",
                 "mt_object_id","gene_id","median_tpm","rsid"))
  dt <- dt[gene_id == ensg]
  if (nrow(dt) == 0) return(NULL)
  probe_used <- NA_character_
  if (length(unique(dt$molecular_trait_id)) > 1) {  # microarray: best probe by min p
    probe_used <- dt[, .(minp = min(pvalue)), by = molecular_trait_id][order(minp)][1, molecular_trait_id]
    dt <- dt[molecular_trait_id == probe_used]
  }
  setorder(dt, pvalue); dt <- dt[!duplicated(variant)]
  dt[, key := paste(chromosome, position, sep = ":")]
  attr(dt, "probe") <- probe_used
  dt
}

build_merged <- function(eq, gw) {
  mg <- merge(eq, gw, by = "key", suffixes = c("", "_g"))
  ex <- mg$ref == mg$ref_g & mg$alt == mg$alt_g
  sw <- mg$ref == mg$alt_g & mg$alt == mg$ref_g
  mg <- mg[ex | sw]; s2 <- sw[ex | sw]
  mg[, beta_gwas := ifelse(s2, -beta_g, beta_g)]
  mg[, maf_use := pmin(maf, 1 - maf)]
  mg[, snp := paste(chromosome, position, ref, alt, sep = ":")]
  setorder(mg, pvalue); mg <- mg[!duplicated(snp)]  # drop multiallelic/dup-key duplicates
  mg
}

run_abf <- function(mg, gw_meta, N1, p12 = 1e-5) {
  d1 <- list(beta = mg$beta, varbeta = mg$se^2, N = N1, type = "quant",
             MAF = mg$maf_use, snp = mg$snp)
  d2 <- list(beta = mg$beta_gwas, varbeta = mg$sebeta^2, N = gw_meta$total, type = "cc",
             s = gw_meta$cases / gw_meta$total, snp = mg$snp)
  coloc.abf(d1, d2, p12 = p12)$summary
}

results <- list(); k <- 1
for (gi in seq_len(nrow(genes))) {
  sym <- genes$symbol[gi]; ensg <- genes$ensg[gi]
  for (di in seq_len(nrow(meta))) {
    eq <- read_eqtl_gene(meta$file[di], ensg)
    if (is.null(eq)) {
      results[[k]] <- data.table(gene = sym, study = meta$study[di], cell = meta$cell[di],
                                 platform = meta$platform[di], endpoint = NA, nsnps = 0,
                                 note = "gene not tested (failed expression QC)"); k <- k + 1; next
    }
    tpm <- if (meta$platform[di] == "RNAseq") unique(eq$median_tpm)[1] else NA_real_
    probe <- attr(eq, "probe")
    for (ep_name in names(gwas_list)) {
      mg <- build_merged(eq, gwas_list[[ep_name]]$dt)
      if (nrow(mg) < 50) {
        results[[k]] <- data.table(gene = sym, study = meta$study[di], cell = meta$cell[di],
                                   platform = meta$platform[di], endpoint = ep_name,
                                   nsnps = nrow(mg), note = "<50 shared SNPs"); k <- k + 1; next
      }
      pp <- run_abf(mg, gwas_list[[ep_name]], meta$N[di])
      sens <- sapply(c(1e-6, 5e-6, 1e-5, 5e-5, 1e-4),
                     function(p12) run_abf(mg, gwas_list[[ep_name]], meta$N[di], p12)["PP.H4.abf"])
      ws <- sapply(c(250000, 125000), function(w) {
        msub <- mg[abs(position - 41487062) <= w]
        if (nrow(msub) < 50) return(NA_real_)
        run_abf(msub, gwas_list[[ep_name]], meta$N[di])["PP.H4.abf"]
      })
      lead_e <- mg[which.min(pvalue)]; lead_g <- mg[which.min(pval)]
      # rs11635906-specific row (signal A SNP)
      r1 <- mg[rsid == "rs11635906"]
      results[[k]] <- data.table(
        gene = sym, study = meta$study[di], cell = meta$cell[di], platform = meta$platform[di],
        endpoint = ep_name, nsnps = nrow(mg), N1 = meta$N[di],
        median_tpm = tpm, probe = probe,
        PP.H0 = pp["PP.H0.abf"], PP.H1 = pp["PP.H1.abf"], PP.H2 = pp["PP.H2.abf"],
        PP.H3 = pp["PP.H3.abf"], PP.H4 = pp["PP.H4.abf"],
        PP.H4_p12_1e6 = sens[1], PP.H4_p12_5e6 = sens[2], PP.H4_p12_1e5 = sens[3],
        PP.H4_p12_5e5 = sens[4], PP.H4_p12_1e4 = sens[5],
        PP.H4_w250 = ws[1], PP.H4_w125 = ws[2],
        min_p_eqtl = min(mg$pvalue), min_p_gwas = min(mg$pval),
        lead_eqtl_rsid = lead_e$rsid, lead_eqtl_p = lead_e$pvalue,
        lead_gwas_rsid = lead_g$rsid, lead_gwas_p = lead_g$pval,
        rs11635906_in = nrow(r1) > 0,
        rs11635906_eqtl_p = if (nrow(r1)) r1$pvalue else NA_real_,
        rs11635906_eqtl_beta = if (nrow(r1)) r1$beta else NA_real_,
        note = "ok")
      k <- k + 1
    }
  }
  cat(sym, "done\n")
}
res <- rbindlist(results, fill = TRUE)
fwrite(res, file.path(DIR, "coloc_tcell_results.csv"))
ok <- res[note == "ok"]
print(ok[, .(gene, study, cell, endpoint, nsnps, min_p_eqtl = signif(min_p_eqtl, 2),
             PP.H0 = round(PP.H0, 3), PP.H3 = round(PP.H3, 3), PP.H4 = round(PP.H4, 3),
             p12_lo = round(PP.H4_p12_1e6, 3), p12_hi = round(PP.H4_p12_1e4, 3))])
