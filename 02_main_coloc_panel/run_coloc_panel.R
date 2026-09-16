#!/usr/bin/env Rscript
# P2-02 位点内多基因 coloc 鉴别分析
# 问题: rs11635906 位点内多个基因被调控, 哪个基因的 eQTL 与湿疹/AD GWAS 共定位?
# 设计: 同一 SNP 宇宙 (变异±500kb), 同一 GWAS, 逐基因逐数据集 coloc.abf + p12 敏感性
suppressMessages({library(data.table); library(coloc)})

DIR <- "/workspace/coloc_panel"
WIN_CHR <- 15; WIN_LO <- 40987062; WIN_HI <- 41987062  # rs11635906 (41,487,062) ±500kb, hg38

# ---------- 基因面板 (数据驱动: GTEx v8 中 rs11635906 p<1e-4 的全部基因) ----------
panel <- data.table(
  gene_id = c("ENSG00000137825","ENSG00000137806","ENSG00000187446","ENSG00000103932","ENSG00000137815","ENSG00000247556"),
  symbol  = c("ITPKA","NDUFAF1","CHP1","RPAP1","RTF1","OIP5-AS1")
)

# ---------- 读取 eQTL 窗口文件 ----------
# GTEx v8 import 格式列
cols_v8 <- c("variant","r2","pvalue","mto_id","molecular_trait_id","maf","gene_id",
             "median_tpm","beta","se","an","ac","chromosome","position","ref","alt","type","rsid")
# catalogue sumstats (QTD) 格式列
cols_qtd <- c("molecular_trait_id","chromosome","position","ref","alt","variant","ma_samples","maf",
              "pvalue","beta","se","type","ac","an","r2","mto_id","gene_id","median_tpm","rsid")

read_eqtl <- function(file, fmt) {
  dt <- fread(file, header = FALSE)
  setnames(dt, if (fmt == "v8") cols_v8 else cols_qtd)
  dt[, `:=`(position = as.integer(position), beta = as.numeric(beta), se = as.numeric(se),
            maf = as.numeric(maf), pvalue = as.numeric(pvalue), an = as.numeric(an))]
  dt
}

datasets <- list(
  GTExv8_WholeBlood = list(file = "eqtl_gtexv8_whole_blood.tsv",  fmt = "v8"),
  GTExv8_SkinSE     = list(file = "eqtl_gtexv8_skin_se.tsv",      fmt = "v8"),
  GTExv8_SkinNSE    = list(file = "eqtl_gtexv8_skin_nse.tsv",     fmt = "v8"),
  BLUEPRINT_Neutro  = list(file = "eqtl_blueprint_neutrophil.tsv",fmt = "qtd"),
  BLUEPRINT_Mono    = list(file = "eqtl_blueprint_monocyte.tsv",  fmt = "qtd")
)

# ---------- 读取 FinnGen GWAS 窗口 ----------
read_gwas <- function(ep) {
  dt <- fread(file.path(DIR, paste0("gwas_finngen_", ep, ".tsv")), header = FALSE)
  setnames(dt, c("chromosome","position","ref","alt","rsids","nearest_genes","pval","mlogp",
                 "beta","sebeta","af_alt","af_alt_cases","af_alt_controls"))
  dt[, `:=`(position = as.integer(position), beta = as.numeric(beta), sebeta = as.numeric(sebeta),
            af_alt = as.numeric(af_alt))]
  dt
}
gwas_list <- list(
  Eczema = list(dt = read_gwas("L12_DERMATITISECZEMA"), cases = 67474, total = 500348),
  AD     = list(dt = read_gwas("L12_ATOPIC"),           cases = 31245, total = 500348)
)

# ---------- 等位基因调和 ----------
harmonise <- function(eq, gw) {
  # 键: chromosome:position; eQTL beta 相对 alt, GWAS beta 相对 alt
  eq[, key := paste(chromosome, position, sep = ":")]
  gw[, key := paste(chromosome, position, sep = ":")]
  m <- merge(eq, gw, by = "key", suffixes = c("_e", "_g"))
  exact <- m$ref_e == m$ref_g & m$alt_e == m$alt_g
  swap  <- m$ref_e == m$alt_g & m$alt_e == m$ref_g
  keep  <- exact | swap
  m <- m[keep]
  sw <- swap[keep]
  m[, beta_g := ifelse(sw, -beta_g, beta_g)]              # 翻转到 eQTL alt 等位
  m[, af_g   := ifelse(sw, 1 - af_alt, af_alt)]
  # 去重 (同一 key 多行: 取 eQTL p 最小行)
  setorder(m, pvalue)
  m <- m[!duplicated(paste(key, molecular_trait_id))]
  m
}

# ---------- 主循环 ----------
results <- list(); k <- 1
for (ds_name in names(datasets)) {
  eq_all <- read_eqtl(file.path(DIR, datasets[[ds_name]]$file), datasets[[ds_name]]$fmt)
  n1 <- round(unique(eq_all$an)[1] / 2)   # 样本量 = 等位基因数/2
  for (i in seq_len(nrow(panel))) {
    gid <- panel$gene_id[i]; sym <- panel$symbol[i]
    eq <- eq_all[gene_id == gid]
    if (nrow(eq) < 50) {
      results[[k]] <- data.table(dataset = ds_name, gene = sym, endpoint = NA, nsnps = nrow(eq),
                                 note = "gene absent or <50 SNPs in dataset"); k <- k + 1; next
    }
    for (ep_name in names(gwas_list)) {
      gw <- gwas_list[[ep_name]]$dt
      m <- harmonise(eq, gw)
      if (nrow(m) < 50) {
        results[[k]] <- data.table(dataset = ds_name, gene = sym, endpoint = ep_name, nsnps = nrow(m),
                                   note = "<50 shared SNPs after harmonisation"); k <- k + 1; next
      }
      m[, snp := paste(chromosome_e, position_e, ref_e, alt_e, sep = ":")]
      d1 <- list(beta = m$beta_e, varbeta = m$se^2, N = n1, type = "quant",
                 MAF = pmin(m$maf, 1 - m$maf), snp = m$snp)
      d2 <- list(beta = m$beta_g, varbeta = m$sebeta^2, N = gwas_list[[ep_name]]$total,
                 type = "cc", s = gwas_list[[ep_name]]$cases / gwas_list[[ep_name]]$total, snp = m$snp)
      # 主分析: 默认先验; 敏感性: p12 扫描
      pp_main <- coloc.abf(d1, d2)$summary
      sens <- sapply(c(1e-6, 5e-6, 1e-5, 5e-5, 1e-4), function(p12)
        coloc.abf(d1, d2, p12 = p12)$summary["PP.H4.abf"])
      lead_e <- m[which.min(pvalue), snp]; lead_g <- m[which.min(pval), snp]
      results[[k]] <- data.table(
        dataset = ds_name, gene = sym, endpoint = ep_name, nsnps = nrow(m), N1 = n1,
        PP.H0 = pp_main["PP.H0.abf"], PP.H1 = pp_main["PP.H1.abf"], PP.H2 = pp_main["PP.H2.abf"],
        PP.H3 = pp_main["PP.H3.abf"], PP.H4 = pp_main["PP.H4.abf"],
        PP.H4_p12_1e6 = sens[1], PP.H4_p12_5e6 = sens[2], PP.H4_p12_1e5 = sens[3],
        PP.H4_p12_5e5 = sens[4], PP.H4_p12_1e4 = sens[5],
        lead_eqtl_snp = lead_e, lead_gwas_snp = lead_g,
        min_p_eqtl = min(m$pvalue), min_p_gwas = min(m$pval),
        beta_lead_eqtl = m[which.min(pvalue), beta_e], beta_lead_gwas_at_eqtl_lead = m[which.min(pvalue), beta_g],
        note = "ok")
      k <- k + 1
    }
  }
  cat(ds_name, "done\n")
}
res <- rbindlist(results, fill = TRUE)
fwrite(res, file.path(DIR, "coloc_panel_results.csv"))
print(res[, .(dataset, gene, endpoint, nsnps, PP.H3 = round(PP.H3,3), PP.H4 = round(PP.H4,3), note)])
