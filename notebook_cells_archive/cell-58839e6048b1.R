suppressMessages(library(DESeq2))

# ---- GSE121212 样本表 ----
mat1 <- as.matrix(c1[, -1]); rownames(mat1) <- genes1
parse1 <- function(x) {
  if (grepl("^CTRL", x)) return(c(patient=x, cond="CTRL"))
  p <- strsplit(x, "_")[[1]]
  patient <- paste(p[1], p[2], sep="_")
  cond <- ifelse(length(p)==3 && p[3]=="lesional", "aLS",
          ifelse(length(p)==3 && p[3]=="non-lesional", "NL", "cLS"))
  c(patient=patient, cond=cond)
}
pd1 <- as.data.frame(t(sapply(colnames(mat1), parse1)))
pd1$cond <- factor(pd1$cond, levels=c("CTRL","NL","aLS","cLS"))
print(table(pd1$cond))

# 配对 aLS vs NL（同一患者两者都有）
both <- intersect(pd1$patient[pd1$cond=="aLS"], pd1$patient[pd1$cond=="NL"])
cat("aLS-NL 完整配对患者数:", length(both), "\n")
sub_p <- rownames(pd1)[pd1$patient %in% both & pd1$cond %in% c("aLS","NL")]
cd_p <- droplevels(pd1[sub_p, ]); cd_p$patient <- factor(cd_p$patient)

dds_p <- DESeqDataSetFromMatrix(round(mat1[, sub_p]), cd_p, ~ patient + cond)
dds_p <- DESeq(dds_p, quiet=TRUE)
res_p <- results(dds_p, contrast=c("cond","aLS","NL"))
it_p <- res_p["ITPKA", ]
cat("\n[GSE121212 配对 aLS vs NL, n=", length(both), "对] ITPKA: log2FC=", round(it_p$log2FoldChange,3),
    " p=", signif(it_p$pvalue,3), " padj=", signif(it_p$padj,3), "\n", sep="")

# 非配对 aLS vs CTRL
sub_u <- rownames(pd1)[pd1$cond %in% c("aLS","CTRL")]
cd_u <- droplevels(pd1[sub_u, ])
dds_u <- DESeqDataSetFromMatrix(round(mat1[, sub_u]), cd_u, ~ cond)
dds_u <- DESeq(dds_u, quiet=TRUE)
res_u <- results(dds_u, contrast=c("cond","aLS","CTRL"))
it_u <- res_u["ITPKA", ]
cat("[GSE121212 aLS vs CTRL, ", sum(cd_u$cond=='aLS'), " vs ", sum(cd_u$cond=='CTRL'),
    "] ITPKA: log2FC=", round(it_u$log2FoldChange,3), " p=", signif(it_u$pvalue,3),
    " padj=", signif(it_u$padj,3), "\n", sep="")

# 非参数交叉验证：配对 Wilcoxon on log2(CPM+1)
lib <- colSums(mat1[, sub_p])
cpm <- t(t(mat1[, sub_p]) / lib * 1e6)
lcpm <- log2(cpm + 1)
i_a <- sub_p[cd_p$cond=="aLS"]; i_n <- sub_p[cd_p$cond=="NL"]
i_a <- i_a[order(cd_p$patient[cd_p$cond=="aLS"])]; i_n <- i_n[order(cd_p$patient[cd_p$cond=="NL"])]
wt <- wilcox.test(lcpm["ITPKA", i_a], lcpm["ITPKA", i_n], paired=TRUE)
cat("[配对Wilcoxon交叉验证] ITPKA log2CPM: aLS中位=", round(median(lcpm['ITPKA', i_a]),2),
    " NL中位=", round(median(lcpm['ITPKA', i_n]),2), " 配对p=", signif(wt$p.value,3), "\n", sep="")

# 全基因组背景：ITPKA在标准过滤后是否进入FDR集合
keep_std <- rowSums(counts(dds_u) >= 10) >= 10
cat("\n标准过滤(>=10计数 in >=10样本)后基因数:", sum(keep_std), "/", nrow(mat1),
    "；ITPKA被过滤掉:", !keep_std["ITPKA"], "\n")
cat("方向判定（预注册：MR称表达↑→风险↑，皮损应↑）: aLS vs NL log2FC", 
    ifelse(it_p$log2FoldChange>0, ">0 方向一致", "<0 方向不一致"), "\n")