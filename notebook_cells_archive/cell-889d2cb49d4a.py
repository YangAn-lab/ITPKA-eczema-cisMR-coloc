# ---- GSE157194 ----
mat2 <- as.matrix(c2[, -1]); rownames(mat2) <- genes2

# 治疗分组：ch1.3 中 "therapy: xxx" 行（m3样本）映射到患者
th <- ph[["characteristics_ch1.3"]]
pid_all <- sub("patient id: ", "", ph[["characteristics_ch1"]])
therapy <- setNames(rep(NA_character_, length(pid_all)), pid_all)
is_th <- grepl("therapy:", th)
therapy[pid_all[is_th]] <- sub("therapy: ", "", th[is_th])
therapy <- therapy[!is.na(therapy)]
print(table(therapy))

meta2$therapy <- therapy[meta2$patient]

# m0 配对 AL vs AN
m0 <- meta2[meta2$month=="m0", ]
both2 <- intersect(m0$patient[m0$site=="AL"], m0$patient[m0$site=="AN"])
cat("m0 完整配对患者数:", length(both2), "\n")
sub2 <- m0$sample[m0$patient %in% both2]
cd2 <- droplevels(m0[m0$sample %in% sub2, ]); cd2$patient <- factor(cd2$patient); cd2$site <- factor(cd2$site, levels=c("AN","AL"))
rownames(cd2) <- cd2$sample
dds2 <- DESeqDataSetFromMatrix(round(mat2[, sub2]), cd2, ~ patient + site)
dds2 <- DESeq(dds2, quiet=TRUE)
res2 <- results(dds2, contrast=c("site","AL","AN"))
it2r <- res2["ENSG00000137825", ]
cat("[GSE157194 m0 配对 AL vs AN, n=", length(both2), "对] ITPKA: log2FC=", round(it2r$log2FoldChange,3),
    " p=", signif(it2r$pvalue,3), " padj=", signif(it2r$padj,3), "\n", sep="")

# 治疗响应：AL m3 vs AL m0，分治疗组（配对）
for (arm in c("dupilumab","cyclosporine")) {
  pa <- names(therapy)[therapy==arm]
  al <- meta2[meta2$site=="AL" & meta2$patient %in% pa, ]
  both3 <- intersect(al$patient[al$month=="m0"], al$patient[al$month=="m3"])
  if (length(both3) < 3) { cat(arm, ": 配对不足 (", length(both3), ")，跳过\n"); next }
  ss <- al$sample[al$patient %in% both3]
  cd3 <- droplevels(al[al$sample %in% ss, ]); cd3$patient <- factor(cd3$patient)
  cd3$month <- factor(cd3$month, levels=c("m0","m3")); rownames(cd3) <- cd3$sample
  dds3 <- DESeqDataSetFromMatrix(round(mat2[, ss]), cd3, ~ patient + month)
  dds3 <- DESeq(dds3, quiet=TRUE)
  r3 <- results(dds3, contrast=c("month","m3","m0"))["ENSG00000137825", ]
  cat("[", arm, " AL m3 vs m0, n=", length(both3), "对] ITPKA: log2FC=", round(r3$log2FoldChange,3),
      " p=", signif(r3$pvalue,3), "\n", sep="")
}