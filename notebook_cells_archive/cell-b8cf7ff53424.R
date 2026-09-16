# ---- GSE157194 ----
c2 <- fread("/workspace/geo/GSE157194_Raw_gene_counts_matrix.txt.gz")
cn2 <- colnames(c2)[-1]
# 解析: Patient_X_AL/AN_m0/m3
parts <- do.call(rbind, strsplit(cn2, "_"))
colnames(parts) <- c("patient0","pid","st","month")
meta2 <- data.frame(sample=cn2, patient=paste(parts[,1], parts[,2], sep="_"),
                    site=parts[,3], month=parts[,4])
print(table(meta2$site, meta2$month))
cat("患者数:", length(unique(meta2$patient)), "\n")

# ITPKA = ENSG00000137825
genes2 <- c2[[1]]
it2 <- c2[genes2 == "ENSG00000137825", ]
cat("ITPKA(ENSG00000137825) 存在:", nrow(it2)>0, "\n")
v2 <- as.numeric(it2[1, -1])
cat("ITPKA 原始计数: median=", median(v2), " mean=", round(mean(v2),1), " max=", max(v2), " 零计数:", sum(v2==0), "/", length(v2), "\n")

# 治疗分组元数据（GEO series matrix）
suppressMessages(library(GEOquery))
gse2 <- getGEO("GSE157194", GSEMatrix = TRUE, getGPL = FALSE)
ph <- pData(gse2[[1]])
keep <- ph[, c("title", intersect(c("characteristics_ch1","characteristics_ch1.1","characteristics_ch1.2","characteristics_ch1.3"), colnames(ph)))]
print(head(keep, 8))
cat("\ncharacteristics 唯一值:\n")
for (cc in grep("characteristics", colnames(ph), value=TRUE)[1:4]) {
  cat(cc, ":", paste(head(unique(ph[[cc]]), 10), collapse=" | "), "\n")
}