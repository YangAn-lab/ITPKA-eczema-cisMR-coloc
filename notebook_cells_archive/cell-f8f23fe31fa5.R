library(data.table)

# ---- GSE121212 ----
c1 <- fread("/workspace/geo/GSE121212_readcount.txt.gz")
genes1 <- c1[[1]]
cn1 <- colnames(c1)[-1]
# 解析样本类型
grp1 <- sub(".*_(lesional|non-lesional)$", "\\1", cn1)
grp1[!grepl("lesional", grp1)] <- sapply(cn1[!grepl("lesional", grp1)], function(x) paste(strsplit(x, "_")[[1]][1:2], collapse="_"))
print(table(grp1))
cat("样本名示例:", head(cn1, 4), "...", tail(cn1, 6), "\n")
cat("总样本数:", length(cn1), " 基因数:", length(genes1), "\n")

# ITPKA 行
itpka1 <- c1[genes1 == "ITPKA", ]
cat("\nITPKA 在GSE121212中存在:", nrow(itpka1) > 0, "\n")
if (nrow(itpka1) > 0) {
  v <- as.numeric(itpka1[1, -1])
  cat("ITPKA 原始计数: median=", median(v), " mean=", round(mean(v),1), " max=", max(v), " 零计数样本:", sum(v==0), "/", length(v), "\n")
}