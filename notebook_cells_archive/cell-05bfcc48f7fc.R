suppressMessages(library(ggplot2))
theme_set(theme_classic(base_size=11, base_family="Liberation Sans"))

# ---- 结果汇总表 ----
res_tab <- data.frame(
  数据集=c("GSE121212","GSE121212","GSE121212","GSE157194","GSE157194","GSE157194"),
  对比=c("皮损aLS vs 非皮损NL (配对48对)","皮损aLS vs 健康CTRL (49 vs 38)","配对Wilcoxon交叉验证",
         "皮损AL vs 非皮损AN m0 (配对54对)","dupilumab: AL m3 vs m0 (配对21对)","环孢素: AL m3 vs m0 (配对8对)"),
  log2FC=c(0.013, 0.196, NA, -0.068, -0.177, 0.13),
  p值=c(0.958, 0.449, 0.471, 0.648, 0.483, 0.768),
  判定=rep("不显著（零结果）", 6))
write.csv(res_tab, "/mnt/results/P2-02_GEO皮损表达验证_结果表.csv", row.names=FALSE, fileEncoding="UTF-8")
print(res_tab)

# ---- 配对图 ----
# GSE121212: NL -> aLS
o <- order(cd_p$patient[cd_p$cond=="NL"])
dfA <- data.frame(
  patient=rep(cd_p$patient[cd_p$cond=="NL"][o], 2),
  cond=rep(c("NL","aLS"), each=length(o)),
  val=c(lcpm["ITPKA", i_n], lcpm["ITPKA", i_a]))
dfA$cond <- factor(dfA$cond, levels=c("NL","aLS"))

# GSE157194: AN -> AL (m0)
lib2 <- colSums(mat2[, sub2]); lcpm2 <- log2(t(t(mat2[, sub2])/lib2*1e6)+1)
j_n <- cd2$sample[cd2$site=="AN"]; j_a <- cd2$sample[cd2$site=="AL"]
ord <- order(cd2$patient[cd2$site=="AN"])
j_n <- j_n[ord]; j_a <- j_a[ord]
dfB <- data.frame(patient=rep(cd2$patient[cd2$site=="AN"][ord], 2),
                  cond=rep(c("AN","AL"), each=length(ord)),
                  val=c(lcpm2["ENSG00000137825", j_n], lcpm2["ENSG00000137825", j_a]))
dfB$cond <- factor(dfB$cond, levels=c("AN","AL"))

mk <- function(df, title, labels) {
  ggplot(df, aes(cond, val)) +
    geom_line(aes(group=patient), color="grey80", linewidth=0.3) +
    geom_boxplot(aes(fill=cond), width=0.35, outlier.shape=NA, alpha=0.85) +
    geom_point(position=position_jitter(width=0.06), size=0.7, alpha=0.6) +
    scale_fill_manual(values=c("#0279EE","#FF9400"), guide="none") +
    scale_x_discrete(labels=labels) +
    labs(title=title, x=NULL, y="ITPKA log2(CPM+1)")
}
pA <- mk(dfA, "GSE121212 (48 pairs)", c("Non-lesional","Lesional AD"))
pB <- mk(dfB, "GSE157194 (54 pairs)", c("Non-lesional","Lesional AD"))
comb <- cowplot::plot_grid(pA, pB, ncol=2)
ggsave("/workspace/fig_itpka_geo.png", comb, width=8, height=4, dpi=300)
ggsave("/workspace/fig_itpka_geo.svg", comb, width=8, height=4)
system("cp /workspace/fig_itpka_geo.png /mnt/results/P2-02_fig_ITPKA_GEO皮损表达.png && cp /workspace/fig_itpka_geo.svg /mnt/results/P2-02_fig_ITPKA_GEO皮损表达.svg")
cat("图表已保存\n")