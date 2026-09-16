# 分析G：固定效应逆方差加权 meta 分析（审稿人 v2 意见一.2）
# 主分析：FinnGen AD + Oliva EUR + EAGLE（均为 AD 表型、欧洲裔）
# 敏感性：FinnGen 湿疹 + Oliva 多祖先 + EAGLE
library(data.table)

dt <- fread("/mnt/results/P2-02_跨队列验证_rs11635906.csv")

meta_fe <- function(beta, se, label) {
  w <- 1 / se^2
  b_meta <- sum(w * beta) / sum(w)
  se_meta <- sqrt(1 / sum(w))
  z <- b_meta / se_meta
  p <- 2 * pnorm(-abs(z))
  Q <- sum(w * (beta - b_meta)^2)
  df <- length(beta) - 1
  p_Q <- pchisq(Q, df, lower.tail = FALSE)
  I2 <- max(0, (Q - df) / Q) * 100
  data.table(analysis = label, k = length(beta),
             beta_meta = round(b_meta, 4), se_meta = round(se_meta, 5),
             z = round(z, 2), p_meta = signif(p, 3),
             OR = round(exp(b_meta), 3),
             OR_lo = round(exp(b_meta - 1.96 * se_meta), 3),
             OR_hi = round(exp(b_meta + 1.96 * se_meta), 3),
             Q = round(Q, 2), Q_df = df, p_het = signif(p_Q, 3), I2_pct = round(I2, 1))
}

# 主分析：AD 表型、欧洲裔
main <- dt[cohort %in% c("FinnGen R12 L12_ATOPIC", "Oliva 2025 EUR", "EAGLE 2015 (European analysis)")]
res_main <- meta_fe(main$beta_G, main$se, "Primary: FinnGen AD + Oliva EUR + EAGLE (AD, EUR)")

# 敏感性 1：FinnGen 湿疹替代 FinnGen AD
s1 <- dt[cohort %in% c("FinnGen R12 L12_DERMATITISECZEMA", "Oliva 2025 EUR", "EAGLE 2015 (European analysis)")]
res_s1 <- meta_fe(s1$beta_G, s1$se, "Sensitivity 1: FinnGen eczema + Oliva EUR + EAGLE")

# 敏感性 2：Oliva 多祖先替代 Oliva EUR
s2 <- dt[cohort %in% c("FinnGen R12 L12_ATOPIC", "Oliva 2025 multi-ancestry", "EAGLE 2015 (European analysis)")]
res_s2 <- meta_fe(s2$beta_G, s2$se, "Sensitivity 2: FinnGen AD + Oliva multi-ancestry + EAGLE")

# 敏感性 3：全部 7 队列（含欠功效的 ASN/AFR）
res_s3 <- meta_fe(dt$beta_G, dt$se, "Sensitivity 3: all 7 cohorts/strata")

res <- rbind(res_main, res_s1, res_s2, res_s3)
print(res)
fwrite(res, "/workspace/coloc_panel/meta_analysis_rs11635906.csv")
