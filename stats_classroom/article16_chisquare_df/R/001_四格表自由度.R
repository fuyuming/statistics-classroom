# 统计课堂16：四格表为什么自由度是1 ----
# 按节选运行；对象名直接查看结果。

# ① 看数据与立假设：数据与边际合计 ----
O <- matrix(c(225, 275, 128, 272), nrow=2, byrow=TRUE,
            dimnames=list(c("女性", "男性"), c("收看", "未收看")))
row_total <- rowSums(O)
col_total <- colSums(O)
N <- sum(O)
O
c(row_total=row_total[1], row_total_male=row_total[2],
  col_watch=col_total[1], col_not_watch=col_total[2], N=N)

# ② 选方法与查条件：自由度与期望数 ----
df <- (nrow(O)-1)*(ncol(O)-1)
E <- outer(row_total, col_total) / N
E
c(df=df, minimum_expected=min(E))

# ③ 做检验与看效应：未校正 Pearson χ²、Yates 校正与两率 ----
cell_contribution <- (O-E)^2/E
chi_uncorrected <- sum(cell_contribution)
p_uncorrected <- pchisq(chi_uncorrected, df=df, lower.tail=FALSE)
c(chi_uncorrected=chi_uncorrected, df=df, p_uncorrected=p_uncorrected)
cell_contribution

# R 官方函数：Yates 默认校正的关联检验与两率检验 ----
fit_association <- chisq.test(O)
fit_rate <- prop.test(c(O[1,1], O[2,1]), row_total)
fit_association
fit_rate

# pooled Z、Z² 与率差；Yates 校正用连续性修正作对照 ----
x <- c(O[1,1], O[2,1])
p_hat <- x/row_total
p_pool <- sum(x)/N
rate_difference <- p_hat[1]-p_hat[2]
se_null <- sqrt(p_pool*(1-p_pool)*sum(1/row_total))
z_pooled <- rate_difference/se_null
z_squared <- z_pooled^2
se_estimate <- sqrt(sum(p_hat*(1-p_hat)/row_total))
rate_difference_ci_uncorrected <- rate_difference + c(-1,1)*qnorm(.975)*se_estimate
# 与R prop.test默认校正一致的本题近似区间；不是精确区间。
half_width_corrected <- qnorm(.975)*se_estimate + 0.5*sum(1/row_total)
ci_corrected <- rate_difference + c(-1,1)*half_width_corrected
corrected_chi <- unname(fit_association$statistic)
corrected_p <- fit_association$p.value
result <- data.frame(rate_difference=rate_difference,
                     rate_difference_pp=100*rate_difference,
                     ci_low=unname(fit_rate$conf.int[1]),
                     ci_high=unname(fit_rate$conf.int[2]),
                     ci_corrected_low=ci_corrected[1],
                     ci_corrected_high=ci_corrected[2],
                     z_pooled=z_pooled, z_squared=z_squared,
                     chi_uncorrected=chi_uncorrected,
                     p_uncorrected=p_uncorrected,
                     chi_yates=corrected_chi, p_yates=corrected_p,
                     df=df)
result

# ④ 讲结果：写出可复核结果（可选保存） ----
dir.create("output/R", recursive=TRUE, showWarnings=FALSE)
write.csv(result, "output/R/001_result.csv", row.names=FALSE)
