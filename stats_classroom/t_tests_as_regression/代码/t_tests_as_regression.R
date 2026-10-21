# -*- coding: UTF-8 -*-
# 统计课堂：t检验，原来也是回归
# 目的：用普通最小二乘模型重写单样本、两独立样本和配对t检验。
# 运行：打开项目Rproj后在RStudio Source本文件；也可用Rscript运行。
# 依赖：R自带的 stats；共享输入在 数据/，结果在 运行结果/R/。
# 约定：A组编码为1，B组编码为0；所有观测值使用原始单位g或分。
# 结果对象在每一节末尾保留，RStudio中可直接点开查看；数值同时写入CSV便于核对。

# %% 定位项目、创建输出目录、读取共享数据
args <- commandArgs(FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
source_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (length(file_arg) > 0) {
  root <- dirname(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
} else if (!is.null(source_file)) {
  root <- dirname(dirname(normalizePath(source_file)))
} else {
  root <- getwd()
}

out <- file.path(root, "运行结果", "R")
fig_dir <- file.path(root, "文章配图")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

rat <- read.csv(file.path(root, "数据", "rat_gain_long.csv"),
                stringsAsFactors = FALSE, fileEncoding = "UTF-8")
paired <- read.csv(file.path(root, "数据", "paired_drinks.csv"),
                   stringsAsFactors = FALSE, fileEncoding = "UTF-8")

# 一行是一只大鼠；先直接查看字段与前几行，确认没有把长表读成矩阵。
rat_head <- head(rat)
paired_head <- head(paired)
rat_head
paired_head

# %% 1. 单样本t：只有截距的模型
# lm(y ~ 1)中没有分组变量，只有截距；默认的原假设是截距=0。
# 这等价于单样本t检验均数是否为0。B组均数约101 g，比较0通常没有教学意义。
b_data <- rat[rat$group == "B", ]
fit_b_default <- lm(gain_g ~ 1, data = b_data)

# I(y - mu0)把明确的参考值移到0，再用同一个截距模型检验。
# 这里的100 g是本项目人为提出的课堂参照，不是教材例4-8的原问题。
mu0 <- 100
fit_b_reference <- lm(I(gain_g - mu0) ~ 1, data = b_data)

extract_fit <- function(model, example, term_index = 1) {
  coefficient_table <- summary(model)$coefficients
  interval <- confint(model, level = 0.95)[term_index, ]
  data.frame(
    example = example,
    term = rownames(coefficient_table)[term_index],
    estimate = coefficient_table[term_index, "Estimate"],
    standard_error = coefficient_table[term_index, "Std. Error"],
    df = df.residual(model),
    t_value = coefficient_table[term_index, "t value"],
    p_value = coefficient_table[term_index, "Pr(>|t|)"],
    ci_low = interval[1],
    ci_high = interval[2],
    row.names = NULL
  )
}

single_results <- rbind(
  extract_fit(fit_b_default, "single_B_vs_0"),
  extract_fit(fit_b_reference, "single_B_vs_100")
)
write.csv(single_results, file.path(out, "01_single_sample.csv"), row.names = FALSE)

# %% 2. 两独立样本：0/1组变量的普通OLS
# model: gain_g = beta0 + beta1 * group01 + error。
# group01=0时预测值是B组均数，所以截距=mean(B)；group01=1时
# 预测值为截距+beta1=mean(A)，因此beta1=mean(A)-mean(B)。
fit_groups <- lm(gain_g ~ group01, data = rat)
group_coefficient_table <- summary(fit_groups)$coefficients
group_interval <- confint(fit_groups, level = 0.95)
group_results <- data.frame(
  example = "two_independent_groups_pooled_OLS",
  term = rownames(group_coefficient_table),
  estimate = group_coefficient_table[, "Estimate"],
  standard_error = group_coefficient_table[, "Std. Error"],
  df = rep(df.residual(fit_groups), nrow(group_coefficient_table)),
  t_value = group_coefficient_table[, "t value"],
  p_value = group_coefficient_table[, "Pr(>|t|)"],
  ci_low = group_interval[, 1],
  ci_high = group_interval[, 2],
  row.names = NULL
)
write.csv(group_results, file.path(out, "02_group_coefficients.csv"), row.names = FALSE)

# 两组只有一个非截距系数，因此整体F与A-B系数的t平方相等。
f_statistic <- summary(fit_groups)$fstatistic
overall_F <- unname(f_statistic["value"])
overall_F_p <- pf(overall_F, unname(f_statistic["numdf"]),
                  unname(f_statistic["dendf"]), lower.tail = FALSE)
overall_F_result <- data.frame(
  model = "group_OLS",
  F_value = overall_F,
  numerator_df = unname(f_statistic["numdf"]),
  denominator_df = unname(f_statistic["dendf"]),
  p_value = overall_F_p
)
write.csv(overall_F_result, file.path(out, "03_group_overall_F.csv"), row.names = FALSE)

# 作为核验，R的t.test(var.equal=TRUE)使用同一个合并方差与自由度17。
a <- rat$gain_g[rat$group == "A"]
b <- rat$gain_g[rat$group == "B"]
pooled_t <- t.test(a, b, var.equal = TRUE)
welch_t <- t.test(a, b, var.equal = FALSE)
welch_result <- data.frame(
  method = c("ordinary_OLS_equal_variance", "Welch_unequal_variance"),
  estimate_A_minus_B = c(mean(a) - mean(b), unname(welch_t$estimate[1] - welch_t$estimate[2])),
  t_value = c(unname(pooled_t$statistic), unname(welch_t$statistic)),
  df = c(unname(pooled_t$parameter), unname(welch_t$parameter)),
  p_value = c(pooled_t$p.value, welch_t$p.value),
  ci_low = c(pooled_t$conf.int[1], welch_t$conf.int[1]),
  ci_high = c(pooled_t$conf.int[2], welch_t$conf.int[2])
)
write.csv(welch_result, file.path(out, "04_pooled_vs_welch.csv"), row.names = FALSE)

# %% 3. 配对t：先按ID作差，再拟合只有截距的模型
# 同一个人的两次评分不是两组独立对象。先计算每个人的delta=new-old，
# 再检验平均delta是否为0；这与配对t检验相同。
paired$delta <- paired$new - paired$old
fit_paired <- lm(delta ~ 1, data = paired)
paired_result <- extract_fit(fit_paired, "paired_new_minus_old")
paired_check <- t.test(paired$delta, mu = 0)
paired_result$check_t_value <- unname(paired_check$statistic)
paired_result$check_p_value <- paired_check$p.value
write.csv(paired_result, file.path(out, "05_paired_difference.csv"), row.names = FALSE)
write.csv(paired, file.path(out, "05_paired_with_delta.csv"), row.names = FALSE)

# %% 4. 用图看见“截距+组差”和新旧款配对变化
# 两张图使用固定绘图参数；英文轴标签可避免不同电脑缺少中文字体，正文图注提供中文解释。
mean_a <- mean(a)
mean_b <- mean(b)
x_a <- rep(1, length(a)) + seq(-0.11, 0.11, length.out = length(a))
x_b <- rep(0, length(b)) + seq(-0.11, 0.11, length.out = length(b))
png(file.path(fig_dir, "01_组均数与OLS截距.png"), width = 1600, height = 1000, res = 160)
par(mar = c(4.5, 4.5, 3.2, 1.2))
plot(c(-0.25, 1.25), range(rat$gain_g) + c(-8, 8), type = "n",
     xaxt = "n", xlab = "Group code (0=B, 1=A)", ylab = "Gain (g)",
     main = "Two-group OLS: intercept and group difference")
axis(1, at = c(0, 1), labels = c("B (0)", "A (1)"))
points(x_b, b, pch = 19, col = "#1769AA")
points(x_a, a, pch = 19, col = "#D65A4A")
segments(-0.16, mean_b, 0.16, mean_b, lwd = 4, col = "#1769AA")
segments(0.84, mean_a, 1.16, mean_a, lwd = 4, col = "#D65A4A")
legend("topleft", legend = c("B observations", "A observations", "group mean"),
       pch = c(19, 19, NA), lty = c(NA, NA, 1), lwd = c(NA, NA, 4),
       col = c("#1769AA", "#D65A4A", "#333333"), bty = "n")
text(0, mean_b + 8, sprintf("intercept = %.2f", mean_b), col = "#1769AA")
text(1, mean_a + 8, sprintf("A-B = %.2f", mean_a - mean_b), col = "#D65A4A")
dev.off()
file.copy(file.path(fig_dir, "01_组均数与OLS截距.png"),
          file.path(out, "01_组均数与OLS截距.png"), overwrite = TRUE)

x_old <- rep(1, nrow(paired))
x_new <- rep(2, nrow(paired))
png(file.path(fig_dir, "02_新旧款配对与差值.png"), width = 1600, height = 1000, res = 160)
par(mar = c(4.5, 4.5, 3.2, 1.2))
plot(c(0.75, 2.25), range(c(paired$old, paired$new)) + c(-0.5, 0.5), type = "n",
     xaxt = "n", xlab = "Condition", ylab = "Score",
     main = "Paired data: keep each person's connection")
axis(1, at = c(1, 2), labels = c("Old", "New"))
for (i in seq_len(nrow(paired))) {
  segments(1, paired$old[i], 2, paired$new[i], col = "#9AA6B2")
}
points(x_old, paired$old, pch = 19, col = "#1769AA")
points(x_new, paired$new, pch = 19, col = "#D65A4A")
points(c(1, 2), c(mean(paired$old), mean(paired$new)), pch = 18, cex = 1.6, col = "#222222")
legend("topleft", legend = c("Old", "New", "mean"), pch = c(19, 19, 18),
       col = c("#1769AA", "#D65A4A", "#222222"), bty = "n")
text(1.5, max(c(paired$old, paired$new)) + 0.25,
     sprintf("mean delta (new-old) = %.2f", mean(paired$delta)), col = "#222222")
dev.off()
file.copy(file.path(fig_dir, "02_新旧款配对与差值.png"),
          file.path(out, "02_新旧款配对与差值.png"), overwrite = TRUE)

# %% 5. 汇总可核对结果
key_results <- rbind(
  single_results[single_results$example == "single_B_vs_100",
                 c("example", "term", "estimate", "standard_error", "df",
                   "t_value", "p_value", "ci_low", "ci_high")],
  group_results[group_results$term == "group01",
                c("example", "term", "estimate", "standard_error", "df",
                  "t_value", "p_value", "ci_low", "ci_high")],
  paired_result[, c("example", "term", "estimate", "standard_error", "df",
                    "t_value", "p_value", "ci_low", "ci_high")]
)
write.csv(key_results, file.path(out, "all_key_results.csv"), row.names = FALSE)

# RStudio中直接查看这些结构化结果；没有自动把结果写成结论。
results <- list(
  single_sample = single_results,
  group_coefficients = group_results,
  overall_F = overall_F_result,
  pooled_vs_welch = welch_result,
  paired_difference = paired_result,
  key_results = key_results
)
results
