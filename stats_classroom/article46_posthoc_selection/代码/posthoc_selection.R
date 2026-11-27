# 统计课堂46｜四种均数比较方法，UTF-8
# 数据：李春喜等《生物统计学》6版例7-1数字；改写支架情境，非真实工程实验。
# 依赖 multcomp；在项目根打开Rproj，依次运行各节，或Source整文件。
# 1. 设置与共享数据 ----
options(stringsAsFactors = FALSE)
library(multcomp)
alpha <- .05
SAVE_OUTPUTS <- FALSE  # 改为TRUE保存CSV和图片；默认只显示。
UNIT_SCALE <- 1        # 改为1000练习MPa转kPa。
if (Sys.getenv("SAVE_OUTPUTS") == "1") SAVE_OUTPUTS <- TRUE
if (nzchar(Sys.getenv("UNIT_SCALE"))) {
  UNIT_SCALE <- as.numeric(Sys.getenv("UNIT_SCALE"))
}
unit <- if (UNIT_SCALE == 1000) "kPa" else "MPa"
root <- if (dir.exists("data")) "." else ".."
dat <- read.csv(file.path(root, "data", "scaffold.csv"))
dat$modulus <- dat$modulus * UNIT_SCALE
dat$group <- factor(dat$group, levels = c("A1", "A2", "A3", "A4"))
print(head(dat))  # group为结构名；modulus为每个独立试件的一次测量。
print(str(dat))

# 2. 整体模型与共同MSE ----
fit <- aov(modulus ~ group, data = dat)
anova_table <- summary(fit)[[1]]
df_error <- df.residual(fit)
mse <- deviance(fit) / df_error  # 所有四组的组内SS，除以N-k。
means <- tapply(dat$modulus, dat$group, mean)
n <- table(dat$group)
k <- length(means)
print(means)
summary_table <- data.frame(MSE = mse, df = df_error,
  F = anova_table[1, "F value"], F_p = anova_table[1, "Pr(>F)"])
print(summary_table)
protected_gate <- summary_table$F_p < alpha
cat("整体F通过保护性LSD的前置关：", protected_gate, "\n")
cat("下列LSD表仍列普通P，不能称作全部组对的强FWER调整。\n")

# 3. LSD计算形式：共同MSE的普通双侧t比较 ----
# combn把四个组名组合成六个不同组对，方向为前组减后组。
pairs <- combn(names(means), 2)
lsd_rows <- list()
for (row in seq_len(ncol(pairs))) {
  g1 <- pairs[1, row]
  g2 <- pairs[2, row]
  difference <- unname(means[g1] - means[g2])
  se <- sqrt(mse * (1 / n[g1] + 1 / n[g2]))
  t_value <- difference / se
  # pt负尾乘2给出双侧P；qt(.975,df)给出双侧5%临界值。
  p <- 2 * pt(-abs(t_value), df_error)
  half_width <- qt(1 - alpha / 2, df_error) * se
  lsd_rows[[row]] <- data.frame(group1 = g1, group2 = g2,
    difference = difference, SE = unname(se), t = unname(t_value),
    df = df_error, p = unname(p), CI_low = difference - half_width,
    CI_high = difference + half_width, reject = unname(p < alpha))
}
lsd <- do.call(rbind, lsd_rows)
rownames(lsd) <- NULL
print(lsd, digits = 7)
cat("普通P小于alpha的组对数：", sum(lsd$reject), "；差值/区间单位：", unit, "\n")

# 4. Tukey：所有六对的同时推断 ----
# TukeyHSD返回diff/lwr/upr/p adj；R方向为后组减前组，保留原行名。
tukey_result <- TukeyHSD(fit, "group", conf.level = 1 - alpha)$group
tukey <- data.frame(comparison = rownames(tukey_result),
  difference = tukey_result[, "diff"], CI_low = tukey_result[, "lwr"],
  CI_high = tukey_result[, "upr"], adjusted_p = tukey_result[, "p adj"])
tukey$reject <- tukey$adjusted_p < alpha
rownames(tukey) <- NULL
q_critical <- qtukey(1 - alpha, k, df_error)
cat("Tukey q临界值：", q_critical, "；同时推断拒绝数：", sum(tukey$reject), "\n")
print(tukey, digits = 7)

# 5. Dunnett：只比较A2/A3/A4减A1 ----
# factor第一个水平A1为共同参照。glht读取完整模型的共同MSE。
set.seed(20261130)
dun_model <- glht(fit, linfct = mcp(group = "Dunnett"))
# single-step按三项的联合t分布调整；积分精度显式提高。
dun_summary <- summary(dun_model,
  test = adjusted("single-step", maxpts = 1000000, abseps = 1e-5))
dun_ci <- confint(dun_model, level = 1 - alpha,
  calpha = adjusted_calpha(maxpts = 1000000, abseps = 1e-5))$confint
dunnett <- data.frame(comparison = names(dun_summary$test$coefficients),
  difference = unname(dun_summary$test$coefficients),
  adjusted_p = unname(dun_summary$test$pvalues),
  CI_low = dun_ci[, "lwr"], CI_high = dun_ci[, "upr"])
dunnett$reject <- dunnett$adjusted_p < alpha
rownames(dunnett) <- NULL
print(dunnett, digits = 7)
cat("共同对照族拒绝数：", sum(dunnett$reject), "；只保护预定三项。\n")
# 共同A1参照使三个差值相关；这里用共同MSE和样本量打印可核对的尺度。
v_control <- mse / as.numeric(n[1])
v_difference <- mse * (1 / as.numeric(n[-1]) + 1 / as.numeric(n[1]))
cov_shared <- mse / as.numeric(n[1])
cor_shared <- cov_shared / sqrt(v_difference[1] * v_difference[2])
dunnett_covariance <- data.frame(v_MSE_over_n_A1 = v_control,
  var_each_treatment_minus_A1 = v_difference[1], cov_A2A1_A3A1 = cov_shared,
  cor_A2A1_A3A1 = cor_shared)
print(dunnett_covariance, digits = 7)
cat("说明：共同参照联合t推断使用该相关结构；MSE/n仅演示共同尺度。\n")

# 6. SNK：等组量的教学自编实现 ----
stopifnot(length(unique(as.numeric(n))) == 1)
# order得到样本均数降序对应的组名；i/j是排序位置，r包含两端。
ordered_names <- names(sort(means, decreasing = TRUE))
snk_rows <- list()
row <- 0
for (i in seq_len(k - 1)) {
  for (j in (i + 1):k) {
    row <- row + 1
    span <- j - i + 1
    difference <- means[ordered_names[i]] - means[ordered_names[j]]
    q <- unname(difference / sqrt(mse / as.numeric(n[1])))
    snk_rows[[row]] <- data.frame(high = ordered_names[i], low = ordered_names[j],
      i = i, j = j, r = span, q = q, qcrit = qtukey(1 - alpha, span, df_error),
      raw_range_P = ptukey(q, span, df_error, lower.tail = FALSE))
  }
}
snk <- do.call(rbind, snk_rows)
# 一个大区间未拒绝，内部组对便不能拒绝。
# 因此每对必须让所有包含它的大区间都通过：取这些局部P的最大值。
for (row in seq_len(nrow(snk))) {
  contains_pair <- snk$i <= snk$i[row] & snk$j >= snk$j[row]
  snk$step_P[row] <- max(snk$raw_range_P[contains_pair])
}
snk$reject <- snk$step_P < alpha
# A1-A3的所有包含区间：按排序位置动态筛选，便于逐步规则核对。
a1_pos <- match("A1", ordered_names)
a3_pos <- match("A3", ordered_names)
contains_a1_a3 <- snk$i <= min(a1_pos, a3_pos) & snk$j >= max(a1_pos, a3_pos)
snk_a1_a3_ranges <- snk[contains_a1_a3, c("high", "low", "r", "raw_range_P")]
print(snk_a1_a3_ranges, digits = 7)
cat("A1-A3包含区间raw_range_P最大值step_P：", max(snk_a1_a3_ranges$raw_range_P), "\n")
snk <- snk[order(-snk$r, snk$i), ]
print(snk, digits = 7)
cat("SNK逐步拒绝数：", sum(snk$reject), "；step_P不是强FWER调整P。\n")

# 7. 实际观测点与组均数；保存开关 ----
# stripchart画出每个观测点；固定偏移便于重现，钻石为组均数。
plot_points <- function() {
  plot(NA, xlim = c(.5, k + .5), ylim = range(dat$modulus),
    xaxt = "n", xlab = "Scaffold group", ylab = paste("Compressive modulus", unit))
  axis(1, at = seq_len(k), labels = names(means))
  for (i in seq_len(k)) {
    values <- dat$modulus[dat$group == names(means)[i]]
    points(i + seq(-.1, .1, length.out = length(values)), values,
      pch = 19, col = "#4D8EC5")
  }
  points(seq_len(k), means, pch = 18, col = "#0000F5", cex = 1.5)
  title("Individual specimens and group means")
}
plot_points()
if (SAVE_OUTPUTS) {
  out <- file.path(root, "运行结果", "R")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  tables <- list(summary = summary_table, LSD = lsd, Tukey = tukey,
    Dunnett = dunnett, SNK = snk)
  for (name in names(tables)) {
    write.csv(tables[[name]], file.path(out, paste0(name, ".csv")), row.names = FALSE)
  }
  png(file.path(out, "means_points.png"), width = 1400, height = 950, res = 180)
  plot_points()
  dev.off()
}
