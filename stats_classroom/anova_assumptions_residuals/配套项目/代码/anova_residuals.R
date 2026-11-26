# 统计课堂43：方差分析的正态性、方差齐性：指什么，一定要检验吗？
# 输入是可复现的模拟教学数据，不是实测数据；运行目录为项目根目录。

# ---- 参数与读取
options(stringsAsFactors = FALSE)
project_root <- normalizePath(getwd(), mustWork = FALSE)
if (!dir.exists(file.path(project_root, "数据")) && basename(project_root) == "代码") {
  project_root <- dirname(project_root)
}
data_path <- file.path(project_root, "数据", "residual_cases.csv")
source(file.path(project_root, "代码", "residplot.R"))
dat <- read.csv(data_path, check.names = FALSE)
cat("字段：", paste(names(dat), collapse = ", "), "\n")
print(head(dat, 6))
SAVE_RESULTS <- identical(Sys.getenv("SAVE_RESULTS", "false"), "true")

# ---- 场景A的描述统计与拟合
make_results <- function(d, label) {
  d$true_error <- d$value - d$mu
  fit <- lm(value ~ group, data = d)
  d$residual <- residuals(fit)
  d$fitted <- fitted(fit)
  summary_tab <- aggregate(value ~ group, d, function(x) c(n = length(x), mean = mean(x), sd = sd(x)))
  summary_tab <- data.frame(scenario = label, group = summary_tab$group, n = summary_tab$value[, 1], mean = summary_tab$value[, 2], sd = summary_tab$value[, 3])
  grand <- mean(d$value)
  groups <- split(d$value, d$group)
  means <- sapply(groups, mean)
  ssb <- sum(sapply(seq_along(groups), function(i) length(groups[[i]]) * (means[i] - grand)^2))
  sse <- sum((d$value - ave(d$value, d$group, FUN = mean))^2)
  dfb <- length(groups) - 1
  dfe <- nrow(d) - length(groups)
  msb <- ssb / dfb
  mse <- sse / dfe
  Fvalue <- msb / mse
  pvalue <- pf(Fvalue, dfb, dfe, lower.tail = FALSE)
  a <- data.frame(scenario = label, SSB = ssb, SSE = sse, df_between = dfb, df_error = dfe, MSB = msb, MSE = mse, F = Fvalue, p = pvalue)
  list(data = d, summary = summary_tab, anova = a, fit = fit)
}
stopifnot(identical(unique(dat$scenario), "A"))
res <- list(A = make_results(dat, "A"))
group_summary <- do.call(rbind, lapply(res, `[[`, "summary"))
anova_tab <- do.call(rbind, lapply(res, `[[`, "anova"))
residuals_tab <- do.call(rbind, lapply(res, function(x) x$data[, c("scenario", "group", "replicate", "value", "mu", "z", "true_error", "fitted", "residual")]))
cat("\n主例题A：每组 n、均数、SD\n")
print(group_summary)
cat("\n手算 ANOVA 与 lm+anova 核对；残差按组和应为0\n")
print(anova_tab)
print(aggregate(residual ~ scenario + group, residuals_tab, sum))
print(lapply(res, function(x) anova(x$fit)))

# ---- 场景A原始value的方差齐性检验（中位数中心的Brown–Forsythe修正版）
# H0：各组总体方差相等；P>=.05时仅表述为未拒绝H0，不等同于证明方差齐。
if (!requireNamespace("car", quietly = TRUE)) stop("场景A Levene检验需要car包；请先安装后再运行。")
a <- subset(dat, scenario == "A")
a$group <- factor(a$group, levels = c("A1", "A2", "A3"))
lev <- car::leveneTest(value ~ group, data = a, center = median)
bar <- bartlett.test(value ~ group, data = a)
variance_tests <- data.frame(
  method = c("Levene_median", "Bartlett"),
  statistic = c(unname(lev[1, "F value"]), unname(bar$statistic)),
  df1 = c(unname(lev[1, "Df"]), unname(bar$parameter)),
  df2 = c(unname(lev[2, "Df"]), NA_real_),
  p = c(unname(lev[1, "Pr(>F)"]), unname(bar$p.value))
)
cat("\n场景A原始value方差齐性检验（H0：各组总体方差相等）：\n")
cat(sprintf("Levene（中位数中心，Brown–Forsythe修正版）F(%d,%d) = %.6f, p = %.6f\n", variance_tests$df1[1], variance_tests$df2[1], variance_tests$statistic[1], variance_tests$p[1]))
cat(if (variance_tests$p[1] >= .05) "解释：P>=.05，未拒绝各组总体方差相等的H0；不表示已证明方差齐。\n" else "解释：P<.05，拒绝各组总体方差相等的H0。\n")
cat(sprintf("Bartlett χ²(%d) = %.6f, p = %.6f\n", variance_tests$df1[2], variance_tests$statistic[2], variance_tests$p[2]))
cat(if (variance_tests$p[2] >= .05) "解释：P>=.05，未拒绝各组总体方差相等的H0；不表示已证明方差齐。\n" else "解释：P<.05，拒绝各组总体方差相等的H0。\n")

# ---- 场景A各组原始value的正态性检验
# 每组n=5：Shapiro-Wilk实际运行；Anderson-Darling仅在n>=8时运行。
a_groups <- split(a$value, a$group)
normality_rows <- lapply(names(a_groups), function(g) {
  x <- a_groups[[g]]
  sw <- shapiro.test(x)
  ad <- if (length(x) >= 8L) nortest::ad.test(x) else NULL
  rbind(
    data.frame(scenario = "A", group = g, n = length(x), method = "Shapiro_Wilk",
               statistic = unname(sw$statistic), p = sw$p.value,
               status = if (sw$p.value >= .05) "not_rejected" else "rejected"),
    data.frame(scenario = "A", group = g, n = length(x), method = "Anderson_Darling",
               statistic = if (is.null(ad)) NA_real_ else unname(ad$statistic),
               p = if (is.null(ad)) NA_real_ else ad$p.value,
               status = if (is.null(ad)) "not_run_n_lt_8" else if (ad$p.value >= .05) "not_rejected" else "rejected")
  )
})
normality_tests <- do.call(rbind, normality_rows)
cat("\n场景A各组原始value正态性检验（每组独立样本；不对合并残差作iid正态性P断言）：\n")
print(normality_tests)

# ---- 场景A排序的Q-Q数值表
qq_parts <- lapply(split(residuals_tab, residuals_tab$scenario), function(d) {
  d <- d[order(d$residual), c("scenario", "group", "replicate", "residual")]
  d$rank <- seq_len(nrow(d))
  d$theoretical <- qnorm((d$rank - 0.5) / nrow(d))
  d[, c("scenario", "group", "replicate", "rank", "residual", "theoretical")]
})
qq <- do.call(rbind, qq_parts)
cat("\nQ-Q表：场景A的N=15排序；不作独立同分布Shapiro P断言\n")
print(qq)

# ---- 场景A的三张诊断图（RStudio默认显示，批量验证可SHOW_PLOTS=true）
# 此节先保留组别看尺度，再把15个残差合并按正态理论位置作图。
SHOW_PLOTS <- interactive() || identical(tolower(Sys.getenv("SHOW_PLOTS", "false")), "true")
if (SHOW_PLOTS) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("绘图需要安装ggplot2包。")
  if (!requireNamespace("car", quietly = TRUE)) stop("Q-Q图需要安装car包。")
  a <- subset(residuals_tab, scenario == "A")
  qdat <- a[order(a$residual), ]
  qdat$theory <- qnorm((seq_len(nrow(qdat)) - 0.5) / nrow(qdat))
  yq <- quantile(a$residual, c(.25, .75), type = 7)
  zq <- qnorm(c(.25, .75))
  slope <- unname(diff(yq) / diff(zq))
  intercept <- unname(yq[1] - slope * zq[1])
  group_cols <- c(A1 = "#3155ad", A2 = "#008e91", A3 = "#ae6373")
  print(ggplot2::ggplot(a, ggplot2::aes(fitted, residual, colour = group)) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::scale_colour_manual(values = group_cols) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::labs(x = "Fitted group mean (MPa)", y = "Residual (MPa)") +
    ggplot2::theme_minimal())
  # car::qqPlot实际接收已排序的原始残差向量；其p=ppoints(n)与qq.csv的
  # (rank-0.5)/n一致，line="quartiles"也使用同一type=7四分位数直线。
  car::qqPlot(qdat$residual, distribution = "norm", line = "quartiles",
              envelope = FALSE, id = FALSE, col = unname(group_cols[as.character(qdat$group)]),
              xlab = "Normal theoretical quantile", ylab = "Residual (MPa)",
              main = "Raw residual Q-Q plot")
  legend("topleft", legend = names(group_cols), col = unname(group_cols), pch = 1, bty = "n")
  # 与第42篇同一函数：输入普通lm，输出外部学生化残差及直方图/正态曲线/KDE。
  fit_a <- lm(value ~ group, data = subset(dat, scenario == "A"))
  residplot(fit_a, nbreaks = 5L, draw = TRUE)
}

# ---- 按需保存统一结果
if (SAVE_RESULTS) {
  out <- file.path(project_root, "结果", "R")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  write.csv(group_summary, file.path(out, "group_summary.csv"), row.names = FALSE)
  write.csv(residuals_tab, file.path(out, "residuals.csv"), row.names = FALSE)
  write.csv(anova_tab, file.path(out, "anova.csv"), row.names = FALSE)
  write.csv(qq, file.path(out, "qq.csv"), row.names = FALSE)
  write.csv(variance_tests, file.path(out, "variance_tests.csv"), row.names = FALSE, na = "")
  write.csv(normality_tests, file.path(out, "normality_tests.csv"), row.names = FALSE, na = "")
}
