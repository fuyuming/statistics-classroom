# 统计课堂07：正态、卡方、t、F。UTF-8；只需要基础 R。
# 命令行：Rscript 代码/normal_chisq_t_f.R
# RStudio：打开本脚本，Source 完整运行。不要只运行中间片段。
# 如逐行运行，请先将工作目录设为 article07_normal_chisq_t_f。
args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep('^--file=', args, value = TRUE)
script_path <- if (length(script_arg)) sub('^--file=', '', script_arg[1]) else NULL
if (is.null(script_path)) {
  for (fr in sys.frames()) if (!is.null(fr$ofile)) script_path <- fr$ofile
}
root <- if (!is.null(script_path)) dirname(dirname(normalizePath(script_path))) else getwd()
stopifnot(file.exists(file.path(root, '数据', 'teaching_data.csv')))
out <- file.path(root, '运行结果'); figs <- file.path(root, '文章配图')
dir.create(out, showWarnings = FALSE); dir.create(figs, showWarnings = FALSE)

# 一、精确分位点与理论错误率；不是模拟结果。
dfs <- c(1, 5, 10, 30, 100)
critical <- data.frame(df = dfs, n_one_sample = dfs + 1,
                      t_975 = qt(.975, dfs),
                      t_975_squared = qt(.975, dfs)^2,
                      f_95 = qf(.95, 1, dfs),
                      rejection_at_1_96 = 2 * pt(1.96, dfs, lower.tail = FALSE))
stopifnot(max(abs(critical$t_975_squared - critical$f_95)) < 1e-8)
write.csv(critical, file.path(out, 'critical_values.csv'), row.names = FALSE)

# 二、导出配图所用的全部曲线坐标，供其他软件复画。
x <- seq(-5, 5, length.out = 1001)
curves <- data.frame(x = x, normal = dnorm(x), t_df1 = dt(x, 1),
                     t_df5 = dt(x, 5), t_df30 = dt(x, 30))
write.csv(curves, file.path(out, 't_normal_curves.csv'), row.names = FALSE)
png(file.path(figs, '02-t分布与标准正态比较.png'), width = 1600, height = 1000, res = 160)
par(mar = c(4.5, 4.5, 3, 1), las = 1)
matplot(x, curves[, -1], type = 'l', lty = c(2, 1, 1, 1),
        col = c('#222222', '#178C88', '#4144A5', '#E36854'), lwd = 2,
        xlab = 'Statistic value', ylab = 'Density',
        main = 'Student t and standard normal distributions')
legend('topright', c('N(0, 1)', 't: df = 1', 't: df = 5', 't: df = 30'),
       col = c('#222222', '#178C88', '#4144A5', '#E36854'),
       lty = c(2, 1, 1, 1), lwd = 2, bty = 'n', cex = .8)
dev.off()

# 二补充、卡方曲线：分面避免 df=1 在零附近的高密度压扁其他曲线。
chi_x <- seq(.001, 55, length.out = 3001)
chi_df <- c(1, 2, 3, 5, 10, 30)
chi_curves <- data.frame(x = chi_x)
for (k in chi_df) chi_curves[[paste0('df_', k)]] <- dchisq(chi_x, k)
write.csv(chi_curves, file.path(out, 'chisq_curves.csv'), row.names = FALSE)
png(file.path(figs, '03-卡方分布是一族曲线.png'), width = 1600, height = 1500, res = 160)
par(mfrow = c(2, 1), mar = c(4, 4.5, 3, 1), las = 1)
plot(chi_x, dchisq(chi_x, 1), type = 'l', col = '#178C88', lwd = 2,
     xlim = c(0, 8), ylim = c(0, 1.5), xlab = 'Chi-square value', ylab = 'Density',
     main = 'Chi-square: df = 1 and 2')
lines(chi_x, dchisq(chi_x, 2), col = '#4144A5', lwd = 2, lty = 2)
legend('topright', c('df = 1', 'df = 2'), col = c('#178C88', '#4144A5'), lty = c(1, 2), lwd = 2, bty = 'n')
mtext('df = 1: density tends to infinity as x approaches 0; upper part clipped', side = 3, cex = .7, line = .2)
cols <- c('#178C88', '#4144A5', '#E36854', '#333333')
matplot(chi_x, chi_curves[, c('df_3', 'df_5', 'df_10', 'df_30')], type = 'l',
        col = cols, lty = 1:4, lwd = 2, xlab = 'Chi-square value', ylab = 'Density',
        main = 'Chi-square: increasing degrees of freedom')
legend('topright', paste('df =', c(3, 5, 10, 30)), col = cols, lty = 1:4, lwd = 2, bty = 'n')
dev.off()

# F 曲线：分别固定分母和分子自由度，观察另一个自由度的影响。
f_x <- seq(.001, 6, length.out = 2001)
pairs <- data.frame(d1 = c(3, 5, 10, 5, 5), d2 = c(10, 10, 10, 5, 30))
f_curves <- data.frame(x = f_x)
for (j in seq_len(nrow(pairs))) {
  name <- paste0('F_', pairs$d1[j], '_', pairs$d2[j])
  f_curves[[name]] <- df(f_x, pairs$d1[j], pairs$d2[j])
}
write.csv(f_curves, file.path(out, 'f_curves.csv'), row.names = FALSE)
png(file.path(figs, '04-F分布的两个自由度.png'), width = 1600, height = 1500, res = 160)
par(mfrow = c(2, 1), mar = c(4, 4.5, 3, 1), las = 1)
for (panel in 1:2) {
  keys <- if (panel == 1) c('F_3_10', 'F_5_10', 'F_10_10') else c('F_5_5', 'F_5_10', 'F_5_30')
  labs <- if (panel == 1) c('F(3, 10)', 'F(5, 10)', 'F(10, 10)') else c('F(5, 5)', 'F(5, 10)', 'F(5, 30)')
  matplot(f_x, f_curves[, keys], type = 'l', col = cols[1:3], lty = 1:3, lwd = 2,
          xlab = 'F value', ylab = 'Density',
          main = if (panel == 1) 'Fixed denominator df = 10' else 'Fixed numerator df = 5')
  abline(v = 1, col = '#999999', lty = 3)
  legend('topright', labs, col = cols[1:3], lty = 1:3, lwd = 2, bty = 'n')
}
dev.off()

# 三、人工教学数据：演示 t²=F 与回归等价，不代表真实实验。
d <- read.csv(file.path(root, '数据', 'teaching_data.csv'))
d$group <- factor(d$group, levels = c('A', 'B'))
two_t <- t.test(response ~ group, data = d, var.equal = TRUE)
fit <- lm(response ~ group, data = d)
a <- anova(fit)
reg_t <- coef(summary(fit))['groupB', 't value']
comparison <- data.frame(method = c('pooled_t_squared', 'one_way_ANOVA_F', 'regression_group_t_squared'),
                         statistic = c(unname(two_t$statistic)^2, a$`F value`[1], reg_t^2),
                         p_value = c(two_t$p.value, a$`Pr(>F)`[1], coef(summary(fit))['groupB', 'Pr(>|t|)']))
stopifnot(diff(range(comparison$statistic)) < 1e-10,
          diff(range(comparison$p_value)) < 1e-10)
write.csv(comparison, file.path(out, 't_anova_regression_equivalence.csv'), row.names = FALSE)

# 四、一元连续自变量回归：斜率 t² 与模型整体 F。
reg <- lm(response ~ x, data = d)
rs <- summary(reg)
regression <- data.frame(slope_t = coef(rs)['x', 't value'],
                         slope_t_squared = coef(rs)['x', 't value']^2,
                         overall_F = unname(rs$fstatistic['value']),
                         residual_df = df.residual(reg))
stopifnot(abs(regression$slope_t_squared - regression$overall_F) < 1e-10)
write.csv(regression, file.path(out, 'simple_regression_equivalence.csv'), row.names = FALSE)

# 五、平方/比值构造的随机演示。有限模拟不用于证明定理。
set.seed(20260927)
B <- 50000L; nu <- 10L
z <- rnorm(B)
u <- rowSums(matrix(rnorm(B * nu), nrow = B)^2)
t_constructed <- z / sqrt(u / nu)
f_constructed <- z^2 / (u / nu)
stopifnot(max(abs(t_constructed^2 - f_constructed)) < 1e-10)
q_demo <- data.frame(probability = c(.025, .5, .975),
                    empirical_t = as.numeric(quantile(t_constructed, c(.025, .5, .975))),
                    theoretical_t = qt(c(.025, .5, .975), nu))
write.csv(q_demo, file.path(out, 'construction_simulation_quantiles.csv'), row.names = FALSE)
# 保存全部模拟变量，便于逐行检查构造；z 与 u 由独立随机抽样生成。
write.csv(data.frame(z = z, u_chisq10 = u, t10 = t_constructed, f1_10 = f_constructed),
          gzfile(file.path(out, 'construction_draws.csv.gz')), row.names = FALSE)

report <- capture.output({
  cat('统计课堂07：逐项核对\n\n1. 理论分位点和实际错误率\n'); print(critical, digits = 9)
  cat('\n2. 两独立组：合并方差 t、方差分析、组别回归\n'); print(comparison, digits = 9)
  cat('\n3. 连续自变量一元回归\n'); print(regression, digits = 9)
  cat('\n4. 独立正态与卡方构造 t；模拟分位点允许抽样误差\n'); print(q_demo)
  cat('\n教学数据是人为构造，不用本例的显著性推断真实生物学结论。\n')
  cat('所有代数恒等式核对已通过。\n\n'); print(sessionInfo())
})
writeLines(report, file.path(out, '逐项计算结果.txt'), useBytes = TRUE)
cat(paste(report, collapse = '\n'), '\n')
