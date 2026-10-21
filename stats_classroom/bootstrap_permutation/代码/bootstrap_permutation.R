# 明哥的微生物世界 · 统计课堂
# Bootstrap与置换检验：反复抽样，究竟在模拟什么？
#
# 运行方式：打开上一级的“统计课堂-Bootstrap与置换检验.Rproj”，Source本文件。
# 只用基础R；所有脚本与文字文件均为UTF-8。固定索引让三种语言逐值核对。
# 主要结果对象直接写在代码中，逐段运行时可在RStudio查看；文字解释放在相邻注释。

# %% 0. 项目路径、输入和输出
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
source_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (length(file_arg)) {
  project_dir <- dirname(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
} else if (!is.null(source_file)) {
  project_dir <- dirname(dirname(normalizePath(source_file)))
} else {
  project_dir <- getwd()
}
if (!file.exists(file.path(project_dir, "data", "rats.csv"))) {
  stop("请从项目目录Source完整脚本，或将project_dir设为配套项目目录。")
}

out_dir <- file.path(project_dir, "运行结果", "R")
figure_dir <- file.path(out_dir, "图")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

rats <- read.csv(file.path(project_dir, "data", "rats.csv"),
                 fileEncoding = "UTF-8", stringsAsFactors = FALSE)
index_A <- as.matrix(read.csv(file.path(project_dir, "data", "bootstrap_indices_A.csv"),
                              header = FALSE))
index_B <- as.matrix(read.csv(file.path(project_dir, "data", "bootstrap_indices_B.csv"),
                              header = FALSE))
A <- rats$gain_g[rats$group == "A"]
B <- rats$gain_g[rats$group == "B"]
n_A <- length(A)
n_B <- length(B)
B_rep <- nrow(index_A)

data_summary <- data.frame(
  group = c("A", "B"),
  n = c(n_A, n_B),
  mean_g = c(mean(A), mean(B)),
  sd_g = c(sd(A), sd(B))
)
data_summary
write.csv(data_summary, file.path(out_dir, "data_summary.csv"), row.names = FALSE)

stopifnot(nrow(index_B) == B_rep,
          ncol(index_A) == n_A,
          ncol(index_B) == n_B,
          all(index_A >= 1 & index_A <= n_A),
          all(index_B >= 1 & index_B <= n_B))

# %% 0.1 实际分析先调用专用包；后面的循环用于理解原理
# 首次使用：install.packages("boot.pval")
# 先预测A−B方向；运行后同时读效应量、CI与P，而不只看P。
if (!requireNamespace("boot.pval", quietly = TRUE)) {
  stop("请先运行 install.packages('boot.pval')，再重新运行本脚本。")
}
set.seed(20261023)
fit_package <- boot.pval::boot_t_test(A, B, var.equal = FALSE,
                                     R = 10000, type = "stud")
fit_package
package_result <- data.frame(
  difference_g = unname(fit_package$estimate[1] - fit_package$estimate[2]),
  T = unname(fit_package$statistic),
  CI_low_g = unname(fit_package$conf.int[1]),
  CI_high_g = unname(fit_package$conf.int[2]),
  P = fit_package$p.value,
  package_version = as.character(utils::packageVersion("boot.pval"))
)
package_result
# 0.7.0版本：差19.1667g；CI约[-2.319064,40.927598]；P=.0730。
# 输出R=10000是重抽样次数，不是t分布自由度。
# 包使用自己的抽样，后文固定索引用于三语言逐轮核对，数值略有差异。
write.csv(package_result, file.path(out_dir, "boot_pval_package_result.csv"), row.names = FALSE)

# %% 1. 观察到的均数差和学生化统计量
# A−B是效应量；SE使用两组各自的样本方差和n，不把两组混成一组。
delta_obs <- mean(A) - mean(B)
SE_obs <- sqrt(var(A) / n_A + var(B) / n_B)
T_obs <- delta_obs / SE_obs

observed_result <- data.frame(
  delta_g = delta_obs,
  SE_g = SE_obs,
  T = T_obs
)
observed_result
write.csv(observed_result, file.path(out_dir, "observed_result.csv"), row.names = FALSE)

# %% 2. 普通按组Bootstrap：估计均数差的波动
# 每行索引是一轮；组内有放回，样本量固定为12和7。
# 这一步以原始经验分布为中心，回答“均数差有多不稳定”，不是H0的零分布。
delta_star <- numeric(B_rep)
for (i in seq_len(B_rep)) {
  A_star <- A[index_A[i, ]]
  B_star <- B[index_B[i, ]]
  delta_star[i] <- mean(A_star) - mean(B_star)
}

bootstrap_SE <- sd(delta_star)
CI_percentile <- quantile(delta_star, probs = c(.025, .975), type = 7)
bootstrap_estimate <- data.frame(
  bootstrap_SE_g = bootstrap_SE,
  percentile_low_g = unname(CI_percentile[1]),
  percentile_high_g = unname(CI_percentile[2])
)
bootstrap_estimate
write.csv(bootstrap_estimate, file.path(out_dir, "bootstrap_estimate.csv"), row.names = FALSE)

# %% 3. 检验的H0参照：先中心化，再对每轮差值学生化
# 中心化 A0=A−mean(A)、B0=B−mean(B)，使两组经验均数都为0；这是构造H0: mu_A=mu_B
# 的参照模型，不是说原始数据被“改成了真实数据”。方差和组内形状仍来自各组经验分布。
# 学生化用每一轮自己的SE，把差值放到“几个估计标准误”的尺度上。
A0 <- A - mean(A)
B0 <- B - mean(B)
Tfun <- function(a, b) {
  d <- mean(a) - mean(b)
  se <- sqrt(var(a) / length(a) + var(b) / length(b))
  d / se
}

T_null <- numeric(B_rep)
SE_star <- numeric(B_rep)
for (i in seq_len(B_rep)) {
  A_star0 <- A0[index_A[i, ]]
  B_star0 <- B0[index_B[i, ]]
  SE_star[i] <- sqrt(var(A_star0) / n_A + var(B_star0) / n_B)
  T_null[i] <- Tfun(A_star0, B_star0)
}

# 这条等价式可帮助检查“中心化”也可以写成 delta_star−delta_obs，再除以本轮SE。
stopifnot(max(abs(T_null - (delta_star - delta_obs) / SE_star)) < 1e-10,
          abs(Tfun(A, B) - T_obs) < 1e-12)

# 双侧等尾规则：较小尾部概率乘2；本例左尾9636、右尾364。
left_count <- sum(T_null <= T_obs)
right_count <- sum(T_null >= T_obs)
p_boot <- min(1, 2 * min(left_count, right_count) / B_rep)
q_T <- quantile(T_null, c(.025, .975), type = 7)
CI_t <- delta_obs - rev(q_T) * SE_obs

bootstrap_test <- data.frame(
  left = left_count,
  right = right_count,
  B = B_rep,
  P_two_sided_equal_tail = p_boot,
  bootstrap_t_low_g = unname(CI_t[1]),
  bootstrap_t_high_g = unname(CI_t[2])
)
bootstrap_test
write.csv(bootstrap_test, file.path(out_dir, "bootstrap_test.csv"), row.names = FALSE)

# %% 4. 置换检验：混合19个观测，再按12/7重分
# H0若是“组标签可交换”，固定19个观测后，枚举C(19,12)=50388种分组。
# 这与组内Bootstrap不同：这里不从A、B各自经验分布重复抽取，原有19个值每轮各用一次。
all_values <- c(A, B)
all_index <- seq_along(all_values)
perm_groups <- t(combn(all_index, n_A))
perm_rep <- nrow(perm_groups)
perm_delta <- numeric(perm_rep)
for (i in seq_len(perm_rep)) {
  perm_A_index <- perm_groups[i, ]
  perm_B_index <- all_index[!all_index %in% perm_A_index]
  perm_delta[i] <- mean(all_values[perm_A_index]) - mean(all_values[perm_B_index])
}
perm_extreme <- sum(abs(perm_delta) >= abs(delta_obs))
perm_p <- perm_extreme / perm_rep
permutation_result <- data.frame(
  total_partitions = perm_rep,
  extreme_abs_delta = perm_extreme,
  P_two_sided_exact = perm_p
)
permutation_result
write.csv(permutation_result, file.path(out_dir, "permutation_result.csv"), row.names = FALSE)

# %% 5. 率的边界例：普通百分位Bootstrap与CP区间的差别
# 对观察到x/n，理想无限Bootstrap分布是K*~Binomial(n, x/n)。
# 这里用离散分布的“最小k使CDF>=q”定义分位点，避免用有限B模拟掩盖边界。
binom_quantile_index <- function(n, p, prob) {
  k <- 0:n
  cdf <- pbinom(k, size = n, prob = p)
  k[which(cdf >= prob - 1e-14)[1]]
}
percentile_rate_ci <- function(x, n, alpha = .05) {
  p_hat <- x / n
  lower_k <- binom_quantile_index(n, p_hat, alpha / 2)
  upper_k <- binom_quantile_index(n, p_hat, 1 - alpha / 2)
  c(lower_k / n, upper_k / n)
}

ci_9_of_10 <- percentile_rate_ci(9, 10)
ci_0_of_10 <- percentile_rate_ci(0, 10)
cp_0_upper <- 1 - (.025)^(1 / 10)
rate_result <- data.frame(
  case = c("9/10", "0/10"),
  percentile_low = c(ci_9_of_10[1], ci_0_of_10[1]),
  percentile_high = c(ci_9_of_10[2], ci_0_of_10[2]),
  CP_low = c(NA, 0),
  CP_high = c(NA, cp_0_upper)
)
rate_result
write.csv(rate_result, file.path(out_dir, "rate_boundary_results.csv"), row.names = FALSE)

# 在给定教学模型p=.9、n=10下，枚举所有观察x；每个x使用上面的无限Bootstrap
# 百分位算法，最后按Binomial(10,.9)的真实概率加权。这是精确枚举，不是模拟。
true_p <- .9
x_values <- 0:10
coverage_rows <- lapply(x_values, function(x) {
  interval <- percentile_rate_ci(x, 10)
  data.frame(
    x = x,
    interval_low = interval[1],
    interval_high = interval[2],
    true_p_covered = interval[1] <= true_p && true_p <= interval[2],
    probability_under_true_p = dbinom(x, 10, true_p)
  )
})
coverage_table <- do.call(rbind, coverage_rows)
coverage_table
coverage_exact <- sum(coverage_table$probability_under_true_p[coverage_table$true_p_covered])
coverage_summary <- data.frame(
  true_p = true_p,
  all_success_probability = true_p^10,
  exact_enumerated_coverage = coverage_exact
)
coverage_summary
write.csv(coverage_table, file.path(out_dir, "rate_coverage_enumeration.csv"), row.names = FALSE)
write.csv(coverage_summary, file.path(out_dir, "rate_coverage_summary.csv"), row.names = FALSE)

# %% 6. 两张解释图
# 图中使用英文标签，避免要求学生额外下载中文字体；正文负责中文解释。
png(file.path(figure_dir, "01_两种重抽样与统计量.png"), width = 1400, height = 700, res = 150)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
hist(delta_star, breaks = 40, col = "#9CCBEA", border = "white",
     main = "Groupwise bootstrap: estimate", xlab = "A mean - B mean (g)")
abline(v = delta_obs, col = "#145A86", lwd = 2, lty = 2)
abline(v = CI_percentile, col = "#E67E22", lwd = 2)
legend("topright", legend = c("observed", "percentile limits"),
       col = c("#145A86", "#E67E22"), lty = c(2, 1), lwd = 2, bty = "n", cex = .8)
hist(T_null, breaks = 40, col = "#B9D9C6", border = "white",
     main = "Centered studentized null", xlab = "T* under H0")
abline(v = T_obs, col = "#E67E22", lwd = 2)
dev.off()

png(file.path(figure_dir, "02_置换分布与观察差值.png"), width = 900, height = 650, res = 150)
par(mar = c(4, 4, 3, 1))
hist(perm_delta, breaks = 40, col = "#C9D7F0", border = "white",
     main = "Permutation reallocation of 19 observations",
     xlab = "A mean - B mean (g)")
abline(v = c(-abs(delta_obs), abs(delta_obs)), col = "#E67E22", lwd = 2, lty = 2)
abline(v = delta_obs, col = "#145A86", lwd = 2)
legend("topright", legend = c("observed", "two-sided cutoffs"),
       col = c("#145A86", "#E67E22"), lty = c(1, 2), lwd = 2, bty = "n", cex = .8)
dev.off()

# %% 7. 最后核对并写入简短记录
stopifnot(abs(delta_obs - 19.1666666666667) < 1e-10,
          abs(T_obs - 1.93193081703987) < 1e-10,
          left_count == 9636,
          right_count == 364,
          abs(p_boot - .0728) < 1e-12,
          abs(CI_t[1] - (-2.28522435160976)) < 1e-8,
          abs(CI_t[2] - 40.92199759211) < 1e-8,
          abs(CI_percentile[1] - 1.04761904761904) < 1e-8,
          abs(CI_percentile[2] - 37.3693452380952) < 1e-8,
          abs(coverage_exact - .6496866225) < 1e-12)
run_notes <- c(
  "R脚本完成：固定索引、中心化学生化Bootstrap、置换重分与率边界枚举。",
  sprintf("delta=%.10f, T=%.10f, Bootstrap P=%.4f", delta_obs, T_obs, p_boot),
  sprintf("Bootstrap-t CI=[%.10f, %.10f]", CI_t[1], CI_t[2]),
  sprintf("permutation exact partitions=%d, P=%.10f", perm_rep, perm_p),
  sprintf("rate coverage exact enumeration=%.10f", coverage_exact),
  capture.output(sessionInfo())
)
writeLines(run_notes, file.path(out_dir, "运行记录.txt"), useBytes = TRUE)
