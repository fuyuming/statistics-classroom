# 单样本率：从样本重抽、在H0下模拟、只重排数据有什么不同？
# 假设独立抽取180人，117人达标。p为目标人群的总体达标率。
# 原题：H0:p=0.45，H1:p>0.45，alpha=0.05。
# 另一个估计问题：用双侧95%百分位区间描述p的估计不确定性。
# 117个1与63个0按汇总数展开；它们不是带原始顺序的逐人记录。
# 用RStudio打开项目Rproj，先逐行运行第一节，再依次运行后续各节。

# 1. 读取题目；事件数除以总人数得到样本率 ----
project_root <- getwd()
script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- script_args[grepl("^--file=", script_args)]
if (length(file_arg) > 0) {
  script_path <- sub("^--file=", "", file_arg[1])
  project_root <- normalizePath(file.path(dirname(script_path), ".."))
}
data_path <- file.path(project_root, "data", "cases.csv")
if (!file.exists(data_path)) {
  stop("请先打开配套项目根目录中的Rproj，再运行本节。")
}
cases <- read.csv(data_path, stringsAsFactors = FALSE)
print(cases)
main_case <- cases[cases$case == "main", ]
n <- main_case$n
x <- main_case$x
p_hat <- x / n
p0 <- main_case$p0
stopifnot(abs(p_hat - main_case$p_hat) < 1e-12)
B <- 10000
seed <- 20261110
result_dir <- file.path(project_root, "results")
dir.create(result_dir, showWarnings = FALSE)

# 2. 非参数Bootstrap：每次从当前样本有放回抽n人 ----
# 1代表达标，0代表未达标；mean(y)就是达标人数/n。
y <- c(rep(1, x), rep(0, n - x))
print(head(y))
print(table(y))
set.seed(seed)
bootstrap_rates <- numeric(B)
for (b in seq_len(B)) {
  sampled_y <- sample(y, size = n, replace = TRUE)
  bootstrap_rates[b] <- mean(sampled_y)
}
# 0.025与0.975各在一端留下2.5%，构成双侧95%百分位区间。
bootstrap_ci <- quantile(bootstrap_rates, probs = c(.025, .975), names = FALSE)
print(bootstrap_ci)
# 给定当前样本，每次抽到1的概率为p_hat，所以计数条件分布为Bin(n,p_hat)。
# qbinom返回离散计数的分位数，再除以n才是率的分位数。
conditional_counts <- qbinom(c(.025, .975), size = n, prob = p_hat)
print(conditional_counts)
print(conditional_counts / n)

# 3. 原假设模型模拟：从Bin(n,p0)生成计数 ----
# 按H1方向定义更极端结果（仍在H0下计算）：本题是达标人数至少x人。
set.seed(seed)
null_counts <- rbinom(B, size = n, prob = p0)
null_hits <- sum(null_counts >= x)
null_hit_rate <- null_hits / B
print(c(hits = null_hits, repetitions = B, hit_rate = null_hit_rate))
# 0次命中只说明本次有限模拟没遇到该尾部，不能把精确P写成0。

# 4. 精确二项检验：直接计算整条右尾 ----
# greater保留原题的单侧备择；函数同时给相应的单侧95%率区间。
exact_test <- binom.test(x, n, p = p0, alternative = "greater", conf.level = .95)
print(exact_test)
# pbinom的lower.tail=FALSE计算P(X>x-1)，即P(X>=x)。
manual_p <- pbinom(x - 1, size = n, prob = p0, lower.tail = FALSE)
stopifnot(abs(manual_p - exact_test$p.value) < 1e-15)
# 这里另求双侧95%CP区间，便于与前面的双侧Bootstrap区间比较。
# 这不是上一个单侧检验返回的区间；调用的alternative不同。
exact_two_sided <- binom.test(x, n, alternative = "two.sided", conf.level = .95)
cp_ci <- unname(exact_two_sided$conf.int)
print(cp_ci)

# 5. 只重排一个向量：成功总数保持不变 ----
# replace=FALSE只改变顺序。重复100次，样本率应始终为p_hat。
set.seed(seed)
reordered_rates <- numeric(100)
for (b in seq_len(100)) {
  reordered_y <- sample(y, size = n, replace = FALSE)
  stopifnot(sum(reordered_y) == x)
  reordered_rates[b] <- sum(reordered_y) / n
}
print(unique(reordered_rates))
stopifnot(all(reordered_rates == p_hat))

# 6. 边界例：20人均无事件，及20人均有事件 ----
# 从全0向量重抽仍是全0；真实执行重抽来核对这个结论。
set.seed(seed)
y_zero <- rep(0, 20)
y_all <- rep(1, 20)
bootstrap_zero <- numeric(B)
bootstrap_all <- numeric(B)
for (b in seq_len(B)) {
  bootstrap_zero[b] <- mean(sample(y_zero, size = 20, replace = TRUE))
  bootstrap_all[b] <- mean(sample(y_all, size = 20, replace = TRUE))
}
zero_ci <- quantile(bootstrap_zero, c(.025, .975), names = FALSE)
all_ci <- quantile(bootstrap_all, c(.025, .975), names = FALSE)
cp_zero <- unname(binom.test(0, 20, conf.level = .95)$conf.int)
cp_all <- unname(binom.test(20, 20, conf.level = .95)$conf.int)
# 双侧95%CP上限可用Beta(1,20)的97.5%分位数核对。
zero_upper_beta <- qbeta(.975, shape1 = 1, shape2 = 20)
stopifnot(abs(cp_zero[2] - zero_upper_beta) < 1e-12)
print(rbind(bootstrap_zero = zero_ci, CP_zero = cp_zero,
            bootstrap_all = all_ci, CP_all = cp_all))

# 7. 保存同名字段，供三种语言核对确定性结果 ----
result <- data.frame(
  language = "R", case = "main", n = n, x = x, p_hat = p_hat, p0 = p0,
  B = B, seed = seed,
  bootstrap_ci_low = bootstrap_ci[1], bootstrap_ci_high = bootstrap_ci[2],
  conditional_count_q025 = conditional_counts[1],
  conditional_count_q975 = conditional_counts[2],
  h0_hits = null_hits, h0_hit_rate = null_hit_rate,
  exact_right_tail_p = exact_test$p.value,
  cp_ci_low = cp_ci[1], cp_ci_high = cp_ci[2],
  boundary_zero_bootstrap_low = zero_ci[1], boundary_zero_bootstrap_high = zero_ci[2],
  boundary_zero_cp_low = cp_zero[1], boundary_zero_cp_high = cp_zero[2],
  boundary_all_bootstrap_low = all_ci[1], boundary_all_bootstrap_high = all_ci[2],
  boundary_all_cp_low = cp_all[1], boundary_all_cp_high = cp_all[2]
)
write.csv(result, file.path(result_dir, "R_results.csv"), row.names = FALSE)
report <- c(
  R.version.string,
  capture.output(print(result)),
  capture.output(print(exact_test)),
  "重排100次的率均为0.65；模拟命中比例不是精确P。",
  "CP为模型下的精确区间；Bootstrap百分位区间不具有同样的精确覆盖保证。"
)
writeLines(report, file.path(result_dir, "R_console.txt"))
cat(paste(report, collapse = "\n"), "\n")

# 8. Bootstrap区间构造与覆盖率精确枚举（v3） ----
# Percentile/Basic这里使用条件二项分布的无限次Bootstrap端点，
# 不是把10000次模拟端点冒充为理想无限内层结果。
z975 <- qnorm(.975)
percentile_counts <- qbinom(c(.025, .975), size = n, prob = p_hat)
percentile_ci_v3 <- percentile_counts / n
basic_ci_v3 <- c(2 * p_hat - percentile_ci_v3[2],
                 2 * p_hat - percentile_ci_v3[1])
wald_half <- z975 * sqrt(p_hat * (1 - p_hat) / n)
wald_ci_v3 <- c(p_hat - wald_half, p_hat + wald_half)
wilson_den <- 1 + z975^2 / n
wilson_center <- (p_hat + z975^2 / (2 * n)) / wilson_den
wilson_half <- z975 / wilson_den * sqrt(p_hat * (1 - p_hat) / n + z975^2 / (4 * n^2))
wilson_ci_v3 <- c(wilson_center - wilson_half, wilson_center + wilson_half)
cp_ci_v3 <- unname(binom.test(x, n, alternative = "two.sided")$conf.int)
observed_intervals_v3 <- data.frame(
  method = c("Percentile", "Basic", "Wald", "Wilson", "CP"),
  lower = c(percentile_ci_v3[1], basic_ci_v3[1], wald_ci_v3[1], wilson_ci_v3[1], cp_ci_v3[1]),
  upper = c(percentile_ci_v3[2], basic_ci_v3[2], wald_ci_v3[2], wilson_ci_v3[2], cp_ci_v3[2])
)
write.csv(observed_intervals_v3, file.path(result_dir, "R_intervals_v3.csv"), row.names = FALSE)

# 覆盖率按每个真实p对x=0,...,n逐点枚举，并以二项PMF加权求和。
p_grid <- seq(.01, .99, by = .01)
interval_for_x <- function(x_value) {
  ph <- x_value / n
  perc <- qbinom(c(.025, .975), n, ph) / n
  basic <- c(2 * ph - perc[2], 2 * ph - perc[1])
  wh <- z975 * sqrt(ph * (1 - ph) / n)
  wald <- c(ph - wh, ph + wh)
  wd <- 1 + z975^2 / n
  wc <- (ph + z975^2 / (2 * n)) / wd
  wterm <- z975 / wd * sqrt(ph * (1 - ph) / n + z975^2 / (4 * n^2))
  wilson <- c(wc - wterm, wc + wterm)
  cp <- if (x_value == 0) c(0, 1 - .025^(1 / n)) else if (x_value == n) {
    c(.025^(1 / n), 1)
  } else c(qbeta(.025, x_value, n - x_value + 1), qbeta(.975, x_value + 1, n - x_value))
  rbind(Percentile = perc, Basic = basic, Wald = wald, Wilson = wilson, CP = cp)
}
coverage_rows <- vector("list", length(p_grid) * 5)
row_i <- 0L
for (true_p in p_grid) {
  for (method_name in rownames(interval_for_x(0))) {
    covered <- logical(n + 1L)
    for (x_value in 0:n) {
      bounds <- interval_for_x(x_value)[method_name, ]
      covered[x_value + 1L] <- true_p >= bounds[1] - 1e-12 && true_p <= bounds[2] + 1e-12
    }
    row_i <- row_i + 1L
    coverage_rows[[row_i]] <- data.frame(
      p = true_p, method = method_name,
      coverage = sum(dbinom(0:n, size = n, prob = true_p) * covered)
    )
  }
}
coverage_v3 <- do.call(rbind, coverage_rows)
write.csv(coverage_v3, file.path(result_dir, "R_coverage_v3.csv"), row.names = FALSE)

depth_console <- file(file.path(result_dir, "R_depth_v3_console.txt"), open = "wt")
sink(depth_console, split = TRUE)
cat("R:", R.version.string, "\n")
cat("v3: Percentile/Basic use infinite conditional Binomial endpoints; coverage is exact enumeration.\n")
cat("Observed intervals (Percentile/Basic endpoints directly from conditional distribution):\n")
print(observed_intervals_v3, row.names = FALSE)
cat("基础主脚本循环重抽的SD（bootstrap_rates）：", sd(bootstrap_rates), "\n")
cat("Analytic SE:", sqrt(.65 * .35 / 180), "\n")
cat("Coverage at p=.65:\n")
print(coverage_v3[coverage_v3$p == .65, ], row.names = FALSE)
cat("Coverage at p=.05:\n")
print(coverage_v3[coverage_v3$p == .05, ], row.names = FALSE)
sink()
close(depth_console)
