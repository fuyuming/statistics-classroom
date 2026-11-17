# 统计课堂30：课件原例为新方案50/100、对照30/100的4周达标率。
# UTF-8；可从项目根或 R/ 运行；按第1至5节顺序执行，默认不保存结果。
options(stringsAsFactors = FALSE)
SAVE_OUTPUTS <- identical(tolower(Sys.getenv("SAVE_OUTPUTS", "false")), "true")

# 入口兼容项目根和 R/ 目录；只读取共享课件CSV。
data_candidates <- c("data/lesson10_independent.csv", "../data/lesson10_independent.csv")
data_path <- data_candidates[file.exists(data_candidates)][1]
if (is.na(data_path)) stop("请从项目根或 R/ 目录运行")
project_root <- dirname(dirname(data_path))
result_dir <- file.path(project_root, "results", "R")
if (SAVE_OUTPUTS) dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
save_if_requested <- function(x, name) {
  if (SAVE_OUTPUTS) write.csv(x, file.path(result_dir, name), row.names = FALSE)
}

# 1. 从组名找到两行，并明确方向为新方案减对照。
teaching_data <- read.csv(data_path, check.names = FALSE)
new_row <- teaching_data[teaching_data$group == "new", ][1, ]
control_row <- teaching_data[teaching_data$group == "control", ][1, ]
n1 <- new_row$n
y1 <- new_row$success
n0 <- control_row$n
y0 <- control_row$success
p1 <- y1 / n1
p0 <- y0 / n0
RD <- p1 - p0
tab <- matrix(c(y1, n1 - y1, y0, n0 - y0), nrow = 2, byrow = TRUE,
              dimnames = list(group = c("new", "control"),
                              outcome = c("success", "failure")))
cat("1 读表：p1=", p1, " p0=", p0, " RD=", RD, "\n", sep = "")
print(tab)

# 2. 角变换、SE和双侧正态95% CI；qnorm(.975)是1.96附近的临界值。
h <- 2 * asin(sqrt(p1)) - 2 * asin(sqrt(p0))
hSE <- sqrt(1 / n1 + 1 / n0)
z95 <- qnorm(.975)
hCI <- h + c(-1, 1) * z95 * hSE
baseline <- data.frame(RD = RD, h = h, hSE = hSE,
                       hCI_low = hCI[1], hCI_high = hCI[2])
print(baseline, digits = 12)
save_if_requested(baseline, "baseline_h.csv")

# effectsize要求列为组、行为结局，所以对组×结局表转置核对。
if (!requireNamespace("effectsize", quietly = TRUE)) stop("需要 effectsize")
h_fit <- effectsize::cohens_h(t(tab), ci = .95)
print(h_fit)

# 3. 固定RD=.2；h的符号保留新方案−对照方向，abs_h便于对照分级。
fixed <- data.frame(p0 = c(.1, .3, .4, .7),
                    p1 = c(.3, .5, .6, .9))
fixed$RD <- fixed$p1 - fixed$p0
fixed$h <- 2 * asin(sqrt(fixed$p1)) - 2 * asin(sqrt(fixed$p0))
fixed$abs_h <- abs(fixed$h)
print(fixed, digits = 12)
save_if_requested(fixed[c("p0", "p1", "RD", "h")], "fixed_RD_h.csv")

# 4. 改变每组n；逐元素计算期望数、Q、h CI，并调用无连续性校正Pearson检验。
sample_sizes <- c(20, 100, 1000)
hh <- 2 * asin(sqrt(.5)) - 2 * asin(sqrt(.3))
sample_rows <- vector("list", length(sample_sizes))
for (i in seq_along(sample_sizes)) {
  n <- sample_sizes[i]
  sample_table <- matrix(c(.5 * n, .5 * n, .3 * n, .7 * n), nrow = 2, byrow = TRUE)
  expected <- outer(rowSums(sample_table), colSums(sample_table)) / sum(sample_table)
  Q_manual <- sum((sample_table - expected)^2 / expected)
  pearson_fit <- chisq.test(sample_table, correct = FALSE)
  h_ci <- hh + c(-1, 1) * z95 * sqrt(2 / n)
  sample_rows[[i]] <- data.frame(n = n, h = hh, hCI_low = h_ci[1],
                                 hCI_high = h_ci[2], Q_manual = Q_manual,
                                 Q_native = unname(pearson_fit$statistic),
                                 P_pearson = pearson_fit$p.value)
}
sample_results <- do.call(rbind, sample_rows)
print(sample_results, digits = 12)
save_if_requested(sample_results, "sample_size_h.csv")

# 5. g(p)=2asin(sqrt(p))的斜率；逐元素比较精确增量和线性近似。
variance_ps <- c(.1, .3, .5)
variance_n <- 100
epsilon <- .01
variance_table <- data.frame(p = variance_ps, n = variance_n)
variance_table$raw_SE <- sqrt(variance_table$p * (1 - variance_table$p) / variance_table$n)
variance_table$g_slope <- 1 / sqrt(variance_table$p * (1 - variance_table$p))
variance_table$transformed_SE <- variance_table$g_slope * variance_table$raw_SE
variance_table$approx_variance <- variance_table$transformed_SE^2
variance_table$exact_increment <- 2 * asin(sqrt(variance_table$p + epsilon)) - 2 * asin(sqrt(variance_table$p))
variance_table$linear_increment <- variance_table$g_slope * epsilon
print(variance_table, digits = 12)
save_if_requested(variance_table, "h_variance_stability.csv")
