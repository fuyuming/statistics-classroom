# 统计课堂26：用 boot::boot 与 rbinom 演示单样本率重抽样。
# 数据：180人中117人达标；观察率为 117/180 = 65%。
# 第一部分从当前样本重抽（ordinary bootstrap），第二部分在 H0:p=.45 下模拟。

# 1. 定位项目根目录并准备输出 ----
project_root <- getwd()
args <- commandArgs(trailingOnly = FALSE)
file_arg <- args[grepl("^--file=", args)]
if (length(file_arg) > 0) {
  script_path <- sub("^--file=", "", file_arg[1])
  project_root <- normalizePath(file.path(dirname(script_path), ".."))
}
if (!file.exists(file.path(project_root, "data", "cases.csv"))) {
  stop("请从本项目运行；未找到 data/cases.csv。")
}
out_dir <- file.path(project_root, "results")
data_dir <- file.path(project_root, "data")
dir.create(out_dir, showWarnings = FALSE)

# 将版本、种子、每一步打印内容完整保存到控制台证据文件。
sink(file.path(out_dir, "boot_v2_console.txt"), split = TRUE)
cat("R:", R.version.string, "\n")
cat("boot:", as.character(packageVersion("boot")), "\n")
cat("Seed: 20261110; n: 180; observed count: 117; observed rate: 0.65\n")

if (!requireNamespace("boot", quietly = TRUE)) {
  stop("需要安装 boot 包。")
}

# 2. 从当前样本重抽：ordinary Bootstrap ----
# y 用 1/0 表示达标/未达标；每轮仍抽取180个观测，但允许重复。
y <- c(rep(1, 117), rep(0, 63))
rate_stat <- function(data, indices) {
  mean(data[indices])
}

set.seed(20261110)
result <- boot::boot(data = y, statistic = rate_stat, R = 10000, sim = "ordinary")
boot_rates <- result$t[, 1]
boot_counts <- round(180 * boot_rates)
cat("\n[1] boot::boot ordinary Bootstrap\n")
cat("R:", nrow(result$t), "\n")
cat("First six counts (round(180 * boot_rates)):\n")
print(head(boot_counts))
cat("Percentile 95% CI for rate:\n")
print(quantile(boot_rates, c(.025, .975)))

boot_reference <- data.frame(
  replicate = seq_along(boot_rates),
  stat = as.numeric(boot_rates),
  count = as.integer(boot_counts)
)
write.csv(boot_reference, file.path(data_dir, "boot_v2_reference.csv"), row.names = FALSE)

# 3. H0 模拟：每轮直接生成180人的达标计数 ----
# 这里不再从当前样本重抽，而是假定 H0:p=.45，并计算达到117人或以上的比例。
set.seed(20261110)
null_counts <- rbinom(n = 10000, size = 180, prob = .45)
null_hits <- sum(null_counts >= 117)
null_hit_proportion <- null_hits / length(null_counts)
cat("\n[2] H0 simulation with rbinom\n")
cat("H0: p = 0.45; size = 180; R = 10000\n")
cat("First six simulated counts:\n")
print(head(null_counts))
cat("Null hits (count >= 117):", null_hits, "\n")
cat("Simulated hit proportion:", null_hit_proportion, "\n")

rbinom_reference <- data.frame(
  replicate = seq_along(null_counts),
  count = as.integer(null_counts),
  stat = as.numeric(null_counts / 180)
)
write.csv(rbinom_reference, file.path(data_dir, "rbinom_v2_reference.csv"), row.names = FALSE)

cat("\nWrote: data/boot_v2_reference.csv\n")
cat("Wrote: data/rbinom_v2_reference.csv\n")
sink()
