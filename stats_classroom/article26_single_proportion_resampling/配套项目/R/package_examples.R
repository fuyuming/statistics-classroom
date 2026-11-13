# 统计课堂26：用infer与confintr核对相同的单样本率问题。
# 假设180人中117人达标；H0:p=.45，H1:p>.45。
# Bootstrap与CP比较用双侧95%；原题精确检验保留greater方向。
# 先打开本项目Rproj，再按节运行；基础三语脚本不依赖这些扩展包。

# 1. 检查依赖与项目根 ----
project_root <- getwd()
args <- commandArgs(trailingOnly = FALSE)
file_arg <- args[grepl("^--file=", args)]
if (length(file_arg) > 0) {
  script_path <- sub("^--file=", "", file_arg[1])
  project_root <- normalizePath(file.path(dirname(script_path), ".."))
}
if (!file.exists(file.path(project_root, "data", "cases.csv"))) {
  stop("请先打开本项目根目录的Rproj。")
}
if (!requireNamespace("infer", quietly = TRUE) || !requireNamespace("confintr", quietly = TRUE)) {
  stop("请先在Console运行 install.packages(c('infer', 'confintr'))，再运行本节。")
}
out <- file.path(project_root, "results")
dir.create(out, showWarnings = FALSE)
sink(file.path(out, "R_packages_console.txt"), split = TRUE)
cat("infer:", as.character(packageVersion("infer")), "\n")
cat("confintr:", as.character(packageVersion("confintr")), "\n")
cat("R:", paste(R.version$major, R.version$minor, sep = "."), "\n")

# 2. 数据、Bootstrap率与双侧百分位区间 ----
sample_data <- data.frame(
  outcome = factor(c(rep("yes", 117), rep("no", 63)), levels = c("no", "yes"))
)

set.seed(20261110)
bootstrap_samples <- sample_data |>
  infer::specify(response = outcome, success = "yes") |>
  infer::generate(reps = 10000, type = "bootstrap") |>
  infer::calculate(stat = "prop")
bootstrap_ci <- infer::get_confidence_interval(
  bootstrap_samples, level = .95, type = "percentile"
)
cat("\ninfer 非参数bootstrap百分位区间\n")
print(bootstrap_ci)


# 3. 在H0下模拟：这里的率由p=.45给定 ----
set.seed(20261110)
null_samples <- sample_data |>
  infer::specify(response = outcome, success = "yes") |>
  infer::hypothesize(null = "point", p = .45) |>
  infer::generate(reps = 10000, type = "draw") |>
  infer::calculate(stat = "prop")
cat("\ninfer 原假设模拟；原题单侧 greater\n")
cat("极端命中数:", sum(null_samples$stat >= .65), "\n")
print(infer::get_p_value(null_samples, obs_stat = .65, direction = "greater"))
cat("注意：上述0是有限次模拟的命中比例。原题精确二项P值如下。\n")
print(stats::binom.test(117, 180, p = .45, alternative = "greater"))


# 4. 用另一个包计算同一类型的百分位区间 ----
# 本篇只比较双侧95%百分位区间，明确覆盖默认BCa设置。
y <- c(rep(1, 117), rep(0, 63))
cat("\nconfintr 百分位bootstrap区间；显式指定perc\n")
confintr_ci <- confintr::ci_proportion(
  x = y, probs = c(.025, .975), type = "bootstrap",
  boot_type = "perc", R = 10000, seed = 20261110
)
print(confintr_ci)
str(confintr_ci)


# 5. 确定性基准与边界例 ----
cat("\n无限次bootstrap的条件二项分位区间（离散分位定义）\n")
print(qbinom(c(.025, .975), size = 180, prob = .65) / 180)
cat("\n边界0/20的双侧95%CP区间\n")
print(stats::binom.test(0, 20, conf.level = .95)$conf.int)


# 6. 保存供配图使用的实际运行分布 ----
write.csv(bootstrap_samples, file.path(project_root, "data", "infer_bootstrap_reference.csv"), row.names = FALSE)
write.csv(null_samples, file.path(project_root, "data", "infer_null_reference.csv"), row.names = FALSE)
write.csv(bootstrap_ci, file.path(out, "infer_bootstrap_ci.csv"), row.names = FALSE)
sink()
