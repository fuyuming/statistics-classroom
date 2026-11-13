# 统计课堂26新增操作：固定参考率的Bootstrap率差区间与单组换序。
# 从配套项目根目录在RStudio按节运行，或用Rscript运行此文件。
project_root <- getwd()
args <- commandArgs(trailingOnly = FALSE)
f <- args[grepl("^--file=", args)]
if (length(f)) project_root <- normalizePath(file.path(dirname(sub("^--file=", "", f[1])), ".."))
stopifnot(file.exists(file.path(project_root, "data", "cases.csv")))
local({
  sink(file.path(project_root, "results", "article26_v8_additions_console.txt"), split=TRUE)
  on.exit(sink(), add=TRUE)
  ci <- read.csv(file.path(project_root, "data", "infer_v4_bootstrap_ci.csv"))
  obs_stat <- data.frame(stat=117/180)
  dat <- data.frame(outcome=factor(c(rep("yes",117),rep("no",63)),levels=c("no","yes")))
  cat("Bootstrap率差：点估计、双侧95%区间\n")
effect <- obs_stat$stat - .45
effect_ci <- c(ci$lower_ci, ci$upper_ci) - .45
print(effect)
print(effect_ci)

  cat("三轮单组换序得到的率\n")
set.seed(20261110)
shuffle_rates <- replicate(
  3, mean(sample(dat$outcome, replace = FALSE) == "yes")
)
print(shuffle_rates)

  stopifnot(abs(effect-.2)<1e-12, all(abs(shuffle_rates-.65)<1e-12))
  print(R.version.string)
})
