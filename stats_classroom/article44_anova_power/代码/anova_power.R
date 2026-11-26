# 统计课堂44：ANOVA功效与样本量（R入口，UTF-8）
# 独立、正态、同方差四组总体的研究前规划；均数和SD不是实际观测。
# 依赖：基础R和pwr。RStudio打开项目后可按01到05逐节运行。

# 01 先看输入：读取规划数据、设置与依赖 ----
if (!requireNamespace("pwr", quietly = TRUE)) stop("请先安装pwr包。")
root <- getwd()
if (basename(root) == "代码") root <- dirname(root)
data_dir <- file.path(root, "数据")
result_dir <- file.path(root, "运行结果", "R")
if (!file.exists(file.path(data_dir, "planning_scenarios.csv"))) {
  script_arg <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(script_arg) == 1) {
    root <- dirname(dirname(normalizePath(sub("^--file=", "", script_arg))))
    data_dir <- file.path(root, "数据")
    result_dir <- file.path(root, "运行结果", "R")
  }
}
scenarios <- read.csv(file.path(data_dir, "planning_scenarios.csv"), stringsAsFactors = FALSE)
settings <- read.csv(file.path(data_dir, "settings.csv"), stringsAsFactors = FALSE)
get_setting <- function(name) as.numeric(settings$value[settings$setting == name])
alpha <- get_setting("alpha")
target_power <- get_setting("target_power")
high_target_power <- get_setting("high_target_power")
n_check <- get_setting("n_check")
k <- length(unique(scenarios$group))
save_results <- identical(Sys.getenv("ARTICLE44_SAVE_RESULTS"), "1")
cat("R", as.character(getRversion()), "；pwr", as.character(packageVersion("pwr")), "\n")

cat("规划数据字段：", paste(names(scenarios), collapse = ", "), "\n")
print(head(scenarios, 8))
cat("设置：\n")
print(settings)
scenario_means <- function(scenario) {
  rows <- scenarios[scenarios$scenario == scenario, ]
  rows <- rows[order(rows$group), ]
  rows$mean_MPa
}
scenario_sd <- function(scenario) {
  values <- unique(scenarios$sd_MPa[scenarios$scenario == scenario])
  if (length(values) != 1 || !is.finite(values) || values <= 0)
    stop("每个scenario必须有一个正的共同SD：", scenario)
  values
}
large_means <- scenario_means("large")
small_means <- scenario_means("small")
large_sd <- scenario_sd("large")
small_sd <- scenario_sd("small")
if (length(large_means) < 2 || length(small_means) < 2)
  stop("每个scenario至少需要两个组均数。")
cat("large均数：", paste(large_means, collapse = ", "), "MPa；small均数：",
    paste(small_means, collapse = ", "), "MPa；large SD=", large_sd,
    "，small SD=", small_sd, " MPa。\n")
cat("均数和SD是研究前设定的总体参数，不是一次小样本观察结果。\n")

# 02 Cohen f与80/90%规划 ----
cohen_f <- function(means, sd) {
  if (!is.finite(sd) || sd <= 0) stop("共同SD必须为正数。")
  if (length(unique(means)) == 1) stop("四组均数完全相同时，f=0，无法用当前规划反解样本量；请先给出有实际意义的差异。")
  sqrt(mean((means - mean(means))^2)) / sd
}
power_at_n <- function(n, means, sd, level = alpha) {
  pwr::pwr.anova.test(k = length(means), n = n, f = cohen_f(means, sd),
                      sig.level = level, power = NULL)$power
}
plan_one <- function(scenario, target) {
  means <- scenario_means(scenario)
  sd <- scenario_sd(scenario)[1]
  f_value <- cohen_f(means, sd)
  power_n2 <- power_at_n(2, means, sd)
  # 大效应可能在每组2个时已达目标，因此不强行求不存在的连续根。
  continuous_n <- if (power_n2 >= target) NA_real_ else
    pwr::pwr.anova.test(k = k, f = f_value, sig.level = alpha, power = target)$n
  n_integer <- if (is.na(continuous_n)) 2 else ceiling(continuous_n)
  previous <- if (n_integer > 2) power_at_n(n_integer - 1, means, sd) else NA_real_
  data.frame(scenario = scenario, target = target, f = f_value,
             continuous_n = continuous_n, n_per_group = n_integer,
             total_n = k * n_integer, power = power_at_n(n_integer, means, sd),
             previous_power = previous)
}
planning <- rbind(plan_one("large", target_power), plan_one("large", high_target_power),
                  plan_one("small", target_power), plan_one("small", high_target_power))
cat("Cohen f及目标功效：\n")
print(planning)
cat("continuous_n向上取整为n_per_group，再回代power和前一整数的功效。\n")

# 03 固定n_check：pwr返回这次样本量下的功效 ----
fixed_n <- data.frame(
  scenario = c("large", "small"), n_per_group = n_check,
  f = c(cohen_f(large_means, large_sd), cohen_f(small_means, small_sd)),
  power = c(power_at_n(n_check, large_means, large_sd),
            power_at_n(n_check, small_means, small_sd)))
cat("固定每组", n_check, "个：\n")
print(fixed_n)
cat("功效是重复研究中拒绝整体H0的长期比例，不是一次实验的P值。\n")

# 04 SD敏感性：small均数组合固定，改变共同SD ----
sd_values <- c(2.5, 5, 10)
sensitivity <- do.call(rbind, lapply(sd_values, function(sd) {
  f_value <- cohen_f(small_means, sd)
  continuous_n <- pwr::pwr.anova.test(k = k, f = f_value,
                                      sig.level = alpha, power = target_power)$n
  n_integer <- ceiling(continuous_n)
  data.frame(sigma = sd, f = f_value, continuous_n = continuous_n,
             n_per_group = n_integer, power = power_at_n(n_integer, small_means, sd),
             previous_power = power_at_n(n_integer - 1, small_means, sd))
}))
cat("SD敏感性（small均数固定）：\n")
print(sensitivity)
# 05 可选保存：三张核心结果表 ----
result_tables <- list(planning = planning, fixed_n = fixed_n,
                      sensitivity = sensitivity)
if (save_results) {
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
  for (name in names(result_tables)) write.csv(result_tables[[name]],
      file.path(result_dir, paste0(name, ".csv")), row.names = FALSE, fileEncoding = "UTF-8")
  cat("已写入3个R结果表。\n")
}

# IDE默认画出small情形的功效曲线；自动运行可设ARTICLE44_PLOT=0。
if (!identical(Sys.getenv("ARTICLE44_PLOT"), "0")) {
  sizes <- 2:max(100, ceiling(1.5 * planning$n_per_group[
    planning$scenario == "small" & planning$target == target_power]))
  curve <- vapply(sizes, function(n) power_at_n(n, small_means, small_sd), numeric(1))
  plot(sizes, curve, type = "l", lwd = 2, col = "#366FAD", ylim = c(0, 1),
       xlab = "Independent units per group", ylab = "Power")
  abline(h = target_power, lty = 2, col = "gray50")
}
