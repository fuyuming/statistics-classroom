# 统计课堂27：配对样本率差的分布、矩和正态近似（R）
# 问题：个人Zi=后测−前测，+1、0、−1的概率分别为q+=.05、q0=.70、q−=.25。
# 研究方向固定为“后测减前测”，所以θ=q+−q−=-.20；主例20人数据为8双阳、6双阴、1阴转阳、5阳转阴。
# 本文件输出4个distribution_*、moments和paired_examples，并自动调用独立绘图文件。

# 问题/条件/路径：定位项目根、准备输出 ----
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else ""
if (nzchar(script_path) && !file.exists(script_path)) script_path <- gsub("~+~", " ", script_path, fixed = TRUE)
script_dir <- if (nzchar(script_path) && file.exists(script_path)) dirname(normalizePath(script_path, mustWork = TRUE)) else getwd()
candidates <- unique(c(getwd(), script_dir, dirname(script_dir)))
valid <- vapply(candidates, function(path) file.exists(file.path(path, "数据", "paired_20.csv")), logical(1))
if (!any(valid)) stop("找不到项目根目录；请从项目根或代码目录运行。")
root <- normalizePath(candidates[which(valid)[1]], mustWork = TRUE)
data_dir <- file.path(root, "数据")
out_dir <- file.path(root, "运行结果", "R")

# 问题/条件/数据：先查看共享字段和前5行，再确认差值方向 ----
people <- read.csv(file.path(data_dir, "paired_20.csv"), check.names = FALSE)
names(people)
head(people, 5)
changes <- people$after - people$before
stopifnot(all(people$before %in% c(0, 1)))
stopifnot(all(people$after %in% c(0, 1)))
stopifnot(all(changes == people$z))

# 问题/条件/主例：观察比例与模型参数是两类不同对象 ----
# 主例的观察比例来自20名教学构造对象；q参数是用于精确分布的给定模型，不由这20人反推。
main_plus <- sum(changes == 1)
main_zero <- sum(changes == 0)
main_minus <- sum(changes == -1)
main_both_positive <- sum(people$before == 1 & people$after == 1)
main_n <- nrow(people)
main_D <- (main_plus - main_minus) / main_n
before_rate <- (main_both_positive + main_minus) / main_n
after_rate <- (main_both_positive + main_plus) / main_n
main_result <- data.frame(scenario = "main_20", n = main_n, nplus = main_plus,
  nzero = main_zero, nminus = main_minus, D = main_D,
  before_rate = before_rate, after_rate = after_rate, n_changed = main_plus + main_minus)
main_result

# 条件/参数：每次卷积加入1名独立受试者；exact_probability是精确PMF ----
q_plus <- 0.05
q_minus <- 0.25
q_zero <- 1 - q_plus - q_minus

# 手动计算：Zi的均值、方差、D的方差和SE（不同人独立） ----
# 先算E(Zi)、E(Zi²)，再用Var(Zi)=E(Zi²)-E(Zi)^2；独立性使Var(D)=Var(Zi)/n。
E_Z <- q_plus - q_minus
E_Z2 <- q_plus + q_minus
Var_Z <- E_Z2 - E_Z^2
Var_D <- Var_Z / main_n
SE_D <- sqrt(Var_D)
moment_check <- data.frame(E_Z = E_Z, E_Z2 = E_Z2, Var_Z = Var_Z, Var_D = Var_D, SE_D = SE_D)
moment_check

# 结果核对：主例D是后测率减前测率；矩描述给定模型的抽样分布 ----
stopifnot(abs(main_D - (after_rate - before_rate)) < 1e-12)

# 手动计算：逐人卷积的自编算法 ----
# probability索引对应总和k=-person,...,person；一次循环把现有每个概率分别
# 向k-1、k、k+1移动并乘q−、q0、q+。这不是抽样模拟，也不是现成检验函数。
distribution <- function(n, q_plus, q_minus) {
  q_zero <- 1 - q_plus - q_minus
  probability <- 1
  for (person in seq_len(n)) {
    next_probability <- numeric(length(probability) + 2L)
    index <- seq_along(probability)
    next_probability[index] <- next_probability[index] + q_minus * probability
    next_probability[index + 1L] <- next_probability[index + 1L] + q_zero * probability
    next_probability[index + 2L] <- next_probability[index + 2L] + q_plus * probability
    probability <- next_probability
  }
  k <- -n:n
  D <- k / n
  theta <- q_plus - q_minus
  variance_Z <- q_plus + q_minus - theta^2
  SE_D <- sqrt(variance_Z / n)
  edges <- (D[-length(D)] + D[-1]) / 2
  normal_cell_probability <- diff(c(0, pnorm(edges, mean = theta, sd = SE_D), 1))
  data.frame(k = k, D = D, exact_probability = probability,
             normal_cell_probability = normal_cell_probability)
}

# 手动计算：四个场景的精确分布与正态单格概率 ----
# 正态近似把每个离散D值周围的中点区间概率作为单格概率；它不是连续密度值。
scenario_spec <- data.frame(
  scenario = c("normal_n2", "normal_n20", "normal_n1000", "rare_n1000"),
  n = c(2L, 20L, 1000L, 1000L),
  q_plus = c(.05, .05, .05, .001),
  q_minus = c(.25, .25, .25, .002))
moment_rows <- vector("list", nrow(scenario_spec))
for (i in seq_len(nrow(scenario_spec))) {
  setting <- scenario_spec[i, ]
  result <- distribution(setting$n, setting$q_plus, setting$q_minus)
  theta <- setting$q_plus - setting$q_minus
  variance_Z <- setting$q_plus + setting$q_minus - theta^2
  cdf_error <- max(abs(cumsum(result$exact_probability) - cumsum(result$normal_cell_probability)))
  moment_rows[[i]] <- data.frame(scenario = setting$scenario, n = setting$n,
    qp = setting$q_plus, q0 = 1 - setting$q_plus - setting$q_minus, qm = setting$q_minus,
    theta = theta, mean_D = theta, var_D = variance_Z / setting$n,
    max_abs_cdf_error = cdf_error)
}
moments <- do.call(rbind, moment_rows)
moments
# 手动计算：按场景表逐行计算教学比较例 ----
# scenarios.csv给出模型人数；这里保留同一字段和行序，便于查看n=20与n=1000。
scenario_table <- read.csv(file.path(data_dir, "scenarios.csv"), check.names = FALSE)
example_rows <- list(main_result)
for (i in seq_len(nrow(scenario_table))) {
  setting <- scenario_table[i, ]
  if (setting$scenario == "main_20") next
  example_rows[[length(example_rows) + 1L]] <- data.frame(
    scenario = setting$scenario, n = setting$n, nplus = setting$nplus,
    nzero = setting$nzero, nminus = setting$nminus,
    D = (setting$nplus - setting$nminus) / setting$n,
    before_rate = (setting$both_positive + setting$nminus) / setting$n,
    after_rate = (setting$both_positive + setting$nplus) / setting$n,
    n_changed = setting$nplus + setting$nminus)
}
paired_examples <- do.call(rbind, example_rows)
paired_examples
# 可选保存：默认不写文件；改TRUE后一次性更新7份CSV ----
save_results <- FALSE
if (save_results) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  for (i in seq_len(nrow(scenario_spec))) {
    setting_i <- scenario_spec[i, ]
    distribution_i <- distribution(setting_i$n, setting_i$q_plus, setting_i$q_minus)
    write.csv(distribution_i,
              file.path(out_dir, paste0("distribution_", scenario_spec$scenario[i], ".csv")),
              row.names = FALSE, fileEncoding = "UTF-8")
  }
  write.csv(moments, file.path(out_dir, "moments.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(paired_examples, file.path(out_dir, "paired_examples.csv"), row.names = FALSE, fileEncoding = "UTF-8")
}

# 主入口：自动调用绘图文件；绘图文件独立生成绘图数据 ----
source(file.path(root, "代码", "plot_figures.R"))

# 思考/可选练习：改q+、q−或n，比较精确PMF和正态单格概率；正态近似不等于精确概率。
