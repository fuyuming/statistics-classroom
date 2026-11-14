# 统计课堂27：配对样本率差的基础分布（R）
# 课堂问题：若每个人的变化 Zi=后测−前测只可能为 +1、0、−1，如何得到D=平均变化的精确概率分布？
# 教学构造：q+=0.05、q0=0.70、q−=0.25；20人主例为8双阳、6双阴、1阴转阳、5阳转阴。

# 问题/条件/数据身份：定位项目、先查看共享数据字段和前5行 ----
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else ""
if (nzchar(script_path) && !file.exists(script_path)) script_path <- gsub("~+~", " ", script_path, fixed = TRUE)
script_dir <- if (nzchar(script_path) && file.exists(script_path)) dirname(normalizePath(script_path, mustWork = TRUE)) else getwd()
candidates <- unique(c(getwd(), script_dir, dirname(script_dir)))
valid <- vapply(candidates, function(path) file.exists(file.path(path, "数据", "paired_20.csv")), logical(1))
if (!any(valid)) stop("找不到项目根目录；请从项目根或代码目录运行。")
project <- normalizePath(candidates[which(valid)[1]], mustWork = TRUE)
data_dir <- file.path(project, "数据")
out_dir <- file.path(project, "运行结果", "R")
people <- read.csv(file.path(data_dir, "paired_20.csv"), check.names = FALSE)
names(people)
head(people, 5)

# 条件/参数：Zi=后测−前测，+1、0、−1分别对应阴转阳、不变、阳转阴 ----
n <- 20L
q_plus <- 0.05
q_zero <- 0.70
q_minus <- 0.25
stopifnot(abs(q_plus + q_zero + q_minus - 1) < 1e-12)
changes <- people$after - people$before
stopifnot(all(changes == people$z))

# 手动计算：两人9种有序组合及其概率 ----
# expand.grid每行是一种(z1,z2)组合；独立条件下，两人的安排概率相乘。
z <- expand.grid(z2 = c(-1L, 0L, 1L), z1 = c(-1L, 0L, 1L))
prob_by_z <- c(`-1` = q_minus, `0` = q_zero, `1` = q_plus)
two <- data.frame(z1 = z$z1, z2 = z$z2, D = (z$z1 + z$z2) / 2,
                  p1 = unname(prob_by_z[as.character(z$z1)]),
                  p2 = unname(prob_by_z[as.character(z$z2)]))
two$probability <- two$p1 * two$p2
two

# 手动计算：按相同D合并两人组合概率 ----
two_rate <- aggregate(probability ~ D, two, sum)
two_rate

# 手动计算：主例人数安排数×单个安排概率 ----
# 先选n_plus人，再从剩余人中选n_zero人；剩余者为n_minus人。
n_plus_obs <- 1L
n_zero_obs <- 14L
n_minus_obs <- 5L
observed_arrangements <- choose(n, n_plus_obs) * choose(n - n_plus_obs, n_zero_obs)
single_arrangement_probability <- q_plus^n_plus_obs * q_zero^n_zero_obs * q_minus^n_minus_obs
observed_probability <- observed_arrangements * single_arrangement_probability
observed_result <- data.frame(n_plus = n_plus_obs, n_zero = n_zero_obs, n_minus = n_minus_obs,
                              arrangements = observed_arrangements,
                              single_arrangement_probability = single_arrangement_probability,
                              probability = observed_probability)
observed_result

# 调用函数：stats::dmultinom核对一个人数安排 ----
# dmultinom直接给出三类人数的多项分布概率；它不是上面的手算过程。
function_probability <- dmultinom(c(n_plus_obs, n_zero_obs, n_minus_obs), size = n,
                                  prob = c(q_plus, q_zero, q_minus))
function_probability

# 结果核对：手算与dmultinom应一致 ----
stopifnot(abs(function_probability - observed_probability) < 1e-12)

# 手动计算：枚举20人的231组三类人数 ----
# 外层固定n_plus，内层枚举n_zero；n_minus由总人数约束得到。
rows <- vector("list", (n + 1L) * (n + 2L) / 2L)
k <- 0L
for (n_plus in 0:n) {
  for (n_zero in 0:(n - n_plus)) {
    n_minus <- n - n_plus - n_zero
    arrangements <- choose(n, n_plus) * choose(n - n_plus, n_zero)
    single <- q_plus^n_plus * q_zero^n_zero * q_minus^n_minus
    k <- k + 1L
    rows[[k]] <- data.frame(n_plus = n_plus, n_zero = n_zero, n_minus = n_minus,
                            arrangements = arrangements,
                            single_arrangement_probability = single,
                            probability = arrangements * single,
                            D = (n_plus - n_minus) / n)
  }
}
table <- do.call(rbind, rows)
table

# 调用函数：stats::dmultinom逐行核对231组概率 ----
# 这里仅核对概率列，枚举结构和D仍由前面的步骤展开。
function_probability_all <- vapply(seq_len(nrow(table)), function(i) {
  dmultinom(as.numeric(table[i, c("n_plus", "n_zero", "n_minus")]), n,
            c(q_plus, q_zero, q_minus))
}, numeric(1))
stopifnot(max(abs(table$probability - function_probability_all)) < 1e-12)

# 手动计算：按D合并并求全概率、均值、方差 ----
rate_diff <- aggregate(probability ~ D, table, sum)
same_difference <- subset(table, n_plus - n_minus == -4L)
mean_D <- sum(table$D * table$probability)
var_D <- sum((table$D - mean_D)^2 * table$probability)
rate_diff
same_difference

# 结果核对：精确分布汇总 ----
summary <- data.frame(n = n, q_plus = q_plus, q_zero = q_zero, q_minus = q_minus,
                      combinations = nrow(table), all_probability = sum(table$probability),
                      mean_D = mean_D, var_D = var_D,
                      observed_arrangements = observed_arrangements,
                      observed_probability = observed_probability)
summary

# 可选保存：默认不写文件；改TRUE后一次性更新六份CSV ----
save_results <- FALSE
if (save_results) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(two, file.path(out_dir, "基础分布_两人9种组合.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(two_rate, file.path(out_dir, "基础分布_两人率差概率.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(table, file.path(out_dir, "基础分布_231组概率.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(rate_diff, file.path(out_dir, "基础分布_率差概率.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(same_difference, file.path(out_dir, "基础分布_同率差9组概率.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(summary, file.path(out_dir, "基础分布_汇总.csv"), row.names = FALSE, fileEncoding = "UTF-8")
}

# 思考/可选练习：改q_plus、q_zero、q_minus或n，先预测D的方向和离散程度再重跑。
