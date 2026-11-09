# 基础分布：固定已知的三类概率，展示率差D的精确分布。
#---- 定位项目目录，支持RStudio项目根或代码目录Source
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
candidates <- unique(c(getwd(), script_dir, dirname(script_dir)))
ok <- vapply(candidates, function(x) dir.exists(file.path(x, "运行结果")) && dir.exists(file.path(x, "数据")), logical(1))
project <- normalizePath(candidates[which(ok)[1]], mustWork = TRUE)
out_dir <- file.path(project, "运行结果", "R")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

#---- 固定模型参数与n=2的D=0小演示
n <- 20L
q_plus <- 0.05
q_zero <- 0.70
q_minus <- 0.25
n2_demo <- data.frame(n_plus=c(0L, 1L), n_zero=c(2L, 0L), n_minus=c(0L, 1L), probability=c(q_zero^2, 2*q_plus*q_minus))
print(n2_demo)
print(data.frame(D=c(-1, -0.5, 0, 0.5, 1),
                 probability=c(q_minus^2, 2*q_minus*q_zero, q_zero^2+2*q_plus*q_minus, 2*q_plus*q_zero, q_plus^2)))

#---- 观察人数：安排数、指定人员安排概率和总概率
n_plus_obs <- 1L
n_zero_obs <- 14L
n_minus_obs <- 5L
observed_arrangements <- choose(n, n_plus_obs) * choose(n-n_plus_obs, n_zero_obs)
single_arrangement_probability <- q_plus^n_plus_obs * q_zero^n_zero_obs * q_minus^n_minus_obs
observed_probability <- observed_arrangements * single_arrangement_probability
print(data.frame(n_plus=n_plus_obs, n_zero=n_zero_obs, n_minus=n_minus_obs, arrangements=observed_arrangements, single_arrangement_probability=single_arrangement_probability, probability=observed_probability))

#---- 双循环枚举231组三类人数，并以dmultinom核对拆解公式
rows <- vector("list", (n+1L)*(n+2L)/2L); k <- 0L
for (n_plus in 0:n) for (n_zero in 0:(n-n_plus)) {
  n_minus <- n - n_plus - n_zero
  arrangements <- choose(n, n_plus) * choose(n-n_plus, n_zero)
  single <- q_plus^n_plus * q_zero^n_zero * q_minus^n_minus
  probability <- dmultinom(c(n_plus, n_zero, n_minus), n, c(q_plus, q_zero, q_minus))
  stopifnot(isTRUE(all.equal(probability, arrangements*single, tolerance=1e-12)))
  k <- k + 1L
  rows[[k]] <- data.frame(n_plus=n_plus, n_zero=n_zero, n_minus=n_minus, arrangements=arrangements, single_arrangement_probability=single, probability=probability, D=(n_plus-n_minus)/n)
}
table <- do.call(rbind, rows)
write.csv(table, file.path(out_dir, "基础分布_231组概率.csv"), row.names=FALSE, fileEncoding="UTF-8")

#---- 按相同D相加，并汇总全概率、均值和方差
rate_diff <- aggregate(probability ~ D, table, sum)
write.csv(rate_diff, file.path(out_dir, "基础分布_率差概率.csv"), row.names=FALSE, fileEncoding="UTF-8")
mean_D <- sum(table$D*table$probability)
summary <- data.frame(n=n, q_plus=q_plus, q_zero=q_zero, q_minus=q_minus, combinations=nrow(table), all_probability=sum(table$probability), mean_D=mean_D, var_D=sum((table$D-mean_D)^2*table$probability), observed_arrangements=observed_arrangements, observed_probability=observed_probability)
write.csv(summary, file.path(out_dir, "基础分布_汇总.csv"), row.names=FALSE, fileEncoding="UTF-8")
print(summary)
