# 三类人数的完整概率表：20人、+1/0/-1三类（UTF-8）。
#---- 读取路径与教学参数
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
here <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
candidates <- unique(c(here, normalizePath(file.path(here, ".."), mustWork = FALSE), getwd()))
has_project_dirs <- vapply(candidates, function(x) dir.exists(file.path(x, "数据")) && dir.exists(file.path(x, "运行结果")), logical(1))
project <- normalizePath(candidates[which(has_project_dirs)[1]], mustWork = TRUE)
data_dir <- file.path(project, "数据")
out_dir <- file.path(project, "运行结果", "R")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
n <- 20L
theta <- -0.45
p_t <- 0.44
q_plus <- (1 - p_t + theta) / 2
q_zero <- p_t
q_minus <- (1 - p_t - theta) / 2

#---- 观察组合：一种安排概率与安排种数
r_obs <- 1L
t_obs <- 14L
s_obs <- 5L
observed_arrangements <- choose(n, r_obs) * choose(n-r_obs, t_obs)
single_observation <- q_plus^r_obs * q_zero^t_obs * q_minus^s_obs
observed_probability <- observed_arrangements * single_observation
observed <- data.frame(n_plus=r_obs, n_zero=t_obs, n_minus=s_obs,
  arrangements=observed_arrangements, single_arrangement_probability=single_observation,
  probability=observed_probability)
print(observed)

#---- 双循环列出全部231个三类人数，并用dmultinom作对照
tail <- read.csv(file.path(data_dir, "wang_lower_tail_n20.csv"), check.names = FALSE)
tail_key <- paste(tail$n_plus, tail$n_zero, tail$n_minus, sep = ",")
rows <- vector("list", (n + 1L) * (n + 2L) / 2L)
k <- 0L
for (r in 0:n) {
  for (t in 0:(n-r)) {
    s <- n-r-t
    k <- k + 1L
    arrangements <- choose(n, r) * choose(n-r, t)
    single <- q_plus^r * q_zero^t * q_minus^s
    probability <- dmultinom(c(r, t, s), size = n, prob = c(q_plus, q_zero, q_minus))
    if (!isTRUE(all.equal(probability, arrangements * single, tolerance = 1e-12))) stop("dmultinom与拆解公式不一致")
    rows[[k]] <- data.frame(n_plus=r, n_zero=t, n_minus=s, arrangements=arrangements,
      single_arrangement_probability=single, probability=probability,
      in_lower_tail=paste(r,t,s,sep=",") %in% tail_key)
  }
}
table <- do.call(rbind, rows)
write.csv(table, file.path(out_dir, "三类人数的231组概率.csv"), row.names = FALSE, fileEncoding = "UTF-8")
tail_probability <- sum(table$probability[table$in_lower_tail])
all_probability <- sum(table$probability)
summary <- data.frame(n=n, theta=theta, p_t=p_t, q_plus=q_plus, q_zero=q_zero, q_minus=q_minus,
  combinations=nrow(table), tail_combinations=sum(table$in_lower_tail), other_combinations=sum(!table$in_lower_tail),
  all_probability=all_probability, tail_probability=tail_probability, other_probability=all_probability-tail_probability,
  observed_probability=observed_probability, observed_arrangements=observed_arrangements)
write.csv(summary, file.path(out_dir, "三类人数概率汇总.csv"), row.names = FALSE, fileEncoding = "UTF-8")
print(summary)

#---- 固定Wang下限尾部的候选theta核对：4001点后在峰值邻域局部优化
candidate_theta <- c(-0.45, -0.44923, -0.44)
cutoff <- 0.025
if (nrow(tail) != 141L || any(!is.finite(as.matrix(tail))) ||
    any(tail$n_plus < 0 | tail$n_zero < 0 | tail$n_minus < 0) ||
    any(tail$n_plus != floor(tail$n_plus) | tail$n_zero != floor(tail$n_zero) |
        tail$n_minus != floor(tail$n_minus)) ||
    any(tail$n_plus + tail$n_zero + tail$n_minus != n)) {
  stop("wang_lower_tail_n20.csv必须恰好有141行，且每行是非负整数并且总和为20")
}
tail_probability <- function(theta_value, p_value) {
  q_plus_value <- (1 - p_value + theta_value) / 2
  q_minus_value <- (1 - p_value - theta_value) / 2
  if (q_plus_value < -1e-12 || p_value < -1e-12 || q_minus_value < -1e-12) return(NA_real_)
  q_plus_value <- max(0, q_plus_value)
  p_value <- max(0, p_value)
  q_minus_value <- max(0, q_minus_value)
  coefficient <- factorial(n) / (factorial(tail$n_plus) * factorial(tail$n_zero) * factorial(tail$n_minus))
  sum(coefficient * q_plus_value^tail$n_plus * p_value^tail$n_zero * q_minus_value^tail$n_minus)
}
candidate_rows <- lapply(candidate_theta, function(theta_value) {
  grid <- seq(0, 1 - abs(theta_value), length.out = 4001L)
  values <- vapply(grid, function(p) tail_probability(theta_value, p), numeric(1))
  peak <- which.max(values)
  left <- if (peak > 1L) grid[peak - 1L] else grid[peak]
  right <- if (peak < length(grid)) grid[peak + 1L] else grid[peak]
  refined <- optimize(function(p) -tail_probability(theta_value, p), interval = c(left, right), tol = 1e-14)
  data.frame(theta = theta_value, p_t_at_max = refined$minimum,
             max_tail_probability = -refined$objective, cutoff = cutoff)
})
candidate_out <- do.call(rbind, candidate_rows)
write.csv(candidate_out, file.path(out_dir, "候选率差的尾部概率.csv"), row.names = FALSE, fileEncoding = "UTF-8")
print(candidate_out)
