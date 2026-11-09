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
