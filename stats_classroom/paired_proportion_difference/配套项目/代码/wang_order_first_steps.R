# Wang排序前三步：仅演示第三位候选的选择（R基础版）
#---- 相对路径与固定教学参数
args <- commandArgs(trailingOnly = FALSE)
file_hits <- grep("^--file=", args, value = TRUE)
file_arg <- if (length(file_hits)) sub("^--file=", "", file_hits[1]) else ""
script_dir <- if (nzchar(file_arg)) dirname(normalizePath(file_arg)) else getwd()
project_candidates <- unique(c(script_dir, dirname(script_dir), getwd(), dirname(getwd())))
is_project <- vapply(project_candidates, function(x) {
  dir.exists(file.path(x, "代码")) && dir.exists(file.path(x, "数据"))
}, logical(1))
if (!any(is_project)) stop("请在本篇项目根目录或代码目录运行。")
project <- project_candidates[which(is_project)[1]]
out <- file.path(project, "运行结果", "R")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
n <- 20
alpha_tail <- 0.025
candidates <- list(A = c(19L, 0L, 1L), B = c(18L, 2L, 0L))

#---- 给定theta、pT时，临时加入一个候选后的尾概率
tail_probability <- function(theta, p_t, candidate) {
  a <- (1 - p_t + theta) / 2
  b <- (1 - p_t - theta) / 2
  base <- a^20 + 20 * a^19 * p_t
  if (candidate == "A") base + 20 * a^19 * b else base + 190 * a^18 * p_t^2
}

#---- 全域网格后，在网格峰附近作optimize，并比较两个端点
peak_over_pt <- function(theta, candidate) {
  right <- 1 - abs(theta)
  grid <- seq(0, right, length.out = 401)
  values <- vapply(grid, function(p) tail_probability(theta, p, candidate), numeric(1))
  k <- which.max(values)
  best_p <- grid[k]
  best_q <- values[k]
  left <- grid[max(1, k - 1)]
  right_local <- grid[min(401, k + 1)]
  opt <- optimize(function(p) -tail_probability(theta, p, candidate), c(left, right_local), tol = 1e-14)
  for (p in c(0, right, opt$minimum)) {
    q <- tail_probability(theta, p, candidate)
    if (q > best_q) {
      best_p <- p
      best_q <- q
    }
  }
  c(peak_p_t = best_p, tail_probability = best_q)
}

#---- 对每个候选单独求maxQ(theta)=.025的最小根
rows <- lapply(names(candidates), function(candidate) {
  lower <- 0
  upper <- 1
  for (i in seq_len(70)) {
    middle <- (lower + upper) / 2
    if (peak_over_pt(middle, candidate)["tail_probability"] >= alpha_tail) {
      upper <- middle
    } else {
      lower <- middle
    }
  }
  peak <- peak_over_pt(upper, candidate)
  counts <- candidates[[candidate]]
  data.frame(candidate = candidate, n_plus = counts[1], n_zero = counts[2], n_minus = counts[3],
             D = (counts[1] - counts[3]) / n, trial_lower = upper, peak_p_t = peak["peak_p_t"],
             tail_probability = peak["tail_probability"], selected_for_third = FALSE)
})
summary <- do.call(rbind, rows)
rownames(summary) <- NULL
# 计算后比较临时下限；相同的最大下限同时选入。
summary$selected_for_third <- abs(summary$trial_lower - max(summary$trial_lower)) <= 1e-10

#---- A候选在pT=0时的解析式复核
a_at_a <- (1 + summary$trial_lower[summary$candidate == "A"]) / 2
print(c(A_analytic_pT0 = a_at_a^20 + 20 * a_at_a^19 * (1 - a_at_a)))

#---- 保存两类教学结果
write.csv(summary, file.path(out, "Wang排序第三位试算.csv"), row.names = FALSE, fileEncoding = "UTF-8")
tail_rows <- do.call(rbind, lapply(names(candidates), function(candidate) {
  x <- rbind(c(20, 0, 0), c(19, 1, 0), candidates[[candidate]])
  data.frame(candidate = candidate, n_plus = x[, 1], n_zero = x[, 2], n_minus = x[, 3])
}))
write.csv(tail_rows, file.path(out, "Wang排序候选尾部.csv"), row.names = FALSE, fileEncoding = "UTF-8")
print(summary)
stopifnot(all(abs(summary$tail_probability - alpha_tail) < 1e-10))
stopifnot(summary$selected_for_third[summary$candidate == "B"])
stopifnot(all(rowSums(tail_rows[, c("n_plus", "n_zero", "n_minus")]) == n))
