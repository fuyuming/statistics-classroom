# 自编数值辅助：固定q/df，求CDF达到目标的非负lambda。
lambda_limit <- function(q, df, target) {
  at_zero <- stats::pchisq(q, df, ncp = 0)
  if (q == 0 || at_zero <= target) {
    return(list(value = 0, boundary = TRUE))
  }
  upper <- max(1, q)
  while (stats::pchisq(q, df, ncp = upper) > target) {
    upper <- upper * 2
  }
  objective <- function(value) stats::pchisq(q, df, ncp = value) - target
  answer <- stats::uniroot(objective, c(0, upper), tol = 1e-12)$root
  list(value = answer, boundary = FALSE)
}

# 选读的自编完整函数：与Python/MATLAB配套采用相同的非中心χ²区间定义。
cramers_v_ci <- function(counts, level = .95) {
  stopifnot(is.matrix(counts), min(dim(counts)) >= 2)
  stopifnot(all(is.finite(counts)), all(counts >= 0))
  stopifnot(all(rowSums(counts) > 0), all(colSums(counts) > 0))
  stopifnot(level > 0, level < 1)
  n <- sum(counts)
  expected <- outer(rowSums(counts), colSums(counts)) / n
  q <- sum((counts - expected)^2 / expected)
  df <- (nrow(counts) - 1) * (ncol(counts) - 1)
  m <- min(nrow(counts) - 1, ncol(counts) - 1)
  lower <- lambda_limit(q, df, (1 + level) / 2)
  upper <- lambda_limit(q, df, (1 - level) / 2)
  data.frame(n = n, Q = q, df = df, m = m, V = sqrt(q/(n*m)),
             lambda_L = lower$value, lambda_U = upper$value,
             V_L = min(1, sqrt(lower$value/(n*m))),
             V_U = min(1, sqrt(upper$value/(n*m))),
             boundary_L = lower$boundary, boundary_U = upper$boundary)
}
