# 本文自编辅助函数：供主脚本调用；先理解主脚本，再阅读递归枚举。
ordinal_stats <- function(tab, row_scores = 1:4, col_scores = 1:4) {
  N <- sum(tab)
  # R按列展开；每一列内部行评分依次变化，列评分保持相同。
  count <- as.vector(tab)
  x <- rep(row_scores, times = ncol(tab))
  y <- rep(col_scores, each = nrow(tab))
  mean_x <- sum(count * x) / N
  mean_y <- sum(count * y) / N
  T <- sum(count * (x - mean_x) * (y - mean_y))
  Lxx <- sum(count * (x - mean_x)^2)
  Lyy <- sum(count * (y - mean_y)^2)
  r <- T / sqrt(Lxx * Lyy)
  Q <- (N - 1) * r^2
  # lower.tail=FALSE读取χ²分布的右尾，df=1是一个评分关联分量。
  c(N = N, T = T, Lxx = Lxx, Lyy = Lyy, r = r,
    VarT = Lxx * Lyy / (N - 1), Q = Q,
    chi2_p = stats::pchisq(Q, df = 1, lower.tail = FALSE))
}

# 给当前行分配整数，每一列不超过剩余名额。
compositions <- function(total, caps) {
  if (length(caps) == 1) {
    if (total <= caps[1]) return(list(c(total)))
    return(list())
  }
  choices <- list()
  for (first in 0:min(total, caps[1])) {
    tails <- compositions(total - first, caps[-1])
    for (tail in tails) {
      choices[[length(choices) + 1]] <- c(first, tail)
    }
  }
  choices
}

enumerate_tables <- function(rows, remain, i = 1, acc = list()) {
  if (i > length(rows)) {
    if (all(remain == 0)) return(list(do.call(rbind, acc)))
    return(list())
  }
  result <- list()
  for (row in compositions(rows[i], remain)) {
    # 当前行分完，从各列名额中扣除，继续分下一行。
    next_tables <- enumerate_tables(rows, remain - row, i + 1, c(acc, list(row)))
    result <- c(result, next_tables)
  }
  result
}

# 仅用于本篇N8小表；大N的直接阶乘会溢出，不对主例完整枚举。
table_probability <- function(tab) {
  numerator <- prod(factorial(rowSums(tab))) * prod(factorial(colSums(tab)))
  denominator <- factorial(sum(tab)) * prod(factorial(tab))
  numerator / denominator
}

conditional_mc <- function(tab, B, seed, row_scores, col_scores) {
  set.seed(seed)
  obs <- ordinal_stats(tab, row_scores, col_scores)
  # r2dtable按正确条件分布抽表，而非均匀抽所有不同表。
  sampled <- stats::r2dtable(B, rowSums(tab), colSums(tab))
  first_table <- sampled[[1]]
  values <- vapply(sampled, function(z) {
    ordinal_stats(z, row_scores, col_scores)["T"]
  }, numeric(1))
  # 浮点容差只用于避免数学上相等的端点因舍入而漏计。
  k <- sum(abs(values) >= abs(obs["T"]) - 1e-10)
  p_hat <- k / B
  interval <- stats::binom.test(k, B, conf.level = 0.95)$conf.int
  summary <- c(B = B, seed = seed, k = k, k_over_B = p_hat,
    plus1_over_plus1 = (k + 1) / (B + 1),
    MCSE = sqrt(p_hat * (1 - p_hat) / B),
    CP95_low = interval[1], CP95_high = interval[2],
    zero_exceedance_one_sided_95_upper = if (k == 0) 1 - 0.05^(1 / B) else NA,
    Tobs = unname(obs["T"]))
  list(summary = summary, T_sim = values, first_table = first_table)
}

plot_distribution <- function(distribution, sim, Tobs) {
  # 相同T的各张表先按概率相加，不能对282张表等权画直方图。
  exact <- aggregate(probability ~ T, distribution, sum)
  support <- exact$T
  frequency <- vapply(support, function(t) mean(abs(sim - t) < 1e-9), numeric(1))
  old <- par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))
  on.exit(par(old))
  colors <- ifelse(abs(support) >= abs(Tobs) - 1e-10, "#F58518", "#4C78A8")
  plot(support, exact$probability, type = "h", lwd = 5, col = colors,
    xlab = "Centered product sum T", ylab = "Exact probability",
    main = "Exact conditional T distribution")
  abline(v = c(-abs(Tobs), abs(Tobs)), lty = 2)
  plot(support, exact$probability, type = "o", col = "#F58518",
    xlab = "T value", ylab = "Probability / frequency",
    main = "Exact probability and Monte Carlo frequency")
  lines(support, frequency, type = "o", col = "#4C78A8", pch = 3)
  legend("topright", c("Exact", "Monte Carlo"), col = c("#F58518", "#4C78A8"), lty = 1)
}
