# 小样本率的置信区间：从6/13找到上下限
# 目标：固定数据文件中的观察成功数x和试验次数n（默认x=6、n=13），直接用二项尾概率反演95%区间。
# 来源：颜艳《医学统计学》第5版例6-2、例6-5，第92—94页；这是教材案例参数，不是个人原始数据。
# 运行：打开小样本率区间估计.Rproj后Source本文件；只使用R基础包。
# 编码：UTF-8。主要结果对象保留在工作区，便于初学者直接查看。

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
root <- if (length(file_arg) > 0) {
  dirname(dirname(normalizePath(file_arg[1])))
} else {
  getwd()
}
out <- file.path(root, "运行结果", "R")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

parameters <- read.csv(
  file.path(root, "数据", "案例参数.csv"),
  stringsAsFactors = FALSE,
  fileEncoding = "UTF-8",
  check.names = FALSE
)
x <- as.integer(parameters$x[1])
n <- as.integer(parameters$n[1])
confidence <- parameters$confidence[1]
alpha <- parameters$alpha[1]
scan_step <- parameters$scan_step[1]
display_rounds <- as.integer(parameters$bisection_rounds[1])
tail_target <- alpha / 2

# 1. 先看样本率和Wald近似区间。
sample_proportion <- x / n
wald_se <- sqrt(sample_proportion * (1 - sample_proportion) / n)
# 双侧95%：左边累计概率=1-0.025=0.975，因此查97.5%分位点。
z_975 <- qnorm(1 - tail_target)
z_975
# 95%时约1.959964，不是从本次13人算出的值。
normal_central_area <- pnorm(z_975) - pnorm(-z_975)
normal_central_area
# [1] 0.95：曲线下的面积，不能用曲线高度替代。
wald_low <- sample_proportion - z_975 * wald_se
wald_high <- sample_proportion + z_975 * wald_se
wald_result <- data.frame(
  quantity = c("sample_proportion", "wald_se", "wald_low", "wald_high"),
  value = c(sample_proportion, wald_se, wald_low, wald_high)
)
write.csv(wald_result, file.path(out, "Wald近似.csv"), row.names = FALSE, fileEncoding = "UTF-8")
# 参考：sample_proportion=0.46153846，Wald区间约为0.1905457—0.7325312。

# 1.1 Wilson区间（课堂外补充），不作连续性校正。
# 与Wald的区别：score检验的分母采用候选p0，而不是固定的样本率。
# 将 n*(p_hat-p0)^2 <= z^2*p0*(1-p0) 解成二次不等式。
wilson_denominator <- 1 + z_975^2 / n
wilson_center <- (sample_proportion + z_975^2 / (2*n)) / wilson_denominator
wilson_half_width <- z_975 * sqrt(sample_proportion*(1-sample_proportion)/n + z_975^2/(4*n^2)) / wilson_denominator
wilson_ci <- wilson_center + c(-1, 1) * wilson_half_width
wilson_check <- prop.test(x, n, correct = FALSE, conf.level = confidence)$conf.int
wilson_ci
# 6/13的95%区间：[0.2320607, 0.7085620]。这是近似方法，不是CP。
stopifnot(max(abs(wilson_ci - wilson_check)) < 1e-10)
# R默认correct=TRUE，实际使用时明确报告是否作连续性校正。
# 它处理用连续正态近似离散二项计数时的半格边界；CP不需要这一校正。
wilson_corrected <- prop.test(x, n, correct = TRUE, conf.level = confidence)$conf.int
wilson_corrected
# 默认6/13：校正区间[0.2040175,0.7387967]，比未校正版本宽。
# 有校正不等于所有情况下都更准确；本节先对照软件默认值，再回到CP主线。
wilson_correction_comparison <- data.frame(
  method=c("Wilson_without_continuity_correction", "Wilson_with_continuity_correction"),
  lower=c(wilson_check[1],wilson_corrected[1]),
  upper=c(wilson_check[2],wilson_corrected[2]))
write.csv(wilson_correction_comparison, file.path(out, "Wilson连续性校正对照.csv"), row.names=FALSE)

normal_wilson <- data.frame(
  quantity=c("z_critical", "normal_central_area", "wilson_low", "wilson_high"),
  value=c(z_975, normal_central_area, wilson_ci))
write.csv(normal_wilson, file.path(out, "正态分位点与Wilson.csv"), row.names=FALSE)
# 思考：Wald在x=0时的估计标准误是0；Wilson的上限仍大于0，为什么？

# 2. 固定输入文件中的x，扫描候选总体率p。
# 1-pbinom(x-1, n, p)计算P_p(X>=x)，包含观察到的x。
# pbinom(x, n, p)计算P_p(X<=x)，也包含观察到的x。
p_grid <- seq(0, 1, by = scan_step)
right_tail <- 1 - pbinom(x - 1, size = n, prob = p_grid)
upper_tail <- pbinom(x, size = n, prob = p_grid)
scan <- data.frame(
  candidate_p = p_grid,
  right_tail_P_X_ge_x = right_tail,
  upper_tail_P_X_le_x = upper_tail
)
write.csv(scan, file.path(out, "扫描表.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# 下限右尾随p增大；上限左尾随p增大而减小。
# x=0时，Pr_p(X>=0)=1，不可能等于0.025；下限直接是参数空间边界p=0。
# x=n时，Pr_p(X<=n)=1，不可能等于0.025；上限直接是参数空间边界p=1。
scan_brackets <- data.frame(
  bound = character(), status = character(), candidate_low = numeric(),
  tail_at_low = numeric(), candidate_high = numeric(),
  tail_at_high = numeric(), stringsAsFactors = FALSE
)
if (x == 0) {
  lower_below_i <- NA_integer_
  lower_above_i <- NA_integer_
  scan_brackets <- rbind(
    scan_brackets,
    data.frame(
      bound = "lower",
      status = "boundary fixed at p=0; Pr_p(X>=0)=1; no interior root",
      candidate_low = NA_real_, tail_at_low = 1,
      candidate_high = NA_real_, tail_at_high = 1,
      stringsAsFactors = FALSE
    )
  )
} else {
  lower_below_i <- max(which(right_tail <= tail_target))
  lower_above_i <- min(which(right_tail >= tail_target))
  scan_brackets <- rbind(
    scan_brackets,
    data.frame(
      bound = "lower", status = "interior root",
      candidate_low = p_grid[lower_below_i],
      tail_at_low = right_tail[lower_below_i],
      candidate_high = p_grid[lower_above_i],
      tail_at_high = right_tail[lower_above_i],
      stringsAsFactors = FALSE
    )
  )
}
if (x == n) {
  upper_low_i <- NA_integer_
  upper_high_i <- NA_integer_
  scan_brackets <- rbind(
    scan_brackets,
    data.frame(
      bound = "upper",
      status = "boundary fixed at p=1; Pr_p(X<=n)=1; no interior root",
      candidate_low = NA_real_, tail_at_low = 1,
      candidate_high = NA_real_, tail_at_high = 1,
      stringsAsFactors = FALSE
    )
  )
} else {
  upper_low_i <- max(which(upper_tail >= tail_target))
  upper_high_i <- min(which(upper_tail <= tail_target))
  scan_brackets <- rbind(
    scan_brackets,
    data.frame(
      bound = "upper", status = "interior root",
      candidate_low = p_grid[upper_low_i],
      tail_at_low = upper_tail[upper_low_i],
      candidate_high = p_grid[upper_high_i],
      tail_at_high = upper_tail[upper_high_i],
      stringsAsFactors = FALSE
    )
  )
}
write.csv(scan_brackets, file.path(out, "扫描夹逼.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# 3. 从扫描夹逼区间开始二分；前20轮保存为教学历史，继续计算到60轮得到高精度端点。
bisect_tail <- function(lower, upper, tail_fun, target, side, shown_rounds, total_rounds = 60) {
  history <- data.frame(
    bound = character(), status = character(), round = integer(), candidate_p = numeric(),
    tail_probability = numeric(), bracket_low_after = numeric(),
    bracket_high_after = numeric(), stringsAsFactors = FALSE
  )
  for (round_number in seq_len(total_rounds)) {
    candidate <- (lower + upper) / 2
    tail <- tail_fun(candidate)
    if (side == "lower") {
      # P_p(X>=x)随p增大；tail较小，端点在右侧。
      if (tail <= target) lower <- candidate else upper <- candidate
    } else {
      # P_p(X<=x)随p增大而减小；tail较大，端点在右侧。
      if (tail >= target) lower <- candidate else upper <- candidate
    }
    if (round_number <= shown_rounds) {
      history <- rbind(
        history,
        data.frame(
          bound = side, status = "interior root", round = round_number, candidate_p = candidate,
          tail_probability = tail, bracket_low_after = lower,
          bracket_high_after = upper, stringsAsFactors = FALSE
        )
      )
    }
  }
  list(endpoint = (lower + upper) / 2, history = history)
}

if (x == 0) {
  lower_endpoint <- 0
  lower_history <- data.frame(
    bound = "lower", status = "boundary fixed at p=0; no interior root",
    round = 0, candidate_p = 0, tail_probability = 1,
    bracket_low_after = NA_real_, bracket_high_after = NA_real_,
    stringsAsFactors = FALSE
  )
} else {
  lower_search <- bisect_tail(
    p_grid[lower_below_i], p_grid[lower_above_i],
    function(candidate) 1 - pbinom(x - 1, size = n, prob = candidate),
    tail_target, "lower", display_rounds
  )
  lower_endpoint <- lower_search$endpoint
  lower_history <- lower_search$history
}
if (x == n) {
  upper_endpoint <- 1
  upper_history <- data.frame(
    bound = "upper", status = "boundary fixed at p=1; no interior root",
    round = 0, candidate_p = 1, tail_probability = 1,
    bracket_low_after = NA_real_, bracket_high_after = NA_real_,
    stringsAsFactors = FALSE
  )
} else {
  upper_search <- bisect_tail(
    p_grid[upper_low_i], p_grid[upper_high_i],
    function(candidate) pbinom(x, size = n, prob = candidate),
    tail_target, "upper", display_rounds
  )
  upper_endpoint <- upper_search$endpoint
  upper_history <- upper_search$history
}
bisection_history <- rbind(lower_history, upper_history)
write.csv(bisection_history, file.path(out, "二分历史.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# beta分位点是Clopper–Pearson端点的另一种数值表达，用来核对二分结果。
cp_lower <- if (x == 0) 0 else qbeta(tail_target, shape1 = x, shape2 = n - x + 1)
cp_upper <- if (x == n) 1 else qbeta(1 - tail_target, shape1 = x + 1, shape2 = n - x)

# 4. 文章中的两个补充核对：9/10 vs p=0.7，以及0/10边界。
reference_x <- 9
reference_n <- 10
reference_p <- 0.7
observed_mass <- dbinom(reference_x, size = reference_n, prob = reference_p)
all_masses <- dbinom(0:reference_n, size = reference_n, prob = reference_p)
# R binom.test的双侧P值按不大于观察结果概率的结果求和。
reference_p_value <- sum(all_masses[all_masses <= observed_mass + 1e-14])
reference_ci_low <- qbeta(tail_target, reference_x, reference_n - reference_x + 1)
reference_ci_high <- qbeta(1 - tail_target, reference_x + 1, reference_n - reference_x)

zero_n <- 10
zero_two_sided_upper <- 1 - tail_target^(1 / zero_n)
zero_one_sided_upper <- 1 - alpha^(1 / zero_n)

key_results <- data.frame(
  quantity = c(
    "sample_proportion", "wald_low", "wald_high",
    "cp_lower_bisection_60", "cp_upper_bisection_60",
    "cp_lower_beta_check", "cp_upper_beta_check",
    "reference_9_of_10_p_value_at_0.7", "reference_9_of_10_cp_low",
    "reference_9_of_10_cp_high", "zero_of_10_two_sided_upper",
    "zero_of_10_one_sided_upper"
  ),
  value = c(
    sample_proportion, wald_low, wald_high, lower_endpoint, upper_endpoint,
    cp_lower, cp_upper, reference_p_value, reference_ci_low,
    reference_ci_high, zero_two_sided_upper, zero_one_sided_upper
  )
)
write.csv(key_results, file.path(out, "关键结果.csv"), row.names = FALSE, fileEncoding = "UTF-8")

reference_9_of_10 <- data.frame(
  quantity = c("two_sided_p_value_at_p0_0.7", "cp95_low", "cp95_high"),
  value = c(reference_p_value, reference_ci_low, reference_ci_high)
)
write.csv(reference_9_of_10, file.path(out, "9_10参考.csv"), row.names = FALSE, fileEncoding = "UTF-8")
zero_reference <- data.frame(
  quantity = c("two_sided_95_upper", "one_sided_95_upper"),
  value = c(zero_two_sided_upper, zero_one_sided_upper)
)
write.csv(zero_reference, file.path(out, "0_10边界.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# 主要结果对象直接留在工作区，便于在RStudio逐段运行查看；不使用cat/print。
validation_objects <- list(
  cp_interval = c(cp_lower, cp_upper),
  reference_p_value = reference_p_value,
  zero_of_10_upper = c(zero_two_sided_upper, zero_one_sided_upper)
)
if (x == 6 && n == 13) {
  stopifnot(abs(cp_lower - 0.1922324418) < 1e-8)
  stopifnot(abs(cp_upper - 0.7486545177) < 1e-8)
  stopifnot(abs(reference_p_value - 0.2995766785) < 1e-12)
}
