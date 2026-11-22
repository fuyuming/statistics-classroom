# 统计课堂38：多次二分类、两次多分类（UTF-8） ----
# 两例都是本篇假设教学计数，不是实测资料：例1为n=60的三法二分类，例2为n=144的前后多分类。
# 对应问题：例1问三种检测的总体阳性率是否相同及配对差异；例2问前后边际是否同质、反向路径是否对称。
# RStudio打开本项目Rproj；从本节起依次运行。首次安装：install.packages("DescTools")。
# 默认只显示结果；要保存CSV，改SAVE_RESULTS为TRUE。每个语言用独立结果目录。
options(stringsAsFactors = FALSE)
SAVE_RESULTS <- identical(Sys.getenv("ARTICLE38_SAVE"), "1")
if (!file.exists("data/binary_records.csv")) {
  stop("未找到data/binary_records.csv；请在RStudio打开Rproj并把工作目录设为项目根目录。")
}
if (!file.exists("data/transition_counts.csv")) {
  stop("未找到data/transition_counts.csv；请在RStudio打开Rproj并把工作目录设为项目根目录。")
}
if (!requireNamespace("DescTools", quietly = TRUE)) {
  stop("缺少DescTools；请在项目根目录的Console运行 install.packages(\"DescTools\")。")
}

# 1 读取逐人记录，汇总模式与转移计数 ----
records <- read.csv("data/binary_records.csv")
stopifnot(identical(names(records), c("id", "A", "B", "C")),
          length(unique(records$id)) == nrow(records),
          !anyNA(records), all(as.matrix(records[c("A", "B", "C")]) %in% c(0, 1)))
# id标识受检者；Y只取三种方案的0/1结果，每行保留同一个人的对应关系。
Y <- as.matrix(records[c("A", "B", "C")])
# 每人贡献1，按A/B/C相同的组合求和，得到频数表。
patterns <- aggregate(list(count = rep(1L, nrow(records))),
                      by = records[c("A", "B", "C")], FUN = sum)
patterns$pattern <- paste0(patterns$A, patterns$B, patterns$C)
patterns <- patterns[c("pattern", "A", "B", "C", "count")]
patterns <- patterns[order(patterns$pattern), ]
rownames(patterns) <- NULL
transitions <- read.csv("data/transition_counts.csv")
cat("逐人记录前10条（id及A/B/C）：\n")
print(head(records, 10), row.names = FALSE)
cat("由逐人记录重汇总的三检测模式：\n")
print(patterns, row.names = FALSE)
print(transitions)
levels_abc <- c("A", "B", "C")
# Y只取三结果列；每人一行、每法60条对应记录。
stopifnot(!anyNA(Y), all(Y %in% c(0, 1)))
# 行为前，列为后，两轴采用完全相同的类别顺序。
M <- matrix(0, 3, 3, dimnames = list(levels_abc, levels_abc))
for (i in seq_len(nrow(transitions))) {
  M[transitions$before[i], transitions$after[i]] <- transitions$count[i]
}
print(head(Y))
print(M)

# 2 本文自编辅助函数：配对率差、SM、Bowker、方向区间 ----
# 函数定义放在第一次调用前；H0下的二项概率为1/2。
paired_summary <- function(b, c_count, n) {
  discordant <- b + c_count
  p_raw <- 1
  if (discordant > 0) {
    p_raw <- binom.test(c_count, discordant, p = 0.5)$p.value
  }
  rd <- (c_count - b) / n
  se <- sqrt((discordant - n * rd^2) / (n * (n - 1)))
  z95 <- qnorm(0.975)            # 双侧95%，每端留2.5%
  z_bonf <- qnorm(1 - 0.05 / 6) # 三个区间，每端留0.05/(2*3)
  data.frame(b = b, c = c_count, rate_diff = rd, se = se, p_raw = p_raw,
             lo95 = rd - z95 * se, hi95 = rd + z95 * se,
             lo_bonf = rd - z_bonf * se, hi_bonf = rd + z_bonf * se)
}
sm_formula <- function(tab) {
  r <- nrow(tab)
  row_total <- rowSums(tab)
  col_total <- colSums(tab)
  difference <- row_total - col_total
  # 所有类别差合计0，只取前r-1个。S由当前数据计算。
  S <- -(tab + t(tab))
  diag(S) <- row_total + col_total - 2 * diag(tab)
  keep <- seq_len(r - 1)
  S_reduced <- S[keep, keep, drop = FALSE]
  d_reduced <- difference[keep]
  statistic <- drop(t(d_reduced) %*% solve(S_reduced, d_reduced))
  list(statistic = statistic, S = S_reduced, d = d_reduced)
}
bowker_formula <- function(tab) {
  use <- upper.tri(tab)
  totals <- (tab + t(tab))[use]
  differences <- (tab - t(tab))[use]
  stopifnot(all(totals > 0)) # 本篇仅实现非零反向格对；零对另审有效df。
  sum(differences^2 / totals)
}

# 3 Cochran Q及Friedman：公式与真实包调用 ----
n <- nrow(Y)
k <- ncol(Y)
C <- colSums(Y)
R <- rowSums(Y)
T <- sum(R)
D <- sum(R * (k - R))
Q <- k * (k - 1) * sum((C - T / k)^2) / D
q_fit <- DescTools::CochranQTest(Y)
f_fit <- stats::friedman.test(Y)
cat("Q节教学输入：n=", n, "，k=", k, "；列和C=", paste(C, collapse = ","),
    "；T=", T, "，D=", D, "；df=", k - 1,
    "；渐近P=", format.pval(pchisq(Q, k - 1, lower.tail = FALSE)), "\n", sep = "")
print(q_fit)
print(f_fit)
# 用平均秩和并列项再核对，不能把未校正的Friedman量当Q。
ranks <- t(apply(Y, 1, rank, ties.method = "average"))
V <- colSums(ranks) - n * (k + 1) / 2
Ties <- sum(R^3 - R + (k - R)^3 - (k - R))
f_manual <- 12 * sum(V^2) / (n * k * (k + 1) - Ties / (k - 1))
q_function_comparison <- data.frame(
  test = c("Cochran_Q", "Friedman_tie_corrected"),
  formula = c(Q, f_manual),
  `function` = c(unname(q_fit$statistic), unname(f_fit$statistic)),
  abs_error = c(abs(Q - q_fit$statistic), abs(f_manual - f_fit$statistic)))
print(q_function_comparison)
stopifnot(abs(Q - q_fit$statistic) < 1e-10)
stopifnot(abs(Q - f_fit$statistic) < 1e-10, abs(Q - f_manual) < 1e-10)
cat("Q回答三法总体阳性概率是否相同；不说明准确性。\n")

# 4 给定每人行和的Cochran Q条件精确枚举 ----
# 本节专用于三法：U数唯一阳性的位置，W数唯一阴性的位置。
# 给定行和后，两组来自不同对象，分别服从多项分布；配置权重相乘。
n1 <- sum(R == 1)
n2 <- sum(R == 2)
n3 <- sum(R == 3)
observed_sum_squares <- sum(C^2)
probability_sum <- 0
p_exact <- 0
configurations <- 0
for (u1 in 0:n1) {
  for (u2 in 0:(n1 - u1)) {
    U <- c(u1, u2, n1 - u1 - u2)
    probability_U <- dmultinom(U, size = n1, prob = rep(1 / 3, k))
    for (v1 in 0:n2) {
      for (v2 in 0:(n2 - v1)) {
        W <- c(v1, v2, n2 - v1 - v2)
        probability_W <- dmultinom(W, size = n2, prob = rep(1 / 3, k))
        conditional_columns <- n3 + n2 + U - W
        probability <- probability_U * probability_W
        probability_sum <- probability_sum + probability
        configurations <- configurations + 1
        if (sum(conditional_columns^2) >= observed_sum_squares) {
          p_exact <- p_exact + probability
        }
      }
    }
  }
}
q_conditional <- data.frame(configurations, probability_sum, p_exact)
print(q_conditional)
stopifnot(abs(probability_sum - 1) < 1e-10)
cat("条件于每人行和；U/W分别是行和1/2对象的位置计数。\n")

# 5 全部三对方法：精确McNemar + Holm；效应分母n=60 ----
pair_index <- combn(seq_len(k), 2)
pair_rows <- list()
for (j in seq_len(ncol(pair_index))) {
  first <- pair_index[1, j]
  second <- pair_index[2, j]
  x <- Y[, first]
  y <- Y[, second]
  b <- sum(x == 1 & y == 0)
  c_count <- sum(x == 0 & y == 1)
  answer <- paired_summary(b, c_count, n)
  answer$comparison <- paste0(colnames(Y)[second], "-", colnames(Y)[first])
  pair_rows[[j]] <- answer
}
pairwise <- do.call(rbind, pair_rows)
pairwise$holm <- p.adjust(pairwise$p_raw, method = "holm")
print(pairwise)
pairwise_teaching_check <- data.frame(comparison = pairwise$comparison,
                                      b = pairwise$b, c = pairwise$c,
                                      m = pairwise$b + pairwise$c)
print(pairwise_teaching_check)
stopifnot(sum(pairwise_teaching_check$m) == D)
cat("例1三对不一致数m=b+c之和=", sum(pairwise_teaching_check$m),
    "，与Q节D相同。\n", sep = "")
cat("rate_diff乘100是百分点；lo95/hi95为逐项区间，lo_bonf/hi_bonf为三对同时区间。\n")

# 6 同一张前后方表：SM问边际，Bowker问反向路径 ----
sm <- sm_formula(M)
SM <- sm$statistic
Bowker <- bowker_formula(M)
print(sm)
sm_fit <- DescTools::StuartMaxwellTest(M)
bow_fit <- stats::mcnemar.test(M, correct = FALSE)
print(sm_fit)
print(bow_fit)
sm_function_comparison <- data.frame(
  test = c("Stuart_Maxwell", "Bowker"),
  formula = c(SM, Bowker),
  `function` = c(unname(sm_fit$statistic), unname(bow_fit$statistic)),
  abs_error = c(abs(SM - sm_fit$statistic), abs(Bowker - bow_fit$statistic)))
print(sm_function_comparison)
stopifnot(abs(SM - sm_fit$statistic) < 1e-10)
stopifnot(abs(Bowker - bow_fit$statistic) < 1e-10)
key_results <- data.frame(
  method = c("Cochran_Q", "Friedman", "Stuart_Maxwell", "Bowker"),
  statistic = c(Q, f_manual, SM, Bowker),
  df = c(k - 1, k - 1, nrow(M) - 1, nrow(M) * (nrow(M) - 1) / 2))
# lower.tail=FALSE计算卡方分布右尾；均为渐近P。
key_results$p_value <- pchisq(key_results$statistic, key_results$df, lower.tail = FALSE)
print(key_results)

# 7 逐类别“此类/其余类”：保持全体144人为分母 ----
category_rows <- list()
for (i in seq_len(nrow(M))) {
  leave <- sum(M[i, ]) - M[i, i]
  enter <- sum(M[, i]) - M[i, i]
  answer <- paired_summary(leave, enter, sum(M))
  answer$category <- levels_abc[i]
  category_rows[[i]] <- answer
}
categories <- do.call(rbind, category_rows)
categories$holm <- p.adjust(categories$p_raw, method = "holm")
print(categories)
cat("此处率差为后减前，逐项CI仍用全体完整配对人数。\n")

# 8 每条反向路径：条件方向概率与Clopper-Pearson区间 ----
path_rows <- list()
for (j in seq_len(ncol(pair_index))) {
  i1 <- pair_index[1, j]
  i2 <- pair_index[2, j]
  forward <- M[i1, i2]
  reverse <- M[i2, i1]
  fit <- binom.test(forward, forward + reverse, p = 0.5, conf.level = 0.95)
  path_rows[[j]] <- data.frame(
    path = paste0(levels_abc[i1], "->", levels_abc[i2]),
    forward = forward, reverse = reverse, total = forward + reverse,
    theta = forward / (forward + reverse), p_raw = fit$p.value,
    lo95 = fit$conf.int[1], hi95 = fit$conf.int[2])
}
paths <- do.call(rbind, path_rows)
paths$holm <- p.adjust(paths$p_raw, method = "holm")
print(paths)
cat("theta的分母只含该路径的互转者；CP区间为逐项95%，不是总体率差区间。\n")

# 9 自编条件枚举：固定反向格对总数，采用Bowker统计量右尾 ----
# 本例恰有三条路径；各方向在对称H0下条件独立二项(总数,1/2)。
sizes <- paths$total
probability_sum <- 0
p_exact <- 0
configurations <- 0
for (a in 0:sizes[1]) {
  for (b in 0:sizes[2]) {
    for (cc in 0:sizes[3]) {
      counts <- c(a, b, cc)
      statistic <- sum((2 * counts - sizes)^2 / sizes)
      probability <- prod(dbinom(counts, sizes, prob = 0.5))
      probability_sum <- probability_sum + probability
      configurations <- configurations + 1
      if (statistic >= Bowker - 1e-12) {
        p_exact <- p_exact + probability
      }
    }
  }
}
enumeration <- data.frame(configurations, probability_sum, p_exact)
print(enumeration)
stopifnot(abs(probability_sum - 1) < 1e-10)
cat("这是自编的统计量右尾枚举；与配置概率排序的exact输出不必相同。\n")

# 10 循环转移练习：保留原M，另作M2 ----
M2 <- M
M2[1, 2] <- 15
M2[3, 3] <- 30
exercise <- data.frame(SM = sm_formula(M2)$statistic, Bowker = bowker_formula(M2))
print(rowSums(M2))
print(colSums(M2))
print(exercise)

# 11 可选保存：不覆盖其他语言结果 ----
if (SAVE_RESULTS) {
  output_dir <- "results/R"
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  result_tables <- list(key_results = key_results, q_conditional = q_conditional,
                        pairwise = pairwise,
                        categories = categories, paths = paths,
                        enumeration = enumeration, exercise = exercise)
  for (name in names(result_tables)) {
    write.csv(result_tables[[name]], file.path(output_dir, paste0(name, ".csv")), row.names = FALSE)
  }
}
