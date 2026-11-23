# 统计课堂38：多次二分类、两次多分类（UTF-8） ----
# 两例都是本篇假设教学计数，不是实测资料：例1、例2均为n=8。
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
if (!file.exists("data/transition_records.csv")) {
  stop("未找到data/transition_records.csv；请在RStudio打开Rproj并把工作目录设为本项目根目录。")
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
transition_records <- read.csv("data/transition_records.csv")
stopifnot(identical(names(transition_records), c("id", "before", "after")),
          length(unique(transition_records$id)) == nrow(transition_records),
          nrow(transition_records) == nrow(records),
          all(as.matrix(transition_records[c("before", "after")]) %in% c("A", "B", "C")))
transitions <- read.csv("data/transition_counts.csv")
cat("逐人记录前10条（id及A/B/C）：\n")
print(head(records, 10), row.names = FALSE)
cat("由逐人记录重汇总的三检测模式：\n")
print(patterns, row.names = FALSE)
print(transitions)
levels_abc <- c("A", "B", "C")
# Y只取三结果列；每人一行、每法8条对应记录。
stopifnot(!anyNA(Y), all(Y %in% c(0, 1)))
# 行为前，列为后，两轴采用完全相同的类别顺序。
M <- table(factor(transition_records$before, levels = levels_abc),
           factor(transition_records$after, levels = levels_abc))
M <- unclass(M)
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
  data.frame(b = b, c = c_count, rate_diff = rd, p_raw = p_raw)
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
# 仅适用于3×3，且三对反向格合计均>0；按固定反向合计的二项条件分布，
# 对Bowker统计量做右尾精确枚举（无蒙特卡罗）。
bowker_exact3 <- function(M) {
  stopifnot(identical(dim(M), c(3L, 3L)))
  use <- upper.tri(M)
  m <- (M + t(M))[use]
  stopifnot(length(m) == 3L, all(m > 0))
  observed <- sum((M - t(M))[use]^2 / m)
  configurations <- expand.grid(lapply(m, function(n) 0:n))
  statistic <- apply(configurations, 1, function(x) sum((2 * x - m)^2 / m))
  probability <- apply(configurations, 1, function(x) prod(dbinom(x, m, prob = 0.5)))
  list(statistic = observed, configurations = nrow(configurations),
       probability_sum = sum(probability),
       p.value = sum(probability[statistic >= observed - 1e-12]))
}

# 3 Cochran Q：公式与真实包调用 ----
n <- nrow(Y)
k <- ncol(Y)
C <- colSums(Y)
R <- rowSums(Y)
T <- sum(R)
D <- sum(R * (k - R))
Q <- k * (k - 1) * sum((C - T / k)^2) / D
q_fit <- DescTools::CochranQTest(Y)
cat("Q节教学输入：n=", n, "，k=", k, "；列和C=", paste(C, collapse = ","),
    "；T=", T, "，D=", D, "；df=", k - 1,
    "；渐近P=", format.pval(pchisq(Q, k - 1, lower.tail = FALSE)), "\n", sep = "")
print(q_fit)
q_function_comparison <- data.frame(
  test = "Cochran_Q", formula = Q,
  `function` = unname(q_fit$statistic),
  abs_error = abs(Q - q_fit$statistic))
print(q_function_comparison)
stopifnot(abs(Q - q_fit$statistic) < 1e-10)
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
person_assignments <- 3^n1
q_conditional <- data.frame(configurations, person_assignments, probability_sum, p_exact)
print(q_conditional)
stopifnot(abs(probability_sum - 1) < 1e-10)
cat("条件于每人行和；U/W分别是行和1/2对象的位置计数。\n")

# 5 全部三对方法：精确McNemar + Holm；效应分母n=8 ----
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
pairwise <- pairwise[c("comparison", "b", "c", "rate_diff", "p_raw", "holm")]
print(pairwise)
pairwise_teaching_check <- data.frame(comparison = pairwise$comparison,
                                      b = pairwise$b, c = pairwise$c,
                                      m = pairwise$b + pairwise$c)
print(pairwise_teaching_check)
stopifnot(sum(pairwise_teaching_check$m) == D)
cat("例1三对不一致数m=b+c之和=", sum(pairwise_teaching_check$m),
    "，与Q节D相同。\n", sep = "")
cat("rate_diff乘100是百分点；两例均以8人为分母。\n")

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
  method = c("Cochran_Q", "Stuart_Maxwell", "Bowker"),
  statistic = c(Q, SM, Bowker),
  df = c(k - 1, nrow(M) - 1, nrow(M) * (nrow(M) - 1) / 2))
# lower.tail=FALSE计算卡方分布右尾；均为渐近P。
key_results$p_value <- pchisq(key_results$statistic, key_results$df, lower.tail = FALSE)
print(key_results)

# 7 逐类别“此类/其余类”：分母为8人 ----
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
cat("此处率差为后减前。\n")

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
bowker_exact <- bowker_exact3(M)
enumeration <- data.frame(configurations = bowker_exact$configurations,
                          probability_sum = bowker_exact$probability_sum,
                          p_exact = bowker_exact$p.value)
print(enumeration)
stopifnot(abs(bowker_exact$probability_sum - 1) < 1e-10,
          abs(bowker_exact$statistic - Bowker) < 1e-10)
cat("这是自编的统计量右尾精确枚举；无蒙特卡罗。\n")

# 10 循环转移练习：另一组8人的M2 ----
M2 <- matrix(c(1,2,0, 0,1,2, 2,0,0), 3, 3, byrow = TRUE,
             dimnames = list(levels_abc, levels_abc))
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

# 可选：首次手动安装 install.packages("ANSM5") 后，用ANSM5穷举核对例1Q。
# 默认上限为100000；本例8×3需6^7=279936，故显式设为300000以保证穷举。
# 超过上限时可能输出pval.mc，不应视作穷举精确值。
# fit_exact <- ANSM5::cochran.q(
#   y = as.vector(Y),
#   groups = factor(rep(colnames(Y), each = nrow(Y))),
#   blocks = factor(rep(seq_len(nrow(Y)), times = ncol(Y))),
#   do.exact = TRUE, do.asymp = FALSE, max.exact.perms = 300000)
# fit_exact$pval.exact  # 例1实跑：0.00411522633744865
