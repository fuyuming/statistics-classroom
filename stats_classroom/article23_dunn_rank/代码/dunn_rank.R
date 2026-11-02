# 统计课堂23：Dunn检验，从统一秩算到Holm
# 第12课33只小鼠完整存活天数；RStudio打开本项目后Source本文件。
# 需要rstatix；输出写入运行结果/R。先完整运行，再修改末尾练习。

# 1. 读数据，使用标准函数 ----
if (!file.exists("数据/小鼠存活天数.csv")) {
  stop("请先打开统计课堂23.Rproj，使工作目录为配套项目根目录。")
}
out <- "运行结果/R"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
save_table <- function(x, name) {
  write.csv(x, file.path(out, name), row.names = FALSE, fileEncoding = "UTF-8")
}
dat <- read.csv("数据/小鼠存活天数.csv", fileEncoding = "UTF-8")
dat$group <- factor(dat$group, levels = c("A", "B", "C"))
stopifnot(nrow(dat) == 33, !anyNA(dat), !anyDuplicated(dat$id))
print(dat)
print(table(dat$group))
# 先检查各组原始分布，再结合研究问题选择KW。
source("代码/check_distribution.R")
distribution_check <- check_distribution(dat)
print(distribution_check$normality)
kw_fit <- kruskal.test(days ~ group, data = dat)
print(kw_fit)
if (!requireNamespace("rstatix", quietly = TRUE)) {
  stop("请先安装rstatix：install.packages('rstatix')")
}
standard_dunn <- rstatix::dunn_test(
  dat, days ~ group, p.adjust.method = "holm", detailed = TRUE
)
print(standard_dunn)

# 2. 手算统一秩、平均秩差和Dunn ----
# 写成函数，以便后面重新计算“所有天数乘10”的结果。
calculate <- function(values, grouping) {
  ranks <- rank(values, ties.method = "average")
  labels <- c("A", "B", "C")
  n <- as.numeric(table(factor(grouping, levels = labels)))
  mean_rank <- as.numeric(tapply(ranks, grouping, mean)[labels])
  rank_sum <- n * mean_rank
  rank_variance <- var(ranks)  # 全体秩，分母N−1
  N <- length(values)
  counts <- as.numeric(table(values))
  tie_total <- sum(counts^3 - counts)
  between <- sum(n * (mean_rank - mean(ranks))^2)
  H <- between / rank_variance  # 平均并列秩使此式已含并列修正
  epsilon2 <- H / (N - 1)
  first <- c(1, 1, 2)
  second <- c(2, 3, 3)
  difference <- mean_rank[second] - mean_rank[first]
  se <- sqrt(rank_variance * (1 / n[first] + 1 / n[second]))
  z <- difference / se
  p <- 2 * pnorm(-abs(z))  # 对称正态两尾相加
  # Holm：排序、乘3/2/1、累计最大值、限制上限、还原顺序。
  ordering <- order(p)
  products <- (length(p):1) * p[ordering]
  adjusted_sorted <- pmin(1, cummax(products))
  adjusted <- numeric(length(p))
  adjusted[ordering] <- adjusted_sorted
  stopifnot(max(abs(adjusted - p.adjust(p, "holm"))) < 1e-12)
  pair <- c("B-A", "C-A", "C-B")
  list(
    ranks = ranks,
    groups = data.frame(group = labels, n, rank_sum, mean_rank),
    kw = data.frame(N, H, P = pchisq(H, 2, lower.tail = FALSE),
                    sr2 = rank_variance, tie_total, epsilon2),
    dunn = data.frame(pair, rank_diff = difference, SE = se, Z = z,
                      P_raw = p, P_holm = adjusted),
    holm = data.frame(pair = pair[ordering], P_raw = p[ordering],
                      multiplier = length(p):1, product = products,
                      cumulative_max = cummax(products), P_holm = adjusted_sorted)
  )
}
a <- calculate(dat$days, dat$group)
stopifnot(abs(a$kw$H - unname(kw_fit$statistic)) < 1e-10)
stopifnot(max(abs(a$dunn$P_raw - standard_dunn$p)) < 1e-10)
stopifnot(max(abs(a$dunn$P_holm - standard_dunn$p.adj)) < 1e-10)
stopifnot(max(abs(a$dunn$Z - standard_dunn$statistic)) < 1e-10)
dat$rank <- a$ranks
save_table(dat, "ranks.csv")
save_table(a$groups, "groups.csv")
save_table(a$kw, "kw.csv")
save_table(a$dunn, "dunn.csv")
save_table(a$holm, "holm_steps.csv")
counts <- table(dat$days)
save_table(data.frame(days = as.numeric(names(counts)), count = as.numeric(counts)), "ties.csv")
print(a$dunn)
cat("均秩差用秩位表示；P_holm已对三对比较调整。\n")

# 2.1 秩ε²对应秩上的单因素ANOVA η² ----
# 总平方和：所有秩相对全体平均秩的平方离差和。
rank_total_ss <- sum((dat$rank - mean(dat$rank))^2)
# 组间平方和：各组平均秩的偏离平方，乘该组只数后相加。
rank_between_ss <- sum(a$groups$n *
  (a$groups$mean_rank - mean(dat$rank))^2)
rank_eta2 <- rank_between_ss / rank_total_ss
# aov只在这里核对平方和比例；KW的P仍采用前面的卡方近似。
rank_aov <- aov(rank ~ group, data = dat)
rank_aov_ss <- summary(rank_aov)[[1]][["Sum Sq"]]
eta2_from_rank_aov <- rank_aov_ss[1] / sum(rank_aov_ss)
stopifnot(abs(rank_eta2 - a$kw$epsilon2) < 1e-12,
          abs(rank_eta2 - eta2_from_rank_aov) < 1e-12)
# 明确借用Field(2013)的η²分界作粗略量级参照。
# 极小/小/中等/大的分界为.01/.06/.14；并非通用的生物学意义阈值。
size_label <- as.character(cut(rank_eta2,
  breaks = c(-Inf, .01, .06, .14, Inf), right = FALSE,
  labels = c("very small", "small", "medium", "large")))
effect_definition <- data.frame(
  rank_SS_between = rank_between_ss, rank_SS_total = rank_total_ss,
  rank_epsilon2 = a$kw$epsilon2, eta2_on_ranks = eta2_from_rank_aov,
  sample_rank_variation_percent = 100 * rank_eta2,
  interpretation_field2013 = size_label)
save_table(effect_definition, "rank_effect_definition.csv")
print(effect_definition)

# 3. 两两Wilcoxon：每对重新编秩，显式关闭连续性校正 ----
first <- c("A", "A", "B")
second <- c("B", "C", "C")
w_results <- data.frame(pair = a$dunn$pair, P_no_continuity = NA_real_)
for (i in 1:3) {
  x <- dat$days[dat$group == first[i]]
  y <- dat$days[dat$group == second[i]]
  test <- wilcox.test(x, y, exact = FALSE, correct = FALSE)
  w_results$P_no_continuity[i] <- test$p.value
}
w_results$P_holm <- p.adjust(w_results$P_no_continuity, "holm")
save_table(w_results, "pairwise_wilcoxon.csv")
ab <- dat[dat$group %in% c("A", "B"), ]
ab$rank_ab <- rank(ab$days)
save_table(ab, "ab_rerank.csv")
stopifnot(all(ab$rank_ab[ab$days == 4] == 7))
print(w_results)

# 4. 固定全体秩，随机分为A11、B10、C12，观察均秩差的波动 ----
set.seed(20261103)
B_permutation <- 10000
difference_sim <- numeric(B_permutation)
for (b in seq_len(B_permutation)) {
  shuffled <- sample(a$ranks)  # 不放回打乱33个秩
  mean_A <- mean(shuffled[1:11])
  mean_B <- mean(shuffled[12:21])
  difference_sim[b] <- mean_B - mean_A
}
save_table(data.frame(difference = difference_sim), "permutation_differences.csv")
permutation_summary <- data.frame(
  B = B_permutation, seed = 20261103,
  mean = mean(difference_sim), SD = sd(difference_sim),
  theoretical_SD = a$dunn$SE[1]
)
save_table(permutation_summary, "permutation_summary.csv")
print(permutation_summary)

# 5. 组内有放回bootstrap：估计整体秩效应的不确定性 ----
# 每轮保持11/10/12只，重新合并编秩；区间不能用来替代KW检验。
set.seed(20261102)
B_bootstrap <- 20000
epsilon_sim <- numeric(B_bootstrap)
labels <- c("A", "B", "C")
bootstrap_group <- factor(rep(labels, c(11, 10, 12)), levels = labels)
for (b in seq_len(B_bootstrap)) {
  draw <- numeric(0)
  for (label in labels) {
    original <- dat$days[dat$group == label]
    draw <- c(draw, sample(original, length(original), replace = TRUE))
  }
  H_draw <- kruskal.test(draw, bootstrap_group)$statistic
  epsilon_sim[b] <- unname(H_draw) / 32
}
interval <- quantile(epsilon_sim, c(0.025, 0.975), type = 7)
bootstrap_summary <- data.frame(
  B = B_bootstrap, seed = 20261102, epsilon2_observed = a$kw$epsilon2,
  percentile_2_5 = unname(interval[1]), percentile_97_5 = unname(interval[2])
)
save_table(data.frame(epsilon2 = epsilon_sim), "bootstrap_epsilon.csv")
save_table(bootstrap_summary, "bootstrap_summary.csv")
print(bootstrap_summary)

# 6. 自我练习：乘10后重新计算全部结果 ----
changed <- calculate(dat$days * 10, dat$group)
stopifnot(identical(a$ranks, changed$ranks))
stopifnot(max(abs(a$dunn$P_holm - changed$dunn$P_holm)) < 1e-12)
stopifnot(abs(a$kw$epsilon2 - changed$kw$epsilon2) < 1e-12)
save_table(data.frame(pair = changed$dunn$pair,
  P_raw_x10 = changed$dunn$P_raw, P_holm_x10 = changed$dunn$P_holm,
  H_x10 = changed$kw$H, epsilon2_x10 = changed$kw$epsilon2), "unit_check.csv")
cat("乘10后秩、P和秩效应不变；原始天数的描述尺度变为10倍。\n")
writeLines(c(R.version.string, paste("rstatix", packageVersion("rstatix"))), file.path(out, "versions.txt"))
