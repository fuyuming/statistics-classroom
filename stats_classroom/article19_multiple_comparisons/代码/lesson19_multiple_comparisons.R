# 统计课堂19：从三组有效率到Bonferroni与Holm
# 数据为课堂/教材例题汇总计数。基础R即可；先打开本目录Rproj。
# 输入：数据/therapies.csv；输出：运行结果/R。全文按节顺序运行。

# 一、读题与数据 -------------------------------------------------
root <- getwd()
if (basename(root) == "代码") root <- dirname(root)
input <- file.path(root, "数据", "therapies.csv")
if (!file.exists(input)) stop("请先打开本篇Rproj，或将工作目录切到配套根目录。")
out <- file.path(root, "运行结果", "R")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
d <- read.csv(input)
print(d)  # physical=物理，drug=药物，topical=膏药。
a <- as.matrix(d[, c("effective", "ineffective")])
rownames(a) <- d$group
alpha <- 0.05

# 二、整体检验及三次两两检验 -----------------------------------
# 三组由不同患者组成；独立性由设计支持，不能从P值判断。
overall <- chisq.test(a, correct = FALSE)
print(overall$expected)
pairs <- list(c(1, 2), c(1, 3), c(2, 3))
names_pair <- c("physical-drug", "physical-topical", "drug-topical")
raw <- numeric(3)
expected_min <- numeric(3)
for (i in seq_along(pairs)) {
  subset <- a[pairs[[i]], , drop = FALSE]
  fit <- chisq.test(subset, correct = TRUE)
  raw[i] <- fit$p.value  # 尚未作多重校正，已经作Yates连续性校正。
  expected_min[i] <- min(fit$expected)
}

# 三、手算调整，再与R内置函数核对 -------------------------------
m <- length(raw)
bonf <- pmin(m * raw, 1)
ord <- order(raw)
sorted_p <- raw[ord]
factor <- m:1
middle <- sorted_p * factor
# cummax保留本行及之前所有乘积的最大值；pmin使结果不超过1。
sorted_holm <- pmin(cummax(middle), 1)
holm <- numeric(m)
holm[ord] <- sorted_holm  # 还原到原来的疗法对顺序。
stopifnot(isTRUE(all.equal(holm, p.adjust(raw, "holm"))))
stopifnot(isTRUE(all.equal(bonf, p.adjust(raw, "bonferroni"))))
result <- data.frame(pair = names_pair, raw_p_yates = raw,
                     bonferroni_p = bonf, holm_p = holm,
                     expected_min = expected_min)
print(result)
cat("0.05标准下：Bonferroni拒绝", sum(bonf <= alpha), "项；Holm拒绝", sum(holm <= alpha), "项。\n")
cat("这些P值不改变原来的样本有效率。\n")

# 四、累计最大值演示与自我练习 ---------------------------------
demo <- c(.010, .011, .012)  # 人工算法例子，不是疗法数据。
demo_middle <- sort(demo) * (3:1)
demo_final <- pmin(cummax(demo_middle), 1)
exercise <- data.frame(pair = names_pair,
                      bonf_at_001 = bonf <= .01, holm_at_001 = holm <= .01)
print(exercise)
# 调换输入顺序，按对应位置放回，结果应不变。
perm <- c(3, 1, 2)
restored <- numeric(m)
restored[perm] <- p.adjust(raw[perm], "holm")
stopifnot(isTRUE(all.equal(restored, holm)))
stopifnot(all(abs(demo_final - .03) < 1e-14))

# 五、保存各步骤，让手算可以对照 --------------------------------
write.csv(data.frame(statistic = unname(overall$statistic),
                     df = unname(overall$parameter), p_value = overall$p.value,
                     expected_min = min(overall$expected)),
          file.path(out, "overall_result.csv"), row.names = FALSE)
write.csv(result, file.path(out, "pair_results.csv"), row.names = FALSE)
write.csv(data.frame(rank = 1:m, pair_sorted = names_pair[ord],
                     raw_p_sorted = sorted_p, factor = factor, middle = middle,
                     holm_cummax = sorted_holm),
          file.path(out, "holm_intermediate.csv"), row.names = FALSE)
write.csv(data.frame(p_input = demo, factor = 3:1,
                     middle = demo_middle, final = demo_final),
          file.path(out, "holm_demo.csv"), row.names = FALSE)
write.csv(exercise, file.path(out, "alpha001_exercise.csv"), row.names = FALSE)
cat("完成：结果在运行结果/R。输入顺序置换核验通过。\n")

# 六、补充：率差和单项95%区间（未作多重比较调整） ----------------
# 物理—药物；报告效应大小与估计精度，不用普通区间替代整族判断。
x_pair <- a[1:2, 1]
n_pair <- rowSums(a[1:2, , drop = FALSE])
fit_ci <- prop.test(x_pair, n_pair, correct = TRUE, conf.level = .95)
delta <- x_pair[1] / n_pair[1] - x_pair[2] / n_pair[2]
ci_result <- data.frame(pair = "physical-drug", estimate_pp = 100 * unname(delta),
                        lower_pp = 100 * fit_ci$conf.int[1],
                        upper_pp = 100 * fit_ci$conf.int[2])
print(ci_result)
write.csv(ci_result, file.path(out, "effect_ci.csv"), row.names = FALSE)
