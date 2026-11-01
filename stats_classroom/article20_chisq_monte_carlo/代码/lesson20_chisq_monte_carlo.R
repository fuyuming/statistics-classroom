# 统计课堂20：一张模拟表怎样变成Monte Carlo P值？
# 输入为课堂249人汇总表。先在RStudio打开本篇Rproj，再Source完整脚本。
# 基础R即可；结果保存到运行结果/R。按一至六节顺序运行。

# 一、读原题，计算每格期望和Pearson统计量 ----------------------
root <- getwd()
if (basename(root) == "代码") root <- dirname(root)
input <- file.path(root, "数据", "ace_genotype.csv")
if (!file.exists(input)) stop("请先打开本篇Rproj，工作目录应为配套根目录。")
out <- file.path(root, "运行结果", "R")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
data <- read.csv(input)
print(data)
O <- as.matrix(data[, c("DD", "ID", "II")])
rownames(O) <- data$group
rows <- rowSums(O)
cols <- colSums(O)
N <- sum(O)
# 标准函数先完成观察表检验；$取出期望表和χ²。
fit_obs <- chisq.test(O)
print(fit_obs)
E <- fit_obs$expected
Q <- unname(fit_obs$statistic)  # 去掉X-squared名称，保存数值。
# 逐格计算用于对照前面的数学公式，六格之和应与函数相同。
contribution <- (O - E)^2 / E
stopifnot(max(abs(E - outer(rows, cols) / N)) < 1e-12)
stopifnot(abs(sum(contribution) - Q) < 1e-12)
print(E)
print(contribution)
print(Q)
# 本表为2×3，不使用四格表Yates校正；右尾概率不乘2。
asymptoticP <- pchisq(Q, df = 2, lower.tail = FALSE)
threshold <- Q * (1 - 64 * .Machine$double.eps)
# 先转置再展开，按DN三格、noDN三格保存。
cells <- data.frame(group = rep(data$group, each = 3),
                    genotype = rep(colnames(O), 2),
                    observed = c(t(O)), expected = c(t(E)),
                    pearson_contribution = c(t(contribution)))
write.csv(cells, file.path(out, "cells.csv"), row.names = FALSE)

# 二、无放回抽111张卡片，生成一张模拟表 -------------------------
labels <- rep(1:3, times = cols)
set.seed(123)
picked <- sample(seq_along(labels), size = rows[1], replace = FALSE)
first_row <- tabulate(labels[picked], nbins = 3)
first <- rbind(first_row, cols - first_row)
print(first)
cat("这一张模拟表的χ²：", sum((first - E)^2 / E), "\n")
write.csv(data.frame(group = data$group, DD = first[, 1],
                     ID = first[, 2], II = first[, 3]),
          file.path(out, "one_table.csv"), row.names = FALSE)

# 三、批量生成B张固定边际的表，每张算同一个χ² -------------------
B <- 100000
set.seed(123)  # 从相同起点重新开始本批模拟。
sim <- r2dtable(B, rows, cols)
qs <- numeric(B)
for (i in seq_len(B)) {
  current <- sim[[i]]
  stopifnot(all(rowSums(current) == rows), all(colSums(current) == cols))
  stopifnot(all(current >= 0), all(current == floor(current)))
  qs[i] <- sum((current - E)^2 / E)
}
b <- sum(qs >= threshold)
P <- (b + 1) / (B + 1)
SE <- sqrt(P * (1 - P) / B)
# 必须重设种子，才是在同一起点与R内置函数核对。
set.seed(123)
builtinP <- chisq.test(O, simulate.p.value = TRUE, B = B)$p.value
stopifnot(abs(P - builtinP) < 1e-14)
write.csv(data.frame(Q = qs), file.path(out, "simulation.csv"), row.names = FALSE)
cat("极端次数：", b, "；模拟P：", P, "；模拟标准误：", SE, "\n")

# 四、穷举条件表，为模拟提供不含随机误差的校验基准 --------------
# 第一行前两格确定后，第三格由111人约束确定，第二行由列合计确定。
# 各表概率=三种基因型选法的组合数乘积/总选法数；用log避免大数溢出。
records <- list()
index <- 1
for (dd in 0:cols[1]) {
  for (id in 0:cols[2]) {
    ii <- rows[1] - dd - id
    if (ii >= 0 && ii <= cols[3]) {
      first_counts <- c(dd, id, ii)
      current <- rbind(first_counts, cols - first_counts)
      current_q <- sum((current - E)^2 / E)
      log_prob <- sum(lchoose(cols, first_counts)) - lchoose(N, rows[1])
      records[[index]] <- c(dd, id, ii, current_q, exp(log_prob))
      index <- index + 1
    }
  }
}
en <- as.data.frame(do.call(rbind, records))
names(en) <- c("DD", "ID", "II", "Q", "probability")
en$pearson_tail <- en$Q >= threshold
exact <- sum(en$probability[en$pearson_tail])
stopifnot(nrow(en) == 4063, abs(sum(en$probability) - 1) < 1e-10)
stopifnot(abs(P - exact) < 5 * SE)
write.csv(en, file.path(out, "enumeration.csv"), row.names = FALSE)
cat("固定边际、按Pearson χ²定义的精确条件尾概率：", exact, "\n")

# 五、效应量及近似95%区间：描述原始249人的关联 ------------------
# 2×3表的Cramér's V=sqrt(Q/N)。区间反演非中心χ²分布。
# 求非中心参数使观察Q分别落在累计概率.975与.025，再换成V尺度。
effect_v <- sqrt(Q / N)
low_ncp <- uniroot(function(x) pchisq(Q, 2, ncp = x) - .975,
                   c(0, 100), tol = 1e-12)$root
high_ncp <- uniroot(function(x) pchisq(Q, 2, ncp = x) - .025,
                    c(0, 100), tol = 1e-12)$root
ci_low <- sqrt(low_ncp / N)
ci_high <- sqrt(high_ncp / N)
summary <- data.frame(N = N, Qobs = Q, df = 2, asymptoticP = asymptoticP,
                       B = B, b = b, P = P, SE = SE, exact_pearson_P = exact,
                       effect_v = effect_v, ci_low = ci_low, ci_high = ci_high,
                       builtinP = builtinP)
print(summary)
write.csv(summary, file.path(out, "summary.csv"), row.names = FALSE)

# 六、取同一模拟序列的不同前缀，观察计算精度 --------------------
sizes <- sort(unique(c(2000, 10000, B)))
sizes <- sizes[sizes <= B]
exercise <- data.frame(B = sizes, b = 0, P = 0, SE = 0)
for (i in seq_along(sizes)) {
  count <- sum(qs[seq_len(sizes[i])] >= threshold)
  estimate <- (count + 1) / (sizes[i] + 1)
  exercise$b[i] <- count
  exercise$P[i] <- estimate
  exercise$SE[i] <- sqrt(estimate * (1 - estimate) / sizes[i])
}
write.csv(exercise, file.path(out, "precision_exercise.csv"), row.names = FALSE)
print(exercise)
cat("练习：改变种子再运行，比较P波动；增大B降低计算误差，N仍为249。\n")
cat("R_CHECKS_PASSED：边际、整数、概率总和、内置核对及Monte Carlo误差。\n")
