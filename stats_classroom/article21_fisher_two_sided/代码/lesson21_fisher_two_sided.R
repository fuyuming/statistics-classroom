# 统计课堂21：先用标准函数，再拆解Fisher双侧P值。
# 课堂例7-4：33名HBsAg阳性孕妇，两组新生儿感染结局的汇总计数。
# RStudio中打开本篇Rproj，Source全文；基础R即可，结果存于运行结果/R。

# 一、读题并检查计数，先用fisher.test检验 -----------------------
root <- getwd()
if (basename(root) == "代码") root <- dirname(root)
input <- file.path(root, "数据", "hbv_prevention.csv")
if (!file.exists(input)) stop("请打开本篇Rproj，把工作目录设为配套根目录。")
data <- read.csv(input)
print(data)
O <- as.matrix(data[, c("positive", "negative")])
rownames(O) <- data$group
if (!is.numeric(O) || !all(dim(O) == c(2, 2))) stop("需要2×2数值计数表。")
if (any(!is.finite(O)) || any(O < 0) || any(O != floor(O))) {
  stop("计数必须是有限的非负整数，缺失值不能作为零。")
}
if (any(rowSums(O) == 0) || any(colSums(O) == 0)) stop("每行和每列合计须大于0。")
# two.sided指定双侧；默认按各表概率与观察表概率比较。
fit <- fisher.test(O, alternative = "two.sided", conf.level = 0.95)
print(fit)
rows <- rowSums(O)
cols <- colSums(O)
N <- sum(O)
a <- O[1, 1]

# 二、逐表计算：一格决定其余三格 -------------------------------
lower <- max(0, rows[1] - cols[2])
upper <- min(rows[1], cols[1])
x <- lower:upper
# dhyper：x是第一组阳性数，m总阳性、n总阴性、k第一组人数。
q <- dhyper(x, m = cols[1], n = cols[2], k = rows[1])
q_obs <- q[x == a]
# choose给出组合数，与标准函数逐项核对。
q_formula <- choose(cols[1], x) * choose(cols[2], rows[1] - x)
q_formula <- q_formula / choose(N, rows[1])
stopifnot(abs(sum(q) - 1) < 1e-12)
stopifnot(max(abs(q - q_formula)) < 1e-12)
print(data.frame(X = x, probability = q))
cat("10种表的概率和（本题）=", sum(q), "\n")

# 三、分别计算两种单侧与两种双侧规则 ---------------------------
included <- q <= q_obs * (1 + 1e-7)
p_manual <- sum(q[included])
p_left <- sum(q[x <= a])
p_right <- sum(q[x >= a])
p_double <- min(1, 2 * min(p_left, p_right))
# 等概率表要包含；1e-7是浮点误差容差，与R概率排序实现一致。
stopifnot(abs(p_manual - fit$p.value) < 1e-12)
stopifnot(abs(p_left - fisher.test(O, alternative = "less")$p.value) < 1e-12)
stopifnot(abs(p_right - fisher.test(O, alternative = "greater")$p.value) < 1e-12)
print(x[included])
print(c(left = p_left, right = p_right,
        double_min = p_double, Fisher = p_manual))

# 四、效应量和区间：样本OR与条件OR分别保存 ----------------------
# 样本OR是交叉乘积之比；零格可能产生0或Inf。
OR_sample <- O[1, 1] * O[2, 2] / (O[1, 2] * O[2, 1])
OR_conditional <- unname(fit$estimate)
OR_interval <- unname(fit$conf.int)
rates <- O[, 1] / rows
print(rates)
print(c(sample_OR = OR_sample, conditional_OR = OR_conditional,
        CI_low = OR_interval[1], CI_high = OR_interval[2]))
# 默认区间是等尾区间，对应central双侧规则；P值是概率排序定义。
# 报告时注明P与区间方法，不把区间误写为率差区间。

# 五、交换行或列练习：复用标准函数 -----------------------------
analyze_variant <- function(current, label) {
  test <- fisher.test(current, alternative = "two.sided")
  left <- fisher.test(current, alternative = "less")$p.value
  right <- fisher.test(current, alternative = "greater")$p.value
  sample_OR <- current[1, 1] * current[2, 2]
  sample_OR <- sample_OR / (current[1, 2] * current[2, 1])
  data.frame(variant = label, fisher_P = test$p.value,
             left_P = left, right_P = right,
             OR_sample = sample_OR, OR_conditional = unname(test$estimate),
             OR_low = test$conf.int[1], OR_high = test$conf.int[2])
}
original <- analyze_variant(O, "original")
swap_rows <- analyze_variant(O[2:1, ], "swap_rows")
swap_columns <- analyze_variant(O[, 2:1], "swap_columns")
exercise <- rbind(original, swap_rows, swap_columns)
print(exercise)
stopifnot(max(abs(exercise$fisher_P - fit$p.value)) < 1e-12)
cat("交换行或列：双侧P不变，左右单侧互换；OR和区间改变比较方向。\n")

# 六、保存可供核对的三张结果表 ---------------------------------
out <- file.path(root, "运行结果", "R")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
enumeration <- data.frame(
  X = x, prevention_positive = x, prevention_negative = rows[1] - x,
  control_positive = cols[1] - x, control_negative = rows[2] - cols[1] + x,
  probability = q, included = included
)
summary <- data.frame(
  N = N, row1 = rows[1], row2 = rows[2], col1 = cols[1], col2 = cols[2],
  observed_X = a, q_observed = q_obs,
  left_P = p_left, right_P = p_right, double_min_P = p_double,
  fisher_P = p_manual, OR_sample = OR_sample, OR_conditional = OR_conditional,
  OR_low = OR_interval[1], OR_high = OR_interval[2], rate1 = rates[1], rate2 = rates[2]
)
write.csv(enumeration, file.path(out, "enumeration.csv"), row.names = FALSE)
write.csv(summary, file.path(out, "summary.csv"), row.names = FALSE)
write.csv(exercise, file.path(out, "swap_exercise.csv"), row.names = FALSE)
