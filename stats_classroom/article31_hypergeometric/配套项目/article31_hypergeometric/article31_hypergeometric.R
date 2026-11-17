# 统计课堂31：超几何分布从哪里来？（UTF-8）
# 题目：33名HBsAg阳性孕妇随机分组，比较新生儿HBV感染率。
# 预防组感染4、未感染18；非预防组感染5、未感染6。
# 来源：颜艳、王彤《医学统计学》第5版例7-4，p111。
# 在项目根打开Rproj后，从第1节按顺序运行；绘图需ggplot2。

# 1. 读取数据与定义随机对象 ----
library(ggplot2)
SAVE_OUTPUT <- FALSE
# Rscript可从任意目录启动；RStudio项目的工作目录为项目根。
args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- args_all[grepl("^--file=", args_all)]
project_root <- getwd()
if (length(file_arg) > 0) {
  project_root <- dirname(normalizePath(sub("^--file=", "", file_arg[1])))
}
hbv <- read.csv(file.path(project_root, "data", "hbv_table.csv"))
hbv
# total是每组人数，infected是感染人数；按组名取值以保持方向。
N <- sum(hbv$total)
K <- sum(hbv$infected)
n <- hbv$total[hbv$group == "prevention"]
x_observed <- hbv$infected[hbv$group == "prevention"]
# 给定行列合计后，X表示预防组感染人数；观察值为4。

# 2. 手算观察表：先数分子，再数分母 ----
positive_choices <- choose(K, x_observed)
negative_choices <- choose(N - K, n - x_observed)
all_choices <- choose(N, n)
observed_probability <- positive_choices * negative_choices / all_choices
combination_steps <- data.frame(positive_choices, negative_choices,
                               all_choices, observed_probability)
combination_steps
# 分子126×134596，分母193536720，观察表概率约.08762728。

# 3. 手算完整分布、期望、方差和尾部 ----
# 下限保证阴性够抽；上限同时受阳性总数与抽取人数约束。
x <- max(0, n - (N - K)):min(n, K)
probability <- choose(K, x) * choose(N - K, n - x) / all_choices
probability_sum <- sum(probability)
mean_x <- sum(x * probability)
variance_x <- sum((x - mean_x)^2 * probability)
mean_formula <- n * K / N
variance_formula <- n * (K / N) * (1 - K / N) * (N - n) / (N - 1)
moments <- data.frame(probability_sum, mean_x, mean_formula,
                      variance_x, variance_formula)
moments
# 总概率1，均值6，方差1.5；方差含有限总体修正因子11/32。
left_tail <- sum(probability[x <= x_observed])
right_tail <- sum(probability[x >= x_observed])
tails <- data.frame(left_tail, right_tail)
tails
# 两尾均包含观察值，故不能相加得到1；此处只是两个事件概率。

# 4. 调用分布函数：核对手算 ----
# dhyper的m/n为总体阳性/阴性数；k才是抽取人数。
probability_function <- dhyper(x, m = K, n = N - K, k = n)
probability_table <- data.frame(x, probability, probability_function)
probability_table
left_function <- phyper(x_observed, m = K, n = N - K, k = n)
# lower.tail=FALSE返回P(X>q)，要包含4，所以输入4-1。
right_function <- phyper(x_observed - 1, m = K, n = N - K,
                         k = n, lower.tail = FALSE)
tail_check <- data.frame(left_tail, left_function, right_tail, right_function)
tail_check

# 5. 手算Fisher双侧P：按表概率筛选后相加 ----
# H0为两组总体感染率相同。1e-7为浮点相等容差。
included <- probability <= observed_probability * (1 + 1e-7)
probability_table$included <- included
probability_table
included_x <- x[included]
included_x
fisher_manual <- sum(probability[included])
# 纳入0、1、2、3、4、9，双侧P约.12104475。

# 6. 调用完整检验函数：核对同一张表 ----
prevention <- hbv[hbv$group == "prevention", ]
control <- hbv[hbv$group == "control", ]
tab <- rbind(Prevention = c(prevention$infected, prevention$uninfected),
             Control = c(control$infected, control$uninfected))
colnames(tab) <- c("Infected", "Uninfected")
fisher_fit <- fisher.test(tab, alternative = "two.sided", conf.int = FALSE)
fisher_check <- data.frame(manual = fisher_manual, function_P = fisher_fit$p.value)
fisher_check
# alpha=.05时不拒绝H0，尚缺足够差异证据；观察率差仍为-27.27个百分点。

# 7. 绘图：横轴人数，纵轴为每个计数的概率 ----
plot_x <- ggplot(probability_table, aes(x = x, y = probability, fill = included)) +
  geom_col(width = .72) +
  geom_vline(xintercept = mean_x, linetype = "dashed") +
  scale_fill_manual(values = c("FALSE" = "#6192C3", "TRUE" = "#D76B26")) +
  scale_x_continuous(breaks = x) +
  labs(x = "X: infected in prevention group", y = "Probability P(X = x)",
       fill = "In two-sided P") +
  theme_minimal(base_size = 13)
plot_x
# 橙柱计入双侧P，虚线是期望6；绘图不重新估计概率。

# 8. 练习：抽11人，比较留下的22人 ----
draws_exercise <- 11
probability_exercise <- dhyper(x, m = K, n = N - K, k = draws_exercise)
remaining_positive <- K - x
exercise_mean <- sum(x * probability_exercise)
exercise_variance <- sum((x - exercise_mean)^2 * probability_exercise)
remaining_mean <- sum(remaining_positive * probability_exercise)
remaining_variance <- sum((remaining_positive - remaining_mean)^2 * probability_exercise)
exercise <- data.frame(exercise_mean, exercise_variance, remaining_mean, remaining_variance)
exercise
# 两个阳性数之和恒为9，期望3和6，方差均为1.5。

# 9. 可选保存：将第1节SAVE_OUTPUT改为TRUE ----
if (SAVE_OUTPUT) {
  out <- file.path(project_root, "results", "r")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  write.csv(probability_table, file.path(out, "probability_table.csv"), row.names = FALSE)
  write.csv(moments, file.path(out, "moments.csv"), row.names = FALSE)
  write.csv(fisher_check, file.path(out, "fisher_check.csv"), row.names = FALSE)
  ggsave(file.path(out, "hypergeometric_distribution.png"), plot_x, width = 7, height = 5, dpi = 180)
}
