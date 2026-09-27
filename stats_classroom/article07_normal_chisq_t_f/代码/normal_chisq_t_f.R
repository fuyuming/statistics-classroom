# 明哥的微生物世界 · 统计课堂07 · R / RStudio
# 学习顺序：标准误 → 临界值 → 两组比较 → 回归 → 模拟 → 配图。
# 所有数据都是教学设定，不是真实实验。脚本与注释保存为UTF-8。
# 双击项目的“统计课堂07.Rproj”，打开本文件，点击Source。
# 只需基础R，不安装额外包。绘图细节单独放在plot_distributions.R。

# ---------- 0. 准备项目目录 ----------
# 完整Source或Rscript运行时自动定位；逐行运行时以Rproj目录为起点。
project_dir <- ""  # 自动定位失败时，在引号内填写项目目录的完整路径。
args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
source_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (!nzchar(project_dir)) {
  if (length(file_arg)) {
    project_dir <- dirname(dirname(normalizePath(file_arg[1])))
  } else if (!is.null(source_file)) {
    project_dir <- dirname(dirname(normalizePath(source_file)))
  } else {
    project_dir <- getwd()
  }
}
if (!file.exists(file.path(project_dir, "数据", "teaching_data.csv"))) {
  stop("请先打开Rproj项目，再Source完整脚本；或填写project_dir。")
}
out <- file.path(project_dir, "运行结果")
figs <- file.path(project_dir, "文章配图")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(figs, recursive = TRUE, showWarnings = FALSE)

# 小工具：把课堂说明同时显示并留存，最后写成运行记录。
messages <- character()
explain <- function(message) {
  cat(message, "\n")
  messages <<- c(messages, message)
}

# ---------- 1. 样本均数的波动，与个体波动一样吗？ ----------
sigma <- 2  # 题设给定总体标准差，不是从这批数据估计的。
n <- 10     # 10个独立观测。
standard_error <- sigma / sqrt(n)
explain("问题1：总体标准差为2，独立样本量10，均数的标准误是多少？")
explain(sprintf("答案：2 / sqrt(10) = %.6f。标准误描述均数的抽样波动。", standard_error))

# ---------- 2. t界值与F界值，为什么能对上？ ----------
# q开头的函数计算分位点；p开头计算累计概率。
# 双侧alpha=0.05，右边留0.025，所以取左侧累计概率0.975的分位点。
nu <- 10
alpha <- 0.05
t_cutoff <- qt(1 - alpha / 2, df = nu)
f_cutoff <- qf(1 - alpha, df1 = 1, df2 = nu)
# lower.tail=FALSE算右尾P(T>1.96)；对称分布双侧概率再乘2。
wrong_error <- 2 * pt(1.96, df = nu, lower.tail = FALSE)
explain("\n问题2：t(10)双侧5%界值平方，与F(1,10)上侧5%界值相同吗？")
explain(sprintf("t=%.7f；t²=%.7f；F界值=%.7f。", t_cutoff, t_cutoff^2, f_cutoff))
explain(sprintf("错用±1.96的第一类错误率=%.2f%%。这是零假设下的理论概率。", 100 * wrong_error))

# 再看多种自由度；R函数会按向量逐个计算，不需要手工重复五次。
dfs <- c(1, 5, 10, 30, 100)
q <- qt(0.975, df = dfs)
critical <- data.frame(
  df = dfs,
  n_one_sample = dfs + 1,
  t_975 = q,
  t_975_squared = q^2,
  f_95 = qf(0.95, df1 = 1, df2 = dfs),
  rejection_at_1_96 = 2 * pt(1.96, df = dfs, lower.tail = FALSE)
)
explain(paste(capture.output(print(critical)), collapse = "\n"))
write.csv(critical, file.path(out, "critical_values.csv"), row.names = FALSE)
# stopifnot核对程序结果；1e-8是允许的数值误差，不是显著性水平。
stopifnot(max(abs(critical$t_975_squared - critical$f_95)) < 1e-8)

# ---------- 3. 同一份数据：两组t检验、方差分析、组别回归 ----------
data <- read.csv(file.path(project_dir, "数据", "teaching_data.csv"),
                 fileEncoding = "UTF-8", stringsAsFactors = FALSE)
# sample_id是虚构编号；group是A/B；x为教学自变量；response为教学响应。
# factor指定组别顺序，A作参考组，回归系数groupB表示B-A。
data$group <- factor(data$group, levels = c("A", "B"))
explain("\n问题3：两组均数比较，t平方是否等于方差分析F？")
explain(paste(capture.output(head(data)), collapse = "\n"))
stopifnot(!anyNA(data), all(is.finite(data$response)))

# response ~ group表示用组别解释响应；var.equal=TRUE选择合并方差t检验。
# t.test默认是双侧；这里不能换成Welch后还照搬同一等价关系。
t_result <- t.test(response ~ group, data = data, var.equal = TRUE)
group_model <- lm(response ~ group, data = data)
anova_result <- anova(group_model)
coefficient_table <- coef(summary(group_model))
group_t <- coefficient_table["groupB", "t value"]
comparison <- data.frame(
  method = c("pooled_t_squared", "one_way_ANOVA_F", "regression_group_t_squared"),
  statistic = c(unname(t_result$statistic)^2, anova_result$`F value`[1], group_t^2),
  p_value = c(t_result$p.value, anova_result$`Pr(>F)`[1], coefficient_table["groupB", "Pr(>|t|)"])
)
explain(paste(capture.output(print(comparison)), collapse = "\n"))
explain("解释：三个统计量约为19.7388，P值相同。t符号随相减顺序改变，平方后不受影响。")
write.csv(comparison, file.path(out, "t_anova_regression_equivalence.csv"), row.names = FALSE)
stopifnot(diff(range(comparison$statistic)) < 1e-10)
stopifnot(diff(range(comparison$p_value)) < 1e-10)

# ---------- 4. 一元回归：系数除以自己的标准误 ----------
# 另拟合response~x，含截距与一个连续自变量，不是同时控制group。
regression_model <- lm(response ~ x, data = data)
model_summary <- summary(regression_model)
coefficient_table <- coef(model_summary)
slope <- coefficient_table["x", "Estimate"]
slope_se <- coefficient_table["x", "Std. Error"]
slope_t <- slope / slope_se  # 零假设：总体斜率=0。
overall_f <- unname(model_summary$fstatistic["value"])
regression_table <- data.frame(
  slope_t = slope_t,
  slope_t_squared = slope_t^2,
  overall_F = overall_f,
  residual_df = df.residual(regression_model)
)
explain("\n问题4：一元回归的斜率t²和模型整体F怎样对应？")
explain(paste(capture.output(print(regression_table)), collapse = "\n"))
explain("两者约为119.3944；多元回归中单个系数t²不能直接对应整体F。")
write.csv(regression_table, file.path(out, "simple_regression_equivalence.csv"), row.names = FALSE)
stopifnot(abs(slope_t^2 - overall_f) < 1e-10)

# ---------- 5. 从独立标准正态变量“造”出χ²、t、F ----------
set.seed(20260927)  # 本语言中固定随机结果，跨语言不保证生成相同序列。
repetitions <- 50000
nu <- 10
# 每行10个独立标准正态变量；rowSums对每行求和。
normal_values <- matrix(rnorm(repetitions * nu), nrow = repetitions, ncol = nu)
chi_values <- rowSums(normal_values^2)
# 再独立抽一批Z。不能取上一行中的一个变量当作独立分子。
z_values <- rnorm(repetitions)
t_values <- z_values / sqrt(chi_values / nu)
f_values <- z_values^2 / (chi_values / nu)
stopifnot(max(abs(t_values^2 - f_values)) < 1e-10)
probabilities <- c(0.025, 0.5, 0.975)
simulation_summary <- data.frame(
  probability = probabilities,
  empirical_t = as.numeric(quantile(t_values, probabilities)),
  theoretical_t = qt(probabilities, df = nu)
)
explain("\n问题5：模拟得到的t分位点，是否接近理论t(10)？")
explain(paste(capture.output(print(simulation_summary)), collapse = "\n"))
explain("经验分位点允许抽样误差。模拟帮助理解，不能代替数学证明。")
write.csv(simulation_summary, file.path(out, "construction_simulation_quantiles.csv"), row.names = FALSE)
simulation_data <- data.frame(z = z_values, u_chisq10 = chi_values, t10 = t_values, f1_10 = f_values)
connection <- gzfile(file.path(out, "construction_draws.csv.gz"), open = "wt")
write.csv(simulation_data, connection, row.names = FALSE)
close(connection)

# ---------- 6. 生成五张图，并保存运行记录 ----------
# 绘图文件使用上面设置的out、figs、dfs、q；先运行主脚本，不单独Source绘图文件。
source(file.path(project_dir, "代码", "plot_distributions.R"), encoding = "UTF-8")
explain("\n读图：01构造关系；02 t与正态；03卡方；04 F；05错用界值的错误率。")
explain("练习：把自由度10改为5或30，看临界值如何变化，再用软件核对。")
explain("全部自动核对通过；图像与曲线坐标已保存。不要把密度高度当作概率。")
writeLines(c(messages, capture.output(sessionInfo())), file.path(out, "逐项计算结果.txt"), useBytes = TRUE)
