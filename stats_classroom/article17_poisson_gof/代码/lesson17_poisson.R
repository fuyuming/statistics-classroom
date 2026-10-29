# 统计课堂17：泊松分布拟合优度；UTF-8。
# 四部分顺序：看数据，查条件，计算检验，解释结果。课堂 CLABSI 例题为教学演示数据。
# 在项目根目录打开 Rproj 后 Source；也可从第一节起逐段运行。

# ---- 第一部分：看数据 → 明问题 → 立假设 ----
# 基础R即可。工作目录为本项目根；$按字段名提取一列。
project_dir <- normalizePath(getwd(), mustWork = TRUE)
data_path <- file.path(project_dir, "数据", "poisson_weekly_counts.csv")
if (!file.exists(data_path)) stop("请先打开本项目Rproj，再从第一部分运行。")
weekly_data <- read.csv(data_path, fileEncoding = "UTF-8")
weekly_data
weekly_count <- weekly_data$weekly_count
frequency <- weekly_data$frequency
N <- sum(frequency)
total_events <- sum(weekly_count * frequency)
raw_mean <- total_events / N
c(N = N, total_events = total_events, raw_mean = raw_mean)
# raw_mean是完整计数的均数，只作引入和初查；正式检验用分组MLE。
# 假定各周独立、暴露机会可比；H0：各周来自同一Poisson，lambda未知。

# ---- 第二部分：选方法 → 查条件 ----
category <- c("0", "1", "2", ">=3")
observed <- c(frequency[1], frequency[2], frequency[3], sum(frequency[4:5]))
# R的索引从1开始；最后一类包含原表中3次、4次的15+4周。
# dpois为恰好某次数的概率；ppois(...,lower.tail=FALSE)为右尾概率。
# 未观察到>=5不表示其理论概率为0；不能仅保留0至4后归一化。
expected_initial <- N * c(dpois(0:4, raw_mean),
                          ppois(4, raw_mean, lower.tail = FALSE))
initial_check <- data.frame(count=c("0","1","2","3","4",">=5"),
                            O=c(frequency,0),E=round(expected_initial,4))
initial_check
# 0/1/2的E约216.22/149.65/51.79，保留；4与>=5合并仍不足5，所以尾部合为>=3。

# ---- 第三部分：做检验 → 判显著 → 看效应 ----
negative_log_likelihood <- function(lambda) {
  probability <- c(dpois(0, lambda), dpois(1, lambda), dpois(2, lambda),
                   ppois(2, lambda, lower.tail = FALSE))
  -sum(observed * log(probability))
}
# 按四组频数最大化似然；负号把最大化变成最小化。
# 搜索区间0.001至5针对本题，不是置信区间。
estimate <- optimize(negative_log_likelihood, interval = c(0.001, 5), tol = 1e-10)
lambda_mle <- estimate$minimum
probability_mle <- c(dpois(0, lambda_mle), dpois(1, lambda_mle), dpois(2, lambda_mle),
                     ppois(2, lambda_mle, lower.tail = FALSE))
expected_mle <- N * probability_mle
all(expected_mle >= 5)
# 本例各类期望>=5，通过常用经验检查；这不能检验各周独立性。
contribution_mle <- (observed - expected_mle)^2 / expected_mle
Q_mle <- sum(contribution_mle)
naive <- chisq.test(observed, p = probability_mle)
naive
# 上述函数把概率视为事先给定，默认df=3，未扣除估计lambda。
df_naive <- unname(naive$parameter)
P_naive <- naive$p.value
df_mle <- length(observed) - 1 - 1
# 两个-1分别对应总数约束与从同一四组数据估计的一个lambda。
P_mle <- pchisq(Q_mle, df = df_mle, lower.tail = FALSE)
# lower.tail=FALSE取卡方右尾；本题不使用2x2列联表Yates校正。
result_mle <- data.frame(N = N, raw_mean = raw_mean, lambda = lambda_mle, Q = Q_mle,
                         df = df_mle, P = P_mle, df_naive = df_naive, P_naive = P_naive)
contribution_table <- data.frame(category, observed, probability = probability_mle,
                                 expected = expected_mle, contribution = contribution_mle)
result_mle
contribution_table

# ---- 第四部分：讲结果 → 做参数练习 ----
# Q约2.401，df=2，P约0.301：未发现四组频数显著偏离拟合模型。
# P>0.05不证明Poisson，也不证明周间独立；尾组内部形状未被检验。
# 练习：若lambda在采样前已指定为0.7，应重新计算期望、Q及df=3的P。
lambda_given <- 0.7
probability_given <- c(dpois(0, lambda_given), dpois(1, lambda_given),
                       dpois(2, lambda_given), ppois(2, lambda_given, lower.tail = FALSE))
expected_given <- N * probability_given
Q_given <- sum((observed - expected_given)^2 / expected_given)
df_given <- length(observed) - 1
P_given <- pchisq(Q_given, df = df_given, lower.tail = FALSE)
given_table <- data.frame(category, expected_given)
given_table
exercise_result <- data.frame(
  model = c("estimated_lambda", "given_lambda"),
  lambda = c(lambda_mle, lambda_given),
  Q = c(Q_mle, Q_given), df = c(df_mle, df_given),
  P = c(P_mle, P_given)
)
exercise_result
# 可选：保存结果；改为FALSE则只保留工作空间中的命名对象。
save_outputs <- TRUE
output_dir <- file.path(project_dir, "运行结果", "R")
if (save_outputs) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(result_mle, file.path(output_dir, "main_result.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(contribution_table, file.path(output_dir, "expected_contributions.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(exercise_result, file.path(output_dir, "lambda_exercise.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(given_table, file.path(output_dir, "given_expected.csv"), row.names = FALSE, fileEncoding = "UTF-8")
}
