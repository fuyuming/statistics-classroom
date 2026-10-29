# 统计课堂17：泊松分布拟合优度；UTF-8。
# 四部分按单元顺序执行；脚本不依赖 __file__，请先把当前工作目录设为项目根目录。

# %% 第一部分：看数据 → 明问题 → 立假设
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.stats import poisson, chi2

project_dir = Path.cwd()
data_path = project_dir / "数据" / "poisson_weekly_counts.csv"
if not data_path.is_file():
    raise FileNotFoundError("请将Spyder工作目录设为配套练习根目录，再从第一单元运行。")
weekly_data = pd.read_csv(data_path)
weekly_data
weekly_count = weekly_data["weekly_count"].to_numpy()
frequency = weekly_data["frequency"].to_numpy()
N = int(frequency.sum())
total_events = int((weekly_count * frequency).sum())
raw_mean = total_events / N
data_summary = {"N": N, "total_events": total_events, "raw_mean": raw_mean}
data_summary
# raw_mean 只作引入和初查；拟合使用四组数据估计 lambda。
# H0：各周来自同一Poisson(lambda)，lambda未知；要求周间独立、暴露机会可比。

# %% 第二部分：选方法 → 查条件
category = np.array(["0", "1", "2", ">=3"])
observed = np.array([frequency[0], frequency[1], frequency[2], frequency[3:].sum()])
# Python索引从0开始，frequency[3:]取第4个元素起直到末尾，即15+4周。
# pmf计算恰好某次数的概率；sf(4)计算大于4次、即>=5次的概率。
expected_initial = N * np.r_[poisson.pmf(np.arange(5), raw_mean), poisson.sf(4, raw_mean)]
initial_check = pd.DataFrame({"count":["0","1","2","3","4",">=5"],
                              "O":np.r_[frequency,0], "E":expected_initial.round(4)})
initial_check
# 0/1/2的E约216.22/149.65/51.79，保留；4与>=5合并仍不足5，所以尾部合为>=3。
# 固定0、1、2、>=3的分组；>=3理论尾部包含所有更高次数。

# %% 第三部分：做检验 → 判显著 → 看效应
def negative_log_likelihood(lambda_value):
    probability = np.array([poisson.pmf(0, lambda_value), poisson.pmf(1, lambda_value),
                            poisson.pmf(2, lambda_value), poisson.sf(2, lambda_value)])
    return -float((observed * np.log(probability)).sum())

# 最小化负对数似然。bounds是本题搜索范围，不是置信区间。
estimate = minimize_scalar(negative_log_likelihood, bounds=(0.001, 5), method="bounded",
                           options={"xatol": 1e-8})
lambda_mle = float(estimate.x)
probability_mle = np.array([poisson.pmf(0, lambda_mle), poisson.pmf(1, lambda_mle),
                            poisson.pmf(2, lambda_mle), poisson.sf(2, lambda_mle)])
expected_mle = N * probability_mle
expected_check = bool(np.all(expected_mle >= 5))
expected_check
# 期望频数>=5是常用经验检查，不能替代独立性等设计条件。
contribution_mle = (observed - expected_mle) ** 2 / expected_mle
Q_mle = float(contribution_mle.sum())
df_naive = len(observed) - 1
P_naive = float(chi2.sf(Q_mle, df_naive))
df_mle = len(observed) - 1 - 1
# 两个-1分别对应固定总数与从同一分组数据估计lambda。
# chi2.sf取右尾P值。df_naive=3模拟R chisq.test把p当成给定的默认行为。
P_mle = float(chi2.sf(Q_mle, df_mle))
result_mle = pd.DataFrame({"N": [N], "raw_mean": [raw_mean], "lambda": [lambda_mle], "Q": [Q_mle],
                           "df": [df_mle], "P": [P_mle], "df_naive": [df_naive], "P_naive": [P_naive]})
contribution_table = pd.DataFrame({"category": category, "observed": observed,
                                   "probability": probability_mle, "expected": expected_mle,
                                   "contribution": contribution_mle})
result_mle
contribution_table

# %% 第四部分：讲结果 → 做参数练习
# Q约2.401、df=2、P约0.301：没有发现四类频数显著偏离模型。
# 不证明Poisson或周间独立；本题也不直接套用2x2列联表Yates校正。
# 练习假定lambda在采样前给定；重新算E/Q，df=3，不从当前数据估计参数。
lambda_given = 0.7
probability_given = np.array([poisson.pmf(0, lambda_given), poisson.pmf(1, lambda_given),
                              poisson.pmf(2, lambda_given), poisson.sf(2, lambda_given)])
expected_given = N * probability_given
Q_given = float(((observed - expected_given) ** 2 / expected_given).sum())
df_given = len(observed) - 1
P_given = float(chi2.sf(Q_given, df_given))
given_table = pd.DataFrame({"category": category, "expected_given": expected_given})
given_table
exercise_result = pd.DataFrame({
    "model": ["estimated_lambda", "given_lambda"],
    "lambda": [lambda_mle, lambda_given], "Q": [Q_mle, Q_given],
    "df": [df_mle, df_given], "P": [P_mle, P_given]
})
exercise_result
# 可选：保存结果；改为False则只保留变量，Spyder变量浏览器可查看表。
save_outputs = True
output_dir = project_dir / "运行结果" / "Python"
if save_outputs:
    output_dir.mkdir(parents=True, exist_ok=True)
    result_mle.to_csv(output_dir / "main_result.csv", index=False, encoding="utf-8")
    contribution_table.to_csv(output_dir / "expected_contributions.csv", index=False, encoding="utf-8")
    exercise_result.to_csv(output_dir / "lambda_exercise.csv", index=False, encoding="utf-8")
    given_table.to_csv(output_dir / "given_expected.csv", index=False, encoding="utf-8")
