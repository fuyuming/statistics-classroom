# 统计课堂36：固定lambda=.5的Storey q估计。
# 输入同一份1000个模拟P值；本篇指定算法得到qhat=pi0hat*BH。
# 在Spyder根目录从02.1开始，无需先运行01。

# %% 02.1 准备输入与本项目自编函数
import os
import sys
from pathlib import Path
import numpy as np
import pandas as pd

if not Path("data/screening_pvalues.csv").is_file():
    raise FileNotFoundError("请将工作目录设为项目根目录，再从02.1开始。")
# helpers是随项目提供的目录；加入模块搜索路径后才能导入自编函数。
sys.path.insert(0, str(Path("helpers").resolve()))
from storey_q_fixed import storey_q_fixed
raw = pd.read_csv("data/screening_pvalues.csv")
print(raw.head())
p = raw["p"].to_numpy(float)
m = len(p)
lambda_value = 0.5
alpha = 0.05

# %% 02.2 手算真实零假设比例
right_tail_count = np.sum(p >= lambda_value)
pi0 = min(1.0, right_tail_count / (m * (1 - lambda_value)))
pi0_result = pd.DataFrame({"lambda": [lambda_value],
                           "tail_count": [right_tail_count], "pi0": [pi0]})
print(pi0_result)

# %% 02.3 手算本篇参数下的q值
order = np.argsort(p)
sorted_p = p[order]
rank_p = np.arange(1, m + 1)
bh_candidate = m * sorted_p / rank_p
bh_sorted = np.minimum(1.0, np.minimum.accumulate(bh_candidate[::-1])[::-1])
bh = np.empty(m)
bh[order] = bh_sorted
q_manual = pi0 * bh

# %% 02.4 调用配套自编函数核对
# 这不是R后端，也不是Python原生qvalue包；helper实现正文的固定lambda算法。
pi0_function, q_function = storey_q_fixed(p, lambda_value)
function_check = pd.DataFrame({
    "pi0_manual": [pi0], "pi0_function": [pi0_function],
    "max_abs_q_error": [np.max(abs(q_manual - q_function))],
})
print(function_check)
assert function_check["max_abs_q_error"].iloc[0] < 1e-10

# %% 02.5 查看结果与发现数
comparison = raw.assign(pi0=pi0_function, q=q_function, BH=bh)
summary_result = pd.DataFrame({
    "method": ["q_Storey_lambda_0.5", "BH"],
    "discoveries": [np.sum(q_function <= alpha), np.sum(bh <= alpha)],
    "pi0": [pi0_function, np.nan],
})
print(summary_result)

# %% 02.6 lambda敏感性选读
lambda_grid = np.array([.3, .5, .7])
pi0_grid = np.empty(len(lambda_grid))
for i, current_lambda in enumerate(lambda_grid):
    pi0_grid[i] = min(1.0, np.mean(p >= current_lambda) / (1 - current_lambda))
lambda_sensitivity = pd.DataFrame({"lambda": lambda_grid, "pi0": pi0_grid})
print(lambda_sensitivity)

# %% 02.7 可选保存
if os.environ.get("ARTICLE36_SAVE_RESULTS") == "1":
    Path("results").mkdir(exist_ok=True)
    comparison.to_csv("results/Python_02_storey_q.csv", index=False)
    summary_result.to_csv("results/Python_02_summary.csv", index=False)
    lambda_sensitivity.to_csv("results/Python_02_lambda_sensitivity.csv", index=False)
