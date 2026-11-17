# %% 初始化：数据路径、依赖与保存开关
"""统计课堂30：Cohen h两个率角变换。UTF-8，按1至5节顺序运行。"""
import os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

# 默认只打印；设置环境变量 SAVE_OUTPUTS=true 才保存到 results/Python/。
SAVE_OUTPUTS = os.getenv("SAVE_OUTPUTS", "false").lower() == "true"
data_candidates = [Path("data/lesson10_independent.csv"), Path("../data/lesson10_independent.csv")]
data_path = next((path for path in data_candidates if path.exists()), None)
if data_path is None:
    raise FileNotFoundError("请从项目根或 Python/ 目录运行")
project_root = data_path.parent.parent
result_dir = project_root / "results" / "Python"
if SAVE_OUTPUTS:
    result_dir.mkdir(parents=True, exist_ok=True)

def show_and_save(frame, name):
    """打印一张表；保存开关打开时再写CSV。"""
    print(frame.to_string(index=False, float_format=lambda value: f"{value:.12g}"))
    if SAVE_OUTPUTS:
        frame.to_csv(result_dir / name, index=False)

# %% 1. 读取课件原例；方向固定为新方案减常规对照。
teaching_data = pd.read_csv(data_path)
new_row = teaching_data[teaching_data["group"] == "new"].iloc[0]
control_row = teaching_data[teaching_data["group"] == "control"].iloc[0]
n1, y1 = int(new_row["n"]), int(new_row["success"])
n0, y0 = int(control_row["n"]), int(control_row["success"])
p1, p0 = y1 / n1, y0 / n0
RD = p1 - p0
table = np.array([[y1, n1 - y1], [y0, n0 - y0]])
print("1 读表：p1=", p1, "p0=", p0, "RD=", RD, "\n", table)

# %% 2. asin/sqrt为弧度角变换；ppf(.975)是双侧95%正态临界值。
h = 2 * np.arcsin(np.sqrt(p1)) - 2 * np.arcsin(np.sqrt(p0))
hSE = np.sqrt(1 / n1 + 1 / n0)
z95 = stats.norm.ppf(.975)
hCI = h + np.array([-1, 1]) * z95 * hSE
baseline = pd.DataFrame({"RD": [RD], "h": [h], "hSE": [hSE],
                         "hCI_low": [hCI[0]], "hCI_high": [hCI[1]]})
show_and_save(baseline, "baseline_h.csv")

# %% 3. 固定RD=.2；保留符号方向，同时输出绝对h供|h|=.2/.5/.8分级练习。
fixed = pd.DataFrame({"p0": [.1, .3, .4, .7], "p1": [.3, .5, .6, .9]})
fixed["RD"] = fixed["p1"] - fixed["p0"]
fixed["h"] = 2 * np.arcsin(np.sqrt(fixed["p1"])) - 2 * np.arcsin(np.sqrt(fixed["p0"]))
fixed["abs_h"] = fixed["h"].abs()
show_and_save(fixed[["p0", "p1", "RD", "h"]], "fixed_RD_h.csv")

# %% 4. 每个n逐行计算期望数、Pearson Q、h CI和无连续性校正双侧P。
hh = 2 * np.arcsin(np.sqrt(.5)) - 2 * np.arcsin(np.sqrt(.3))
sample_rows = []
for n in [20, 100, 1000]:
    sample_table = np.array([[.5 * n, .5 * n], [.3 * n, .7 * n]])
    expected = np.outer(sample_table.sum(axis=1), sample_table.sum(axis=0)) / sample_table.sum()
    Q_manual = ((sample_table - expected) ** 2 / expected).sum()
    native = stats.chi2_contingency(sample_table, correction=False)
    h_ci = hh + np.array([-1, 1]) * z95 * np.sqrt(2 / n)
    sample_rows.append([n, hh, h_ci[0], h_ci[1], Q_manual, native.statistic, native.pvalue])
sample_results = pd.DataFrame(sample_rows, columns=["n", "h", "hCI_low", "hCI_high",
                                                    "Q_manual", "Q_native", "P_pearson"])
show_and_save(sample_results, "sample_size_h.csv")

# %% 5. g斜率逐元素作用于原始SE，比较epsilon=.01的exact/linear增量。
variance_table = pd.DataFrame({"p": [.1, .3, .5], "n": [100, 100, 100]})
epsilon = .01
variance_table["raw_SE"] = np.sqrt(variance_table["p"] * (1 - variance_table["p"]) / variance_table["n"])
variance_table["g_slope"] = 1 / np.sqrt(variance_table["p"] * (1 - variance_table["p"]))
variance_table["transformed_SE"] = variance_table["g_slope"] * variance_table["raw_SE"]
variance_table["approx_variance"] = variance_table["transformed_SE"] ** 2
variance_table["exact_increment"] = 2 * np.arcsin(np.sqrt(variance_table["p"] + epsilon)) - 2 * np.arcsin(np.sqrt(variance_table["p"]))
variance_table["linear_increment"] = variance_table["g_slope"] * epsilon
show_and_save(variance_table, "h_variance_stability.csv")
