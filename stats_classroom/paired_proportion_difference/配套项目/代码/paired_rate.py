"""统计课堂27｜主例：从个人前后差值到样本率差D的分布。

问题：同一批人前后检测，后测阳性率减前测阳性率如何变化？本题固定给定q模型，
不做区间估计。主例20人的共享数据是教学构造的逐人记录；q+=(.05)、q0=.70、q−=.25是演算模型，
虽然与主例三类比例相同，也不能据此认定为真实总体参数。
"""
# %% 问题/条件：导入、项目路径与可选保存开关
from pathlib import Path
import math
import pandas as pd
import numpy as np
from scipy.stats import norm

SAVE_RESULTS = False  # 如需保存CSV，将此开关改为True。


def find_project():
    places = [Path.cwd()]
    if "__file__" in globals():
        places.append(Path(__file__).resolve().parent)
    places += [p.parent for p in list(places)]
    for place in places:
        if (place / "数据" / "paired_20.csv").is_file():
            return place
    raise FileNotFoundError("请从配套项目根目录或代码目录启动")


project = find_project()
data_dir = project / "数据"
output_dir = project / "运行结果" / "Python"

# %% 问题/条件：读取共享数据字段和前5行
# before、after是0/1检测结果，z=after−before；先看字段和前5行，再做配对差值。
people = pd.read_csv(data_dir / "paired_20.csv")
data_columns = list(people.columns)
first_five_people = people.head()
data_columns
first_five_people
assert np.isin(people[["before", "after"]].to_numpy(), [0, 1]).all()
individual_difference = people["after"] - people["before"]
assert np.array_equal(individual_difference.to_numpy(), people["z"].to_numpy())

# %% 手动计算：主例人数、前后率与观察到的D
# +1是阴转阳，−1是阳转阴；所有人都进入D的分母，D=个人Z的平均数。
n_main = len(people)
n_plus_main = int((individual_difference == 1).sum())
n_zero_main = int((individual_difference == 0).sum())
n_minus_main = int((individual_difference == -1).sum())
n_both_positive = int(((people["before"] == 1) & (people["after"] == 1)).sum())
observed_before_rate = (n_both_positive + n_minus_main) / n_main
observed_after_rate = (n_both_positive + n_plus_main) / n_main
observed_D = (n_plus_main - n_minus_main) / n_main
observed_summary = pd.DataFrame([{
    "scenario": "main_20", "n": n_main, "nplus": n_plus_main,
    "nzero": n_zero_main, "nminus": n_minus_main, "D": observed_D,
    "before_rate": observed_before_rate, "after_rate": observed_after_rate,
    "n_changed": n_plus_main + n_minus_main,
}])

observed_summary
# D=-.20：后测率.45减前测率.65，下降20个百分点；6个人发生变化。

# %% 条件：给定模型与四个教学场景
scenarios = [
    ("normal_n2", 2, 0.05, 0.25),
    ("normal_n20", 20, 0.05, 0.25),
    ("normal_n1000", 1000, 0.05, 0.25),
    ("rare_n1000", 1000, 0.001, 0.002),
]

# %% 手动计算：逐人卷积得到D的精确概率
# kernel按−1、0、+1排列；每次convolve加入一个独立人的Z，返回总和的概率数组。
def exact_distribution(n, q_plus, q_minus):
    q_zero = 1 - q_plus - q_minus
    probability = np.array([1.0])
    kernel = np.array([q_minus, q_zero, q_plus])
    for _ in range(n):
        probability = np.convolve(probability, kernel)
    k = np.arange(-n, n + 1)
    d = k / n
    return k, d, probability, q_zero


moment_rows = []
distribution_tables = {}
for name, n, q_plus, q_minus in scenarios:
    k, d, exact_probability, q_zero = exact_distribution(n, q_plus, q_minus)
    theta = q_plus - q_minus
    EZ2 = q_plus + q_minus
    individual_variance = EZ2 - theta ** 2
    variance_D = individual_variance / n
    standard_error_D = math.sqrt(variance_D)
    # 正态近似按相邻D格点的半格边界计算，不把离散质量当连续密度。
    edges = (d[:-1] + d[1:]) / 2
    normal_cell_probability = np.empty(2 * n + 1)
    normal_cell_probability[0] = norm.cdf(edges[0], theta, standard_error_D)
    normal_cell_probability[-1] = 1 - norm.cdf(edges[-1], theta, standard_error_D)
    normal_cell_probability[1:-1] = norm.cdf(edges[1:], theta, standard_error_D) - norm.cdf(edges[:-1], theta, standard_error_D)
    table = pd.DataFrame({"k": k, "D": d, "exact_probability": exact_probability,
                          "normal_cell_probability": normal_cell_probability})
    distribution_tables[name] = table
    cdf_error = np.max(np.abs(table["exact_probability"].cumsum() - table["normal_cell_probability"].cumsum()))
    moment_rows.append({"scenario": name, "n": n, "qp": q_plus, "q0": q_zero,
                        "qm": q_minus, "theta": theta, "mean_D": theta,
                        "var_D": variance_D, "max_abs_cdf_error": cdf_error})
moments = pd.DataFrame(moment_rows)

# %% 手动计算：从场景表计数得到前后率与率差
# scenarios.csv提供四个教学场景的计数；这里逐行计算前后率，避免隐藏在函数内。
scenario_input = pd.read_csv(data_dir / "scenarios.csv")
description_rows = [observed_summary.iloc[0].to_dict()]
for row in scenario_input.itertuples(index=False):
    if row.scenario == "main_20":
        continue
    n = int(row.n)
    n_plus = int(row.nplus)
    n_zero = int(row.nzero)
    n_minus = int(row.nminus)
    before_rate = (int(row.both_positive) + n_minus) / n
    after_rate = (int(row.both_positive) + n_plus) / n
    description_rows.append({"scenario": row.scenario, "n": n, "nplus": n_plus,
                             "nzero": n_zero, "nminus": n_minus,
                             "D": (n_plus - n_minus) / n,
                             "before_rate": before_rate, "after_rate": after_rate,
                             "n_changed": n_plus + n_minus})
paired_examples = pd.DataFrame(description_rows)[["scenario", "n", "nplus", "nzero", "nminus", "D", "before_rate", "after_rate", "n_changed"]]

# %% 结果核对：θ、E(Z²)、个人方差、Var(D)与SE(D)
# E(Z)=q+−q−=θ；E(Z²)=q++q−；独立人的平均数方差为个人方差除以n。
theta_main = 0.05 - 0.25
EZ2_main = 0.05 + 0.25
individual_variance_main = EZ2_main - theta_main ** 2
variance_D_main = individual_variance_main / 20
SE_D_main = math.sqrt(variance_D_main)
assert math.isclose(theta_main, -0.20)
assert math.isclose(individual_variance_main, 0.26)
assert math.isclose(variance_D_main, 0.013)
assert math.isclose(SE_D_main, math.sqrt(0.013))

model_moments = pd.DataFrame([{
    "E_Z": theta_main, "E_Z2": EZ2_main, "Var_Z": individual_variance_main,
    "Var_D": variance_D_main, "SE_D": SE_D_main,
}])
model_moments
# .30−(−.20)^2=.26是个人方差；再除以20得.013，开方约.114。
paired_examples
moments

# %% 可选保存：默认False；将开关改True才保存4个distribution、moments和paired_examples
if SAVE_RESULTS:
    output_dir.mkdir(parents=True, exist_ok=True)
    for name, table in distribution_tables.items():
        table.to_csv(output_dir / f"distribution_{name}.csv", index=False)
    moments.to_csv(output_dir / "moments.csv", index=False)
    paired_examples.to_csv(output_dir / "paired_examples.csv", index=False)

# %% 查看图形：自动运行独立绘图文件；保存图形的开关在该文件末尾。
import runpy
runpy.run_path(str(project / "代码" / "plot_figures.py"), run_name="__main__")

# %% 思考：改变q_plus、q_minus或n，先预测θ、Var(D)和SE(D)，再运行比较精确分布与近似。
