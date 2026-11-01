"""统计课堂19；Python/Spyder F5完整运行。
依赖numpy、pandas、scipy、matplotlib。读取教材例题汇总频数，
输出运行结果/Python，再调用同目录绘图脚本。不是患者个体原始记录。
"""
from pathlib import Path
import runpy
import numpy as np
import pandas as pd
from scipy import stats

# %% 一、读数据：每行是一种疗法，每列是结局人数
root = Path(__file__).resolve().parents[1]
out = root / "运行结果" / "Python"
out.mkdir(parents=True, exist_ok=True)
d = pd.read_csv(root / "数据" / "therapies.csv")
print(d)
a = d[["effective", "ineffective"]].to_numpy()
pair_indices = [(0, 1), (0, 2), (1, 2)]  # Python行号从0开始。
pair_names = ["physical-drug", "physical-topical", "drug-topical"]

# %% 二、整体与两两：correct不等于多重比较调整
chi, p_all, df, expected = stats.chi2_contingency(a, correction=False)
print("整体期望频数：", expected)
raw = np.empty(3)
expected_min = np.empty(3)
for i, rows in enumerate(pair_indices):
    subset = a[list(rows), :]
    fit = stats.chi2_contingency(subset, correction=True)
    raw[i] = fit.pvalue
    expected_min[i] = fit.expected_freq.min()

# %% 三、Bonferroni与Holm逐步计算
m = len(raw)
bonf = np.minimum(m * raw, 1)
order = np.argsort(raw)
sorted_p = raw[order]
factor = np.arange(m, 0, -1)
middle = sorted_p * factor
# maximum.accumulate就是本行和之前各行乘积的累计最大值。
sorted_holm = np.minimum(np.maximum.accumulate(middle), 1)
holm = np.empty(m)
holm[order] = sorted_holm
result = pd.DataFrame({"pair": pair_names, "raw_p_yates": raw,
                       "bonferroni_p": bonf, "holm_p": holm,
                       "expected_min": expected_min})
print(result)
print("0.05标准下Bonferroni、Holm拒绝数：", (bonf <= .05).sum(), (holm <= .05).sum())
print("样本有效率并未因多重比较方法改变。")

# %% 四、人工演示与练习；先计算，再保存答案
example = np.array([.010, .011, .012])
example_middle = np.sort(example) * factor
example_final = np.minimum(np.maximum.accumulate(example_middle), 1)
exercise = pd.DataFrame({"pair": pair_names, "bonf_at_001": bonf <= .01,
                         "holm_at_001": holm <= .01})
print(exercise)
permutation = np.array([2, 0, 1])
permuted_p = raw[permutation]
permuted_order = np.argsort(permuted_p)
permuted_adjusted = np.empty(m)
permuted_adjusted[permuted_order] = np.minimum(
    np.maximum.accumulate(permuted_p[permuted_order] * factor), 1)
restored = np.empty(m)
restored[permutation] = permuted_adjusted
np.testing.assert_allclose(restored, holm, rtol=1e-13)
np.testing.assert_allclose(example_final, [.03, .03, .03])

# %% 五、保存共享格式的结果，调用科学绘图
pd.DataFrame({"statistic": [chi], "df": [df], "p_value": [p_all],
              "expected_min": [expected.min()]}).to_csv(out / "overall_result.csv", index=False)
result.to_csv(out / "pair_results.csv", index=False)
pd.DataFrame({"rank": np.arange(1, m + 1), "pair_sorted": np.array(pair_names)[order],
              "raw_p_sorted": sorted_p, "factor": factor, "middle": middle,
              "holm_cummax": sorted_holm}).to_csv(out / "holm_intermediate.csv", index=False)
pd.DataFrame({"p_input": example, "factor": factor, "middle": example_middle,
              "final": example_final}).to_csv(out / "holm_demo.csv", index=False)
exercise.to_csv(out / "alpha001_exercise.csv", index=False)
runpy.run_path(str(Path(__file__).with_name("plot_multiple_comparisons.py")), run_name="__main__")
print("完成：结果在运行结果/Python，图在文章配图。输入顺序置换核验通过。")

# %% 六、补充率差区间：复现R prop.test两样本的近似区间
# 单项95%区间，没有多重调整；与Bonferroni P所控制的层级不同。
x_pair = a[:2, 0]
n_pair = a[:2, :].sum(axis=1)
rates = x_pair / n_pair
delta = rates[0] - rates[1]
standard_error = np.sqrt(np.sum(rates * (1 - rates) / n_pair))
# R在两样本区间中加连续性修正，并限制不超过观察率差。
correction = min(.5 * np.sum(1 / n_pair), abs(delta))
half_width = stats.norm.ppf(.975) * standard_error + correction
ci_result = pd.DataFrame({"pair": ["physical-drug"], "estimate_pp": [100 * delta],
                          "lower_pp": [100 * max(delta - half_width, -1)],
                          "upper_pp": [100 * min(delta + half_width, 1)]})
print(ci_result)
ci_result.to_csv(out / "effect_ci.csv", index=False)
