"""统计课堂22：Tukey检验与学生化极差q。
数据：第12课例9-1，三个小麦品种各10个独立地块；原题未给单位。
Spyder打开项目文件夹后F5运行。依赖见requirements.txt。
输出在运行结果/Python；本文件末尾调用plot_tukey_q.py生成文章配图。
"""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import scipy
import matplotlib.pyplot as plt
from scipy import stats
from plot_tukey_q import make_figures

# %% 1. 看数据并先做条件诊断
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "运行结果/Python"
OUT.mkdir(parents=True, exist_ok=True)
wide = pd.read_csv(ROOT / "数据/小麦产量.csv", encoding="utf-8-sig")
print(wide)
assert wide.shape == (10, 3) and not wide.isna().any().any()
values = wide.to_numpy(dtype=float)
# 每一列是一个品种的10个独立地块；列与列不配对。
samples = [values[:, i] for i in range(3)]

def save_table(table, name):
    table.to_csv(OUT / name, index=False, encoding="utf-8-sig")

# Shapiro看各组正态性，Levene/Bartlett看组间方差；P>.05只是诊断线索，
# 不能证明前提，也不应机械地据此选择检验。Q-Q图帮助观察形状与离群点。
shapiro_p = [stats.shapiro(x).pvalue for x in samples]
# Levene比较绝对偏差，中位数中心较稳健；Bartlett对非正态较敏感。
levene_median = stats.levene(*samples, center="median")
levene_mean = stats.levene(*samples, center="mean")
bartlett = stats.bartlett(*samples)
variance_tests = pd.DataFrame([
    {"method": "Levene median", "statistic": levene_median.statistic, "P": levene_median.pvalue},
    {"method": "Levene mean", "statistic": levene_mean.statistic, "P": levene_mean.pvalue},
    {"method": "Bartlett", "statistic": bartlett.statistic, "P": bartlett.pvalue},
])
save_table(pd.DataFrame({"group": [1, 2, 3], "shapiro_P": shapiro_p}), "diagnostics.csv")
save_table(variance_tests, "variance_tests.csv")
print(variance_tests)
print(pd.DataFrame({"group": [1, 2, 3], "shapiro_P": shapiro_p}))
fig, axes = plt.subplots(1, 3, figsize=(12, 4), layout="constrained")
for i, ax in enumerate(axes):
    stats.probplot(samples[i], dist="norm", plot=ax)
    ax.set_title(f"Variety {i+1}")
fig.savefig(OUT / "qq_diagnostics_before_anova.png", dpi=150)
if plt.get_backend().lower() != "agg":
    plt.show(block=False)

# %% 2. ANOVA及Tukey主分析（诊断之后）
fit = stats.f_oneway(*samples)
tukey = stats.tukey_hsd(*samples)
intervals = tukey.confidence_interval(confidence_level=0.95)
print(fit)
print(tukey)

# %% 3. 分步计算MSE、自由度、q的尺度
n, k = values.shape
N = n * k
nu = N - k
means = values.mean(axis=0)
# axis=0沿各列计算；广播使每列减去自己的均值。
residuals = values - means
within_ss = (residuals**2).sum(axis=0)
sse = within_ss.sum()
mse = sse / nu
ssa = n * ((means - values.mean())**2).sum()
f_value = (ssa / (k - 1)) / mse
p_anova = stats.f.sf(f_value, k - 1, nu)
scale_q = np.sqrt(mse / n)
group_result = pd.DataFrame({"group": [1, 2, 3], "n": n, "mean": means,
                             "variance": within_ss / (n - 1), "within_ss": within_ss})
anova_result = pd.DataFrame([{"SSA": ssa, "SSE": sse, "df1": k - 1, "df2": nu,
                              "MSE": mse, "F": f_value, "P": p_anova}])
save_table(group_result, "groups.csv")
save_table(anova_result, "anova.csv")
assert np.isclose(f_value, fit.statistic) and np.isclose(p_anova, fit.pvalue)

# %% 4. 三对均值差、调整P和同时区间
# Python从0编号；(1,0)对应品种2减品种1。
pair_index = [(1, 0), (2, 0), (2, 1)]
q_critical = stats.studentized_range.ppf(0.95, k, nu)
half_width = q_critical * scale_q
rows = []
for i, j in pair_index:
    difference = means[i] - means[j]
    q_value = abs(difference) / scale_q
    # sf是右尾概率；k输入整个比较族的组数。
    p_adj = stats.studentized_range.sf(q_value, k, nu)
    lower = difference - half_width
    upper = difference + half_width
    rows.append({"pair": f"{i+1}-{j+1}", "difference": difference, "q": q_value,
                 "p_adj": p_adj, "lower": lower, "upper": upper})
    assert np.isclose(difference, tukey.statistic[i, j])
    assert np.isclose(p_adj, tukey.pvalue[i, j], atol=1e-10)
    assert np.isclose(lower, intervals.low[i, j], atol=1e-7)
    assert np.isclose(upper, intervals.high[i, j], atol=1e-7)
pairs = pd.DataFrame(rows)
save_table(pairs, "pairs.csv")
print(pairs.to_string(index=False))
print("品种2-1的均值差为-90；同时95%区间约[-137.59,-42.41]。")
print("品种3-1区间跨0。结合差异量级和区间解释调整P值。")

# %% 5. k=2时与双侧t核对
q_grid = np.array([0, 1, 3, 6], dtype=float)
p_from_q = stats.studentized_range.sf(q_grid, 2, nu)
p_from_t = 2 * stats.t.sf(q_grid / np.sqrt(2), nu)
save_table(pd.DataFrame({"q": q_grid, "q_tail": p_from_q, "t_two_sided": p_from_t}), "k2_check.csv")
assert np.allclose(p_from_q, p_from_t, atol=1e-7)


# 最大绝对均值差就是极差；至少一对越界等价于最大差越界。
# 这里计算总体均值全相等时的理论FWER，重复抽样时MSE也随样本变化。
q_limits = np.array([np.sqrt(2) * stats.t.ppf(0.975, nu), q_critical])
family_control = pd.DataFrame({
    "method": ["Unadjusted pairs", "Tukey family"],
    "qcrit": q_limits, "diff_threshold": q_limits * scale_q,
    "fwer_global_null": stats.studentized_range.sf(q_limits, k, nu),
})
save_table(family_control, "family_control.csv")
print(family_control)
print("差72：单对P与Tukey调整P",
      2 * stats.t.sf(72 / (np.sqrt(2) * scale_q), nu),
      stats.studentized_range.sf(72 / scale_q, k, nu))

# %% 6. 模拟q分布，比较组数2、3、5
# Q=(max Z-min Z)/sqrt(U/nu)。Z是独立标准正态，U为独立卡方(nu)。
# 正态模型下，它分别表示标准化组均值误差和共同方差估计的波动。
rng = np.random.default_rng(20261101)
B = 100000
z = rng.normal(size=(B, 5))
u_chisq = rng.chisquare(nu, size=B)
denominator = np.sqrt(u_chisq / nu)
k_values = [2, 3, 5]
q_samples = {}
simulation_rows = []
for group_count in k_values:
    selected = z[:, :group_count]
    # axis=1沿每次模拟的多个组求最大、最小。
    sample_range = selected.max(axis=1) - selected.min(axis=1)
    simulated_q = sample_range / denominator
    cutoff = stats.studentized_range.ppf(0.95, group_count, nu)
    q_samples[f"k{group_count}"] = simulated_q
    simulation_rows.append({"k": group_count, "df": nu, "q95": cutoff,
                            "simulated_q95": np.quantile(simulated_q, 0.95),
                            "exceedance": np.mean(simulated_q > cutoff)})
simulation = pd.DataFrame(simulation_rows)
save_table(pd.DataFrame(q_samples), "q_samples.csv")
save_table(simulation, "simulation.csv")
mc_se = np.sqrt(0.05 * 0.95 / B)
assert np.all(np.abs(simulation.exceedance - 0.05) < 5 * mc_se)
print(simulation.to_string(index=False))

# %% 7. 改单位：全部产量乘10
scaled = stats.tukey_hsd(*(sample * 10 for sample in samples))
scaled_ci = scaled.confidence_interval(0.95)
unit_rows = []
for r, (i, j) in enumerate(pair_index):
    unit_rows.append({"pair": pairs.loc[r, "pair"], "difference_x10": scaled.statistic[i, j],
                      "lower_x10": scaled_ci.low[i, j], "upper_x10": scaled_ci.high[i, j],
                      "p_x10": scaled.pvalue[i, j]})
unit_check = pd.DataFrame(unit_rows)
save_table(unit_check, "unit_check.csv")
assert np.allclose(unit_check.difference_x10, pairs.difference * 10)
assert np.allclose(unit_check.lower_x10, pairs.lower * 10)
assert np.allclose(unit_check.p_x10, pairs.p_adj, atol=1e-10)
print("差值、区间随单位乘10；q和P不变。")

# %% 8. 全局H0下的整组误报教学模拟（每轮三组同时判断）
rng_null = np.random.default_rng(20261102)
sim_trials, sim_n, sim_mean, sim_variance = 10000, 10, 720, 1842
sim_df = 3 * (sim_n - 1)
sim_values = rng_null.normal(sim_mean, np.sqrt(sim_variance),
                             size=(sim_trials, 3, sim_n))
sim_means = sim_values.mean(axis=2)
sim_mse = ((sim_values - sim_means[:, :, None]) ** 2).sum(axis=(1, 2)) / sim_df
sim_max_q = (sim_means.max(axis=1) - sim_means.min(axis=1)) / np.sqrt(sim_mse / sim_n)
sim_thresholds = {"Unadjusted pairs": np.sqrt(2) * stats.t.ppf(.975, sim_df),
                  "Tukey family": stats.studentized_range.ppf(.95, 3, sim_df)}
sim_rows = []
for method, threshold in sim_thresholds.items():
    false_count = int(np.count_nonzero(sim_max_q > threshold))
    sim_rows.append({"method": method, "trials": sim_trials,
                     "false_positive_trials": false_count,
                     "rate": false_count / sim_trials})
null_family_simulation = pd.DataFrame(sim_rows)
save_table(null_family_simulation, "null_family_simulation.csv")
print(null_family_simulation.to_string(index=False))

# 教学汇总例子：按不同人数重新计算总均值及SSA。
demo_means = np.array([756, 666, 738], dtype=float)
ssa_rows = []
for label, sizes in [("equal_n_10_10_10", [10, 10, 10]),
                     ("unequal_n_5_10_15", [5, 10, 15])]:
    sizes = np.asarray(sizes)
    grand = np.sum(sizes * demo_means) / sizes.sum()
    ssa_rows.append({"design": label, "grand_mean": grand,
                     "ssa": np.sum(sizes * (demo_means-grand)**2)})
ssa_demo = pd.DataFrame(ssa_rows)
save_table(ssa_demo, "ssa_demo.csv")
assert np.allclose(ssa_demo.grand_mean, [720, 717])
assert np.allclose(ssa_demo.ssa, [45360, 40230])

# %% 9. 文章配图（诊断已在ANOVA前完成）
make_figures(ROOT, wide, pairs, nu)
from plot_tukey_q import make_null_family_figure
make_null_family_figure(ROOT, null_family_simulation)
(OUT / "version.txt").write_text(
    f"Python {sys.version.split()[0]}\nNumPy {np.__version__}\n"
    f"pandas {pd.__version__}\nSciPy {scipy.__version__}\n所有核对通过。\n", encoding="utf-8")
