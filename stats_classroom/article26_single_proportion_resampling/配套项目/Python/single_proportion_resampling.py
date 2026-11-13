# %% 1. 题目、依赖与共享输入
# 假设独立抽取180人，117人达标。p是目标人群的总体达标率。
# 原题H0:p=.45，H1:p>.45，alpha=.05；另估计p的双侧95%区间。
# 练习向量按汇总计数展开，不是带原始顺序的逐人记录。
# Spyder先运行本节，后续用Ctrl+Enter逐节运行；F9运行选中行。
from pathlib import Path
import platform
import numpy as np
import pandas as pd
import scipy
from scipy import stats

if "__file__" in globals():
    project_root = Path(__file__).resolve().parent.parent
else:
    project_root = Path.cwd()
data_path = project_root / "data" / "cases.csv"
if not data_path.exists():
    raise FileNotFoundError("请将Spyder工作目录设为配套项目根，再运行第一节。")
cases = pd.read_csv(data_path)
print(cases)
main_case = cases.loc[cases["case"] == "main"].iloc[0]
n = int(main_case["n"])
x = int(main_case["x"])
p_hat = x / n
p0 = float(main_case["p0"])
assert abs(p_hat - float(main_case["p_hat"])) < 1e-12
B = 10000
seed = 20261110
result_dir = project_root / "results"
result_dir.mkdir(exist_ok=True)

# %% 2. 非参数Bootstrap：有放回抽n个观察值，再计算率
# 数组中的1代表达标，0代表未达标；mean就是达标人数/n。
y = np.concatenate([np.ones(x, dtype=int), np.zeros(n - x, dtype=int)])
print(pd.Series(y).value_counts())
rng = np.random.default_rng(seed)
bootstrap_rates = np.zeros(B)
for b in range(B):
    sampled_y = rng.choice(y, size=n, replace=True)
    bootstrap_rates[b] = np.mean(sampled_y)
# 两端各留下2.5%，得到双侧95%百分位区间。
bootstrap_ci = np.quantile(bootstrap_rates, [.025, .975])
print("Bootstrap双侧95%百分位区间:", bootstrap_ci)
# 给定样本，每次抽到1的概率为p_hat；计数为Bin(n,p_hat)，率再除以n。
conditional_counts = stats.binom.ppf([.025, .975], n=n, p=p_hat)
print("条件二项计数分位:", conditional_counts)
print("条件二项率分位:", conditional_counts / n)

# %% 3. 原假设模拟：按p0生成二项计数，数右尾命中次数
# 按H1方向定义更极端结果（仍在H0下计算）：本题是达标人数至少x人。
null_rng = np.random.default_rng(seed)
null_counts = null_rng.binomial(n=n, p=p0, size=B)
null_hits = int(np.sum(null_counts >= x))
null_hit_rate = null_hits / B
print(f"H0模拟命中{null_hits}/{B}，比例{null_hit_rate:.8f}")
print("0次命中不等于真实P=0；有限模拟在本题难以表达很小的尾概率。")

# %% 4. 精确检验与区间：单侧检验、双侧估计分别调用
exact_test = stats.binomtest(x, n, p=p0, alternative="greater")
print(exact_test)
one_sided_ci = exact_test.proportion_ci(confidence_level=.95, method="exact")
print("原题单侧95%率区间:", one_sided_ci)
# sf(k)=P(X>k)，所以输入x-1得到P(X>=x)，避免1-CDF的消减误差。
manual_p = stats.binom.sf(x - 1, n=n, p=p0)
assert abs(manual_p - exact_test.pvalue) < 1e-15
# 另求双侧95%CP区间，与前面的双侧Bootstrap区间比较。
exact_two_sided = stats.binomtest(x, n, alternative="two-sided")
cp_ci = exact_two_sided.proportion_ci(confidence_level=.95, method="exact")
print("双侧95%CP区间:", cp_ci)

# %% 5. 只重新排序一个0/1向量，总成功数不会改变
permutation_rng = np.random.default_rng(seed)
reordered_rates = np.zeros(100)
for b in range(100):
    reordered_y = permutation_rng.permutation(y)
    reordered_rates[b] = np.mean(reordered_y)
print("重排得到的不同率:", np.unique(reordered_rates))
assert np.all(reordered_rates == p_hat)

# %% 6. 边界：20人全为0、20人全为1，实际重抽核对
boundary_rng = np.random.default_rng(seed)
y_zero = np.zeros(20)
y_all = np.ones(20)
bootstrap_zero = np.zeros(B)
bootstrap_all = np.zeros(B)
for b in range(B):
    bootstrap_zero[b] = np.mean(boundary_rng.choice(y_zero, size=20, replace=True))
    bootstrap_all[b] = np.mean(boundary_rng.choice(y_all, size=20, replace=True))
zero_ci = np.quantile(bootstrap_zero, [.025, .975])
all_ci = np.quantile(bootstrap_all, [.025, .975])
# CP边界由Beta分位数计算；0/20下限0，20/20上限1。
cp_zero = [0.0, stats.beta.ppf(.975, 1, 20)]
cp_all = [stats.beta.ppf(.025, 20, 1), 1.0]
print("0/20: Bootstrap", zero_ci, "CP", cp_zero)
print("20/20: Bootstrap", all_ci, "CP", cp_all)

# %% 7. 同名字段保存；随机端点不要求跨语言逐项相同
result = {
    "language": "Python", "case": "main", "n": n, "x": x, "p_hat": p_hat,
    "p0": p0, "B": B, "seed": seed,
    "bootstrap_ci_low": bootstrap_ci[0], "bootstrap_ci_high": bootstrap_ci[1],
    "conditional_count_q025": conditional_counts[0],
    "conditional_count_q975": conditional_counts[1],
    "h0_hits": null_hits, "h0_hit_rate": null_hit_rate,
    "exact_right_tail_p": exact_test.pvalue,
    "cp_ci_low": cp_ci.low, "cp_ci_high": cp_ci.high,
    "boundary_zero_bootstrap_low": zero_ci[0], "boundary_zero_bootstrap_high": zero_ci[1],
    "boundary_zero_cp_low": cp_zero[0], "boundary_zero_cp_high": cp_zero[1],
    "boundary_all_bootstrap_low": all_ci[0], "boundary_all_bootstrap_high": all_ci[1],
    "boundary_all_cp_low": cp_all[0], "boundary_all_cp_high": cp_all[1],
}
result_table = pd.DataFrame([result])
result_table.to_csv(result_dir / "Python_results.csv", index=False)
report = [
    f"Python {platform.python_version()}, NumPy {np.__version__}, pandas {pd.__version__}, SciPy {scipy.__version__}",
    result_table.to_string(index=False),
    f"单侧精确P={exact_test.pvalue:.12g}；单侧95%率区间={one_sided_ci}",
    "重排100次的率均为0.65；模拟命中比例不是精确P。",
    "CP为模型下的精确区间；Bootstrap百分位区间不具有同样的精确覆盖保证。",
]
(result_dir / "Python_console.txt").write_text("\n".join(report), encoding="utf-8")
print("\n".join(report))

# %% 8. Bootstrap区间构造与覆盖率精确枚举（v3）
# Percentile/Basic使用条件二项分布的无限次Bootstrap端点，
# 不把10000次模拟端点冒充为理想无限内层结果。
z975 = stats.norm.ppf(.975)
percentile_counts = stats.binom.ppf([.025, .975], n=n, p=p_hat)
percentile_ci_v3 = percentile_counts / n
basic_ci_v3 = np.array([2 * p_hat - percentile_ci_v3[1], 2 * p_hat - percentile_ci_v3[0]])
wald_half = z975 * np.sqrt(p_hat * (1 - p_hat) / n)
wald_ci_v3 = np.array([p_hat - wald_half, p_hat + wald_half])
wilson_den = 1 + z975**2 / n
wilson_center = (p_hat + z975**2 / (2 * n)) / wilson_den
wilson_half = z975 / wilson_den * np.sqrt(p_hat * (1 - p_hat) / n + z975**2 / (4 * n**2))
wilson_ci_v3 = np.array([wilson_center - wilson_half, wilson_center + wilson_half])
cp_ci_v3_obj = stats.binomtest(x, n, alternative="two-sided").proportion_ci(confidence_level=.95, method="exact")
cp_ci_v3 = np.array([cp_ci_v3_obj.low, cp_ci_v3_obj.high])
observed_intervals_v3 = pd.DataFrame({
    "method": ["Percentile", "Basic", "Wald", "Wilson", "CP"],
    "lower": [percentile_ci_v3[0], basic_ci_v3[0], wald_ci_v3[0], wilson_ci_v3[0], cp_ci_v3[0]],
    "upper": [percentile_ci_v3[1], basic_ci_v3[1], wald_ci_v3[1], wilson_ci_v3[1], cp_ci_v3[1]],
})
observed_intervals_v3.to_csv(result_dir / "Python_intervals_v3.csv", index=False)

def interval_for_x_v3(x_value):
    ph = x_value / n
    perc = stats.binom.ppf([.025, .975], n=n, p=ph) / n
    basic = np.array([2 * ph - perc[1], 2 * ph - perc[0]])
    wh = z975 * np.sqrt(ph * (1 - ph) / n)
    wald = np.array([ph - wh, ph + wh])
    wd = 1 + z975**2 / n
    wc = (ph + z975**2 / (2 * n)) / wd
    wterm = z975 / wd * np.sqrt(ph * (1 - ph) / n + z975**2 / (4 * n**2))
    wilson = np.array([wc - wterm, wc + wterm])
    if x_value == 0:
        cp = np.array([0.0, 1 - .025 ** (1 / n)])
    elif x_value == n:
        cp = np.array([.025 ** (1 / n), 1.0])
    else:
        cp = np.array([stats.beta.ppf(.025, x_value, n - x_value + 1),
                       stats.beta.ppf(.975, x_value + 1, n - x_value)])
    return {"Percentile": perc, "Basic": basic, "Wald": wald, "Wilson": wilson, "CP": cp}

p_grid = np.arange(.01, 1.00, .01)
coverage_rows_v3 = []
for true_p in p_grid:
    for method_name in ["Percentile", "Basic", "Wald", "Wilson", "CP"]:
        covered = []
        for x_value in range(n + 1):
            lower, upper = interval_for_x_v3(x_value)[method_name]
            covered.append(true_p >= lower - 1e-12 and true_p <= upper + 1e-12)
        pmf = stats.binom.pmf(np.arange(n + 1), n=n, p=true_p)
        coverage_rows_v3.append({"p": true_p, "method": method_name, "coverage": float(np.sum(pmf * covered))})
coverage_v3 = pd.DataFrame(coverage_rows_v3)
coverage_v3.to_csv(result_dir / "Python_coverage_v3.csv", index=False)

depth_report = [
    f"Python {platform.python_version()}, NumPy {np.__version__}, SciPy {scipy.__version__}",
    "v3: Percentile/Basic use infinite conditional Binomial endpoints; coverage is exact enumeration.",
    "Observed intervals (Percentile/Basic endpoints directly from conditional distribution):",
    observed_intervals_v3.to_string(index=False),
    f"基础主脚本循环重抽的SD（bootstrap_rates）：{np.std(bootstrap_rates, ddof=1):.15g}",
    f"Analytic SE: {np.sqrt(.65 * .35 / 180):.15g}",
    "Coverage at p=.65:", coverage_v3[coverage_v3["p"] == .65].to_string(index=False),
    "Coverage at p=.05:", coverage_v3[coverage_v3["p"] == .05].to_string(index=False),
]
(result_dir / "Python_depth_v3_console.txt").write_text("\n".join(depth_report), encoding="utf-8")
print("\n".join(depth_report))
