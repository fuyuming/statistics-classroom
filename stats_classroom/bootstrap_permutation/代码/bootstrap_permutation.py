# %% 明哥的微生物世界 · 统计课堂
# 学习方式：按分节逐段运行，先读注释，再查看结果对象；预期数值用于核对。
# Bootstrap与置换检验：反复抽样，究竟在模拟什么？
#
# Spyder：在项目目录中打开本文件，按F5完整运行；也可按# %%分节运行。
# 依赖见requirements.txt。脚本使用__file__定位data/，输出写入运行结果/Python/。
# 主要结果用单独对象表达，便于Spyder在分节运行时查看；不把print当作主计算接口。

# %% 0. 项目路径、输入和输出
from pathlib import Path
import sys
from itertools import combinations

import numpy as np
import pandas as pd
import scipy
from scipy import stats
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

# 显式寻找支持中文的字体，避免运行成功却生成方框字。
from matplotlib import font_manager
font_names = ["Microsoft YaHei", "SimHei", "Noto Sans CJK SC", "PingFang SC", "Heiti SC", "Arial Unicode MS", "WenQuanYi Zen Hei"]
font_path = None
for font_name in font_names:
    try:
        font_path = font_manager.findfont(font_name, fallback_to_default=False)
        break
    except ValueError:
        continue
chinese_font = font_manager.FontProperties(fname=font_path).get_name() if font_path else "DejaVu Sans"
plt.rcParams["font.family"] = chinese_font
plt.rcParams["axes.unicode_minus"] = False


project_dir = Path(__file__).resolve().parents[1]
out_dir = project_dir / "运行结果" / "Python"
figure_dir = out_dir / "图"
project_figure_dir = project_dir / "文章配图"
article_figure_dir = project_figure_dir  # 输出限定在本项目
out_dir.mkdir(parents=True, exist_ok=True)
figure_dir.mkdir(parents=True, exist_ok=True)
project_figure_dir.mkdir(parents=True, exist_ok=True)
article_figure_dir.mkdir(parents=True, exist_ok=True)

rats = pd.read_csv(project_dir / "data" / "rats.csv", encoding="utf-8-sig")
index_a = pd.read_csv(
    project_dir / "data" / "bootstrap_indices_A.csv", header=None
).to_numpy(dtype=int) - 1
index_b = pd.read_csv(
    project_dir / "data" / "bootstrap_indices_B.csv", header=None
).to_numpy(dtype=int) - 1

group_a = rats.loc[rats["group"] == "A", "gain_g"].to_numpy(dtype=float)
group_b = rats.loc[rats["group"] == "B", "gain_g"].to_numpy(dtype=float)
n_a = len(group_a)
n_b = len(group_b)
b_rep = index_a.shape[0]

data_summary = pd.DataFrame({
    "group": ["A", "B"],
    "n": [n_a, n_b],
    "mean_g": [group_a.mean(), group_b.mean()],
    "sd_g": [group_a.std(ddof=1), group_b.std(ddof=1)],
})
data_summary
data_summary.to_csv(out_dir / "data_summary.csv", index=False, encoding="utf-8-sig")

assert index_b.shape == (b_rep, n_b)
assert index_a.shape == (b_rep, n_a)
assert np.all((index_a >= 0) & (index_a < n_a))
assert np.all((index_b >= 0) & (index_b < n_b))

# %% 1. 观察到的均数差和学生化统计量
# SE使用两组各自的样本方差和n；T是均数差除以估计SE。
delta_obs = group_a.mean() - group_b.mean()
se_obs = np.sqrt(
    group_a.var(ddof=1) / n_a + group_b.var(ddof=1) / n_b
)
t_obs = delta_obs / se_obs
observed_result = pd.DataFrame({
    "delta_g": [delta_obs],
    "SE_g": [se_obs],
    "T": [t_obs],
})
observed_result
observed_result.to_csv(out_dir / "observed_result.csv", index=False, encoding="utf-8-sig")

# %% 2. 普通按组Bootstrap：估计均数差的波动
# 每行索引是一轮；组内有放回且每轮仍是12和7只。
# 这个分布围绕观察差值，回答“估计有多不稳定”，不是H0的零分布。
delta_star = (
    group_a[index_a].mean(axis=1)
    - group_b[index_b].mean(axis=1)
)
bootstrap_se = delta_star.std(ddof=1)
ci_percentile = np.quantile(delta_star, [0.025, 0.975], method="linear")
bootstrap_estimate = pd.DataFrame({
    "bootstrap_SE_g": [bootstrap_se],
    "percentile_low_g": [ci_percentile[0]],
    "percentile_high_g": [ci_percentile[1]],
})
bootstrap_estimate
bootstrap_estimate.to_csv(out_dir / "bootstrap_estimate.csv", index=False, encoding="utf-8-sig")
pd.DataFrame({"delta_star_g": delta_star}).to_csv(
    out_dir / "delta_star.csv", index=False, encoding="utf-8-sig"
)

# %% 3. 检验的H0参照：先中心化，再对每轮差值学生化
# 中心化使两组经验均数都为0，表达“假定H0:两个总体均数相等”的参照模型。
# 学生化用每轮自己的SE，把差值改写为“几个估计标准误”。
group_a0 = group_a - group_a.mean()
group_b0 = group_b - group_b.mean()
a0_star = group_a0[index_a]
b0_star = group_b0[index_b]
delta0_star = a0_star.mean(axis=1) - b0_star.mean(axis=1)
se_star = np.sqrt(
    a0_star.var(axis=1, ddof=1) / n_a
    + b0_star.var(axis=1, ddof=1) / n_b
)
t_null = delta0_star / se_star

assert np.allclose(
    t_null, (delta_star - delta_obs) / se_star, rtol=1e-12, atol=1e-12
)
assert np.isclose(
    (group_a.mean() - group_b.mean())
    / np.sqrt(group_a.var(ddof=1) / n_a + group_b.var(ddof=1) / n_b),
    t_obs,
)

# 等尾双侧规则：较小尾部比例乘2；固定索引下应为左9636、右364。
left_count = int(np.sum(t_null <= t_obs))
right_count = int(np.sum(t_null >= t_obs))
p_boot = min(1.0, 2 * min(left_count, right_count) / b_rep)
q_t = np.quantile(t_null, [0.025, 0.975], method="linear")
ci_t = delta_obs - q_t[::-1] * se_obs
bootstrap_test = pd.DataFrame({
    "left": [left_count],
    "right": [right_count],
    "B": [b_rep],
    "P_two_sided_equal_tail": [p_boot],
    "bootstrap_t_low_g": [ci_t[0]],
    "bootstrap_t_high_g": [ci_t[1]],
})
bootstrap_test
bootstrap_test.to_csv(out_dir / "bootstrap_test.csv", index=False, encoding="utf-8-sig")
pd.DataFrame({"T_null": t_null}).to_csv(
    out_dir / "T_null.csv", index=False, encoding="utf-8-sig"
)

# %% 4. 置换检验：混合19个观测，再按12/7重分
# 若H0包含“组标签可交换”，固定19个观测，枚举C(19,12)=50388种分组。
# 这里不做组内有放回抽样；每个观测在一轮中恰好使用一次。
all_values = np.concatenate([group_a, group_b])
all_index = np.arange(len(all_values))
permutation_groups = list(combinations(all_index, n_a))
permutation_rep = len(permutation_groups)
permutation_delta = np.empty(permutation_rep)
for i, perm_a_tuple in enumerate(permutation_groups):
    perm_a_index = np.asarray(perm_a_tuple)
    perm_b_index = np.setdiff1d(all_index, perm_a_index, assume_unique=False)
    permutation_delta[i] = (
        all_values[perm_a_index].mean() - all_values[perm_b_index].mean()
    )
perm_extreme = int(np.sum(np.abs(permutation_delta) >= abs(delta_obs)))
perm_p = perm_extreme / permutation_rep
permutation_result = pd.DataFrame({
    "total_partitions": [permutation_rep],
    "extreme_abs_delta": [perm_extreme],
    "P_two_sided_exact": [perm_p],
})
permutation_result
permutation_result.to_csv(out_dir / "permutation_result.csv", index=False, encoding="utf-8-sig")
pd.DataFrame({"permutation_delta_g": permutation_delta}).to_csv(
    out_dir / "permutation_delta.csv", index=False, encoding="utf-8-sig"
)

# %% 5. 率的边界例：普通百分位Bootstrap与CP区间
# 对x/n，理想无限Bootstrap分布是K*~Binomial(n,x/n)。
# 用“最小k使CDF>=q”定义离散分位点，直接看边界退化，不靠有限B模拟。
def binom_quantile_index(n, p, probability):
    for k in range(n + 1):
        if stats.binom.cdf(k, n, p) >= probability - 1e-14:
            return k
    raise RuntimeError("二项分布分位点没有找到")


def percentile_rate_ci(x, n, alpha=0.05):
    p_hat = x / n
    lower_k = binom_quantile_index(n, p_hat, alpha / 2)
    upper_k = binom_quantile_index(n, p_hat, 1 - alpha / 2)
    return np.array([lower_k / n, upper_k / n], dtype=float)


ci_9_of_10 = percentile_rate_ci(9, 10)
ci_0_of_10 = percentile_rate_ci(0, 10)
cp_0_upper = 1 - 0.025 ** (1 / 10)
rate_result = pd.DataFrame({
    "case": ["9/10", "0/10"],
    "percentile_low": [ci_9_of_10[0], ci_0_of_10[0]],
    "percentile_high": [ci_9_of_10[1], ci_0_of_10[1]],
    "CP_low": [np.nan, 0.0],
    "CP_high": [np.nan, cp_0_upper],
})
rate_result
rate_result.to_csv(out_dir / "rate_boundary_results.csv", index=False, encoding="utf-8-sig")

# 给定教学模型true_p=.9，精确枚举观察到的x=0,...,10；不是模拟。
true_p = 0.9
coverage_rows = []
for x in range(11):
    interval = percentile_rate_ci(x, 10)
    coverage_rows.append({
        "x": x,
        "interval_low": interval[0],
        "interval_high": interval[1],
        "true_p_covered": bool(interval[0] <= true_p <= interval[1]),
        "probability_under_true_p": stats.binom.pmf(x, 10, true_p),
    })
coverage_table = pd.DataFrame(coverage_rows)
coverage_table
coverage_exact = coverage_table.loc[
    coverage_table["true_p_covered"], "probability_under_true_p"
].sum()
coverage_summary = pd.DataFrame({
    "true_p": [true_p],
    "all_success_probability": [true_p ** 10],
    "exact_enumerated_coverage": [coverage_exact],
})
coverage_summary
coverage_table.to_csv(out_dir / "rate_coverage_enumeration.csv", index=False, encoding="utf-8-sig")
coverage_summary.to_csv(out_dir / "rate_coverage_summary.csv", index=False, encoding="utf-8-sig")

# %% 6. 两张解释图
# 英文坐标标签避免额外字体依赖；文章图注负责中文说明。
def draw_figures(delta_star_values, t_null_values, permutation_values, target_dir, chinese=False):
    target_dir.mkdir(parents=True, exist_ok=True)
    font = None
    chinese = chinese and font_path is not None  # 无中文字体时生成英文图，数值计算照常完成
    if chinese:
        # 文章图使用本机已有字体；不联网下载字体。项目README说明英文图的跨平台版本。
        font = FontProperties(fname=font_path, size=14)
        matplotlib.rcParams["axes.unicode_minus"] = False
        bootstrap_title = "按组 Bootstrap：估计均数差"
        null_title = r"中心化、学生化的 $H_0$ 参照"
        permutation_title = "置换检验：19个观测按12/7重分"
        difference_label = "A组均数 − B组均数（g）"
        count_label = "重抽样次数"
        t_label = r"$T^*$（$H_0$下）"
        observed_label = "观察差值"
        percentile_label = "百分位端点"
        observed_t_label = "观察到的 T"
        cutoff_label = "双侧极端界线"
    else:
        bootstrap_title = "Groupwise bootstrap: estimate"
        null_title = "Centered studentized null"
        permutation_title = "Permutation reallocation of 19 observations"
        difference_label = "A mean - B mean (g)"
        count_label = "Count"
        t_label = "T* under H0"
        observed_label = "observed"
        percentile_label = "percentile limits"
        observed_t_label = "observed T"
        cutoff_label = "two-sided cutoffs"

    fig1, axes = plt.subplots(2, 1, figsize=(6.4, 9.2)) if chinese else plt.subplots(1, 2, figsize=(12, 5))
    axes[0].hist(delta_star_values, bins=40, color="#9CCBEA", edgecolor="white")
    axes[0].axvline(delta_obs, color="#145A86", linewidth=2, linestyle="--", label=observed_label)
    axes[0].axvline(ci_percentile[0], color="#E67E22", linewidth=2, label=percentile_label)
    axes[0].axvline(ci_percentile[1], color="#E67E22", linewidth=2)
    axes[0].set_title(bootstrap_title, fontproperties=font)
    axes[0].set_xlabel(difference_label, fontproperties=font)
    axes[0].set_ylabel(count_label, fontproperties=font)
    axes[0].legend(prop=font, fontsize=9, frameon=False)
    axes[1].hist(t_null_values, bins=40, color="#B9D9C6", edgecolor="white")
    axes[1].axvline(t_obs, color="#E67E22", linewidth=2, label=observed_t_label)
    axes[1].set_title(null_title, fontproperties=font)
    axes[1].set_xlabel(t_label, fontproperties=font)
    axes[1].set_ylabel(count_label, fontproperties=font)
    axes[1].legend(prop=font, fontsize=9, frameon=False)
    for axis in axes:
        axis.tick_params(labelsize=12)
    fig1.tight_layout()
    fig1.savefig(target_dir / "01_两种重抽样与统计量.png", dpi=160)
    plt.close(fig1)

    fig2, ax = plt.subplots(figsize=(6.4, 5.5))
    ax.hist(permutation_values, bins=40, color="#C9D7F0", edgecolor="white")
    ax.axvline(-abs(delta_obs), color="#E67E22", linewidth=2, linestyle="--", label=cutoff_label)
    ax.axvline(abs(delta_obs), color="#E67E22", linewidth=2, linestyle="--")
    ax.axvline(delta_obs, color="#145A86", linewidth=2, label=observed_label)
    ax.set_title(permutation_title, fontproperties=font)
    ax.set_xlabel(difference_label, fontproperties=font)
    ax.set_ylabel("分组方案数" if chinese else "Number of partitions", fontproperties=font)
    ax.legend(prop=font, fontsize=9, frameon=False)
    ax.tick_params(labelsize=12)
    fig2.tight_layout()
    fig2.savefig(target_dir / "02_置换分布与观察差值.png", dpi=160)
    plt.close(fig2)


draw_figures(delta_star, t_null, permutation_delta, figure_dir)
draw_figures(delta_star, t_null, permutation_delta, project_figure_dir, chinese=True)
draw_figures(delta_star, t_null, permutation_delta, article_figure_dir, chinese=True)

# %% 7. 最后核对并写入简短记录
assert np.isclose(delta_obs, 19.1666666666667, atol=1e-10)
assert np.isclose(t_obs, 1.93193081703987, atol=1e-10)
assert left_count == 9636
assert right_count == 364
assert np.isclose(p_boot, 0.0728, atol=1e-12)
assert np.isclose(ci_t[0], -2.28522435160976, atol=1e-8)
assert np.isclose(ci_t[1], 40.92199759211, atol=1e-8)
assert np.isclose(ci_percentile[0], 1.04761904761904, atol=1e-8)
assert np.isclose(ci_percentile[1], 37.3693452380952, atol=1e-8)
assert np.isclose(coverage_exact, 0.6496866225, atol=1e-12)

run_notes = [
    "Python脚本完成：固定索引、中心化学生化Bootstrap、置换重分与率边界枚举。",
    f"delta={delta_obs:.10f}, T={t_obs:.10f}, Bootstrap P={p_boot:.4f}",
    f"Bootstrap-t CI=[{ci_t[0]:.10f}, {ci_t[1]:.10f}]",
    f"permutation exact partitions={permutation_rep}, P={perm_p:.10f}",
    f"rate coverage exact enumeration={coverage_exact:.10f}",
    f"Python {sys.version.split()[0]}; NumPy {np.__version__}; pandas {pd.__version__}; SciPy {scipy.__version__}; Matplotlib {matplotlib.__version__}",
]
(out_dir / "运行记录.txt").write_text("\n".join(run_notes) + "\n", encoding="utf-8")

# %% 8. 公众号公式图：数学排版避免星号被Markdown误认成斜体
# 这里只负责排版，不参与前面的数值计算。
from matplotlib.mathtext import math_to_image
from matplotlib.font_manager import FontProperties
formula_dir = project_figure_dir
formula_dir.mkdir(exist_ok=True)
math_to_image(r"$T^*=\frac{\overline{A_0^*}-\overline{B_0^*}}{\sqrt{s_{A^*}^2/12+s_{B^*}^2/7}}$",
              formula_dir / "03_Bootstrap学生化公式.png",
              prop=FontProperties(size=22), dpi=240, format="png", color="#1769aa")
math_to_image(r"$\mathrm{SE}_{\mathrm{obs}}=\sqrt{\frac{s_A^2}{12}+\frac{s_B^2}{7}}$",
              formula_dir / "04_观察标准误公式.png",
              prop=FontProperties(size=22), dpi=240, format="png", color="#1769aa")

if matplotlib.get_backend().lower() != "agg":
    plt.show()
