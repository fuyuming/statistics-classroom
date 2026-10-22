# %% 小样本率的置信区间：从6/13找到上下限
# 目标：固定数据文件中的观察成功数x和试验次数n（默认x=6、n=13），直接用二项尾概率反演95%区间。
# 来源：颜艳《医学统计学》第5版例6-2、例6-5，第92—94页；这是教材案例参数，不是个人原始数据。
# 运行：Spyder打开本文件后按F5；依赖见项目根目录requirements.txt。
# 编码：UTF-8。主要结果对象保留在工作区，便于初学者直接查看。

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
import pandas as pd
from scipy.stats import beta, binom, norm

root = Path(__file__).resolve().parents[1]
out = root / "运行结果" / "Python"
out.mkdir(parents=True, exist_ok=True)
project_fig_dir = root / "文章配图"
article_fig_dir = project_fig_dir  # 输出限定在本项目，不写到相邻目录
project_fig_dir.mkdir(parents=True, exist_ok=True)
article_fig_dir.mkdir(parents=True, exist_ok=True)

parameters = pd.read_csv(root / "数据" / "案例参数.csv", encoding="utf-8")
x = int(parameters.loc[0, "x"])
n = int(parameters.loc[0, "n"])
confidence = float(parameters.loc[0, "confidence"])
alpha = float(parameters.loc[0, "alpha"])
scan_step = float(parameters.loc[0, "scan_step"])
display_rounds = int(parameters.loc[0, "bisection_rounds"])
tail_target = alpha / 2

# %% 1. 先看样本率与Wald近似区间
sample_proportion = x / n
wald_se = np.sqrt(sample_proportion * (1 - sample_proportion) / n)
# ppf输入累计概率，返回分位点；95%双侧查0.975，得到约1.959964。
z_975 = norm.ppf(1 - tail_target)
normal_central_area = norm.cdf(z_975) - norm.cdf(-z_975)
# 中间的面积=0.95；PDF曲线的高度本身不是概率。
wald_low = sample_proportion - z_975 * wald_se
wald_high = sample_proportion + z_975 * wald_se

wald_result = pd.DataFrame(
    {
        "quantity": ["sample_proportion", "wald_se", "wald_low", "wald_high"],
        "value": [sample_proportion, wald_se, wald_low, wald_high],
    }
)
wald_result.to_csv(out / "Wald近似.csv", index=False, encoding="utf-8-sig")
# 参考：sample_proportion=0.46153846，Wald区间约为0.1905457—0.7325312。
# 这是正态近似结果；下面的端点直接使用二项尾概率。

# %% 1.1 Wilson区间：反演使用候选p0标准误的score检验
# 解 n*(p_hat-p0)^2 <= z^2*p0*(1-p0)，两根之间就是区间。
wilson_denominator = 1 + z_975**2 / n
wilson_center = (sample_proportion + z_975**2 / (2*n)) / wilson_denominator
wilson_half_width = z_975 * np.sqrt(sample_proportion*(1-sample_proportion)/n + z_975**2/(4*n**2)) / wilson_denominator
wilson_low = wilson_center - wilson_half_width
wilson_high = wilson_center + wilson_half_width
# 默认6/13：Wilson区间[0.2320607,0.7085620]；与反正弦变换无关。
normal_wilson = pd.DataFrame({
    "quantity": ["z_critical", "normal_central_area", "wilson_low", "wilson_high"],
    "value": [z_975, normal_central_area, wilson_low, wilson_high]})
normal_wilson.to_csv(out / "正态分位点与Wilson.csv", index=False, encoding="utf-8-sig")
normal_wilson

# %% 2. 固定输入文件中的x，扫描候选总体率p
# sf(x-1)计算P_p(X>=x)，必须包含观察到的x。
# cdf(x)计算P_p(X<=x)，同样包含观察到的x。
p_grid = np.round(np.arange(0, 1 + scan_step / 2, scan_step), 10)
right_tail = binom.sf(x - 1, n, p_grid)
upper_tail = binom.cdf(x, n, p_grid)
scan = pd.DataFrame(
    {
        "candidate_p": p_grid,
        "right_tail_P_X_ge_x": right_tail,
        "upper_tail_P_X_le_x": upper_tail,
    }
)
scan.to_csv(out / "扫描表.csv", index=False, encoding="utf-8-sig")

# 找到离0.025最近的两端夹逼区间，供后面的二分使用。
# x=0时，Pr_p(X>=0)=1，不可能等于0.025；下限直接是参数空间边界p=0。
# x=n时，Pr_p(X<=n)=1，不可能等于0.025；上限直接是参数空间边界p=1。
scan_bracket_rows = []
if x == 0:
    lower_below = None
    lower_above = None
    scan_bracket_rows.append(
        {
            "bound": "lower",
            "status": "boundary fixed at p=0; Pr_p(X>=0)=1; no interior root",
            "candidate_low": np.nan,
            "tail_at_low": 1.0,
            "candidate_high": np.nan,
            "tail_at_high": 1.0,
        }
    )
else:
    lower_below = scan.loc[scan["right_tail_P_X_ge_x"] <= tail_target].iloc[-1]
    lower_above = scan.loc[scan["right_tail_P_X_ge_x"] >= tail_target].iloc[0]
    scan_bracket_rows.append(
        {
            "bound": "lower",
            "status": "interior root",
            "candidate_low": lower_below["candidate_p"],
            "tail_at_low": lower_below["right_tail_P_X_ge_x"],
            "candidate_high": lower_above["candidate_p"],
            "tail_at_high": lower_above["right_tail_P_X_ge_x"],
        }
    )
if x == n:
    upper_below = None
    upper_above = None
    scan_bracket_rows.append(
        {
            "bound": "upper",
            "status": "boundary fixed at p=1; Pr_p(X<=n)=1; no interior root",
            "candidate_low": np.nan,
            "tail_at_low": 1.0,
            "candidate_high": np.nan,
            "tail_at_high": 1.0,
        }
    )
else:
    upper_below = scan.loc[scan["upper_tail_P_X_le_x"] <= tail_target].iloc[0]
    upper_above = scan.loc[scan["upper_tail_P_X_le_x"] >= tail_target].iloc[-1]
    scan_bracket_rows.append(
        {
            "bound": "upper",
            "status": "interior root",
            "candidate_low": upper_above["candidate_p"],
            "tail_at_low": upper_above["upper_tail_P_X_le_x"],
            "candidate_high": upper_below["candidate_p"],
            "tail_at_high": upper_below["upper_tail_P_X_le_x"],
        }
    )
scan_brackets = pd.DataFrame(scan_bracket_rows)
scan_brackets.to_csv(out / "扫描夹逼.csv", index=False, encoding="utf-8-sig")

# %% 3. 从扫描夹逼区间开始二分
def bisect_tail(lower, upper, tail_function, target, side, shown_rounds, total_rounds=60):
    """重复寻找尾概率=target的p；前shown_rounds轮写入教学历史。"""
    history = []
    for round_number in range(1, total_rounds + 1):
        candidate = (lower + upper) / 2
        tail = float(tail_function(candidate))
        if side == "lower":
            # P_p(X>=x)随p增大；tail较小，端点在右侧。
            if tail <= target:
                lower = candidate
            else:
                upper = candidate
        else:
            # P_p(X<=x)随p增大而减小；tail较大，端点在右侧。
            if tail >= target:
                lower = candidate
            else:
                upper = candidate
        if round_number <= shown_rounds:
            history.append(
                {
                    "bound": side,
                    "status": "interior root",
                    "round": round_number,
                    "candidate_p": candidate,
                    "tail_probability": tail,
                    "bracket_low_after": lower,
                    "bracket_high_after": upper,
                }
            )
    endpoint = (lower + upper) / 2
    return endpoint, pd.DataFrame(history)

if x == 0:
    lower_endpoint = 0.0
    lower_history = pd.DataFrame(
        [
            {
                "bound": "lower",
                "status": "boundary fixed at p=0; no interior root",
                "round": 0,
                "candidate_p": 0.0,
                "tail_probability": 1.0,
                "bracket_low_after": np.nan,
                "bracket_high_after": np.nan,
            }
        ]
    )
else:
    lower_endpoint, lower_history = bisect_tail(
        float(lower_below["candidate_p"]),
        float(lower_above["candidate_p"]),
        lambda candidate: binom.sf(x - 1, n, candidate),
        tail_target,
        "lower",
        display_rounds,
    )
if x == n:
    upper_endpoint = 1.0
    upper_history = pd.DataFrame(
        [
            {
                "bound": "upper",
                "status": "boundary fixed at p=1; no interior root",
                "round": 0,
                "candidate_p": 1.0,
                "tail_probability": 1.0,
                "bracket_low_after": np.nan,
                "bracket_high_after": np.nan,
            }
        ]
    )
else:
    upper_endpoint, upper_history = bisect_tail(
        float(upper_above["candidate_p"]),
        float(upper_below["candidate_p"]),
        lambda candidate: binom.cdf(x, n, candidate),
        tail_target,
        "upper",
        display_rounds,
    )
bisection_history = pd.concat([lower_history, upper_history], ignore_index=True)
bisection_history.to_csv(out / "二分历史.csv", index=False, encoding="utf-8-sig")

# beta分位点是Clopper–Pearson端点的另一种数值表达，用于核对60轮二分结果。
cp_lower = 0.0 if x == 0 else beta.ppf(tail_target, x, n - x + 1)
cp_upper = 1.0 if x == n else beta.ppf(1 - tail_target, x + 1, n - x)

# %% 4. 文章中的两个补充核对：9/10 vs p=0.7，以及0/10边界
reference_x = 9
reference_n = 10
reference_p = 0.7
observed_mass = binom.pmf(reference_x, reference_n, reference_p)
all_masses = binom.pmf(np.arange(reference_n + 1), reference_n, reference_p)
# R binom.test的双侧P值按不大于观察结果概率的结果求和。
reference_p_value = float(all_masses[all_masses <= observed_mass + 1e-14].sum())
reference_ci_low = beta.ppf(tail_target, reference_x, reference_n - reference_x + 1)
reference_ci_high = beta.ppf(1 - tail_target, reference_x + 1, reference_n - reference_x)

zero_n = 10
zero_two_sided_upper = 1 - tail_target ** (1 / zero_n)
zero_one_sided_upper = 1 - (alpha ** (1 / zero_n))

key_results = pd.DataFrame(
    {
        "quantity": [
            "sample_proportion",
            "wald_low",
            "wald_high",
            "cp_lower_bisection_60",
            "cp_upper_bisection_60",
            "cp_lower_beta_check",
            "cp_upper_beta_check",
            "reference_9_of_10_p_value_at_0.7",
            "reference_9_of_10_cp_low",
            "reference_9_of_10_cp_high",
            "zero_of_10_two_sided_upper",
            "zero_of_10_one_sided_upper",
        ],
        "value": [
            sample_proportion,
            wald_low,
            wald_high,
            lower_endpoint,
            upper_endpoint,
            cp_lower,
            cp_upper,
            reference_p_value,
            reference_ci_low,
            reference_ci_high,
            zero_two_sided_upper,
            zero_one_sided_upper,
        ],
    }
)
key_results.to_csv(out / "关键结果.csv", index=False, encoding="utf-8-sig")

reference_9_of_10 = pd.DataFrame(
    {
        "quantity": ["two_sided_p_value_at_p0_0.7", "cp95_low", "cp95_high"],
        "value": [reference_p_value, reference_ci_low, reference_ci_high],
    }
)
reference_9_of_10.to_csv(out / "9_10参考.csv", index=False, encoding="utf-8-sig")
zero_reference = pd.DataFrame(
    {
        "quantity": ["two_sided_95_upper", "one_sided_95_upper"],
        "value": [zero_two_sided_upper, zero_one_sided_upper],
    }
)
zero_reference.to_csv(out / "0_10边界.csv", index=False, encoding="utf-8-sig")

# %% 5. 生成两张排除区域图
# 数值CSV已经写出后才检查中文字体；无字体时不阻断前面的计算结果。
from matplotlib import font_manager
font_names = ["Microsoft YaHei", "SimHei", "Noto Sans CJK SC", "PingFang SC", "Heiti SC", "Arial Unicode MS", "WenQuanYi Zen Hei"]
font_path = None
for font_name in font_names:
    try:
        font_path = font_manager.findfont(font_name, fallback_to_default=False)
        break
    except ValueError:
        continue
if font_path is None:
    raise RuntimeError("未找到中文字体。数值CSV已写出；请安装 Noto Sans CJK SC 或使用系统微软雅黑/黑体后重跑配图。")
chinese_font = font_manager.FontProperties(fname=font_path).get_name()
plt.rcParams["font.family"] = chinese_font
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["font.size"] = 14

def save_twice(figure, filename):
    figure.savefig(project_fig_dir / filename, dpi=180, bbox_inches="tight")
    figure.savefig(article_fig_dir / filename, dpi=180, bbox_inches="tight")
    plt.close(figure)

plot_p = np.linspace(0, 1, 1001)
plot_right = binom.sf(x - 1, n, plot_p)
fig_lower, ax_lower = plt.subplots(figsize=(6.4, 4.8))
ax_lower.plot(plot_p, plot_right, color="#1769aa", linewidth=2, label=f"P(X ≥ {x} | p)")
if x == 0:
    ax_lower.axvline(0, color="#333333", linestyle=":", linewidth=1.3, label="pL = 0边界")
    ax_lower.text(
        0.04,
        0.08,
        "pL = 0：边界\nPr_p(X≥0)=1，无内点交点",
        transform=ax_lower.transAxes,
        fontsize=11,
        bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "none"},
    )
else:
    ax_lower.axvspan(0, cp_lower, color="#f2b6a0", alpha=0.35, label="排除区域 p < pL")
    ax_lower.axvline(cp_lower, color="#333333", linestyle=":", linewidth=1.3)
    ax_lower.scatter([cp_lower], [tail_target], color="#c23b22", zorder=3)
    ax_lower.annotate(
        f"pL = {cp_lower:.4f}",
        xy=(cp_lower, tail_target),
        xytext=(0.36, 0.055),
        arrowprops={"arrowstyle": "->", "color": "#333333"},
    )
ax_lower.axhline(tail_target, color="#c23b22", linestyle="--", linewidth=1.5, label="0.025")
ax_lower.set_xlabel("候选总体率 p")
ax_lower.set_ylabel("尾概率")
ax_lower.set_title(
    "下限：pL = 0（边界）" if x == 0 else f"下限：P(X ≥ {x} | p) = {tail_target:g}"
)
ax_lower.set_xlim(0, 1)
# 局部放大0—0.10，便于在手机图上看清0.025交点；高于0.10的曲线段有意裁掉。
ax_lower.set_ylim(0, 1.05 if x == 0 else 0.10)
ax_lower.grid(alpha=0.25)
ax_lower.legend(frameon=False)
save_twice(fig_lower, "01_下限排除尾概率.png")

plot_upper = binom.cdf(x, n, plot_p)
fig_upper, ax_upper = plt.subplots(figsize=(6.4, 4.8))
ax_upper.plot(plot_p, plot_upper, color="#1769aa", linewidth=2, label=f"P(X ≤ {x} | p)")
if x == n:
    ax_upper.axvline(1, color="#333333", linestyle=":", linewidth=1.3, label="pU = 1边界")
    ax_upper.text(
        0.04,
        0.08,
        "pU = 1：边界\nPr_p(X≤n)=1，无内点交点",
        transform=ax_upper.transAxes,
        fontsize=11,
        bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "none"},
    )
else:
    ax_upper.axvspan(cp_upper, 1, color="#f2b6a0", alpha=0.35, label="排除区域 p > pU")
    ax_upper.axvline(cp_upper, color="#333333", linestyle=":", linewidth=1.3)
    ax_upper.scatter([cp_upper], [tail_target], color="#c23b22", zorder=3)
    ax_upper.annotate(
        f"pU = {cp_upper:.4f}",
        xy=(cp_upper, tail_target),
        xytext=(0.24, 0.055),
        arrowprops={"arrowstyle": "->", "color": "#333333"},
    )
ax_upper.axhline(tail_target, color="#c23b22", linestyle="--", linewidth=1.5, label="0.025")
ax_upper.set_xlabel("候选总体率 p")
ax_upper.set_ylabel("尾概率")
ax_upper.set_title(
    "上限：pU = 1（边界）" if x == n else f"上限：P(X ≤ {x} | p) = {tail_target:g}"
)
ax_upper.set_xlim(0, 1)
# 局部放大0—0.10，便于在手机图上看清0.025交点；高于0.10的曲线段有意裁掉。
ax_upper.set_ylim(0, 1.05 if x == n else 0.10)
ax_upper.grid(alpha=0.25)
ax_upper.legend(frameon=False)
save_twice(fig_upper, "02_上限排除尾概率.png")

# Spyder中可直接查看这些对象；邻近注释给出文章核对值。
validation_objects = {
    "cp_interval": (cp_lower, cp_upper),
    "reference_p_value": reference_p_value,
    "zero_of_10_upper": (zero_two_sided_upper, zero_one_sided_upper),
}

# %% 6. 文章补图：1.96来源、二项近似与公式推导
# 排版代码单独保存，主流程仍按统计问题分节。示意图使用文中固定参数。
from plot_interval_explanations import draw_explanations
draw_explanations(project_fig_dir, out)
