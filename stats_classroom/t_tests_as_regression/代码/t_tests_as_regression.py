# -*- coding: utf-8 -*-
# 学习方式：按分节逐段运行，先读注释，再查看结果对象；预期数值用于核对。
"""统计课堂：t检验，原来也是回归。

Spyder：打开本文件后按 F5 完整运行；命令行运行也可以。
依赖 numpy、pandas、scipy、statsmodels、matplotlib；共享数据位于 数据/。
每一行主要步骤都保留为可检查的结果对象，数值另写入 运行结果/Python/。
图保存在项目的 文章配图/，并在 Python 输出目录留一份供核对。
"""

# %% 准备：定位项目、读取数据、先认识长表与配对表
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats

root = Path(__file__).resolve().parents[1]
out = root / "运行结果" / "Python"
fig_dir = root / "文章配图"
out.mkdir(parents=True, exist_ok=True)
fig_dir.mkdir(parents=True, exist_ok=True)

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
if font_path is None:
    raise RuntimeError("未找到中文字体。请安装 Noto Sans CJK SC 或使用系统微软雅黑/黑体，重启Python后运行；数值计算不需要中文字体。")
chinese_font = font_manager.FontProperties(fname=font_path).get_name()
plt.rcParams["font.family"] = chinese_font
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["font.size"] = 14

rat = pd.read_csv(root / "数据" / "rat_gain_long.csv")
paired = pd.read_csv(root / "数据" / "paired_drinks.csv")

# 一行是一只大鼠；group01=1表示A高蛋白组，0表示B低蛋白组。
rat_head = rat.head()
paired_head = paired.head()
rat_head
paired_head


def coefficient_record(fit, example, term_index=0):
    """把模型系数、t检验和95%区间整理成可核对的结构对象。"""
    params = np.asarray(fit.params, dtype=float)
    bse = np.asarray(fit.bse, dtype=float)
    tvalues = np.asarray(fit.tvalues, dtype=float)
    pvalues = np.asarray(fit.pvalues, dtype=float)
    interval = np.asarray(fit.conf_int(alpha=0.05), dtype=float)
    term_names = list(fit.params.index) if hasattr(fit.params, "index") else ["Intercept"]
    return {
        "example": example,
        "term": term_names[term_index],
        "estimate": params[term_index],
        "standard_error": bse[term_index],
        "df": float(fit.df_resid),
        "t_value": tvalues[term_index],
        "p_value": pvalues[term_index],
        "ci_low": interval[term_index, 0],
        "ci_high": interval[term_index, 1],
    }


# %% 1. 单样本t：只有截距的模型
# 公式y ~ 1只放截距，.fit()才执行拟合；与R的lm公式含义相同。
b_data = rat.loc[rat["group"] == "B"].copy()
b_y = b_data["gain_g"].to_numpy(dtype=float)
fit_b_default = smf.ols("gain_g ~ 1", data=b_data).fit()
# 截距101，95%CI[81.926,120.074]，默认检验均数=0。
mean_summary = fit_b_default.summary()
mean_summary

# 默认截距检验的是均数=0。把参考值移到0，就可以检验均数=mu0。
# 100 g是本项目人为提出的教学参照，不是教材例4-8原问题。
mu0 = 100.0
b_data["gain_minus_reference"] = b_data["gain_g"] - mu0
fit_b_reference = smf.ols("gain_minus_reference ~ 1", data=b_data).fit()
# 现在截距1表示均数减100，P约0.902；差值区间包含0。
single_results = pd.DataFrame(
    [
        coefficient_record(fit_b_default, "single_B_vs_0"),
        coefficient_record(fit_b_reference, "single_B_vs_100"),
    ]
)
single_results.to_csv(out / "01_single_sample.csv", index=False)

fit_b_default
fit_b_reference

# %% 2. 两独立样本：0/1组变量的普通OLS
# B组是0，所以截距是B均数；A组是1，所以group01系数是A-B。
# statsmodels默认使用普通同方差OLS推断，这正是合并方差两独立样本t的接口。
fit_groups = smf.ols("gain_g ~ group01", data=rat).fit()
group_rows = [coefficient_record(fit_groups, "two_independent_groups_pooled_OLS", i)
              for i in range(len(fit_groups.params))]
group_results = pd.DataFrame(group_rows)
group_results.to_csv(out / "02_group_coefficients.csv", index=False)

overall_F_result = pd.DataFrame(
    [{
        "model": "group_OLS",
        "F_value": float(fit_groups.fvalue),
        "numerator_df": float(fit_groups.df_model),
        "denominator_df": float(fit_groups.df_resid),
        "p_value": float(fit_groups.f_pvalue),
    }]
)
overall_F_result.to_csv(out / "03_group_overall_F.csv", index=False)

# 用scipy的equal_var=True和False分别核对合并方差t与Welch t。
a = rat.loc[rat["group"] == "A", "gain_g"].to_numpy(dtype=float)
b = rat.loc[rat["group"] == "B", "gain_g"].to_numpy(dtype=float)
pooled_check = stats.ttest_ind(a, b, equal_var=True)
welch_check = stats.ttest_ind(a, b, equal_var=False)
v_a = np.var(a, ddof=1)
v_b = np.var(b, ddof=1)
welch_se = np.sqrt(v_a / len(a) + v_b / len(b))
welch_df = (v_a / len(a) + v_b / len(b)) ** 2 / (
    (v_a / len(a)) ** 2 / (len(a) - 1)
    + (v_b / len(b)) ** 2 / (len(b) - 1)
)
estimate_ab = float(a.mean() - b.mean())
pooled_ci = tuple(np.asarray(pooled_check.confidence_interval(), dtype=float))
welch_ci = (
    estimate_ab - stats.t.ppf(0.975, welch_df) * welch_se,
    estimate_ab + stats.t.ppf(0.975, welch_df) * welch_se,
)
welch_result = pd.DataFrame(
    [
        {
            "method": "ordinary_OLS_equal_variance",
            "estimate_A_minus_B": estimate_ab,
            "t_value": float(pooled_check.statistic),
            "df": float(len(a) + len(b) - 2),
            "p_value": float(pooled_check.pvalue),
            "ci_low": pooled_ci[0],
            "ci_high": pooled_ci[1],
        },
        {
            "method": "Welch_unequal_variance",
            "estimate_A_minus_B": estimate_ab,
            "t_value": float(welch_check.statistic),
            "df": float(welch_df),
            "p_value": float(welch_check.pvalue),
            "ci_low": welch_ci[0],
            "ci_high": welch_ci[1],
        },
    ]
)
welch_result.to_csv(out / "04_pooled_vs_welch.csv", index=False)

fit_groups
group_results

# %% 3. 配对t：先按ID作差，再拟合只有截距的模型
# 同一个人的两次评分不是两个独立组。delta保留了ID配对关系；随后检验平均delta=0。
paired["delta"] = paired["new"] - paired["old"]
fit_paired = smf.ols("delta ~ 1", data=paired).fit()
# 每行一个人的配对差值：平均差1.3分，95%CI[0.231,2.369]。
paired_check = stats.ttest_1samp(paired["delta"], popmean=0.0)
paired_result = pd.DataFrame([coefficient_record(fit_paired, "paired_new_minus_old")])
paired_result["check_t_value"] = float(paired_check.statistic)
paired_result["check_p_value"] = float(paired_check.pvalue)
paired_result.to_csv(out / "05_paired_difference.csv", index=False)
paired.to_csv(out / "05_paired_with_delta.csv", index=False)

fit_paired
paired_result

# %% 4. 图：看见截距+组差，以及新旧款的配对连线
def save_figure(fig, name):
    """同一组参数同时保存给项目和自动核对目录。"""
    fig.savefig(fig_dir / name, dpi=160, bbox_inches="tight")
    fig.savefig(out / name, dpi=160, bbox_inches="tight")


# 用固定的横向抖动只分开点，不改变任何纵坐标。
fig, ax = plt.subplots(figsize=(6.4, 5.2))
x_b = np.zeros(len(b)) + np.linspace(-0.11, 0.11, len(b))
x_a = np.ones(len(a)) + np.linspace(-0.11, 0.11, len(a))
ax.scatter(x_b, b, color="#1769AA", label="B组观测")
ax.scatter(x_a, a, color="#D65A4A", label="A组观测")
ax.hlines(b.mean(), -0.16, 0.16, color="#1769AA", linewidth=4)
ax.hlines(a.mean(), 0.84, 1.16, color="#D65A4A", linewidth=4)
ax.text(0, b.mean() + 8, f"截距 = {b.mean():.2f}", ha="center", color="#1769AA")
ax.text(0.50, 145, f"A−B = {a.mean() - b.mean():.2f} g", ha="center", color="#D65A4A", bbox={"facecolor": "white", "edgecolor": "none"})
ax.set(xlim=(-0.25, 1.25), xlabel="组别编码（0=B，1=A）", ylabel="增重量（g）",
       title="两组 OLS：截距与组间差")
ax.set_xticks([0, 1], ["B（0）", "A（1）"])
ax.legend(frameon=False)
fig.tight_layout()
save_figure(fig, "01_组均数与OLS截距.png")

fig, ax = plt.subplots(figsize=(6.4, 5.2))
for row in paired.itertuples(index=False):
    ax.plot([1, 2], [row.old, row.new], color="#9AA6B2", linewidth=1)
ax.scatter(np.ones(len(paired)), paired["old"], color="#1769AA", label="旧款")
ax.scatter(np.ones(len(paired)) * 2, paired["new"], color="#D65A4A", label="新款")
ax.scatter([1, 2], [paired["old"].mean(), paired["new"].mean()],
           color="#222222", marker="D", s=60, label="均数")
ax.text(1.5, max(paired["old"].max(), paired["new"].max()) - 0.25,
        f"平均差（新款−旧款）= {paired['delta'].mean():.2f}", ha="center",
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.85))
ax.set(xlim=(0.75, 2.25), xlabel="款式", ylabel="评分",
       title="同一消费者的旧款与新款评分")
ax.set_xticks([1, 2], ["旧款", "新款"])
ax.legend(frameon=False)
fig.tight_layout()
save_figure(fig, "02_新旧款配对与差值.png")

# %% 5. 汇总：把结构化结果交给人核验，不自动生成结论
key_results = pd.concat(
    [
        single_results.loc[single_results["example"] == "single_B_vs_100"],
        group_results.loc[group_results["term"].astype(str).str.contains("group01")],
        paired_result.drop(columns=["check_t_value", "check_p_value"]),
    ],
    ignore_index=True,
)
key_results.to_csv(out / "all_key_results.csv", index=False)

results = {
    "single_sample": single_results,
    "group_coefficients": group_results,
    "overall_F": overall_F_result,
    "pooled_vs_welch": welch_result,
    "paired_difference": paired_result,
    "key_results": key_results,
}
results

# Spyder会保留图形窗口；无图形环境下可用 MPLBACKEND=Agg 运行并检查PNG。
if plt.get_backend().lower() != "agg":
    plt.show()
