# %% 统计课堂43：方差分析的正态性、方差齐性：指什么，一定要检验吗？
import os
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.formula.api as smf
from statsmodels.stats.diagnostic import normal_ad

# %% 参数与读取：Spyder工作目录设为项目根目录或代码目录
project_root = Path.cwd()
if project_root.name == "代码":
    project_root = project_root.parent
# 先定位项目，再从代码目录加载配套函数，整文件和分节运行均使用同一入口。
helper_directory = str(project_root / "代码")
if helper_directory not in sys.path:
    sys.path.insert(0, helper_directory)
from residplot import residplot
data_path = project_root / "数据" / "residual_cases.csv"
dat = pd.read_csv(data_path)
print("字段：", list(dat.columns))
print(dat.head(6))
SAVE_RESULTS = os.getenv("SAVE_RESULTS", "false").lower() == "true"


# %% 场景A描述统计、真实误差、残差与手算ANOVA
def make_results(d, label):
    d = d.copy()
    d["true_error"] = d["value"] - d["mu"]
    d["fitted"] = d.groupby("group")["value"].transform("mean")
    d["residual"] = d["value"] - d["fitted"]
    summary = (
        d.groupby("group")["value"]
        .agg(n="size", mean="mean", sd=lambda x: x.std(ddof=1))
        .reset_index()
    )
    summary.insert(0, "scenario", label)
    grand = d["value"].mean()
    means = d.groupby("group")["value"].mean()
    ns = d.groupby("group").size()
    ssb = float((ns * (means - grand) ** 2).sum())
    sse = float((d["residual"] ** 2).sum())
    dfb = means.size - 1
    dfe = len(d) - means.size
    msb = ssb / dfb
    mse = sse / dfe
    Fvalue = msb / mse
    pvalue = stats.f.sf(Fvalue, dfb, dfe)
    anova = pd.DataFrame(
        [
            dict(
                scenario=label,
                SSB=ssb,
                SSE=sse,
                df_between=dfb,
                df_error=dfe,
                MSB=msb,
                MSE=mse,
                F=Fvalue,
                p=pvalue,
            )
        ]
    )
    return d, summary, anova


if dat["scenario"].nunique() != 1 or dat["scenario"].iloc[0] != "A":
    raise ValueError("residual_cases.csv must contain only scenario A")
parts = [make_results(dat, "A")]
group_summary = pd.concat([x[1] for x in parts], ignore_index=True)
anova_tab = pd.concat([x[2] for x in parts], ignore_index=True)
residuals_tab = pd.concat([x[0] for x in parts], ignore_index=True)
print("\n主例题A：")
print(group_summary)
print("\n手算ANOVA：")
print(anova_tab)
print("\n每组残差和：")
print(residuals_tab.groupby(["scenario", "group"])["residual"].sum())
print("\nscipy f_oneway equal_var=True（合并方差ANOVA）：",
      stats.f_oneway(*[g["value"].to_numpy() for _, g in dat.groupby("group")]))

# %% 场景A原始value的方差齐性检验（中位数中心的Brown-Forsythe修正版）
# H0：各组总体方差相等；P>=.05时仅表述为未拒绝H0，不等同于证明方差齐。
a = dat.loc[dat["scenario"] == "A"]
a_groups = [g["value"].to_numpy() for _, g in a.groupby("group", sort=False)]
levene = stats.levene(*a_groups, center="median")
bartlett = stats.bartlett(*a_groups)
variance_tests = pd.DataFrame([
    {"method": "Levene_median", "statistic": levene.statistic, "df1": len(a_groups) - 1,
     "df2": len(a) - len(a_groups), "p": levene.pvalue},
    {"method": "Bartlett", "statistic": bartlett.statistic, "df1": len(a_groups) - 1,
     "df2": np.nan, "p": bartlett.pvalue},
])
print("\n场景A原始value方差齐性检验（H0：各组总体方差相等）：")
for row in variance_tests.itertuples(index=False):
    if row.method == "Levene_median":
        print(f"Levene（中位数中心，Brown-Forsythe修正版）F({int(row.df1)},{int(row.df2)}) = {row.statistic:.6f}, p = {row.p:.6f}")
    else:
        print(f"Bartlett χ²({int(row.df1)}) = {row.statistic:.6f}, p = {row.p:.6f}")
    print("解释：P>=.05，未拒绝各组总体方差相等的H0；不表示已证明方差齐。" if row.p >= .05
          else "解释：P<.05，拒绝各组总体方差相等的H0。")

# %% 场景A各组原始value的正态性检验
# 每组n=5：Shapiro-Wilk实际运行；Anderson-Darling仅在n>=8时运行。
normality_rows = []
for _, g in a.groupby("group", sort=False):
    x = g["value"].to_numpy()
    sw_stat, sw_p = stats.shapiro(x)
    normality_rows.append({"scenario": "A", "group": g["group"].iloc[0], "n": len(x),
                           "method": "Shapiro_Wilk", "statistic": sw_stat, "p": sw_p,
                           "status": "not_rejected" if sw_p >= 0.05 else "rejected"})
    if len(x) >= 8:
        ad_stat, ad_p = normal_ad(x)
        ad_status = "not_rejected" if ad_p >= 0.05 else "rejected"
    else:
        ad_stat, ad_p, ad_status = np.nan, np.nan, "not_run_n_lt_8"
    normality_rows.append({"scenario": "A", "group": g["group"].iloc[0], "n": len(x),
                           "method": "Anderson_Darling", "statistic": ad_stat, "p": ad_p,
                           "status": ad_status})
normality_tests = pd.DataFrame(normality_rows,
    columns=["scenario", "group", "n", "method", "statistic", "p", "status"])
print("\n场景A各组原始value正态性检验（每组独立样本；不对合并残差作iid正态性P断言）：")
print(normality_tests)

# %% 场景A排序的Q-Q数值表
qq_parts = []
for label, frame in residuals_tab.groupby("scenario", sort=False):
    part = frame.sort_values("residual")[
        ["scenario", "group", "replicate", "residual"]
    ].reset_index(drop=True)
    part["rank"] = np.arange(1, len(part) + 1)
    part["theoretical"] = stats.norm.ppf((part["rank"] - 0.5) / len(part))
    qq_parts.append(part)
qq = pd.concat(qq_parts, ignore_index=True)
qq = qq[["scenario", "group", "replicate", "rank", "residual", "theoretical"]]
print("\nQ-Q表：场景A的N=15排序；不作独立同分布Shapiro P断言：")
print(qq)

# %% 场景A的诊断图：组别保留在点的颜色中
# 在Spyder中查看图；批处理可SHOW_PLOTS=false，数学计算仍执行。
SHOW_PLOTS = os.getenv("SHOW_PLOTS", "true").lower() == "true"
if SHOW_PLOTS:
    import matplotlib.pyplot as plt
    import seaborn as sns

    sns.set_theme(style="whitegrid")
    a = residuals_tab.loc[residuals_tab["scenario"] == "A"]
    qa = qq.loc[qq["scenario"] == "A"]
    yq = np.quantile(a["residual"], [0.25, 0.75], method="linear")
    zq = stats.norm.ppf([0.25, 0.75])
    slope = (yq[1] - yq[0]) / (zq[1] - zq[0])
    intercept = yq[0] - slope * zq[0]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    for group in ["A1", "A2", "A3"]:
        g = a.loc[a["group"] == group]
        qg = qa.loc[qa["group"] == group]
        axes[0].scatter(g["fitted"], g["residual"], label=group)
        axes[1].scatter(qg["theoretical"], qg["residual"], label=group)
    axes[0].axhline(0, color="grey", linestyle="--")
    xline = np.array([qa["theoretical"].min(), qa["theoretical"].max()])
    axes[1].plot(xline, intercept + slope * xline, color="grey", linestyle="--")
    axes[0].set(xlabel="Fitted group mean (MPa)", ylabel="Residual (MPa)")
    axes[1].set(xlabel="Normal theoretical quantile", ylabel="Residual (MPa)")
    axes[0].legend()
    axes[1].legend()
    plt.show()
    # 与第42篇同一函数：外部学生化残差直方图、拟合正态曲线和Gaussian KDE。
    fit_a = smf.ols("value ~ C(group)", data=a).fit()
    residplot(fit_a, nbreaks=5, draw=True)
    plt.show()

# %% 按需保存统一结果
if SAVE_RESULTS:
    out = project_root / "结果" / "Python"
    out.mkdir(parents=True, exist_ok=True)
    group_summary.to_csv(out / "group_summary.csv", index=False)
    residuals_tab.to_csv(out / "residuals.csv", index=False)
    anova_tab.to_csv(out / "anova.csv", index=False)
    qq.to_csv(out / "qq.csv", index=False)
    variance_tests.to_csv(out / "variance_tests.csv", index=False)
    normality_tests.to_csv(out / "normality_tests.csv", index=False)
