# 模拟教学数据：同一非致病微生物在两类培养基中的固定终点OD600。
# 每行假设为一个独立培养物；不是同一培养物的技术重复，也不是时间序列。
# OD600已按教学设定扣除空白，是无量纲浊度读数，不直接等于活菌数或生长速率。
# 明哥的微生物世界｜Python + Spyder 快速入门｜2026-09-17
# 自拟教学数据，不是真实实验。完整运行本文件，不依赖控制台已有变量。
# 依赖：pandas、numpy、matplotlib。重复运行只覆盖 output/Python 中的练习结果。
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

# 1. __file__是正在运行的脚本；Path让Windows/macOS都能使用相同路径写法。
ROOT = Path(__file__).resolve().parent
OUT = ROOT / "output" / "Python"
OUT.mkdir(parents=True, exist_ok=True)

# 2. 读取数据。CSV空白自动识别为NaN；不把缺失值改成0。
microbes = pd.read_csv(ROOT / "data" / "microbes_od600_demo.csv", encoding="utf-8-sig")
print(microbes.head())
microbes.info()
assert {"culture_id", "group", "od600"}.issubset(microbes.columns)
assert not microbes["culture_id"].duplicated().any()
assert microbes["group"].isin(["medium_A", "medium_B"]).all()
assert pd.api.types.is_numeric_dtype(microbes["od600"])
assert (np.isfinite(microbes["od600"]) | microbes["od600"].isna()).all()
print(microbes.isna().sum())
print(microbes.loc[microbes["od600"].isna()])
groups = ["medium_A", "medium_B"]
microbes["group"] = pd.Categorical(microbes["group"], categories=groups, ordered=True)

# 3. 分组统计。Python用缩进划分函数/循环，不使用R或MATLAB的花括号/end。
def summarise_one(raw):
    """raw是一组数值；先记录缺失，再用有效观测计算统计量。"""
    x = raw.dropna()
    if len(x) < 2:
        raise ValueError("每组至少需要2个有效观测。")
    # ddof=1明确使用样本SD；NumPy std的默认ddof=0，不能不加区分地照抄。
    # linear分位数等价于R type=7，与MATLAB脚本的显式插值保持一致。
    q = x.quantile([0.25, 0.75], interpolation="linear")
    return dict(n_total=len(raw), n=len(x), n_missing=int(raw.isna().sum()),
                mean=x.mean(), sd=x.std(ddof=1), median=x.median(),
                q1=q.loc[0.25], q3=q.loc[0.75], iqr=q.loc[0.75]-q.loc[0.25],
                min=x.min(), max=x.max())

rows = []
for group in groups:
    raw = microbes.loc[microbes["group"] == group, "od600"]
    rows.append(dict(group=group, **summarise_one(raw)))
summary_table = pd.DataFrame(rows)
print(summary_table.to_string(index=False))
summary_table.to_csv(OUT / "summary.csv", index=False, encoding="utf-8-sig")

# 4. 画散点与均值。横向偏移仅让点更容易辨认，不改变测量值。
colors = ["#267D87", "#D69645"]
fig, ax = plt.subplots(figsize=(8, 5), layout="constrained")
for idx, group in enumerate(groups, start=1):
    x = microbes.loc[microbes["group"] == group, "od600"].dropna().to_numpy()
    ax.scatter(idx + np.linspace(-0.06, 0.06, len(x)), x, color=colors[idx-1], s=60)
    ax.hlines(x.mean(), idx-.13, idx+.13, color=colors[idx-1], linewidth=3)
ax.set(xticks=[1,2], xticklabels=["Medium A", "Medium B"], xlim=(.6,2.4),
       ylim=(.3,.65), xlabel="Group", ylabel="OD600 (blank-corrected)",
       title="Microbial OD600: synthetic teaching data")
ax.spines[["top", "right"]].set_visible(False)
fig.savefig(OUT / "group_plot.png", dpi=200)
# 不调用close：Spyder中可以通过Plots窗格查看；同时保留独立图片文件。

# 5. 同样分箱的分组直方图；小样本下不要仅凭直方图判定分布类型。
fig_hist, axes = plt.subplots(1, 2, figsize=(9,4.5), layout="constrained")
for idx, group in enumerate(groups):
    x = microbes.loc[microbes["group"] == group, "od600"].dropna()
    axes[idx].hist(x, bins=np.linspace(.3,.65,8), color=colors[idx], edgecolor="white")
    axes[idx].set(title=group, xlabel="OD600 (blank-corrected)", ylabel="Count", ylim=(0,4))
fig_hist.savefig(OUT / "histogram.png", dpi=200)
assert summary_table["n_total"].sum() == 16
assert summary_table["n"].sum() == 14
assert summary_table["n_missing"].sum() == 2
(OUT / "versions.txt").write_text(
    f"Python: {sys.version}\n"
    f"pandas: {pd.__version__}\nnumpy: {np.__version__}\nmatplotlib: {matplotlib.__version__}\n",
    encoding="utf-8")
print(f"完成：结果保存在 {OUT}")
plt.show()
