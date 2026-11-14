"""统计课堂27｜独立绘制图1/图2；不依赖已有运行结果文件。"""
# %% 问题/条件：按同一给定q模型生成可画的数据表
from pathlib import Path
import math
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import norm

SAVE_FIGURES = False  # 如需保存图形，将此开关改为True。


def exact_table(n, q_plus, q_minus):
    q_zero = 1 - q_plus - q_minus
    probability = np.array([1.0])
    for _ in range(n):
        probability = np.convolve(probability, [q_minus, q_zero, q_plus])
    k = np.arange(-n, n + 1)
    d = k / n
    theta = q_plus - q_minus
    variance_D = (q_plus + q_minus - theta ** 2) / n
    sd = math.sqrt(variance_D)
    edges = (d[:-1] + d[1:]) / 2
    normal_cell_probability = np.empty(2 * n + 1)
    normal_cell_probability[0] = norm.cdf(edges[0], theta, sd)
    normal_cell_probability[-1] = 1 - norm.cdf(edges[-1], theta, sd)
    normal_cell_probability[1:-1] = norm.cdf(edges[1:], theta, sd) - norm.cdf(edges[:-1], theta, sd)
    return pd.DataFrame({"k": k, "D": d, "exact_probability": probability,
                         "normal_cell_probability": normal_cell_probability})


# %% 手动计算：图层数据的含义
# 每根柱是离散D取值的精确概率；红点是正态近似在该格点半格边界内的概率。
normal_n20 = exact_table(20, 0.05, 0.25)
normal_n1000 = exact_table(1000, 0.05, 0.25)
rare_n1000 = exact_table(1000, 0.001, 0.002)

# %% 调用绘图函数：Matplotlib概率柱与seaborn近似概率点
sns.set_theme(style="whitegrid", context="notebook")
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["font.sans-serif"] = ["SimHei", "Microsoft YaHei", "Noto Sans CJK SC", "PingFang SC", "STHeiti", "Arial Unicode MS", "DejaVu Sans"]


def draw(ax, data, title, limits):
    # between只选定显示范围；完整精确分布仍保留在data中。
    shown = data[data["D"].between(*limits)].copy()
    cell_width = float(np.diff(data["D"]).min())
    bar_width = 0.82 * cell_width
    # Matplotlib按真实D坐标画概率柱，不对概率再求平均或计算区间。
    ax.bar(shown["D"], shown["exact_probability"], width=bar_width,
           color="#4C78A8", alpha=0.8, label="精确概率")
    # seaborn只叠加已经算好的正态单格概率点，纵轴与柱子都是概率。
    sns.scatterplot(data=shown, x="D", y="normal_cell_probability", ax=ax,
                    s=20, color="#E45756", label="正态单格概率", zorder=3)
    ax.set(xlabel="率差 D=k/n（后测减前测）", ylabel="单个取值的概率", title=title)
    ax.set_xlim(limits)
    largest_probability = max(shown["exact_probability"].max(), shown["normal_cell_probability"].max())
    ax.set_ylim(0, largest_probability * 1.25)
    ax.legend(loc="upper right", fontsize=9)

# %% 结果核对：先看命名数据表，再看图层
fig1, axes = plt.subplots(2, 1, figsize=(6.6, 8.5), dpi=190)
draw(axes[0], normal_n20, "给定变化概率，n=20", (-0.65, 0.25))
draw(axes[1], normal_n1000, "同一变化概率，n=1000（横轴局部放大）", (-0.27, -0.13))
fig1.tight_layout(h_pad=2.4)

fig2, ax2 = plt.subplots(figsize=(6.6, 4.8), dpi=190)
draw(ax2, rare_n1000, "n=1000，平均只有3人发生变化", (-0.009, 0.007))
fig2.tight_layout()

# %% 可选保存：默认只显示图；改SAVE_FIGURES=True才导出到运行结果/Python/图像
if SAVE_FIGURES:
    candidates = [Path.cwd(), Path.cwd().parent]
    if "__file__" in globals():
        candidates.append(Path(__file__).resolve().parents[1])
    for candidate in candidates:
        if (candidate / "数据" / "paired_20.csv").is_file():
            project = candidate
            break
    else:
        raise FileNotFoundError("请从配套项目根目录或代码目录运行")
    figure_dir = project / "运行结果" / "Python" / "图像"
    figure_dir.mkdir(parents=True, exist_ok=True)
    fig1.savefig(figure_dir / "图1_配对率差精确与近似.png")
    fig1.savefig(figure_dir / "图1_配对率差精确与近似.pdf")
    fig2.savefig(figure_dir / "图2_变化稀少.png")
    fig2.savefig(figure_dir / "图2_变化稀少.pdf")

if plt.get_backend().lower() == "agg":
    plt.close(fig1)
    plt.close(fig2)
else:
    plt.show()
