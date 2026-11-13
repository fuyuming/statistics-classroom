"""统计课堂26科学配图；由R包实跑分布及显式二项模型生成。
Spyder可运行本文件；只重绘图形，不重新生成模拟样本。
"""
from pathlib import Path
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import font_manager
from scipy.stats import binom, beta

def make_figures(project_root=None):
    if project_root is None:
        project_root = Path(__file__).resolve().parent.parent
    project_root = Path(project_root)
    article_root = project_root.parent
    output = article_root / "文章配图"
    output.mkdir(exist_ok=True)
    sns.set_theme(style="whitegrid", context="talk")
    candidates = [
        os.environ.get("CJK_FONT"),
        "/System/Library/Fonts/PingFang.ttc",
        str(Path.home() / "Library/Fonts/MSYH.TTC"),
        "/System/Library/Fonts/Hiragino Sans GB.ttc",
        "C:/Windows/Fonts/msyh.ttc",
        "C:/Windows/Fonts/simhei.ttf",
        "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc",
    ]
    for candidate in candidates:
        if candidate and Path(candidate).exists():
            font_manager.fontManager.addfont(candidate)
            plt.rcParams["font.family"] = font_manager.FontProperties(fname=candidate).get_name()
            break
    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams["font.size"] = 17
    blue = "#2345AB"
    orange = "#CE6C22"
    boot = pd.read_csv(project_root / "data" / "infer_bootstrap_reference.csv")
    counts = np.rint(boot["stat"] * 180).astype(int)
    k = np.arange(90, 146)
    frequency = pd.Series(counts).value_counts().reindex(k, fill_value=0) / len(boot)
    low, high = np.quantile(boot["stat"], [.025, .975]) * 100
    fig, ax = plt.subplots(figsize=(7, 4.8), dpi=170)
    ax.bar(k / 180 * 100, frequency, width=.46, color="#6B97DE", edgecolor="white", linewidth=.3)
    ax.axvline(low, color=orange, ls="--", lw=2)
    ax.axvline(high, color=orange, ls="--", lw=2)
    ax.axvline(65, color=blue, ls=":", lw=2)
    ax.text(low-1, .077, f"{low:.2f}%", ha="right", color=orange, fontsize=17)
    ax.text(high+1, .077, f"{high:.2f}%", ha="left", color=orange, fontsize=17)
    ax.text(65, .086, "观察率65%", ha="center", color=blue, fontsize=17)
    ax.set(xlim=(48, 82), ylim=(0, .095), xlabel="每轮达标率（%）", ylabel="10000轮中的比例")
    ax.set_title("从当前样本重抽", fontsize=23, loc="left", pad=15)
    ax.tick_params(labelsize=16)
    sns.despine(ax=ax)
    fig.tight_layout()
    fig.savefig(output / "01_Bootstrap区间.png")
    plt.close(fig)

    k = np.arange(0, 181)
    model = pd.DataFrame({"count": k, "rate": k/180, "bootstrap_probability": binom.pmf(k, 180, .65), "null_probability": binom.pmf(k, 180, .45)})
    model.to_csv(project_root / "data" / "绘图模型概率.csv", index=False)
    fig, axes = plt.subplots(2, 1, figsize=(7, 7), dpi=170, sharex=True, sharey=True)
    titles = ["从样本重抽：达标概率0.65", "在H0下模拟：达标概率0.45"]
    for ax, p, color, title in zip(axes, [.65, .45], [blue, orange], titles):
        ax.bar(k/180*100, binom.pmf(k, 180, p), width=.46, color=color, alpha=.75)
        ax.axvline(65, color=blue, ls=":", lw=2)
        ax.set(xlim=(25, 85), ylim=(0, .08), ylabel="每个计数的概率")
        ax.set_title(title, fontsize=20, loc="left", pad=10)
        ax.tick_params(labelsize=16)
        sns.despine(ax=ax)
    axes[1].set_xlabel("计数除以180得到的达标率（%）")
    axes[1].text(65.8, .062, "观察到65%", color=blue, fontsize=16)
    fig.tight_layout(h_pad=1.2)
    fig.savefig(output / "02_两个参考分布.png")
    plt.close(fig)

    upper = beta.ppf(.975, 1, 20) * 100
    fig, ax = plt.subplots(figsize=(7, 4.4), dpi=170)
    ax.text(0, 1.9, "Bootstrap百分位区间", fontsize=19, color=blue)
    ax.scatter([0], [1.5], s=100, color=blue, zorder=4)
    ax.text(1.2, 1.5, "[0, 0]", fontsize=18, va="center")
    ax.text(0, .9, "Clopper–Pearson区间", fontsize=19, color=orange)
    ax.hlines(.5, 0, upper, color=orange, linewidth=5)
    ax.vlines([0, upper], .39, .61, color=orange, linewidth=2)
    ax.text(upper, .72, f"{upper:.2f}%", ha="right", fontsize=18, color=orange)
    ax.set(xlim=(-.8, 20), ylim=(.1, 2.2), yticks=[], xticks=[0, 5, 10, 15, 20], xlabel="总体事件率（%）")
    ax.set_title("20人中0人观察到事件", fontsize=23, loc="left", pad=15)
    ax.grid(axis="y", visible=False)
    ax.tick_params(labelsize=17)
    sns.despine(ax=ax, left=True)
    fig.tight_layout()
    fig.savefig(output / "03_全零样本.png")
    plt.close(fig)

if __name__ == "__main__":
    make_figures()
