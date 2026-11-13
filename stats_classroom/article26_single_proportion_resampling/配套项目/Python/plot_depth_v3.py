"""统计课堂26 v3：Bootstrap区间构造与精确覆盖率图。"""
from pathlib import Path
import os

import matplotlib.pyplot as plt
from matplotlib import font_manager
import pandas as pd
import seaborn as sns


METHOD_COLORS = {
    "Percentile": "#2345AB",
    "Basic": "#2A9D8F",
    "Wald": "#CE6C22",
    "Wilson": "#7A5195",
    "CP": "#222222",
}


def _setup_font():
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
    plt.rcParams["font.size"] = 15


def make_figures(project_root=None):
    if project_root is None:
        project_root = Path(__file__).resolve().parent.parent
    project_root = Path(project_root)
    output = project_root.parent / "文章配图"
    output.mkdir(exist_ok=True)
    sns.set_theme(style="whitegrid", context="talk")
    _setup_font()

    intervals = pd.read_csv(project_root / "results" / "Python_intervals_v3.csv")
    coverage = pd.read_csv(project_root / "results" / "Python_coverage_v3.csv")

    # 图3：同一样本的五种区间，Percentile/Basic端点来自条件分布直接计算。
    fig, ax = plt.subplots(figsize=(8.2, 4.8), dpi=170)
    y_positions = {method: len(intervals) - 1 - i for i, method in enumerate(intervals["method"])}
    for _, row in intervals.iterrows():
        method = row["method"]
        ypos = y_positions[method]
        color = METHOD_COLORS[method]
        ax.hlines(ypos, row["lower"], row["upper"], color=color, lw=5, alpha=.9)
        ax.plot([row["lower"], row["upper"]], [ypos, ypos], "o", color=color, ms=7)
        ax.text(row["upper"] + .004, ypos, f"[{row['lower']:.3f}, {row['upper']:.3f}]", va="center", fontsize=12)
    ax.axvline(.65, color="#2345AB", ls=":", lw=2)
    ax.text(.65, 4.42, "观察率 0.65", ha="center", color="#2345AB", fontsize=14)
    ax.set_yticks([4, 3, 2, 1, 0], labels=["Percentile", "Basic", "Wald", "Wilson", "CP"])
    ax.set_xlim(.55, .82)
    ax.set_ylim(-.7, 4.7)
    ax.set_xlabel("达标率 p")
    ax.set_title("同一样本的五种95%区间", fontsize=22, loc="left", pad=14)
    ax.text(.55, -0.53, "Percentile/Basic：条件分布直接计算的无限Bootstrap端点", fontsize=11, color="#555555")
    sns.despine(ax=ax, left=True)
    fig.tight_layout()
    fig.savefig(output / "03_区间构造_v3.png")
    plt.close(fig)

    # 图4：每个p点由二项分布枚举x=0,...,180得到精确覆盖率。
    fig, axes = plt.subplots(2, 1, figsize=(7.2, 8.2), dpi=170, sharex=True, sharey=True)
    panels = [(axes[0], ["Percentile", "Basic"]), (axes[1], ["Wald", "Wilson", "CP"])]
    for ax, methods in panels:
        for method in methods:
            sub = coverage[coverage["method"] == method]
            ax.plot(sub["p"], sub["coverage"], marker="o", ms=2.5, lw=1.6,
                    label=method, color=METHOD_COLORS[method])
        ax.axhline(.95, color="#777777", ls="--", lw=1.4, label="95%目标")
        ax.set_xlim(.01, .99)
        ax.set_ylim(.80, 1.005)
        ax.set_xticks([.1, .3, .5, .7, .9])
        ax.legend(frameon=False, fontsize=16, loc="lower center", ncol=2)
        sns.despine(ax=ax)
    axes[0].set_ylabel("覆盖率")
    axes[1].set_ylabel("覆盖率")
    axes[1].set_xlabel("真实总体率 p")
    fig.suptitle("95%区间的精确覆盖率", fontsize=22, x=.08, ha="left", y=1.03)
    fig.text(.08, -.005, "n=180；Bootstrap内层无限重抽；覆盖率由二项枚举（p=.01,.02,…,.99）", fontsize=11, color="#555555")
    fig.tight_layout(w_pad=1.4)
    fig.savefig(output / "04_覆盖率_v3.png", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    make_figures()
