"""统计课堂26：只用 boot::boot 实跑 CSV 重绘 Bootstrap 配图 v2。"""
from pathlib import Path
import os

import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
import pandas as pd
import seaborn as sns


def make_figure(project_root=None):
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

    boot = pd.read_csv(project_root / "data" / "boot_v2_reference.csv")
    if len(boot) != 10000 or not {"stat", "count"}.issubset(boot.columns):
        raise ValueError("boot_v2_reference.csv 必须包含10000轮 stat/count 数据")
    counts = boot["count"].astype(int)
    k = np.arange(90, 146)
    frequency = counts.value_counts().reindex(k, fill_value=0) / len(boot)
    low, high = np.quantile(boot["stat"], [.025, .975]) * 100

    blue = "#2345AB"
    orange = "#CE6C22"
    fig, ax = plt.subplots(figsize=(7, 4.8), dpi=170)
    ax.bar(k / 180 * 100, frequency, width=.46, color="#6B97DE", edgecolor="white", linewidth=.3)
    ax.axvline(low, color=orange, ls="--", lw=2)
    ax.axvline(high, color=orange, ls="--", lw=2)
    ax.axvline(65, color=blue, ls=":", lw=2)
    ax.text(low - 1, .077, f"{low:.2f}%", ha="right", color=orange, fontsize=17)
    ax.text(high + 1, .077, f"{high:.2f}%", ha="left", color=orange, fontsize=17)
    ax.text(65, .086, "观察率65%", ha="center", color=blue, fontsize=17)
    ax.set(xlim=(48, 82), ylim=(0, .095), xlabel="每轮达标率（%）", ylabel="10000轮中的比例")
    ax.set_title("从当前样本重抽", fontsize=23, loc="left", pad=15)
    ax.tick_params(labelsize=16)
    sns.despine(ax=ax)
    fig.tight_layout()
    fig.savefig(output / "01_Bootstrap区间_v2.png")
    plt.close(fig)


if __name__ == "__main__":
    make_figure()
