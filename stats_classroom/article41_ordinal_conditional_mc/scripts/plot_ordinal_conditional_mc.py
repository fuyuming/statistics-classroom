# 自编绘图函数：从主脚本传入的枚举概率及本次模拟T值构图。
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

def plot_results(enum_df, observed, mc, output=None):
    threshold = abs(observed["mini"]["T"])
    exact = enum_df.groupby("T", as_index=False)["probability"].sum()
    support = exact["T"].to_numpy()
    simulated = np.asarray(mc["mini"]["T_sim"])
    # 同一T值的多个表先按条件概率归并，模拟则计该T的出现比例。
    frequency = np.array([np.mean(np.abs(simulated - t) < 1e-9) for t in support])
    tail = np.abs(support) >= threshold - 1e-10
    colors = np.where(tail, "#F58518", "#4C78A8")
    fig, axes = plt.subplots(2, 1, figsize=(5.2, 8.0))
    axes[0].bar(support, exact["probability"], width=.78, color=colors)
    axes[0].set(title="Exact conditional T distribution", xlabel="Centered product sum T", ylabel="Exact probability")
    axes[0].axvline(threshold, color="#D62728", ls="--", label=f"Tail: |T| >= {threshold:g}")
    axes[0].axvline(-threshold, color="#D62728", ls="--")
    axes[0].legend(fontsize=12, loc="upper right")
    axes[1].plot(support, exact["probability"], "o-", ms=4, color="#F58518", label="Exact probability")
    axes[1].plot(support, frequency, "x--", ms=5, color="#4C78A8", label="Monte Carlo frequency")
    axes[1].set(title=f"Monte Carlo comparison (B={len(simulated)})", xlabel="T value", ylabel="Probability / frequency")
    axes[1].legend(fontsize=12, loc="upper right")
    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
        # 窄画布配较大字号，缩到390px手机正文后仍可直接读数。
        ax.tick_params(labelsize=12)
        ax.title.set_fontsize(15)
        ax.xaxis.label.set_fontsize(14)
        ax.yaxis.label.set_fontsize(14)
        ax.set_xticks([-10, -5, 0, 5, 10])
        ax.set_ylim(0, max(exact["probability"].max(), frequency.max()) * 1.65)
    fig.tight_layout(pad=1.2, h_pad=2.4)
    if output is not None:
        output = Path(output)
        output.parent.mkdir(exist_ok=True)
        fig.savefig(output, dpi=200, facecolor="white")
    else:
        plt.show()
    return fig
