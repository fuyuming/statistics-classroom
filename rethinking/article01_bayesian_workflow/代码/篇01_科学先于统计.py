# -*- coding: utf-8 -*-
"""
《Statistical Rethinking》精读 第 01 篇 —— Python 配套脚本
用 numpy/scipy/matplotlib 把 R 版脚本走一遍，数字必须与 R、MATLAB 完全一致。
运行： python3 篇01_科学先于统计.py     （需 numpy scipy matplotlib）
输出： 运行结果/01_网格近似的后验汇总_python.csv、运行结果/02_两种过程的近似程度_python.csv
       文章配图/06-Python版-先验似然后验与收窄.png
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import binom, poisson, nbinom

# 中文字体：macOS 用 PingFang SC；Windows 换 Microsoft YaHei；Linux 换 Noto Sans CJK SC
plt.rcParams["font.sans-serif"] = ["PingFang SC", "Heiti SC", "Noto Sans CJK SC"]
plt.rcParams["axes.unicode_minus"] = False

OUT, RES = "文章配图", "运行结果"
os.makedirs(OUT, exist_ok=True)
os.makedirs(RES, exist_ok=True)


def grid_approx(W, N, n_grid=20):
    """网格近似：先验 × 似然 → 归一化后验（与 R 的 grid_approx 一一对应）"""
    p = np.linspace(0, 1, n_grid)
    prior = np.ones(n_grid)                       # 平坦先验
    like = binom.pmf(W, N, p)
    post = like * prior
    post = post / post.sum()
    return p, prior / prior.sum(), like / like.sum(), post


cases = [(2, 3), (6, 9), (20, 30)]
rows, curves = [], []
for W, N in cases:
    p, prior, like, post = grid_approx(W, N)
    cum = np.cumsum(post)
    lo = p[np.argmax(cum >= 0.055)]
    hi = p[np.argmax(cum >= 0.945)]
    rows.append([N, W, round(W / N, 4), round(float((p * post).sum()), 6),
                 round(float(p[post.argmax()]), 6), round(float(lo), 4),
                 round(float(hi), 4), round(float(post.max()), 6)])
    curves.append((f"{N} 次里 {W} 次是水", p, prior, like, post))

header = "样本量N,水上次数W,观测比例,网格后验均值,网格最大值点,区间下89,区间上89,后验密度峰值"
with open(os.path.join(RES, "01_网格近似的后验汇总_python.csv"), "w", encoding="utf-8") as f:
    f.write(header + "\n")
    for r in rows:
        f.write(",".join(str(x) for x in r) + "\n")
print(header)
for r in rows:
    print(r)

# 两种过程有多像：直接比理论概率质量函数的全变差距离（确定性，可跨语言复核）
k = np.arange(0, 41)
pm_pois = poisson.pmf(k, 3)
pm_nb = nbinom.pmf(k, n=20, p=20 / (20 + 3))
tv = 0.5 * np.abs(pm_pois - pm_nb).sum()
maxdiff = np.abs(pm_pois - pm_nb).max()
with open(os.path.join(RES, "02_两种过程的近似程度_python.csv"), "w", encoding="utf-8") as f:
    f.write("指标,数值\n")
    f.write(f"全变差距离,{tv:.6f}\n")
    f.write(f"最大单点概率差,{maxdiff:.6f}\n")
print(f"泊松 vs 负二项(size=20) 全变差距离 = {tv:.6f}，最大单点概率差 = {maxdiff:.6f}")

fig, ax = plt.subplots(1, 2, figsize=(11, 4.2), dpi=240)
_, p, prior, like, post = curves[1]
ax[0].plot(p, prior, color="#8a8a8a", lw=1.6, marker="o", ms=3, label="先验（还没看数据）")
ax[0].plot(p, like, color="#c1462c", lw=1.4, ls="--", label="似然（数据说什么）")
ax[0].plot(p, post, color="#167d80", lw=1.6, marker="o", ms=3, label="后验（看完之后）")
ax[0].set(title="把 20 个候选比例各自算一遍，就得到后验",
          xlabel="水面比例 p", ylabel="相对密度（各自归一化）")
ax[0].legend(loc="upper center", ncol=3, fontsize=8, frameon=False,
             bbox_to_anchor=(0.5, 1.22))
for name, pp, _, _, po in curves:
    ax[1].plot(pp, po, lw=1.8, label=name)
ax[1].set(title="观测比例都是 2/3，样本越多，后验越尖",
          xlabel="水面比例 p", ylabel="后验密度")
ax[1].legend(loc="upper center", ncol=3, fontsize=8, frameon=False,
             bbox_to_anchor=(0.5, 1.22))
for a in ax:
    a.grid(alpha=0.25, lw=0.6)
    for s in ("top", "right"):
        a.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "06-Python版-先验似然后验与收窄.png"), bbox_inches="tight",
            facecolor="white")
print("Python 图已保存；numpy", np.__version__)
