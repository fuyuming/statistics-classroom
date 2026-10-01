# -*- coding: utf-8 -*-
# ============================================================
# 精读 12｜多层模型（Python 版，与 R/MATLAB 逐位一致）
# 对应：McElreath 2023 第 12 讲 Multilevel Models；教材第 13 章
# 用法： python 代码/篇12_多层模型.py
# 模型 y_ij = α_j + b·x_ij + ε，α_j ~ Normal(μ, σ_α)；对 (b, σ_α, σ) 铺三维网格，
# 组截距与 μ 解析积分掉（闭式边际似然）。
# ============================================================
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES, OUT = os.path.join(_ROOT, "运行结果"), os.path.join(_ROOT, "文章配图")
os.makedirs(RES, exist_ok=True)
os.makedirs(OUT, exist_ok=True)
for cand in ["PingFang SC", "Hiragino Sans GB", "Heiti SC", "Arial Unicode MS"]:
    if any(f.name == cand for f in font_manager.fontManager.ttflist):
        plt.rcParams["font.family"] = cand
        break
plt.rcParams["axes.unicode_minus"] = False
TEAL, RED, PINK, GREY = "#167d80", "#c1462c", "#e8b4a8", "#9aa8a8"


def num(x, d=6):
    s = f"{round(float(x), d):.{d}f}".rstrip("0").rstrip(".")
    return "0" if s in ("", "-0") else s


def write_csv(cols, rows, path):
    out = [",".join(f'"{c}"' for c in cols)]
    for r in rows:
        vals = []
        for v in r:
            if v is None:
                vals.append("NA")
            elif isinstance(v, (int, float, np.floating, np.integer)):
                vals.append(num(v))
            else:
                vals.append(str(v))
        out.append(",".join(vals))
    open(path, "w", encoding="utf-8").write("\n".join(out) + "\n")


def r4(x):
    return round(float(x), 4)


def wq(v, w, p):
    o = np.argsort(v)
    v, w = np.asarray(v)[o], np.asarray(w)[o]
    cw = np.cumsum(w) / w.sum()
    k = int(np.searchsorted(cw, p))
    if k == 0:
        return float(v[0])
    return float(v[k - 1] + (v[k] - v[k - 1]) * (p - cw[k - 1]) / (cw[k] - cw[k - 1]))


x = np.tile(np.array([1.0, 2.0, 3.0]), 12)
alpha_true = np.array([4.60, 4.90, 5.15, 5.35, 5.50, 5.60, 5.70, 5.85, 6.05, 6.30, 6.60, 7.00])
resid = np.array([-0.66, 0.11, 0.55,  0.55, -0.66, 0.11, -0.26, 0.78, -0.52,
                  0.39, -0.65, 0.26, -0.66, 0.11, 0.55,  0.55, -0.66, 0.11,
                  -0.26, 0.78, -0.52, 0.39, -0.65, 0.26, -0.66, 0.11, 0.55,
                  0.55, -0.66, 0.11, -0.26, 0.78, -0.52, 0.39, -0.65, 0.26])
gid = np.repeat(np.arange(1, 13), 3)
y = alpha_true[gid - 1] + 1.0 * x + resid
write_csv(["瓶", "剂量", "长势"], [[int(gid[i]), int(x[i]), y[i]] for i in range(36)],
          os.path.join(RES, "01_数据_python.csv"))

n_j, J, N = 3, 12, 36
gmean = np.array([y[gid == j].mean() for j in range(1, 13)])
gxbar, gbar = x.mean(), y.mean()
within_y = y - np.array([gmean[g - 1] for g in gid])
gxbar_per_bottle = np.array([x[gid == j].mean() for j in range(1, 13)])
within_x = x - np.array([gxbar_per_bottle[g - 1] for g in gid])

Xp = np.column_stack([np.ones(N), x])
bp = np.linalg.solve(Xp.T @ Xp, Xp.T @ y)
b_fe = float((within_y * within_x).sum() / (within_x ** 2).sum())
unpooled = gmean - b_fe * gxbar
pooled = np.full(J, unpooled.mean())


def marg_ll(b, sa, s):
    r = gmean - b * gxbar
    v = np.full(J, s ** 2 / n_j + sa ** 2)
    S = float((1 / v).sum())
    mu_hat = float((r / v).sum() / S)
    q_between = float(((r - mu_hat) ** 2 / v).sum())
    w = within_y - b * within_x
    return -0.5 * (math.log(S) + float(np.log(v).sum()) + q_between) - 0.5 * float((w ** 2).sum()) / s ** 2 - (N - J) * math.log(s)


b_seq = np.linspace(0, 2, 60)
sa_seq = np.linspace(0.005, 2, 60)
s_seq = np.linspace(0.05, 2, 60)
B, SA, S_ = np.meshgrid(b_seq, sa_seq, s_seq, indexing="ij")
B, SA, S_ = B.ravel(order="F"), SA.ravel(order="F"), S_.ravel(order="F")
lp = np.array([marg_ll(B[i], SA[i], S_[i]) for i in range(len(B))])
post = np.exp(lp - lp.max())
post = post / post.sum()

shrink_f = lambda b, sa, s: sa ** 2 / (sa ** 2 + s ** 2 / n_j)
mu_draw = gmean.mean() - B * gxbar
alpha_draw = np.array([mu_draw[i] + (gmean - B[i] * gxbar - mu_draw[i]) * shrink_f(B[i], SA[i], S_[i])
                       for i in range(len(B))])
partial = alpha_draw.T @ post
b_post = float((B * post).sum())
mu_post = float((mu_draw * post).sum())
sa_post = float((SA * post).sum())
s_post = float((S_ * post).sum())
shrink_post = float(np.array([shrink_f(B[i], SA[i], S_[i]) for i in range(len(B))]) @ post)

write_csv(["瓶", "瓶内均值", "完全合并", "完全不合并", "部分合并", "收缩量"],
          [[j + 1, r4(unpooled[j]), r4(pooled[j]), r4(unpooled[j]), r4(partial[j]),
            r4(unpooled[j] - partial[j])] for j in range(J)],
          os.path.join(RES, "02_三种算法_python.csv"))

write_csv(["参数", "后验均值", "下界89", "上界89"],
          [["共同斜率 b", r4(b_post), r4(wq(B, post, 0.055)), r4(wq(B, post, 0.945))],
           ["总体水平 μ", r4(mu_post), None, None],
           ["组间标准差 σ_α", r4(sa_post), r4(wq(SA, post, 0.055)), r4(wq(SA, post, 0.945))],
           ["残差标准差 σ", r4(s_post), r4(wq(S_, post, 0.055)), r4(wq(S_, post, 0.945))],
           ["收缩因子（每瓶都是 3 个孔）", r4(shrink_post), None, None]],
          os.path.join(RES, "03_超参数后验_python.csv"))

se_existing = math.sqrt(s_post ** 2 + s_post ** 2 / n_j)
se_new = math.sqrt(s_post ** 2 + sa_post ** 2)
write_csv(["情形", "预测标准差", "区间半宽89", "说明"],
          [["同一个瓶再测一个孔", r4(se_existing), r4(1.598 * se_existing), "只多了观测误差与瓶均值自身的不确定性"],
           ["换一个全新的瓶再测一个孔", r4(se_new), r4(1.598 * se_new), "还要加上瓶与瓶之间的真实差异 σ_α"]],
          os.path.join(RES, "04_两种不确定性_python.csv"))

write_csv(["量", "值"],
          [["完全合并的水平", r4(unpooled.mean())], ["完全不合并的极差", r4(unpooled.max() - unpooled.min())],
           ["部分合并的极差", r4(partial.max() - partial.min())],
           ["收缩量与偏离的相关", r4(float(np.corrcoef(unpooled - unpooled.mean(), unpooled - partial)[0, 1]))],
           ["瓶间真实差异（生成时设定）", r4(alpha_true.std(ddof=1))], ["瓶内噪声（生成时设定）", r4(float(np.sqrt((resid ** 2).mean())))]],
          os.path.join(RES, "05_汇总_python.csv"))

# 两张对照图
fig, ax = plt.subplots(figsize=(7.4, 4.4), dpi=150)
for j in range(1, 13):
    ax.scatter([j] * 3, y[gid == j], color=TEAL, s=25, alpha=0.85, zorder=2)
ax.scatter(np.arange(1, 13), gmean, color=RED, s=55, zorder=3)
ax.axhline(gbar, ls="--", color=GREY)
ax.set_title("12 个培养瓶，每瓶 3 个孔（Python 版对照图）", fontsize=11)
ax.set_xlabel("培养瓶"); ax.set_ylabel("长势")
fig.tight_layout(); fig.savefig(os.path.join(OUT, "06-Python版-数据.png")); plt.close(fig)

fig, ax = plt.subplots(figsize=(7.0, 5.0), dpi=150)
ax.plot([unpooled.min(), unpooled.max()], [unpooled.min(), unpooled.max()], ls="--", color=GREY)
ax.axhline(unpooled.mean(), ls=":", color="#163d48")
ax.scatter(unpooled, partial, color=TEAL, s=45, zorder=3)
for j in range(J):
    ax.annotate("", xy=(unpooled[j], unpooled[j]), xytext=(unpooled[j], partial[j]),
                arrowprops=dict(arrowstyle="-|>", color=RED, lw=1))
ax.set_title("部分合并＝把每个瓶的估计往整体拉一点（Python 版）", fontsize=11)
ax.set_xlabel("完全不合并（瓶内均值）"); ax.set_ylabel("部分合并")
fig.tight_layout(); fig.savefig(os.path.join(OUT, "07-Python版-收缩.png")); plt.close(fig)

print("完成：2 张对照图 + 5 份 CSV 已写入")
