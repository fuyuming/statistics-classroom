# -*- coding: utf-8 -*-
# ============================================================
# 精读 09｜事件数据：0/1 结果、计数、率与有序类别（Python 版，与 R 版、MATLAB 版逐位一致）
# 对应：McElreath 2023 第 09 讲 Modeling Events；教材第 10、11 章
# 用法： python 代码/篇09_事件数据.py
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
        out.append(",".join(num(v) if isinstance(v, (int, float, np.floating, np.integer)) else str(v) for v in r))
    open(path, "w", encoding="utf-8").write("\n".join(out) + "\n")


def r4(x):
    return round(float(x), 4)


def wq(v, w, p):
    o = np.argsort(v)
    v = np.asarray(v, float)[o]; w = np.asarray(w, float)[o]
    cw = np.cumsum(w) / w.sum()
    k = int(np.argmax(cw >= p))
    if k == 0:
        return float(v[0])
    return float(v[k - 1] + (v[k] - v[k - 1]) * (p - cw[k - 1]) / (cw[k] - cw[k - 1]))


def grid2(a_seq, b_seq, loglik):
    """两参数网格后验：与 R 的 expand.grid 同序（a 变得最快）。"""
    A, B = np.meshgrid(a_seq, b_seq, indexing="ij")
    A = A.ravel(order="F"); B = B.ravel(order="F")
    lp = np.array([loglik(A[i], B[i]) for i in range(len(A))])
    w = np.exp(lp - lp.max())
    return A, B, w / w.sum()


# ---------- 1) 0/1 结果：定植成功与否 ----------
dose = np.arange(1, 9, dtype=float)
n_units = np.repeat(10.0, 8)
n_succ = np.array([1, 2, 4, 5, 7, 8, 9, 10], float)
write_csv(["剂量", "单元数", "成功数", "成功比例"],
          [[int(dose[i]), int(n_units[i]), int(n_succ[i]), r4(n_succ[i] / n_units[i])] for i in range(8)],
          os.path.join(RES, "01_定植数据_python.csv"))


def loglik_binom(a, b):
    p = 1 / (1 + np.exp(-(a + b * dose)))
    p = np.clip(p, 1e-12, 1 - 1e-12)
    return float((n_succ * np.log(p) + (n_units - n_succ) * np.log(1 - p)).sum())


a_seq = np.linspace(-6, 2, 300); b_seq = np.linspace(0, 1.2, 300)
A_b, B_b, W_b = grid2(a_seq, b_seq, loglik_binom)
write_csv(["参数", "后验均值", "下界89", "上界89"],
          [["截距 a", r4((A_b * W_b).sum()), r4(wq(A_b, W_b, 0.055)), r4(wq(A_b, W_b, 0.945))],
           ["斜率 b", r4((B_b * W_b).sum()), r4(wq(B_b, W_b, 0.055)), r4(wq(B_b, W_b, 0.945))]],
          os.path.join(RES, "02_二项后验_python.csv"))

dose_fine = np.linspace(0, 9, 60)
rows3 = []
for d in dose_fine:
    p = 1 / (1 + np.exp(-(A_b + B_b * d)))
    rows3.append([round(d, 3), round(float((p * W_b).sum()), 3), round(wq(p, W_b, 0.055), 3), round(wq(p, W_b, 0.945), 3)])
write_csv(["剂量", "均值", "下", "上"], rows3, os.path.join(RES, "03_定植概率曲线_python.csv"))
p5 = (1 / (1 + np.exp(-(A_b + B_b * 5))))
print("  剂量 5 的成功概率后验均值 =", r4((p5 * W_b).sum()))

# ---------- 2) 计数：菌落数（泊松）----------
# 每个剂量 3 个重复：过离散要在"同一剂量内"算
treat = np.repeat(np.arange(0, 8), 3).astype(float)
counts = np.array([1, 2, 3, 1, 3, 5, 2, 5, 8, 2, 6, 10, 4, 9, 14, 4, 11, 18, 5, 14, 23, 6, 18, 30], float)
write_csv(["处理水平", "菌落数"], [[int(treat[i]), int(counts[i])] for i in range(len(treat))],
          os.path.join(RES, "04_计数数据_python.csv"))
mean_cnt = float(counts.mean())
within = np.array([[float(counts[treat == d].var(ddof=1)), float(counts[treat == d].mean())] for d in range(8)])
var_cnt = float(within[:, 0].sum() / within[:, 1].sum())


def loglik_pois(a, b):
    lam = np.exp(a + b * treat)
    return float((counts * np.log(lam) - lam - np.vectorize(math.lgamma)(counts + 1)).sum())


a_p = np.linspace(-1, 3, 300); b_p = np.linspace(0, 0.6, 300)
A_p, B_p, W_p = grid2(a_p, b_p, loglik_pois)
write_csv(["参数", "后验均值", "下界89", "上界89"],
          [["截距 a", r4((A_p * W_p).sum()), r4(wq(A_p, W_p, 0.055)), r4(wq(A_p, W_p, 0.945))],
           ["斜率 b", r4((B_p * W_p).sum()), r4(wq(B_p, W_p, 0.055)), r4(wq(B_p, W_p, 0.945))]],
          os.path.join(RES, "05_泊松后验_python.csv"))
rr = np.exp(B_p * 7)
write_csv(["量", "值"],
          [["每皿平均菌落数（观测）", r4(mean_cnt)], ["方差/均值（同一剂量内合并，过离散检查）", r4(var_cnt)],
           ["率比 exp(7b) 后验均值", r4((rr * W_p).sum())],
           ["率比 89% 区间下", r4(wq(rr, W_p, 0.055))], ["率比 89% 区间上", r4(wq(rr, W_p, 0.945))]],
          os.path.join(RES, "06_率比与过离散_python.csv"))

# ---------- 3) 率：暴露量不同（offset）----------
plate = np.array([10.0, 10, 10, 10]); obs_small = np.array([3.0, 5, 4, 4])
plate_big = np.array([40.0, 40, 40, 40]); obs_big = np.array([12.0, 18, 15, 19])
naive = obs_big.mean() / obs_small.mean()
offset_correct = (obs_big.sum() / plate_big.sum()) / (obs_small.sum() / plate.sum())
write_csv(["比较", "比值", "说明"],
          [["直接比计数（忽略体积）", r4(naive), "大体积皿当然菌落更多，这是体积的功劳"],
           ["按体积折算成率再比", r4(offset_correct), "把体积放进模型（对数偏移），比的是每微升的密度"]],
          os.path.join(RES, "07_offset对照_python.csv"))

# ---------- 4) 有序类别：轻/中/重（切点固定，只估斜率）----------
ord_x = np.array([1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6], float)
ord_y = np.array([0, 0, 0, 1, 0, 1, 1, 1, 1, 2, 2, 2], float)
write_csv(["剂量", "等级"], [[int(ord_x[i]), int(ord_y[i])] for i in range(12)],
          os.path.join(RES, "08_有序数据_python.csv"))
c1_fix, c2_fix = 1.5, 0.0


def loglik_ord1(b):
    p_le0 = np.clip(1 / (1 + np.exp(-(c1_fix - b * ord_x))), 1e-9, 1 - 1e-9)
    p_le1 = np.clip(1 / (1 + np.exp(-(c2_fix - b * ord_x))), 1e-9, 1 - 1e-9)
    ll = 0.0
    for i in range(len(ord_y)):
        if ord_y[i] == 0:
            ll += math.log(p_le0[i])
        elif ord_y[i] == 1:
            ll += math.log(max(p_le1[i] - p_le0[i], 1e-9))
        else:
            ll += math.log(1 - p_le1[i])
    return ll


b_grid = np.linspace(-0.5, 3, 400)
lp_ord = np.array([loglik_ord1(b) for b in b_grid])
w_ord = np.exp(lp_ord - lp_ord.max()); w_ord = w_ord / w_ord.sum()
write_csv(["参数", "后验均值", "下界89", "上界89"],
          [["斜率 b", r4((b_grid * w_ord).sum()), r4(wq(b_grid, w_ord, 0.055)), r4(wq(b_grid, w_ord, 0.945))]],
          os.path.join(RES, "09_有序后验_python.csv"))
print("  有序斜率 b =", r4((b_grid * w_ord).sum()))
cont_fit = np.polyfit(ord_x, ord_y, 1)
rows10 = []
for d in range(1, 7):
    pk = float((1 / (1 + np.exp(-(c1_fix - b_grid * d))) * w_ord).sum())
    pz = float((1 / (1 + np.exp(-(c2_fix - b_grid * d))) * w_ord).sum())
    rows10.append([d, r4(pk), r4(pz), r4(cont_fit[0] * d + cont_fit[1])])
write_csv(["剂量", "轻的累积概率", "中的累积概率", "当连续变量_预测等级"], rows10,
          os.path.join(RES, "10_连续vs有序_python.csv"))

# ---------- 图（两张对照图）----------
fig, ax = plt.subplots(figsize=(7.4, 4.4), dpi=150)
ax.scatter(dose, n_succ / n_units, s=40, color=RED, zorder=3)
mean_curve = [r[1] for r in rows3]
ax.plot(dose_fine, mean_curve, color=TEAL)
ax.fill_between(dose_fine, [r[2] for r in rows3], [r[3] for r in rows3], color=TEAL, alpha=0.18)
ax.set_title("0/1 结果要的是概率曲线，不是直线", fontsize=11)
ax.set_xlabel("剂量"); ax.set_ylabel("成功概率"); ax.set_ylim(0, 1)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "06-Python版-logit.png")); plt.close(fig)

fig, ax = plt.subplots(figsize=(7.4, 4.4), dpi=150)
ax.bar([str(d) for d in range(1, 7)], [r[3] for r in rows10], color=PINK, width=0.5, label="把等级当连续变量")
ax.plot(range(6), [r[1] for r in rows10], "-o", color=TEAL, label="累积 logit：P(等级 ≤ 轻)")
ax.plot(range(6), [r[2] for r in rows10], "-o", color=RED, label="累积 logit：P(等级 ≤ 中)")
ax.set_title("两种读法给出不同的东西", fontsize=11)
ax.set_xlabel("剂量"); ax.set_ylabel("值（概率或折算后的等级）"); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "07-Python版-连续vs有序.png")); plt.close(fig)

print("完成：2 张对照图 + 10 份 CSV 已写入")
