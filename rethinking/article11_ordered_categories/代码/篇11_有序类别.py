# -*- coding: utf-8 -*-
# ============================================================
# 精读 11｜有序类别：切点自由估计（Python 版，与 R/MATLAB 逐位一致）
# 对应：McElreath 2023 第 11 讲 Ordered Categories；教材第 11 章
# 用法： python 代码/篇11_有序类别.py
#
# 参数化：c1 = u，c2 = u + exp(v)，保证 c1 < c2（有序约束由参数化承担）
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


x = np.array([1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8,9,9,9,10,10,10], float)
y = np.array([0,0,1, 0,0,1, 0,1,1, 0,1,1, 1,1,1, 1,1,2, 1,1,2, 1,2,2, 2,2,2, 2,2,2], int)
write_csv(["剂量", "等级"], [[int(x[i]), int(y[i])] for i in range(30)],
          os.path.join(RES, "01_数据_python.csv"))


def loglik(b, u, v):
    c1 = u
    c2 = u + math.exp(v)
    F1 = 1 / (1 + np.exp(-(c1 - b * x)))
    F2 = 1 / (1 + np.exp(-(c2 - b * x)))
    P = np.column_stack([F1, F2 - F1, 1 - F2])
    return float(np.log(np.maximum(P[np.arange(len(y)), y], 1e-300)).sum())


b_seq = np.linspace(-2, 3, 70)
u_seq = np.linspace(-4, 4, 70)
v_seq = np.linspace(math.log(0.05), math.log(5), 70)
B, U, V = np.meshgrid(b_seq, u_seq, v_seq, indexing="ij")
B, U, V = B.ravel(order="F"), U.ravel(order="F"), V.ravel(order="F")   # 与 R 的 expand.grid 同序：b 变得最快
lp = np.array([loglik(B[i], U[i], V[i]) for i in range(len(B))])
post = np.exp(lp - lp.max())
post = post / post.sum()
c1_arr, c2_arr = U, U + np.exp(V)

rows2 = [["斜率 b", r4((B * post).sum()), r4(wq(B, post, 0.055)), r4(wq(B, post, 0.945))],
         ["切点 c1（轻/中）", r4((c1_arr * post).sum()), r4(wq(c1_arr, post, 0.055)), r4(wq(c1_arr, post, 0.945))],
         ["切点 c2（中/重）", r4((c2_arr * post).sum()), r4(wq(c2_arr, post, 0.055)), r4(wq(c2_arr, post, 0.945))]]
write_csv(["参数", "后验均值", "下界89", "上界89"], rows2, os.path.join(RES, "02_后验_python.csv"))

rows3 = []
for dd in range(1, 11):
    F1 = 1 / (1 + np.exp(-(c1_arr - B * dd)))
    F2 = 1 / (1 + np.exp(-(c2_arr - B * dd)))
    pk, pz, pzh = (F1 * post).sum(), ((F2 - F1) * post).sum(), ((1 - F2) * post).sum()
    rpk, rpz, rpzh = r4(pk), r4(pz), r4(pzh)
    rows3.append([dd, r4((F1 * post).sum()), r4((F2 * post).sum()), rpk, rpz, rpzh, r4(rpk + rpz + rpzh)])
write_csv(["剂量", "累积_轻", "累积_中", "概率_轻", "概率_中", "概率_重", "概率合计"], rows3,
          os.path.join(RES, "03_累积概率与等级概率_python.csv"))

obs_prop = [float((y == k).sum()) / len(y) for k in range(3)]
pred_by_dose = []
for dd in range(1, 11):
    F1 = 1 / (1 + np.exp(-(c1_arr - B * dd)))
    F2 = 1 / (1 + np.exp(-(c2_arr - B * dd)))
    pred_by_dose.append([(F1 * post).sum(), ((F2 - F1) * post).sum(), ((1 - F2) * post).sum()])
pred_prop = np.array(pred_by_dose).mean(axis=0)
write_csv(["等级", "观测比例", "预测比例", "差"],
          [["轻", r4(obs_prop[0]), r4(pred_prop[0]), r4(pred_prop[0] - obs_prop[0])],
           ["中", r4(obs_prop[1]), r4(pred_prop[1]), r4(pred_prop[1] - obs_prop[1])],
           ["重", r4(obs_prop[2]), r4(pred_prop[2]), r4(pred_prop[2] - obs_prop[2])]],
          os.path.join(RES, "04_后验预测检查_python.csv"))

# 合并成二分类：等级 ≥ 2 记 1
y_bin = (y >= 2).astype(float)


def loglik_bin(a, b):
    p = 1 / (1 + np.exp(-(a + b * x)))
    p = np.clip(p, 1e-12, 1 - 1e-12)
    return float((y_bin * np.log(p) + (1 - y_bin) * np.log(1 - p)).sum())


a_seq = np.linspace(-8, 3, 300)
b2_seq = np.linspace(-1.5, 2.5, 300)
Aa, Bb = np.meshgrid(a_seq, b2_seq, indexing="ij")
Aa, Bb = Aa.ravel(order="F"), Bb.ravel(order="F")
lpa = np.array([loglik_bin(Aa[i], Bb[i]) for i in range(len(Aa))])
wa = np.exp(lpa - lpa.max())
wa = wa / wa.sum()
se_bin = float(np.sqrt(((Bb - (Bb * wa).sum()) ** 2 * wa).sum()))
se_ord = float(np.sqrt(((B - (B * post).sum()) ** 2 * post).sum()))
write_csv(["写法", "斜率后验均值", "斜率后验标准差", "说明"],
          [["有序模型（三档全用）", r4((B * post).sum()), r4(se_ord), "斜率是「累积几率的共同位移」"],
           ["合并成二分类（重 vs 其余）", r4((Bb * wa).sum()), r4(se_bin),
            "斜率是「重 vs 其余」的几率变化——不是同一个量"]],
          os.path.join(RES, "05_有序vs二分类_python.csv"))

rows6 = []
for dd in (3, 5, 7):
    F1 = 1 / (1 + np.exp(-(c1_arr - B * dd)))
    F2 = 1 / (1 + np.exp(-(c2_arr - B * dd)))
    p_bin = 1 / (1 + np.exp(-((Aa * wa).sum() + (Bb * wa).sum() * dd)))
    rows6.append([dd, r4((F1 * post).sum()), r4(((F2 - F1) * post).sum()), r4(((1 - F2) * post).sum()),
                  r4(p_bin), r4(1 - p_bin)])
write_csv(["剂量", "有序_轻", "有序_中", "有序_重", "二分类_重", "二分类_其余"], rows6,
          os.path.join(RES, "06_信息损失的形态_python.csv"))

# 图（两张对照图）
fig, ax = plt.subplots(figsize=(7.4, 4.4), dpi=150)
doses = np.arange(1, 11)
p_light = np.array([r[3] for r in rows3]); p_mid = np.array([r[4] for r in rows3]); p_heavy = np.array([r[5] for r in rows3])
ax.bar(doses, p_heavy, color=RED, label="重")
ax.bar(doses, p_mid, bottom=p_heavy, color=PINK, label="中")
ax.bar(doses, p_light, bottom=p_heavy + p_mid, color=TEAL, label="轻")
ax.set_title("三档等级的概率怎么随剂量变化（Python 版对照图）", fontsize=11)
ax.set_xlabel("剂量"); ax.set_ylabel("概率"); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "06-Python版-等级概率.png")); plt.close(fig)

fig, ax = plt.subplots(figsize=(7.4, 4.4), dpi=150)
labels = ["斜率 b", "切点 c1", "切点 c2"]
means = [(B * post).sum(), (c1_arr * post).sum(), (c2_arr * post).sum()]
los = [wq(B, post, 0.055), wq(c1_arr, post, 0.055), wq(c2_arr, post, 0.055)]
his = [wq(B, post, 0.945), wq(c1_arr, post, 0.945), wq(c2_arr, post, 0.945)]
ax.errorbar(means, labels, xerr=[np.array(means) - los, np.array(his) - np.array(means)], fmt="o", color=TEAL, ecolor="#163d48", capsize=5)
ax.axvline(0, ls="--", color=GREY)
ax.set_title("三个参数一起估：斜率与两个切点（Python 版对照图）", fontsize=11)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "07-Python版-参数.png")); plt.close(fig)

print("完成：2 张对照图 + 6 份 CSV 已写入")
