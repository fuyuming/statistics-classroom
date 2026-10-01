# -*- coding: utf-8 -*-
# ============================================================
# 精读 05｜混杂的四种基本形状（Python 版，与 R 版、MATLAB 版逐位一致）
# 对应：McElreath 2023 第 05 讲 Elemental Confounds；教材第 5、6 章
# 用法： python 代码/篇05_四种形状.py
# 前置：先跑 R 版生成 运行结果/00_教材数据_华夫饼与离婚率.csv
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


def ols_ls(X, y):
    X = np.asarray(X, float); y = np.asarray(y, float)
    b = np.linalg.solve(X.T @ X, X.T @ y)
    return b, math.sqrt(float((y - X @ b) @ (y - X @ b)) / (len(y) - X.shape[1]))


def pearson(x, y):
    x = np.asarray(x, float) - np.mean(x); y = np.asarray(y, float) - np.mean(y)
    return float(x @ y) / math.sqrt(float(x @ x) * float(y @ y))


def r4(x):
    return round(float(x), 4)


# ---------- 0) 教材对照：读 R 导出的华夫饼数据 ----------
wd = np.genfromtxt(os.path.join(RES, "00_教材数据_华夫饼与离婚率.csv"), delimiter=",", names=True,
                   encoding="utf-8", dtype=None)
state = [s.decode() if isinstance(s, bytes) else str(s) for s in wd["州"]]
div = np.asarray(wd["离婚率"], float); waff = np.asarray(wd["华夫饼店"], float); south = np.asarray(wd["南方"], float)
r_all = r4(pearson(waff, div)); r_s = r4(pearson(waff[south == 1], div[south == 1])); r_n = r4(pearson(waff[south == 0], div[south == 0]))
write_csv(["范围", "华夫饼店与离婚率的相关", "州数"],
          [["全体 50 个州", r_all, len(wd)], ["南方各州", r_s, int(south.sum())], ["非南方各州", r_n, int((south == 0).sum())]],
          os.path.join(RES, "01_教材对照_分层相关_python.csv"))
print("  分层相关：全体", r_all, "南方", r_s, "非南方", r_n)

# ---------- 2) 叉（混杂）----------
treat_c = [1, 1, 1, 1, 0, 1, 0, 0, 0, 0]
batch_c = [0, 0, 0, 0, 0, 1, 1, 1, 1, 1]
resid_c = [0.6, -0.8, 1.1, -0.5, 0.3, 0.9, -1.2, 0.4, -0.6, 0.8]
growth_c = [10 + 3.0 * t + 2.5 * b + e for t, b, e in zip(treat_c, batch_c, resid_c)]
write_csv(["单元", "处理", "批次", "长势"],
          [[i + 1, treat_c[i], "批次2" if batch_c[i] == 1 else "批次1", growth_c[i]] for i in range(10)],
          os.path.join(RES, "02_我们的案例_叉_python.csv"))
fn, _ = ols_ls(np.column_stack([np.ones(10), treat_c]), growth_c)
fa, _ = ols_ls(np.column_stack([np.ones(10), treat_c, batch_c]), growth_c)
write_csv(["写法", "处理效应估计", "真值"],
          [["不调整（长势 ~ 处理）", r4(fn[1]), 3.0], ["把批次放进模型（长势 ~ 处理 + 批次）", r4(fa[1]), 3.0]],
          os.path.join(RES, "03_叉_调整前后_python.csv"))
print("  叉：不调整", r4(fn[1]), "调整", r4(fa[1]), "｜处理与批次相关", r4(pearson(treat_c, batch_c)))

# ---------- 3) 管（中介）----------
resid_d = [0.4, -0.6, 0.9, -0.3, 0.5, -0.8, 0.7, -0.2, 0.6, -0.9]
resid_g = [0.7, 0.5, -0.8, 0.9, -0.4, 0.3, -0.6, 0.8, -0.5, 0.2]
density_d = [5 + 2.0 * t + e for t, e in zip(treat_c, resid_d)]
growth_d = [8 + 1.5 * d + e for d, e in zip(density_d, resid_g)]
write_csv(["单元", "处理", "菌群密度", "长势"],
          [[i + 1, treat_c[i], round(density_d[i], 6), round(growth_d[i], 6)] for i in range(10)],
          os.path.join(RES, "04_我们的案例_管_python.csv"))
pt, _ = ols_ls(np.column_stack([np.ones(10), treat_c]), growth_d)
pd_, _ = ols_ls(np.column_stack([np.ones(10), treat_c, density_d]), growth_d)
write_csv(["写法", "处理效应估计", "说明"],
          [["不控制中介（长势 ~ 处理）", r4(pt[1]), "总效应"],
           ["控制中介（长势 ~ 处理 + 菌群密度）", r4(pd_[1]), "直接效应（总效应的一部分被中介带走）"]],
          os.path.join(RES, "05_管_总效应与直接效应_python.csv"))
print("  管：总效应", r4(pt[1]), "直接效应", r4(pd_[1]))

# ---------- 4) 对撞（选择）----------
base_b = [-1.8, -1.2, -0.7, -0.3, 0.1, 0.4, 0.9, 1.3, 1.7, 2.0]
treat_b = [0, 1, 0, 1, 0, 1, 0, 1, 0, 1]
resid_s = [0.5, -0.4, 0.2, 0.8, -0.9, 0.3, -0.2, 0.6, -0.5, 0.1]
score_s = [0.6 * t + 1.2 * b + e for t, b, e in zip(treat_b, base_b, resid_s)]
med = float(np.median(score_s))
selected = [1 if s > med else 0 for s in score_s]
write_csv(["单元", "处理", "基线值", "选中", "倾向分"],
          [[i + 1, treat_b[i], base_b[i], selected[i], round(score_s[i], 6)] for i in range(10)],
          os.path.join(RES, "06_我们的案例_对撞_python.csv"))
rc_all = r4(pearson(treat_b, base_b))
idx = [i for i in range(10) if selected[i] == 1]
rc_sel = r4(pearson([treat_b[i] for i in idx], [base_b[i] for i in idx]))
write_csv(["范围", "处理与基线值的相关", "单元数"],
          [["全部 10 个单元", rc_all, 10], ["只看被选中的 5 个单元", rc_sel, len(idx)]],
          os.path.join(RES, "07_对撞_条件化前后_python.csv"))
print("  对撞：全体", rc_all, "子集", rc_sel)

# ---------- 5) 后代（治疗后变量）----------
resid_p = [0.3, -0.5, 0.8, -0.2, 0.6, -0.7, 0.4, -0.1, 0.5, -0.8]
resid_y = [0.6, 0.4, -0.7, 0.8, -0.3, 0.2, -0.5, 0.7, -0.4, 0.1]
post_d = [3 + 1.6 * t + e for t, e in zip(treat_c, resid_p)]
growth_p = [9 + 2.8 * t + 0.9 * p + e for t, p, e in zip(treat_c, post_d, resid_y)]
write_csv(["单元", "处理", "处理后指标", "长势"],
          [[i + 1, treat_c[i], round(post_d[i], 6), round(growth_p[i], 6)] for i in range(10)],
          os.path.join(RES, "08_我们的案例_后代_python.csv"))
dt, _ = ols_ls(np.column_stack([np.ones(10), treat_c]), growth_p)
dc, _ = ols_ls(np.column_stack([np.ones(10), treat_c, post_d]), growth_p)
truth = 2.8 + 0.9 * 1.6
write_csv(["写法", "处理效应估计", "真值"],
          [["不控制后代（长势 ~ 处理）", r4(dt[1]), r4(truth)],
           ["把后代当协变量（长势 ~ 处理 + 处理后指标）", r4(dc[1]), r4(truth)]],
          os.path.join(RES, "09_后代_控制前后_python.csv"))
print("  后代：不控制", r4(dt[1]), "控制", r4(dc[1]), "｜真值", r4(truth))

# ---------- 图（两张对照图）----------
fig, ax = plt.subplots(figsize=(7.0, 4.2), dpi=150)
vals = [r4(fa[1]), r4(fn[1])]
ax.bar(["把批次放进模型\n（长势 ~ 处理 + 批次）", "不调整\n（长势 ~ 处理）"], vals, color=[TEAL, PINK], width=0.55)
ax.axhline(3.0, linestyle="--", color=RED)
for i, v in enumerate(vals):
    ax.text(i, v + 0.08, f"{v:.2f}", ha="center", fontsize=10)
ax.set_title("叉（混杂）：不调整批次会把效应估偏", fontsize=11)
ax.set_ylabel("处理效应估计")
fig.tight_layout(); fig.savefig(os.path.join(OUT, "06-Python版-叉.png")); plt.close(fig)

fig, ax = plt.subplots(figsize=(7.4, 4.4), dpi=150)
xs = np.array(base_b)
for grp, col, mk in [("未被选中", TEAL, "o"), ("被选中", RED, "o")]:
    sel = [i for i in range(10) if (selected[i] == 1) == (grp == "被选中")]
    ax.scatter([base_b[i] for i in sel], [treat_b[i] for i in sel], color=col, s=60, label=grp)
xs_line = np.array([-1.9, 2.1])
sel_idx = idx
k = np.polyfit([base_b[i] for i in sel_idx], [treat_b[i] for i in sel_idx], 1)
ax.plot(xs_line, np.polyval(k, xs_line), color=RED, label="被选中子集趋势")
k2 = np.polyfit(base_b, treat_b, 1)
ax.plot(xs_line, np.polyval(k2, xs_line), color=GREY, linestyle="--", label="全体趋势")
ax.set_title("对撞（选择）：只看被选中的样本，就出现了假的负相关", fontsize=11)
ax.set_xlabel("基线值"); ax.set_ylabel("处理（1=接种，0=对照）")
ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "07-Python版-对撞.png")); plt.close(fig)

print("完成：2 张对照图 + 10 份 CSV 已写入")
