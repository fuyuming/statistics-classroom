# -*- coding: utf-8 -*-
# ============================================================
# 精读 06｜好控制与坏控制（Python 版，与 R 版、MATLAB 版逐位一致）
# 对应：McElreath 2023 第 06 讲 Good and Bad Controls；教材第 6 章
# 用法： python 代码/篇06_好控制与坏控制.py
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
TEAL, RED, PINK = "#167d80", "#c1462c", "#e8b4a8"


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
    res = y - X @ b
    return b, math.sqrt(float(res @ res) / (len(y) - X.shape[1])), 1 - float(res @ res) / float((y - y.mean()) @ (y - y.mean()))


def orth_resid(v, X):
    v = np.asarray(v, float); X = np.asarray(X, float)
    return v - X @ np.linalg.solve(X.T @ X, X.T @ v)


def r4(x):
    return round(float(x), 4)


n = 12
treat_x = np.array([1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0], float)
X0 = np.column_stack([np.ones(n), treat_x])
resid_z = orth_resid([0.4, -0.3, 0.2, -0.5, 0.6, -0.2, 0.3, -0.4, 0.5, -0.6, 0.1, 0.4], X0)
z_conf = 0.7 * treat_x + resid_z
resid_y = orth_resid([0.5, -0.4, 0.3, -0.6, 0.2, -0.3, 0.6, -0.5, 0.4, -0.2, 0.1, 0.5],
                     np.column_stack([np.ones(n), treat_x, z_conf]))
y_out = 3.0 * treat_x + 2.0 * z_conf + resid_y
resid_w = orth_resid([0.3, 0.5, -0.4, 0.2, -0.5, 0.4, -0.2, 0.6, -0.3, 0.1, -0.6, 0.4],
                     np.column_stack([np.ones(n), treat_x, y_out]))
w_coll = 1.0 * y_out + 0.5 * treat_x + resid_w
write_csv(["单元", "处理", "混杂", "结果", "对撞"],
          [[i + 1, int(treat_x[i]), round(z_conf[i], 4), round(y_out[i], 4), round(w_coll[i], 4)] for i in range(n)],
          os.path.join(RES, "01_我们的案例_python.csv"))

m_none = ols_ls(np.column_stack([np.ones(n), treat_x]), y_out)
m_good = ols_ls(np.column_stack([np.ones(n), treat_x, z_conf]), y_out)
m_bad = ols_ls(np.column_stack([np.ones(n), treat_x, z_conf, w_coll]), y_out)
rows = [["不控制任何变量", r4(m_none[0][1]), 3.0],
        ["只控制混杂 z（好的控制）", r4(m_good[0][1]), 3.0],
        ["同时控制混杂 z 与对撞 w（坏的控制）", r4(m_bad[0][1]), 3.0]]
write_csv(["写法", "处理效应估计", "真值"], rows, os.path.join(RES, "02_好控制与坏控制_python.csv"))
for r in rows:
    print("  ", r)

# 对撞偏倚
y_zero = orth_resid([0.5, -0.4, 0.3, -0.6, 0.2, -0.3, 0.6, -0.5, 0.4, -0.2, 0.1, 0.5], X0)
resid_w2 = orth_resid([0.2, -0.3, 0.4, -0.2, 0.5, -0.4, 0.3, -0.5, 0.1, 0.4, -0.1, 0.3],
                      np.column_stack([np.ones(n), treat_x, y_zero]))
w_coll2 = 1.2 * treat_x + 1.5 * y_zero + resid_w2
z0 = ols_ls(np.column_stack([np.ones(n), treat_x]), y_zero)
z1 = ols_ls(np.column_stack([np.ones(n), treat_x, w_coll2]), y_zero)
rows3 = [["不控制对撞（结果 ~ 处理）", r4(z0[0][1]), 0.0], ["控制对撞（结果 ~ 处理 + 对撞）", r4(z1[0][1]), 0.0]]
write_csv(["写法", "处理效应估计", "真值"], rows3, os.path.join(RES, "03_对撞偏倚_python.csv"))

# 精度寄生虫
resid_p = orth_resid([0.1, -0.2, 0.3, -0.1, 0.2, -0.3, 0.1, 0.2, -0.2, 0.3, -0.1, 0.2],
                     np.column_stack([np.ones(n), treat_x, y_out]))
parasite = 0.9 * treat_x + resid_p
r2_xz = ols_ls(np.column_stack([np.ones(n), parasite]), treat_x)[2]
se_a = m_none[1] / math.sqrt(float(((treat_x - treat_x.mean()) ** 2).sum()))
se_b = ols_ls(np.column_stack([np.ones(n), treat_x, parasite]), y_out)[1] / \
    math.sqrt(float(((treat_x - treat_x.mean()) ** 2).sum()) * (1 - r2_xz))
rows4 = [["处理与寄生虫的决定系数 R²", r4(r2_xz)], ["不控制时的系数标准误", r4(se_a)],
         ["控制寄生虫后的系数标准误", r4(se_b)], ["标准误放大倍数", r4(se_b / se_a)],
         ["理论值 1/sqrt(1-R²)", r4(1 / math.sqrt(1 - r2_xz))]]
write_csv(["量", "值"], rows4, os.path.join(RES, "04_精度寄生虫_python.csv"))
for r in rows4:
    print("  ", r)

# 偏倚放大（含不可测混杂）
u_hidden = 0.6 * treat_x + orth_resid([0.2, -0.4, 0.5, -0.1, 0.3, -0.5, 0.4, -0.2, 0.1, 0.5, -0.3, 0.2], X0)
resid_amp = orth_resid([0.4, -0.3, 0.2, -0.5, 0.6, -0.1, 0.3, -0.4, 0.5, -0.6, 0.1, 0.3],
                       np.column_stack([np.ones(n), treat_x, z_conf, u_hidden]))
y_amp = 3.0 * treat_x + 2.0 * z_conf + 1.5 * u_hidden + resid_amp
w_amp = 1.0 * y_amp + 0.5 * treat_x + orth_resid([0.3, 0.5, -0.4, 0.2, -0.5, 0.4, -0.2, 0.6, -0.3, 0.1, -0.6, 0.4],
                                                np.column_stack([np.ones(n), treat_x, y_amp]))
a0 = ols_ls(np.column_stack([np.ones(n), treat_x]), y_amp)[0][1]
a1 = ols_ls(np.column_stack([np.ones(n), treat_x, z_conf]), y_amp)[0][1]
a2 = ols_ls(np.column_stack([np.ones(n), treat_x, z_conf, w_amp]), y_amp)[0][1]
rows5 = [[t, r4(v), 3.0, r4(abs(v - 3.0))] for t, v in
         [("不控制任何变量", a0), ("只控制混杂 z", a1), ("控制混杂 z 与后果 w", a2)]]
write_csv(["写法", "处理效应估计", "真值", "偏离真值的量"], rows5, os.path.join(RES, "05_偏倚放大_python.csv"))
for r in rows5:
    print("  ", r)

# 图（两张对照图）
fig, ax = plt.subplots(figsize=(7.4, 4.2), dpi=150)
labels = ["不控制\n任何变量", "只控制混杂 z\n（好的控制）", "控制混杂 z\n与对撞 w（坏控制）"]
vals = [r4(m_none[0][1]), r4(m_good[0][1]), r4(m_bad[0][1])]
ax.bar(labels, vals, color=[PINK, TEAL, RED], width=0.55)
ax.axhline(3.0, linestyle="--", color="#163d48")
for i, v in enumerate(vals):
    ax.text(i, v + 0.08, f"{v:.2f}", ha="center", fontsize=10)
ax.set_title("好控制、不控制、坏控制：同一份数据三种答案", fontsize=11)
ax.set_ylabel("处理效应估计")
fig.tight_layout(); fig.savefig(os.path.join(OUT, "06-Python版-三种答案.png")); plt.close(fig)

fig, ax = plt.subplots(figsize=(7.4, 4.2), dpi=150)
labels5 = ["不控制\n任何变量", "只控制混杂 z", "控制混杂 z\n与后果 w"]
vals5 = [r4(a0), r4(a1), r4(a2)]
ax.bar(labels5, vals5, color=[PINK, TEAL, RED], width=0.55)
ax.axhline(3.0, linestyle="--", color="#163d48")
for i, v in enumerate(vals5):
    ax.text(i, v + 0.08, f"{v:.2f}", ha="center", fontsize=10)
ax.set_title("有不可测混杂时：控制「后果」会把正确的调整毁掉", fontsize=11)
ax.set_ylabel("处理效应估计")
fig.tight_layout(); fig.savefig(os.path.join(OUT, "07-Python版-偏倚放大.png")); plt.close(fig)

print("完成：2 张对照图 + 5 份 CSV 已写入")
