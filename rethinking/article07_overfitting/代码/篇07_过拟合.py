# -*- coding: utf-8 -*-
# ============================================================
# 精读 07｜过拟合与正则化（Python 版，与 R 版、MATLAB 版逐位一致）
# 对应：McElreath 2023 第 07 讲 Overfitting；教材第 7、8 章
# 用法： python 代码/篇07_过拟合.py
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
    e = y - X @ b
    yy = y - y.mean()
    return b, e, math.sqrt(float(e @ e) / (len(y) - X.shape[1])), float(e @ e), 1 - float(e @ e) / float(yy @ yy)


def loo_error(X, y):
    X = np.asarray(X, float); y = np.asarray(y, float)
    XtX_inv = np.linalg.inv(X.T @ X)
    h = np.einsum("ij,ij->i", X @ XtX_inv, X)
    b, e, _, _, _ = ols_ls(X, y)
    loo_i = e / (1 - h)
    return math.sqrt(float((loo_i ** 2).mean()))


def orth_resid(v, X):
    v = np.asarray(v, float); X = np.asarray(X, float)
    return v - X @ np.linalg.solve(X.T @ X, X.T @ v)


def r4(x):
    return round(float(x), 4)


n = 12
t_x = np.arange(1, n + 1, dtype=float)
t_s = (t_x - t_x.mean()) / t_x.std(ddof=1)
ts_of = lambda tv: (np.asarray(tv, float) - t_x.mean()) / t_x.std(ddof=1)
mu_true = 20 + 4 * t_x - 0.35 * t_x ** 2 + 0.010 * t_x ** 3
resid_g = orth_resid([0.8, -1.2, 0.6, -0.9, 1.1, -0.5, 0.7, -1.0, 0.4, 0.9, -0.6, 1.0],
                     np.column_stack([np.ones(n), t_x, t_x ** 2, t_x ** 3]))
y_obs = mu_true + resid_g
write_csv(["时间", "观测", "生成均值"],
          [[int(t_x[i]), r4(y_obs[i]), r4(mu_true[i])] for i in range(n)],
          os.path.join(RES, "01_数据_python.csv"))


def poly_X(k, tv=None):
    base = t_s if tv is None else ts_of(tv)
    out = np.column_stack([np.ones(len(base)), base])
    for j in range(2, k + 1):
        out = np.column_stack([out, base ** j])
    return out


rows = []
for k in range(1, 8):
    X = poly_X(k)
    b, e, s_hat, sse, _ = ols_ls(X, y_obs)
    rows.append([k, k + 1, r4(s_hat), r4(sse), r4(loo_error(X, y_obs))])
write_csv(["阶数", "系数个数", "样本内残差标准差", "样本内平方误差和", "留一交叉验证误差"], rows,
          os.path.join(RES, "02_阶数与误差_python.csv"))
for r in rows:
    print("  ", r)
best_k = min(rows, key=lambda r: r[4])[0]

t_fine = np.linspace(1, 16, 200)
truth_fine = 20 + 4 * t_fine - 0.35 * t_fine ** 2 + 0.010 * t_fine ** 3


def fit_curve(k, y=None):
    yy = y_obs if y is None else np.asarray(y, float)
    b = ols_ls(poly_X(k), yy)[0]
    return poly_X(k, t_fine) @ b


write_csv(["时间", "真值形状", "三阶", "七阶"],
          [[r4(t_fine[i]), r4(truth_fine[i]), r4(fit_curve(3)[i]), r4(fit_curve(7)[i])] for i in range(200)],
          os.path.join(RES, "03_两条拟合曲线_python.csv"))

# 正则化＝先验宽度
X1 = np.column_stack([np.ones(n), t_x])
s1 = ols_ls(X1, y_obs)[2]
a_seq = 15 + np.arange(200) * (30 / 199)
b_seq = 0 + np.arange(200) * (8 / 199)


def grid_post(prior_b):
    A, B = np.meshgrid(a_seq, b_seq, indexing="ij")
    A = A.ravel(order="F"); B = B.ravel(order="F")   # 与 R 的 expand.grid 同序：a 变得最快
    lp = np.zeros_like(A)
    for i in range(len(A)):
        mu = A[i] + B[i] * t_x
        lp[i] = -0.5 * float((((y_obs - mu) / s1) ** 2).sum()) - n * (math.log(s1) + 0.5 * math.log(2 * math.pi))
        if prior_b is not None:
            lp[i] -= 0.5 * ((B[i] - prior_b[0]) / prior_b[1]) ** 2
    w = np.exp(lp - lp.max()); w = w / w.sum()
    return A, B, w


prior_list = [(3.0, 0.2), (3.0, 1.0), None]
prior_names = ["窄先验 N(3, 0.2)", "中等先验 N(3, 1.0)", "平坦先验"]
reg_rows = []
for nm, pb in zip(prior_names, prior_list):
    A, B, w = grid_post(pb)
    bm = float((B * w).sum()); bs = math.sqrt(float(((B - bm) ** 2 * w).sum()))
    reg_rows.append([nm, r4(bm), r4(bs), r4(float(((A + B * 16) * w).sum()))])
write_csv(["先验", "后验均值斜率", "后验标准差", "预测第16期"], reg_rows,
          os.path.join(RES, "04_先验宽度_python.csv"))
for r in reg_rows:
    print("  ", r)

# 厚尾 vs 正态
y_out2 = y_obs.copy(); y_out2[8] = y_out2[8] + 22
X3 = poly_X(3)
f_norm = ols_ls(X3, y_out2)
z = f_norm[1] / f_norm[2]
w_t = (4 + 1) / (4 + z ** 2)
f_t = ols_ls(X3 * np.sqrt(w_t)[:, None], y_out2 * np.sqrt(w_t))
write_csv(["观测", "标准化残差", "厚尾权重"],
          [[i + 1, r4(z[i]), r4(w_t[i])] for i in range(n)],
          os.path.join(RES, "05_厚尾权重_python.csv"))
write_csv(["模型", "三次项系数", "残差标准差"],
          [["正态最小二乘", r4(f_norm[0][3]), r4(f_norm[2])],
           ["厚尾（t, df=4，IRLS 一步）", r4(f_t[0][3]), r4(f_t[2])]],
          os.path.join(RES, "06_正态与厚尾_python.csv"))
print("  正态三次项", r4(f_norm[0][3]), "厚尾三次项", r4(f_t[0][3]))

# 图
fig, ax = plt.subplots(figsize=(7.6, 4.4), dpi=150)
ks = [r[0] for r in rows]
ax.plot(ks, [r[3] for r in rows], "-o", color=TEAL, label="样本内平方误差和")
ax.plot(ks, [r[4] ** 2 * 12 for r in rows], "-o", color=RED, label="留一交叉验证误差（放大到同一尺度）")
ax.axvline(best_k, linestyle="--", color="#9aa8a8")
ax.set_title("样本内误差一路变小，样本外误差却先降后升", fontsize=11)
ax.set_xlabel("多项式阶数"); ax.set_ylabel("误差"); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "06-Python版-阶数与误差.png")); plt.close(fig)

fig, ax = plt.subplots(figsize=(7.6, 4.4), dpi=150)
ax.scatter(t_x, y_obs, s=40, color=RED, zorder=3)
for k, col, lb in [(3, TEAL, "3 阶"), (7, "#c1462c", "7 阶")]:
    ax.plot(t_fine, fit_curve(k), color=col, label=lb)
ax.plot(t_fine, truth_fine, color="#9aa8a8", linestyle="--", label="真值形状")
ax.axvline(12, linestyle="--", color="#9aa8a8")
ax.set_title("3 阶与 7 阶：样本内都贴得住，出了数据范围就分道扬镳", fontsize=11)
ax.set_xlabel("时间"); ax.set_ylabel("观测值"); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "07-Python版-两阶对照.png")); plt.close(fig)

print("完成：2 张对照图 + 6 份 CSV 已写入")
