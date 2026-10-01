# -*- coding: utf-8 -*-
# ============================================================
# 精读 08｜MCMC（Python 版，与 R 版、MATLAB 版逐位一致）
# 对应：McElreath 2023 第 08 讲 Markov Chain Monte Carlo；教材第 9 章
#
# 设计：MCMC 天生要用随机数，而本系列要求三个语言给出同一个答案。
# 这里用确定性低差异序列（Halton）代替伪随机数：提议步长由 qnorm(Halton) 给出，
# 接受判断用 Halton 的均匀序列。算法结构与真实 Metropolis 一致，但结果可复现。
# ============================================================
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from scipy.stats import norm as _norm

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


def halton(n, base):
    out = []
    for i in range(1, n + 1):
        f, r, k = 1.0, 0.0, i
        while k > 0:
            f /= base
            r += f * (k % base)
            k //= base
        out.append(r)
    return np.array(out)


def wq(v, w, p):
    o = np.argsort(v)
    v = np.asarray(v)[o]; w = np.asarray(w)[o]
    cw = np.cumsum(w) / w.sum()
    k = int(np.argmax(cw >= p))
    if k == 0:
        return float(v[0])
    return float(v[k - 1] + (v[k] - v[k - 1]) * (p - cw[k - 1]) / (cw[k] - cw[k - 1]))


def orth_resid(v, X):
    v = np.asarray(v, float); X = np.asarray(X, float)
    return v - X @ np.linalg.solve(X.T @ X, X.T @ v)


n = 16
x = np.arange(1, n + 1, dtype=float)
resid_y = orth_resid([0.9, -1.1, 0.7, -0.8, 1.2, -0.4, 0.6, -1.0, 0.5, 1.1, -0.7, 0.8, -0.9, 0.4, -0.6, 1.0],
                     np.column_stack([np.ones(n), x]))
y = 12 + 2.5 * x + resid_y
write_csv(["单元", "时间", "观测"], [[i + 1, int(x[i]), r4(y[i])] for i in range(n)],
          os.path.join(RES, "01_数据_python.csv"))

X0 = np.column_stack([np.ones(n), x])
b_ols = np.linalg.solve(X0.T @ X0, X0.T @ y)
s_y = float(np.sqrt(((y - X0 @ b_ols) ** 2).sum() / (n - 2)))


def log_lik(a, b):
    return -0.5 * float((((y - (a + b * x)) / s_y) ** 2).sum())


# 网格后验（对照）
a_seq = 0 + np.arange(300) * (24 / 299)
b_seq = 0 + np.arange(300) * (5 / 299)
A, B = np.meshgrid(a_seq, b_seq, indexing="ij")
Ag = A.ravel(order="F"); Bg = B.ravel(order="F")   # 与 R 的 expand.grid 同序
lp = np.array([log_lik(Ag[i], Bg[i]) for i in range(len(Ag))])
gw = np.exp(lp - lp.max()); gw = gw / gw.sum()
grid_row = ["网格（300×300 = 9 万格点）", r4((Ag * gw).sum()), r4(wq(Ag, gw, 0.055)), r4(wq(Ag, gw, 0.945)),
            r4((Bg * gw).sum()), r4(wq(Bg, gw, 0.055)), r4(wq(Bg, gw, 0.945))]


def metropolis(start, n_step, step_a, step_b, offset):
    chain = np.zeros((n_step, 2))
    chain[0] = start
    z1 = _norm.ppf(halton(n_step, 2)); z2 = _norm.ppf(halton(n_step, 3)); u = halton(n_step, 5)
    lp_cur = log_lik(chain[0, 0], chain[0, 1])
    n_acc = 0
    for i in range(1, n_step):
        j = ((i + offset - 1) % n_step)
        prop = chain[i - 1] + np.array([step_a * z1[j], step_b * z2[j]])
        lp_prop = log_lik(prop[0], prop[1])
        if np.log(u[j]) < (lp_prop - lp_cur):
            chain[i] = prop; lp_cur = lp_prop; n_acc += 1
        else:
            chain[i] = chain[i - 1]
    return chain, n_acc / (n_step - 1)


n_step, burn = 20000, 2000
ch1, acc1 = metropolis([3, 0.5], n_step, 0.9, 0.9, 0)
ch2, acc2 = metropolis([3, 0.5], n_step, 0.25, 0.025, 0)
ch3, acc3 = metropolis([3, 0.5], n_step, 0.05, 0.005, 0)
ch2b, _ = metropolis([20, 4.0], n_step, 0.25, 0.025, 777)
keep1, keep2, keep3, keep2b = ch1[burn:], ch2[burn:], ch3[burn:], ch2b[burn:]
all_keep = np.vstack([keep2, keep2b])


def acf_of(v, lag):
    v = np.asarray(v, float) - np.mean(v)
    return float((v[lag:] * v[:-lag]).sum() / (v ** 2).sum())


acf_tab = {lb: [r4(acf_of(k[:, 1], l)) for l in range(1, 21)]
           for lb, k in [("步长太大", keep1), ("合适步长", keep2), ("步长太小", keep3)]}

run_mean1 = np.cumsum(keep2[:, 1]) / np.arange(1, len(keep2) + 1)
run_mean2 = np.cumsum(keep2b[:, 1]) / np.arange(1, len(keep2b) + 1)
w_var = float(np.mean([keep2[:, 1].var(ddof=1), keep2b[:, 1].var(ddof=1)]))
b_var = float(np.var([keep2[:, 1].mean(), keep2b[:, 1].mean()], ddof=1))
r_hat = float(np.sqrt((((len(keep2) - 1) / len(keep2)) * w_var + b_var / len(keep2)) / w_var))
neff = float(len(keep2[:, 1]) * (1 - acf_of(keep2[:, 1], 1)) / (1 + acf_of(keep2[:, 1], 1)))


def wq_eq(v, p):
    return wq(v, np.repeat(1.0 / len(v), len(v)), p)


mcmc_row = ["Metropolis（两条链，预热 2000 步，各留 18000 步）", r4(all_keep[:, 0].mean()),
            r4(wq_eq(all_keep[:, 0], 0.055)), r4(wq_eq(all_keep[:, 0], 0.945)),
            r4(all_keep[:, 1].mean()), r4(wq_eq(all_keep[:, 1], 0.055)), r4(wq_eq(all_keep[:, 1], 0.945))]
write_csv(["方法", "截距均值", "截距89下", "截距89上", "斜率均值", "斜率89下", "斜率89上"],
          [grid_row, mcmc_row], os.path.join(RES, "02_网格与MCMC对照_python.csv"))
for r in (grid_row, mcmc_row):
    print("  ", r)

diag_rows = [["步长太大：接受率", r4(acc1)], ["合适步长：接受率", r4(acc2)], ["步长太小：接受率", r4(acc3)],
             ["步长太大：滞后 1 自相关", r4(acf_tab["步长太大"][0])], ["合适步长：滞后 1 自相关", r4(acf_tab["合适步长"][0])],
             ["步长太小：滞后 1 自相关", r4(acf_tab["步长太小"][0])],
             ["合适步长：近似有效样本量（3.6 万步）", r4(neff)], ["简化 R-hat（斜率，两条合适步长的链）", r4(r_hat)],
             ["斜率后验均值（链 1）", r4(keep2[:, 1].mean())], ["斜率后验均值（链 2）", r4(keep2b[:, 1].mean())]]
write_csv(["量", "值"], diag_rows, os.path.join(RES, "03_诊断_python.csv"))
for r in diag_rows:
    print("  ", r)

write_csv(["参数个数", "格点数"], [[k, r4(100.0 ** k)] for k in range(1, 9)],
          os.path.join(RES, "04_维度诅咒_python.csv"))
steps = np.arange(100, len(keep2), 500)
write_csv(["步数", "链1运行均值", "链2运行均值"],
          [[int(steps[i]), r4(run_mean1[steps[i] - 1]), r4(run_mean2[steps[i] - 1])] for i in range(len(steps))],
          os.path.join(RES, "05_运行均值_python.csv"))

# 图（两张对照图）
fig, ax = plt.subplots(figsize=(7.4, 4.6), dpi=150)
ax.scatter(Ag[::37], Bg[::37], c="#eef5f5", s=1)
ax.scatter(all_keep[::20, 0], all_keep[::20, 1], s=1.5, color=RED, alpha=0.5)
ax.set_title("MCMC 做的事：在后验上走一圈", fontsize=11)
ax.set_xlabel("截距 a"); ax.set_ylabel("斜率 b")
fig.tight_layout(); fig.savefig(os.path.join(OUT, "06-Python版-MCMC与网格.png")); plt.close(fig)

fig, ax = plt.subplots(figsize=(7.4, 4.4), dpi=150)
for lb, col in [("步长太大", GREY), ("合适步长", TEAL), ("步长太小", RED)]:
    ax.plot(range(1, 21), acf_tab[lb], "-o", ms=3, color=col, label=lb)
ax.set_title("自相关：步长太大或太小都不好", fontsize=11)
ax.set_xlabel("滞后"); ax.set_ylabel("自相关"); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "07-Python版-自相关.png")); plt.close(fig)

print("完成：2 张对照图 + 5 份 CSV 已写入")
