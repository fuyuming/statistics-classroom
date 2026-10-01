# -*- coding: utf-8 -*-
# ============================================================
# 精读 04｜分类变量、中心化与曲线（Python 版，与 R 版、MATLAB 版逐位一致）
# 对应：McElreath 2023 第 04 讲 Categories and Curves；教材第 4、5 章
# 用法： python 代码/篇04_分类与曲线.py
# 前置：先跑 R 版生成 运行结果/00_教材数据_Howell1*.csv
# ============================================================
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from scipy.stats import norm

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES, OUT = os.path.join(_ROOT, "运行结果"), os.path.join(_ROOT, "文章配图")
os.makedirs(RES, exist_ok=True)
os.makedirs(OUT, exist_ok=True)
for cand in ["PingFang SC", "Hiragino Sans GB", "Heiti SC", "Arial Unicode MS"]:
    if any(f.name == cand for f in font_manager.fontManager.ttflist):
        plt.rcParams["font.family"] = cand
        break
plt.rcParams["axes.unicode_minus"] = False
TEAL, RED, PINK, GREY, DARK = "#167d80", "#c1462c", "#e8b4a8", "#9aa8a8", "#163d48"


def num(x, d=6):
    s = f"{round(float(x), d):.{d}f}".rstrip("0").rstrip(".")
    return "0" if s in ("", "-0") else s


def write_csv(cols, rows, path):
    out = [",".join(f'"{c}"' for c in cols)]
    for r in rows:
        out.append(",".join(num(v) if isinstance(v, (int, float, np.floating, np.integer)) else str(v) for v in r))
    open(path, "w", encoding="utf-8").write("\n".join(out) + "\n")


def ols_ls(X, y):
    """显式最小二乘：解正规方程（与 R 的 solve(t(X)X, t(X)y) 同一套运算）"""
    X = np.asarray(X, float); y = np.asarray(y, float)
    b = np.linalg.solve(X.T @ X, X.T @ y)
    res = y - X @ b
    return b, math.sqrt(float(res @ res) / (len(y) - X.shape[1]))


# ---------- 0) 教材数据 ----------
hw = np.genfromtxt(os.path.join(RES, "00_教材数据_Howell1成年人.csv"), delimiter=",", names=True, encoding="utf-8")
all_hw = np.genfromtxt(os.path.join(RES, "00_教材数据_Howell1全部.csv"), delimiter=",", names=True, encoding="utf-8")

# ---------- 1) 三组终点厚度 ----------
labels = ["对照", "低剂量", "高剂量"]
vals = {"对照": [41.2, 43.0, 39.8, 42.5, 40.9],
        "低剂量": [46.5, 48.1, 45.2, 47.4, 46.0],
        "高剂量": [52.8, 54.6, 51.9, 53.7, 52.2]}
write_csv(["处理", "厚度"], [[g, v] for g in labels for v in vals[g]],
          os.path.join(RES, "01_我们的案例_三组终点厚度_python.csv"))

g_mean = {g: float(np.mean(vals[g])) for g in labels}
g_n = {g: len(vals[g]) for g in labels}
ss_within = sum(float(np.sum((np.array(vals[g]) - g_mean[g]) ** 2)) for g in labels)
df_within = sum(g_n.values()) - len(labels)
sigma_pool = math.sqrt(ss_within / df_within)
post_sd = {g: sigma_pool / math.sqrt(g_n[g]) for g in labels}
z89 = norm.ppf(0.945)

rows2 = [[g, round(g_mean[g], 4), round(post_sd[g], 4),
          round(g_mean[g] - z89 * post_sd[g], 4), round(g_mean[g] + z89 * post_sd[g], 4)] for g in labels]
write_csv(["组", "样本均值", "后验标准差", "区间89下", "区间89上"], rows2,
          os.path.join(RES, "02_组均值与后验区间_python.csv"))

pairs = [("低剂量 - 对照", "低剂量", "对照"), ("高剂量 - 对照", "高剂量", "对照"), ("高剂量 - 低剂量", "高剂量", "低剂量")]
rows3 = []
for name, A, B in pairs:
    d = round(g_mean[A] - g_mean[B], 4)
    se = round(math.sqrt(post_sd[A] ** 2 + post_sd[B] ** 2), 4)
    rows3.append([name, d, se, round(d - z89 * se, 4), round(d + z89 * se, 4), round(float(norm.cdf(d / se)), 4)])
write_csv(["对比", "均值差", "差的标准误", "区间89下", "区间89上", "P大于0"], rows3,
          os.path.join(RES, "03_组间对比_python.csv"))
for r in rows3:
    print("  对比：", r)

# ---------- 2) 饱和曲线 ----------
t_hours = np.arange(1, 13)
resid_fix = np.array([0.5, -1.1, 0.9, -0.6, 1.2, -0.4, 0.7, -1.0, 0.3, 1.1, -0.8, 0.6])
satu = 58 * (1 - np.exp(-0.32 * t_hours)) + resid_fix
write_csv(["时间", "厚度"], [[int(t_hours[i]), round(satu[i], 6)] for i in range(12)],
          os.path.join(RES, "04_饱和曲线数据_python.csv"))

def design(t, kind):
    t = np.asarray(t, float)
    if kind == "lin":
        return np.column_stack([np.ones_like(t), t])
    if kind == "poly":
        return np.column_stack([np.ones_like(t), t, t ** 2, t ** 3])
    return np.column_stack([np.ones_like(t), t, np.maximum(t - 4, 0), np.maximum(t - 8, 0)])

fits = {k: ols_ls(design(t_hours, k), satu) for k in ["lin", "poly", "spl"]}
fit_rows = [[int(t_hours[i]), round(satu[i], 6),
             round(float(design(t_hours, "lin")[i] @ fits["lin"][0]), 4),
             round(float(design(t_hours, "poly")[i] @ fits["poly"][0]), 4),
             round(float(design(t_hours, "spl")[i] @ fits["spl"][0]), 4)] for i in range(12)]
write_csv(["时间", "观测", "直线", "三次多项式", "线性样条"], fit_rows,
          os.path.join(RES, "05_三种拟合_python.csv"))

fit_summary = [["直线", 2, round(fits["lin"][1], 4)], ["三次多项式", 4, round(fits["poly"][1], 4)],
               ["线性样条（结点 4、8）", 4, round(fits["spl"][1], 4)]]
write_csv(["模型", "参数个数", "残差标准差"], fit_summary, os.path.join(RES, "06_拟合优度_python.csv"))
for r in fit_summary:
    print("  拟合：", r)

ext_rows = []
for t in range(13, 17):
    ext_rows.append([t,
                     round(float(design([t], "lin")[0] @ fits["lin"][0]), 4),
                     round(float(design([t], "poly")[0] @ fits["poly"][0]), 4),
                     round(float(design([t], "spl")[0] @ fits["spl"][0]), 4)])
write_csv(["时间", "直线", "三次多项式", "线性样条"], ext_rows, os.path.join(RES, "07_外推四小时_python.csv"))

# 教材对照：年龄→体重（结点 10/15/25/40）
age = np.asarray(all_hw["age"], float); wt = np.asarray(all_hw["weight"], float)
Xl = np.column_stack([np.ones_like(age), age])
Xp = np.column_stack([np.ones_like(age), age, age ** 2, age ** 3])
Xs = np.column_stack([np.ones_like(age), age, np.maximum(age - 10, 0), np.maximum(age - 15, 0),
                      np.maximum(age - 25, 0), np.maximum(age - 40, 0)])
hw_sum = [["直线", 2, round(ols_ls(Xl, wt)[1], 4)], ["三次多项式", 4, round(ols_ls(Xp, wt)[1], 4)],
          ["线性样条（结点 10/15/25/40）", 6, round(ols_ls(Xs, wt)[1], 4)]]
write_csv(["模型", "参数个数", "残差标准差"], hw_sum, os.path.join(RES, "08_教材对照_拟合优度_python.csv"))

# 中心化
tbar = float(t_hours.mean()); n_obs = len(t_hours)
Sxx = float(np.sum((t_hours - tbar) ** 2))
sd_a_raw = fits["lin"][1] * math.sqrt(1 / n_obs + tbar ** 2 / Sxx)
sd_a_ctr = fits["lin"][1] / math.sqrt(n_obs)
ctr = [["未中心化：截距＝时间 0 处的厚度", round(fits["lin"][0][0], 4), round(sd_a_raw, 4)],
       ["中心化：截距＝平均时间处的厚度", round(fits["lin"][0][0] + fits["lin"][0][1] * tbar, 4), round(sd_a_ctr, 4)]]
write_csv(["写法", "截距后验均值", "截距后验标准差"], ctr, os.path.join(RES, "09_中心化对照_python.csv"))

# ---------- 图 ----------
fig, ax = plt.subplots(figsize=(6.6, 4.0), dpi=150)
rng = np.random.default_rng(1)
for i, g in enumerate(labels):
    xs = i + rng.uniform(-0.08, 0.08, size=g_n[g])
    ax.scatter(xs, vals[g], s=26, color=RED, zorder=5)
ax.errorbar(range(len(labels)), [g_mean[g] for g in labels],
            yerr=[z89 * post_sd[g] for g in labels], fmt="o", color=DARK, capsize=5, markersize=8, zorder=6)
ax.set_xticks(range(len(labels))); ax.set_xticklabels(labels)
ax.set_title("三个处理组的终点厚度：先把「组均值」的后验画出来", fontsize=10)
ax.set_ylabel("终点厚度（微米）")
fig.tight_layout(); fig.savefig(os.path.join(OUT, "06-Python版-三组终点厚度.png")); plt.close(fig)

xs = np.linspace(-12, 14, 400)
fig, ax = plt.subplots(figsize=(7.4, 4.2), dpi=150)
for name, d, se, lo, hi, p in rows3:
    ax.plot(xs, norm.pdf(xs, d, se), linewidth=1.4, label=name)
ax.axvline(0, color=GREY, linestyle="--")
ax.set_title("报系数不如报对比：三组两两差多少", fontsize=11)
ax.set_xlabel("厚度差（微米）"); ax.set_ylabel("后验密度"); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "07-Python版-组间对比后验.png")); plt.close(fig)

print("完成：2 张对照图 + 9 份 CSV 已写入")
