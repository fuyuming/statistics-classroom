# -*- coding: utf-8 -*-
# ============================================================
# 精读 03｜12 个时间点，就敢画一条生长曲线吗（Python 版，与 R 版、MATLAB 版逐位一致）
# 对应：McElreath 2023 第 03 讲 Geocentric Models；教材第 4 章
# 用法： python 代码/篇03_一条直线.py
# 前置：先跑 R 版生成 运行结果/00_教材数据_Howell1成年人.csv
# ============================================================
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from scipy.stats import binom, norm

# 目录自动定位：脚本在 <篇目录>/代码/ 下
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
A_LIM, B_LIM, N_GRID = (10.0, 26.0), (2.0, 6.0), 200


def num(x, d=6):
    s = f"{round(float(x), d):.{d}f}"
    s = s.rstrip("0").rstrip(".")
    return "0" if s in ("", "-0") else s


def write_csv(cols, rows, path):
    out = [",".join(f'"{c}"' for c in cols)]
    for r in rows:
        out.append(",".join(num(v) if isinstance(v, (int, float, np.floating, np.integer)) else str(v) for v in r))
    open(path, "w", encoding="utf-8").write("\n".join(out) + "\n")


def halton(n, base):
    out = np.zeros(n)
    for i in range(1, n + 1):
        f, r, k = 1.0, 0.0, i
        while k > 0:
            f /= base; r += f * (k % base); k //= base
        out[i - 1] = r
    return out


def grid_post(x, y, n_grid=N_GRID, sigma=None, a_lim=A_LIM, b_lim=B_LIM, prior_a=None, prior_b=None):
    step_a = (a_lim[1] - a_lim[0]) / (n_grid - 1)
    step_b = (b_lim[1] - b_lim[0]) / (n_grid - 1)
    a_seq = a_lim[0] + np.arange(n_grid) * step_a
    b_seq = b_lim[0] + np.arange(n_grid) * step_b
    x = np.asarray(x, float); y = np.asarray(y, float)
    M = np.outer(x, b_seq)
    lp = np.zeros(n_grid * n_grid)
    for i in range(n_grid * n_grid):
        a = a_seq[i % n_grid]                     # 与 R 的 expand.grid 同序：a 变得最快
        b = b_seq[i // n_grid]
        e = (y - (a + M[:, i // n_grid])) / sigma
        lp[i] = -0.5 * float(np.sum(e * e)) - len(y) * (math.log(sigma) + 0.5 * math.log(2 * math.pi))
        if prior_a is not None:
            lp[i] += -0.5 * ((a - prior_a[0]) / prior_a[1]) ** 2 - math.log(prior_a[1]) - 0.5 * math.log(2 * math.pi)
        if prior_b is not None:
            lp[i] += -0.5 * ((b - prior_b[0]) / prior_b[1]) ** 2 - math.log(prior_b[1]) - 0.5 * math.log(2 * math.pi)
    post = np.exp(lp - lp.max())
    return a_seq, b_seq, post / post.sum()


def ols(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    b = float(np.sum((x - x.mean()) * (y - y.mean())) / np.sum((x - x.mean()) ** 2))
    a = float(y.mean() - b * x.mean())
    res = y - (a + b * x)
    s2 = float(np.sum(res ** 2) / (len(y) - 2))
    return a, b, math.sqrt(s2), math.sqrt(s2 / float(np.sum((x - x.mean()) ** 2)))


def marginal(a_seq, b_seq, post, which):
    n = len(a_seq)
    idx = (np.arange(len(post)) % n) if which == "a" else (np.arange(len(post)) // n)
    p = np.zeros(n)
    for i, j in enumerate(idx):
        p[j] += post[i]
    return (a_seq if which == "a" else b_seq), p / p.sum()


def interval89(seq_, p):
    c_ = np.cumsum(p)
    return seq_[int(np.argmax(c_ >= 0.055))], seq_[int(np.argmax(c_ >= 0.945))]


def joint_draw(a_seq, b_seq, post, n, base=5):
    u = halton(n, base); c_ = np.cumsum(post)
    rows = [int(np.argmax(c_ >= uu)) for uu in u]
    A = np.array([a_seq[r % len(a_seq)] for r in rows])
    B = np.array([b_seq[r // len(a_seq)] for r in rows])
    return A, B


def wquantile(v, w, p):
    o = np.argsort(v, kind="stable"); v = np.asarray(v)[o]; cw = np.cumsum(np.asarray(w)[o]) / np.sum(w)
    k = int(np.argmax(cw >= p))
    if k == 0:
        return float(v[0])
    return float(v[k - 1] + (v[k] - v[k - 1]) * (p - cw[k - 1]) / (cw[k] - cw[k - 1]))


# ---------- 0) 教材对照数据 ----------
p_hw = os.path.join(RES, "00_教材数据_Howell1成年人.csv")
if not os.path.exists(p_hw):
    raise SystemExit("缺 运行结果/00_教材数据_Howell1成年人.csv —— 请先跑 R 版")
raw = np.genfromtxt(p_hw, delimiter=",", names=True, encoding="utf-8")
height, height_c, weight = raw["height"], raw["height_c"], raw["weight"]

# ---------- 1) 主线案例 ----------
t_hours = np.arange(1, 13)
resid_fix = np.array([0.6, -0.9, 1.2, -0.4, 0.8, -1.3, 0.3, 1.1, -0.7, 0.2, -0.5, 1.4])
thickness = 18 + 3.4 * t_hours + resid_fix
write_csv(["时间", "厚度"], list(zip(t_hours, thickness)), os.path.join(RES, "01_我们的案例_生物膜厚度_python.csv"))

# ---------- 2) 网格后验 ----------
ob_a, ob_b, ob_s, ob_se = ols(t_hours, thickness)
oh_a, oh_b, oh_s, oh_se = ols(height_c, weight)
a_bio, b_bio, post_bio = grid_post(t_hours, thickness, N_GRID, ob_s)
a_hw, b_hw, post_hw = grid_post(height_c, weight, N_GRID, oh_s, a_lim=(30.0, 60.0), b_lim=(0.4, 0.9))

mb_a, mb_p = marginal(a_bio, b_bio, post_bio, "b")
mh_a, mh_p = marginal(a_hw, b_hw, post_hw, "b")
iv_bio, iv_hw = interval89(mb_a, mb_p), interval89(mh_a, mh_p)
mean_a_bio = float(np.sum(a_bio[np.arange(len(post_bio)) % N_GRID] * post_bio))
mean_b_bio = float(np.sum(b_bio[np.arange(len(post_bio)) // N_GRID] * post_bio))
mean_a_hw = float(np.sum(a_hw[np.arange(len(post_hw)) % N_GRID] * post_hw))
mean_b_hw = float(np.sum(b_hw[np.arange(len(post_hw)) // N_GRID] * post_hw))


def num4(x):
    s = f"{round(float(x), 4):.4f}"
    return s.rstrip("0").rstrip(".") or "0"


write_csv(["案例", "最小二乘斜率", "后验均值斜率", "后验均值截距", "斜率89下", "斜率89上", "残差标准差"],
          [["我们的案例：厚度~时间", num4(ob_b), num4(mean_b_bio), num4(mean_a_bio), num4(iv_bio[0]), num4(iv_bio[1]), num4(ob_s)],
           ["教材对照：体重~身高", num4(oh_b), num4(mean_b_hw), num4(mean_a_hw), num4(iv_hw[0]), num4(iv_hw[1]), num4(oh_s)]],
          os.path.join(RES, "02_后验摘要_python.csv"))

# ---------- 3) 二项概率 vs 正态的整点区间概率 ----------
n_bin = 10; ks = np.arange(n_bin + 1); sd_bin = math.sqrt(n_bin * 0.25)
pk = binom.pmf(ks, n_bin, 0.5)
pn = norm.cdf(ks + 0.5, n_bin / 2, sd_bin) - norm.cdf(ks - 0.5, n_bin / 2, sd_bin)
write_csv(["成功数", "概率"], list(zip(ks, pk)), os.path.join(RES, "03_二项概率_python.csv"))
write_csv(["成功数", "区间概率"], list(zip(ks, pn)), os.path.join(RES, "03_正态的整点区间概率_python.csv"))

# ---------- 4) 先验预测（从分析所用先验抽线）----------
n_prior = 60
pa = A_LIM[0] + (A_LIM[1] - A_LIM[0]) * halton(n_prior, 2)
pb = B_LIM[0] + (B_LIM[1] - B_LIM[0]) * halton(n_prior, 3)
write_csv(["线号", "截距", "斜率"],
          [[i + 1, round(pa[i], 6), round(pb[i], 6)] for i in range(n_prior)],
          os.path.join(RES, "04_先验预测的直线_python.csv"))

# ---------- 5) 后验预测 ----------
n_draw = 2000
da, db = joint_draw(a_bio, b_bio, post_bio, n_draw, 5)
eps = ob_s * norm.ppf(halton(n_draw, 11))
w = np.full(n_draw, 1.0 / n_draw)
rows5 = []
for tt in t_hours:
    mu = da + db * tt
    rows5.append([tt, round(float(mu.mean()), 3),
                  round(wquantile(mu, w, 0.055), 3), round(wquantile(mu, w, 0.945), 3),
                  round(wquantile(mu + eps, w, 0.055), 3), round(wquantile(mu + eps, w, 0.945), 3)])
write_csv(["时间", "预测均值", "均值89下", "均值89上", "单次89下", "单次89上"], rows5,
          os.path.join(RES, "05_后验预测_python.csv"))

# ---------- 6) 先验敏感性 ----------
sens = [("平坦先验", None), ("斜率先验 N(3.4, 0.2)", (3.4, 0.2)), ("斜率先验 N(5.0, 1.2)", (5.0, 1.2))]
rows6 = []
for name, pv in sens:
    A_, B_, P_ = grid_post(t_hours, thickness, N_GRID, ob_s, prior_b=pv)
    seq_, p_ = marginal(A_, B_, P_, "b")
    lo_, hi_ = interval89(seq_, p_)
    mean_b_ = float(np.sum(B_[np.arange(len(P_)) // N_GRID] * P_))
    wt = 0.0 if pv is None else ob_se ** 2 / (ob_se ** 2 + pv[1] ** 2)
    rows6.append([name, num4(mean_b_), num4(lo_), num4(hi_), num4(wt)])
write_csv(["斜率先验", "斜率均值", "斜率89下", "斜率89上", "先验权重"], rows6,
          os.path.join(RES, "06_先验敏感性_python.csv"))

# ---------- 7) 网格分辨率 ----------
rows7 = []
for ng in [20, 50, 100, 200]:
    A_, B_, P_ = grid_post(t_hours, thickness, ng, ob_s)
    mean_b_ = float(np.sum(B_[np.arange(len(P_)) // ng] * P_))
    mode_b = float(B_[int(np.argmax(P_)) // ng])
    rows7.append([ng, num(mean_b_, 4), num(mode_b, 4), num(ob_b, 4)])
write_csv(["网格点数", "后验均值斜率", "后验密度最高点斜率", "最小二乘参照"], rows7,
          os.path.join(RES, "07_网格分辨率_python.csv"))

# ---------- 图（与 R 版同内容）----------
fig, ax = plt.subplots(figsize=(7.4, 4.2), dpi=150)
xv = np.array([0.5, 12.5])
for i in range(n_prior):
    ax.plot(xv, pa[i] + pb[i] * xv, color=TEAL, alpha=0.20, linewidth=0.4)
ax.scatter(t_hours, thickness, s=26, color=RED, zorder=5)
ax.plot(xv, mean_a_bio + mean_b_bio * xv, color=DARK, linewidth=1.4)
ax.set_title("先验给出一个「直线框」，数据把它收成一条窄带", fontsize=11)
ax.set_xlabel("培养时间（小时）"); ax.set_ylabel("生物膜厚度（微米）")
fig.tight_layout(); fig.savefig(os.path.join(OUT, "06-Python版-先验与后验直线.png")); plt.close(fig)

fig, ax = plt.subplots(figsize=(7.4, 4.2), dpi=150)
ax.bar(ks, pk, color=TEAL, width=0.55, alpha=0.85)
ax.plot(ks, pn, color=RED, linewidth=1.1, marker="o", markersize=4)
ax.set_title("十个小波动加起来，就近似钟形了", fontsize=11)
ax.set_xlabel("成功次数（0 到 10）"); ax.set_ylabel("概率")
fig.tight_layout(); fig.savefig(os.path.join(OUT, "07-Python版-小波动相加.png")); plt.close(fig)

print("完成：2 张对照图 + 8 份 CSV 已写入")
