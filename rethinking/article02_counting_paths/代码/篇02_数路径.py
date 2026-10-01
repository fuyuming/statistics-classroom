# -*- coding: utf-8 -*-
# ============================================================
# 精读 02｜数路径（Python 版，与 R 版、MATLAB 版逐位一致）
# 对应：McElreath 2023 第 02 讲 Garden of Forking Data；教材第 2 章
# 用法： python 代码/篇02_数路径.py
# ============================================================
import math
import os
from math import lgamma

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import beta as beta_dist

# 中文字体（macOS；Windows 换 Microsoft YaHei，Linux 换 Noto Sans CJK SC）
plt.rcParams["font.sans-serif"] = ["PingFang SC", "Microsoft YaHei", "Noto Sans CJK SC"]
plt.rcParams["axes.unicode_minus"] = False

OUT, RES = "文章配图", "运行结果"
os.makedirs(OUT, exist_ok=True)
os.makedirs(RES, exist_ok=True)


def num(x):
    """按 R 的习惯打印数字：保留 6 位后去掉多余的 0，保证三语言 CSV 逐位一致"""
    s = f"{round(float(x), 6):.6f}".rstrip("0").rstrip(".")
    return s if s not in ("", "-") else "0"


def write_csv(path, cols, rows):
    """模仿 R 的 write.csv：表头与字符列加引号，数字列不加引号"""
    with open(path, "w", encoding="utf-8") as f:
        f.write(",".join(f'"{c}"' for c in cols) + "\n")
        for r in rows:
            cells = [f'"{v}"' if isinstance(v, str) else num(v) for v in r]
            f.write(",".join(cells) + "\n")


# ------------------------------------------------------------
# 1) 四面地球仪：数路径
# ------------------------------------------------------------
n_faces = 4                                   # 练习 1：改成 6 就是六面地球仪（三语言同名）
face_water = list(range(n_faces + 1))
face_land = [n_faces - f for f in face_water]
ways = [w * l * w2 for w, l, w2 in zip(face_water, face_land, face_water)]
total_ways = sum(ways)
post_globe = [w / total_ways for w in ways]

cols1 = ["p", "水面数", "陆面数", "W1取法", "L1取法", "W2取法", "路径数", "后验"]
p_globe = [fw / n_faces for fw in face_water]      # 比例由"水面数 / 面数"推出
rows1 = [[p, fw, fl, w1, l1, w2, wy, po] for p, fw, fl, w1, l1, w2, wy, po in zip(
    p_globe, face_water, face_land, face_water, face_land, face_water, ways, post_globe)]
write_csv(os.path.join(RES, "01_四面地球仪路径计数_python.csv"), cols1, rows1)
print("路径数 ways =", ways, " 总数 =", total_ways)
print("后验 =", [num(x) for x in post_globe])

# ------------------------------------------------------------
# 2) 20 点网格：三种先验
# ------------------------------------------------------------
W, L = 6, 3
N = W + L
p_grid = np.linspace(0, 1, 20)


def binom_pmf(k, n, p):
    """与 R 的 dbinom、MATLAB 的 binopdf 一致：直接按组合数算，避免库之间的尾位差异"""
    out = []
    for pi in p:
        if pi <= 0:
            out.append(1.0 if k == 0 else 0.0)
        elif pi >= 1:
            out.append(1.0 if k == n else 0.0)
        else:
            out.append(math.comb(n, k) * pi ** k * (1 - pi) ** (n - k))
    return np.array(out)


def grid_post(prior):
    like = binom_pmf(W, N, p_grid)
    post = like * prior
    return p_grid, prior / prior.sum(), like / like.sum(), post / post.sum()


priors = [("平坦先验", np.ones(20)),
          ("exp(-5|p-0.5|)", np.exp(-5 * np.abs(p_grid - 0.5))),
          ("p<0.5 时为 0", np.where(p_grid < 0.5, 0.0, 1.0))]

cols2 = ["先验", "后验均值", "最大点", "区间下", "区间上"]
rows2 = []
curves = {}
for name, pr in priors:
    p, prior_n, like_n, post = grid_post(pr)
    curves[name] = (p, prior_n, like_n, post)
    lo = p[np.argmax(np.cumsum(post) >= 0.055)]
    hi = p[np.argmax(np.cumsum(post) >= 0.945)]
    rows2.append([name, round(float(np.sum(p * post)), 6), round(float(p[np.argmax(post)]), 6),
                  round(float(lo), 4), round(float(hi), 4)])
write_csv(os.path.join(RES, "02_三种先验下的后验汇总_python.csv"), cols2, rows2)
for r in rows2:
    print(r)

# ------------------------------------------------------------
# 3) 解析解 Beta(W+1, L+1)
# ------------------------------------------------------------
a, b = W + 1, L + 1
an_rows = [["后验均值", round(a / (a + b), 6)],
           ["后验众数", round((a - 1) / (a + b - 2), 6)],
           ["后验标准差", round(math.sqrt(a * b / ((a + b) ** 2 * (a + b + 1))), 6)],
           ["89% 区间下", round(float(beta_dist.ppf(0.055, a, b)), 6)],
           ["89% 区间上", round(float(beta_dist.ppf(0.945, a, b)), 6)]]
write_csv(os.path.join(RES, "03_解析后验对照_python.csv"), ["指标", "解析值"], an_rows)
print(f"解析后验 Beta({a},{b}):", an_rows)

# ------------------------------------------------------------
# 4) 从后验到预测（beta-二项，确定性）
# ------------------------------------------------------------
def logB(x, y):
    return lgamma(x) + lgamma(y) - lgamma(x + y)


pred_next = a / (a + b)
p_hat = pred_next          # 对照：把后验均值当成固定 p（会丢掉后验宽度）
pred_rows = [[k,
              round(math.exp(math.log(math.comb(2, k)) + logB(a + k, b + 2 - k) - logB(a, b)), 6),
              round(math.comb(2, k) * p_hat ** k * (1 - p_hat) ** (2 - k), 6)]
             for k in range(3)]
write_csv(os.path.join(RES, "04_预测分布_python.csv"),
          ["未来两次里的水数", "概率（整条后验）", "概率（固定均值）"], pred_rows)
print("下一次取到水的概率 =", round(pred_next, 6), pred_rows)

# ------------------------------------------------------------
# 5) 测量误差：路径计数
# ------------------------------------------------------------
true_water, true_land = 3, 1
# 练习 2：判定取法改成 10 种（判对 9、判错 1）→ judge_correct, judge_wrong = 9, 1
judge_correct, judge_wrong = 2, 1
rows5 = [["真实是水", "记录为水", true_water, judge_correct, true_water * judge_correct],
         ["真实是水", "记录为陆", true_water, judge_wrong,   true_water * judge_wrong],
         ["真实是陆", "记录为水", true_land,  judge_wrong,   true_land * judge_wrong],
         ["真实是陆", "记录为陆", true_land,  judge_correct, true_land * judge_correct]]
write_csv(os.path.join(RES, "05_误分类路径计数_python.csv"),
          ["来源", "判定结果", "真样本数", "每次判定的取法", "路径数"], rows5)
ways_obs_water = sum(r[4] for r in rows5 if r[1] == "记录为水")
ways_obs_land = sum(r[4] for r in rows5 if r[1] == "记录为陆")
print(f"观测到水 = 6 + 1 = {ways_obs_water}；观测到陆 = 3 + 2 = {ways_obs_land}；合计 = {ways_obs_water + ways_obs_land}")

# ============================================================
# 配图：01、02 两张对照图（R 版是主图，这里只做同内容对照）
# ============================================================
fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))

# 图 2 对照：先验、似然、后验
p, prior_n, like_n, post = curves["平坦先验"]
ax[0].plot(p, prior_n, color="#b9b9b9", lw=1.4, label="先验（平坦）")
ax[0].plot(p, like_n, color="#c1462c", lw=1.4, ls="--", label="似然（归一化）")
ax[0].plot(p, post, color="#167d80", lw=1.6, marker="o", ms=3, label="后验")
ax[0].set(title="把每个候选比例算一遍，就得到后验",
          xlabel="水面比例 p", ylabel="归一化权重")
ax[0].legend(loc="upper center", ncol=3, fontsize=8, frameon=False, bbox_to_anchor=(0.5, 1.22))

# 图 3 对照：三种先验下的后验
for (name, _), col in zip(priors, ["#b9b9b9", "#167d80", "#c1462c"]):
    pp, _, _, po = curves[name]
    ax[1].plot(pp, po, lw=1.6, color=col, label=name)
ax[1].set(title="换一个先验，后验被拉多远？",
          xlabel="水面比例 p", ylabel="后验概率")
ax[1].legend(loc="upper center", ncol=2, fontsize=8, frameon=False, bbox_to_anchor=(0.5, 1.22))

fig.tight_layout()
fig.savefig(os.path.join(OUT, "06-Python版-网格与先验.png"), dpi=200)
plt.close(fig)
print("Python 图已保存；numpy", np.__version__)
