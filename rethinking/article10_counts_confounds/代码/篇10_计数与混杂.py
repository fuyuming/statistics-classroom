# -*- coding: utf-8 -*-
# ============================================================
# 精读 10｜计数里的混杂：零膨胀、暴露量与调整（Python 版，与 R/MATLAB 逐位一致）
# 对应：McElreath 2023 第 10 讲 Counts and Confounds；教材第 11、12 章
# 用法： python 代码/篇10_计数与混杂.py
#
# 主要拟合工具＝显式泊松 IRLS（固定 12 次迭代、不提前停止），批次系数被真正估出来。
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


treat = np.array([0,0,0,0,0, 0,0,0, 1,1,1, 1,1,1,1,1], float)
batch = np.array([1,1,1,1,1, 0,0,0, 1,1,1, 0,0,0,0,0], float)
plate = np.array([10,20,10,20,10, 20,10,20, 10,20,10, 20,10,20,10,20], float)
ok = np.array([1,1,1,0,1, 1,1,1, 1,1,0, 1,1,1,1,1], float)
y = np.array([6,13,5,0,7, 6,3,5, 8,19,0, 9,4,8,5,7], float)
offs = np.log(plate / 10)

write_csv(["单元", "处理", "批次", "体积微升", "接种成功", "菌落数"],
          [[i + 1, int(treat[i]), int(batch[i]), int(plate[i]), int(ok[i]), int(y[i])] for i in range(16)],
          os.path.join(RES, "01_计数数据_python.csv"))


def pois_irls(X, yy, off, n_iter=12):
    """泊松 IRLS（固定 12 次迭代）。off＝偏移项（log 暴露量），系数固定为 1。"""
    b = np.zeros(X.shape[1])
    for _ in range(n_iter):
        eta = X @ b + off
        mu = np.exp(eta)
        W = mu
        z = eta + (yy - mu) / mu - off
        b = np.linalg.solve(X.T @ (W[:, None] * X), X.T @ (W * z))
    return b, X.T @ (W[:, None] * X)


XA = np.column_stack([np.ones(16), treat])              # 偏移单独传，不进设计矩阵
XB = np.column_stack([np.ones(16) + offs, treat, batch])
XB = np.column_stack([np.ones(16), treat, batch])
bA, infoA = pois_irls(XA, y, offs)
bB, infoB = pois_irls(XB, y, offs)


def se_log(info, idx=1):
    return float(np.sqrt(np.linalg.inv(info)[idx, idx]))


rows2 = []
for name, b, info in [("只放处理（含体积偏移）", bA, infoA), ("再加批次（混杂）", bB, infoB)]:
    rr = math.exp(b[1]); s = se_log(info)
    rows2.append([name, r4(rr), r4(math.exp(math.log(rr) - 1.598 * s)), r4(math.exp(math.log(rr) + 1.598 * s)),
                  "（未估计）" if len(b) <= 2 else num(b[2], 4)])
write_csv(["模型", "率比", "下界89", "上界89", "批次系数"], rows2, os.path.join(RES, "02_率比调整前后_python.csv"))
for r in rows2:
    print("  ", r)

muA = np.exp(XA @ bA + offs)
p_zero_pois = float(np.mean(np.exp(-muA)))
keep = ok == 1
bC, _ = pois_irls(XB[keep], y[keep], offs[keep])
muC = np.exp(XB[keep] @ bC + offs[keep])
write_csv(["量", "值", "说明"],
          [["观测到的零比例", r4(float(np.mean(y == 0))), "16 个单元里有 2 个计数为 0（都是没接种上的）"],
           ["单一泊松预测的零比例（含全部单元）", r4(p_zero_pois), "泊松把零当成「恰好长了 0 个」，于是严重低估零的出现"],
           ["只看接种成功单元的率比", r4(math.exp(bC[1])), "把没接种上的单元排除后重估"],
           ["只看接种成功单元：预测的零比例", r4(float(np.mean(np.exp(-muC)))), "排除结构性零之后，泊松对「真零」的预测就合理了"]],
          os.path.join(RES, "03_零的比例_python.csv"))

muB = np.exp(XB @ bB + offs)
pearson = float(((y - muB) ** 2 / muB).sum())
df_res = 16 - XB.shape[1]
phi = pearson / df_res
se_pois = se_log(infoB)
se_quasi = se_pois * math.sqrt(phi)
write_csv(["量", "值"],
          [["Pearson 卡方", r4(pearson)], ["残差自由度", r4(df_res)], ["离散系数 φ（卡方/自由度）", r4(phi)],
           ["泊松假设下的 log 率比标准误", r4(se_pois)], ["按 φ 校正后的标准误", r4(se_quasi)], ["校正倍数", r4(math.sqrt(phi))]],
          os.path.join(RES, "04_过离散_python.csv"))
pearson_fail = float(((y[~keep] - muB[~keep]) ** 2 / muB[~keep]).sum())
pearson_share = pearson_fail / pearson
phi_succ = float(((y[keep] - muC) ** 2 / muC).sum() / (int(keep.sum()) - XB.shape[1]))
write_csv(["量", "值"],
          [["Pearson 卡方", r4(pearson)], ["残差自由度", r4(df_res)], ["离散系数 φ（卡方/自由度）", r4(phi)],
           ["泊松假设下的 log 率比标准误", r4(se_pois)], ["按 φ 校正后的标准误", r4(se_quasi)], ["校正倍数", r4(math.sqrt(phi))],
           ["两个失败单元贡献的 Pearson 量", r4(pearson_fail)], ["其占 Pearson 总量的比例", r4(pearson_share)],
           ["只看成功子集的离散系数 φ（14 个观测、3 个参数）", r4(phi_succ)]],
          os.path.join(RES, "04_过离散_python.csv"))

write_csv(["单元", "处理", "批次", "体积微升", "接种成功", "菌落数", "模型B预测均值"],
          [[i + 1, int(treat[i]), int(batch[i]), int(plate[i]), int(ok[i]), int(y[i]), r4(muB[i])] for i in range(16)],
          os.path.join(RES, "05_拟合值_python.csv"))

b_adj = math.log(math.exp(bB[1]))
write_csv(["量", "值"],
          [["调整批次后的率比", r4(math.exp(b_adj))],
           ["泊松假设下的 89% 区间下", r4(math.exp(b_adj - 1.598 * se_pois))],
           ["泊松假设下的 89% 区间上", r4(math.exp(b_adj + 1.598 * se_pois))],
           ["按 φ 校正后的 89% 区间下", r4(math.exp(b_adj - 1.598 * se_quasi))],
           ["按 φ 校正后的 89% 区间上", r4(math.exp(b_adj + 1.598 * se_quasi))]],
          os.path.join(RES, "06_校正后的率比_python.csv"))

# 两张对照图
fig, ax = plt.subplots(figsize=(7.4, 4.4), dpi=150)
for bt, col in [(0, TEAL), (1, RED)]:
    m = batch == bt
    ax.scatter(plate[m], y[m], s=45, color=col, label=f"批次 {bt}")
ax.set_title("先看清三件事：体积、批次、以及两个空皿", fontsize=11)
ax.set_xlabel("接种体积（微升）"); ax.set_ylabel("每皿菌落数"); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "06-Python版-数据总览.png")); plt.close(fig)

fig, ax = plt.subplots(figsize=(7.4, 4.4), dpi=150)
ax.bar([r[0] for r in rows2], [r[1] for r in rows2], color=[PINK, TEAL], width=0.5)
ax.errorbar([r[0] for r in rows2], [r[1] for r in rows2],
            yerr=[[r[1] - r[2] for r in rows2], [r[3] - r[1] for r in rows2]], fmt="none", ecolor="#163d48", capsize=6)
ax.axhline(1, ls="--", color="#9aa8a8")
ax.set_title("不调整批次，率比就偏了", fontsize=11)
ax.set_ylabel("率比（处理 vs 对照）")
fig.tight_layout(); fig.savefig(os.path.join(OUT, "07-Python版-率比调整.png")); plt.close(fig)

print("完成：2 张对照图 + 6 份 CSV 已写入")
