# %% 样本率与率差：精确分布和正态近似
# 输入是教学构造数据；所有分布计算均确定性完成，不用模拟。
from pathlib import Path
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import norm, beta, binom

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUT = ROOT / "运行结果" / "Python"
FIG = ROOT.parent / "文章配图"
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)
data = pd.read_csv(ROOT / "数据" / "共享教学数据.csv", skipinitialspace=True)
print("共享数据字段：", list(data.columns))
print(data.head(2).to_string(index=False))
plt.rcParams["font.sans-serif"] = ["STHeiti", "PingFang SC", "Arial Unicode MS", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["font.size"] = 11
plt.rcParams["axes.titlesize"] = 10

def binom_pmf(n, p):
    k = np.arange(n + 1)
    # scipy's stable implementation avoids integer-combination overflow at n=1000.
    values = binom.pmf(k, n, p)
    return k, values

def normal_cell(k, mean, sd, width=1.0):
    if sd == 0:
        return np.where(np.asarray(k) == mean, 1.0, 0.0)
    k = np.asarray(k, dtype=float)
    return norm.cdf(k + width/2, loc=mean, scale=sd) - norm.cdf(k - width/2, loc=mean, scale=sd)

def one_curve(n, p):
    k, pmf = binom_pmf(n, p)
    # At k=0 and n, use tail intervals matching the boundary convention.
    approx = normal_cell(k, n*p, math.sqrt(n*p*(1-p)))
    sd = math.sqrt(n*p*(1-p))
    approx[0] = norm.cdf(0.5, loc=n*p, scale=sd)
    approx[-1] = norm.sf(n-0.5, loc=n*p, scale=sd)
    return k, pmf, approx

def save_one(n, p, name, title, xlim=None):
    k, exact, approx = one_curve(n, p)
    pd.DataFrame({"k": k, "rate": k/n, "exact_pmf": exact, "normal_cell": approx}).to_csv(OUT / f"{name}.csv", index=False)
    shown = np.ones(len(k), dtype=bool) if xlim is None else (k >= xlim[0]) & (k <= xlim[1])
    fig, ax = plt.subplots(figsize=(5.2, 4.0), dpi=180)
    ax.bar(k[shown], exact[shown], color="#4C78A8", alpha=.75, label="精确概率")
    ax.plot(k[shown], approx[shown], "o", ms=3.5, color="#E45756", label="连续性校正正态单格概率")
    ax.set(xlabel="计数 k（X）", ylabel="概率", title=title)
    # 图例上方预留空间，避免遮住概率最高的柱和点。
    highest_probability = max(float(np.max(exact[shown])), float(np.max(approx[shown])))
    ax.set_ylim(0, highest_probability * 1.35)
    ax.grid(axis="y", alpha=.25)
    ax.legend(fontsize=9)
    fig.tight_layout()
    suffix = ".png" if name == "one_p05_100" else ("_低概率.png" if name == "one_p01_100" else f"_{name}.png")
    fig.savefig(FIG / ("01_样本率精确与近似" + suffix))
    plt.close(fig)
    return k, exact, approx

one_results = {}
for _, row in data[data["kind"].eq("one")].iterrows():
    n, p, name = int(row["n"]), float(row["p"]), row["scenario"]
    xlim = (0, 8) if p < 0.1 else None
    title = f"样本率：n={n}，p={p:g}"
    one_results[name] = save_one(n, p, name, title, xlim)
    print(f"{name}: rate_variance(pq/n)={p*(1-p)/n:.8g}, P(X=0)={(1-p)**n:.8g}")

def difference_curve(n, p1, p2):
    k1, a = binom_pmf(n, p1); k2, b = binom_pmf(n, p2)
    diff_counts = np.arange(-n, n+1)
    # coefficient for X1-X2=d: sum over x2 of a[d+x2] b[x2]
    pmf = np.zeros(len(diff_counts))
    for j, d in enumerate(diff_counts):
        lo, hi = max(0, -d), min(n, n-d)
        pmf[j] = sum(a[d+x2] * b[x2] for x2 in range(lo, hi+1))
    sd = math.sqrt(p1*(1-p1)/n + p2*(1-p2)/n)
    rates = diff_counts / n
    mean = p1 - p2
    approx = normal_cell(rates, mean, sd, width=1/n)
    # Merge normal tails into the finite-support endpoint cells.
    approx[0] = norm.cdf(-1 + 0.5/n, loc=mean, scale=sd)
    approx[-1] = norm.sf(1 - 0.5/n, loc=mean, scale=sd)
    return diff_counts, rates, pmf, approx, sd

fig, axes = plt.subplots(2, 1, figsize=(5.2, 8.0), dpi=180)
for ax, (_, row) in zip(axes, data[data["kind"].eq("difference")].iterrows()):
    n, p1, p2, name = int(row["n"]), float(row["p1"]), float(row["p2"]), row["scenario"]
    lim = (-8, 8) if n >= 1000 else (-20, 20)
    title = f"率差 D：n={n}，p1={p1:g}，p2={p2:g}"
    d, rates, exact, approx, sd = difference_curve(n, p1, p2)
    pd.DataFrame({"difference_count": d, "D": rates, "exact_pmf": exact, "normal_cell": approx}).to_csv(OUT / f"{name}.csv", index=False)
    shown = (d >= lim[0]) & (d <= lim[1])
    ax.bar(d[shown], exact[shown], color="#59A14F", alpha=.75, label="精确概率")
    ax.plot(d[shown], approx[shown], "o", ms=3.5, color="#F28E2B", label="连续性校正正态单格概率")
    ax.set(xlabel="差值计数 k=X1−X2（D=k/n）", ylabel="概率", title=title)
    # 图例上方预留空间，避免遮住概率最高的柱和点。
    highest_probability = max(float(np.max(exact[shown])), float(np.max(approx[shown])))
    ax.set_ylim(0, highest_probability * 1.35)
    ax.grid(axis="y", alpha=.25)
    ax.legend(fontsize=9)
fig.tight_layout(); fig.savefig(FIG / "02_率差精确与近似.png"); plt.close(fig)

# Observed comparison: pooled Z/P, Pearson Q, unpooled Wald CI.
x1, n1, x2, n2 = 60, 300, 20, 200
r1, r2 = x1/n1, x2/n2; pooled = (x1+x2)/(n1+n2)
se0 = math.sqrt(pooled*(1-pooled)*(1/n1+1/n2)); z = (r1-r2)/se0
wald_se = math.sqrt(r1*(1-r1)/n1 + r2*(1-r2)/n2)
z_p = 2*norm.sf(abs(z)); q = z*z; ci = (r1-r2-1.96*wald_se, r1-r2+1.96*wald_se)
summary = pd.DataFrame([{"pooled_p": pooled, "z": z, "two_sided_p": z_p, "pearson_Q": q, "unpooled_se": wald_se, "wald95_low": ci[0], "wald95_high": ci[1]}])
summary.to_csv(OUT / "observed_comparison.csv", index=False); print(summary.to_string(index=False))

# Boundary x=0,n=20: Wald versus exact Clopper-Pearson.
n = 20; x = 0; wald = (0.0, 0.0); cp = (0.0, beta.ppf(.975, 1, n))
pd.DataFrame([{"x": x, "n": n, "wald_low": wald[0], "wald_high": wald[1], "clopper_pearson_low": cp[0], "clopper_pearson_high": cp[1]}]).to_csv(OUT / "boundary_zero_20.csv", index=False)
print(f"boundary x=0,n=20: Wald={wald}, Clopper-Pearson={cp}")
