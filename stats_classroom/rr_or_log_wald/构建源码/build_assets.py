"""生成本文的18幅公式和2幅科学示意图；无模拟数据。
运行：python 构建源码/build_assets.py。依赖numpy/scipy/matplotlib。
封面沿用第10课已确认的无姓名课堂封面，不由本脚本改绘。
"""
from pathlib import Path
import os
import json

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".matplotlib_cache"))
import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt

OUT = ROOT / "文章配图"
OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans", "mathtext.fontset": "dejavusans",
                     "font.size": 16, "axes.spines.top": False,
                     "axes.spines.right": False, "savefig.facecolor": "white"})
BLUE = "#0000F5"
INK = "#27343a"
ORANGE = "#ce6e21"

FORMULAS = [
 [r"\hat p_1=\frac{50}{100}=0.50,\quad\hat p_0=\frac{30}{100}=0.30",
  r"\widehat{\mathrm{RD}}=0.50-0.30=0.20",
  r"\widehat{\mathrm{RR}}=\frac{0.50}{0.30}=\frac{5}{3}\approx1.67"],
 [r"\widehat{\mathrm{odds}}_1=\frac{50}{50}=1",
  r"\widehat{\mathrm{odds}}_0=\frac{30}{70}=\frac{3}{7}",
  r"\widehat{\mathrm{OR}}=\frac{50/50}{30/70}=\frac{7}{3}\approx2.33"],
 [r"\log(\mathrm{RR})=\log p_1-\log p_0",
  r"\log(\mathrm{OR})=\mathrm{logit}(p_1)-\mathrm{logit}(p_0)",
  r"\log 1=0"],
 [r"\ell=\log\hat\theta,\quad s=\widehat{\mathrm{SE}}(\log\hat\theta)",
  r"\mathrm{CI}_{\log\theta}=[\ell-zs,\ \ell+zs]",
  r"\mathrm{CI}_{\theta}=[e^{\ell-zs},\ e^{\ell+zs}]"],
 [r"\log\widehat{\mathrm{RR}}\approx0.5108",
  r"\mathrm{CI}_{\log\mathrm{RR}}\approx0.5108\ \pm\ 0.3578",
  r"\approx[0.1530,\ 0.8687]",
  r"\mathrm{CI}_{\mathrm{RR}}\approx[1.1653,\ 2.3837]"],
 [r"\log\widehat{\mathrm{OR}}\approx0.8473",
  r"\mathrm{CI}_{\log\mathrm{OR}}\approx0.8473\ \pm\ 0.5802",
  r"\approx[0.2671,\ 1.4275]",
  r"\mathrm{CI}_{\mathrm{OR}}\approx[1.3062,\ 4.1681]"],
 [r"[\exp(\log\hat\theta-h),\ \exp(\log\hat\theta+h)]",
  r"=[\hat\theta/e^h,\ \hat\theta\,e^h]"],
 [r"\mathrm{Var}(\hat p)=\frac{p(1-p)}{n}"],
 [r"g(p+h)=g(p)+g'(p)h",
  r"\quad+\frac{g''(p)}{2}h^2+R_2(h)",
  r"g(p+h)\approx g(p)+g'(p)h"],
 [r"g(\hat p)\approx g(p)+g'(p)(\hat p-p)",
  r"g(\hat p)-g(p)\approx g'(p)(\hat p-p)"],
 [r"\mathrm{Var}[g(\hat p)]",
  r"\approx\mathrm{Var}[g(p)+g'(p)(\hat p-p)]",
  r"=[g'(p)]^2\mathrm{Var}(\hat p)",
  r"\mathrm{SE}[g(\hat p)]\approx|g'(p)|\,\mathrm{SE}(\hat p)"],
 [r"g(p)=\log p,\quad g'(p)=\frac{1}{p}",
  r"\mathrm{Var}(\log\hat p)\approx\frac{1}{p^2}\mathrm{Var}(\hat p)",
  r"=\frac{1}{p^2}\times\frac{p(1-p)}{n}",
  r"=\frac{p(1-p)}{np^2}",
  r"=\frac{1-p}{np}"],
 [r"\log\widehat{\mathrm{RR}}=\log\hat p_1-\log\hat p_0",
  r"\mathrm{Var}(\log\widehat{\mathrm{RR}})",
  r"\approx\frac{1-p_1}{n_1p_1}+\frac{1-p_0}{n_0p_0}"],
 [r"\widehat V_{\log\mathrm{RR}}=\frac{1}{a}-\frac{1}{n_1}+\frac{1}{c}-\frac{1}{n_0}",
  r"\widehat{\mathrm{SE}}(\log\widehat{\mathrm{RR}})",
  r"=\sqrt{\frac{1}{50}-\frac{1}{100}+\frac{1}{30}-\frac{1}{100}}",
  r"\approx0.1826"],
 [r"g(p)=\mathrm{logit}(p)=\log\frac{p}{1-p}",
  r"g'(p)=\frac{1}{p}+\frac{1}{1-p}=\frac{1}{p(1-p)}",
  r"\mathrm{Var}[\mathrm{logit}(\hat p)]",
  r"\approx\frac{1}{p^2(1-p)^2}\frac{p(1-p)}{n}",
  r"=\frac{1}{np(1-p)}"],
 [r"\frac{1}{n_1\hat p_1(1-\hat p_1)}",
  r"=\frac{n_1}{ab}=\frac{a+b}{ab}=\frac{1}{a}+\frac{1}{b}"],
 [r"\log\widehat{\mathrm{OR}}",
  r"=\mathrm{logit}(\hat p_1)-\mathrm{logit}(\hat p_0)"],
 [r"\widehat V_{\log\mathrm{OR}}=\frac{1}{a}+\frac{1}{b}+\frac{1}{c}+\frac{1}{d}",
  r"\widehat{\mathrm{SE}}(\log\widehat{\mathrm{OR}})",
  r"=\sqrt{\frac{1}{50}+\frac{1}{50}+\frac{1}{30}+\frac{1}{70}}",
  r"\approx0.2960"],
]

formula_bounds = []
for index, lines in enumerate(FORMULAS, start=1):
    # 固定画布宽度及字体，长公式换行，避免手机上因缩放而变小。
    height = max(1.25, 0.85 * len(lines) + 0.25)
    fig = plt.figure(figsize=(7.2, height), dpi=200, facecolor="white")
    texts = []
    for row, line in enumerate(lines):
        y = 1 - (row + 0.5) / len(lines)
        texts.append(fig.text(.5, y, f"${line}$", ha="center", va="center",
                              fontsize=22, color=INK))
    fig.canvas.draw()
    width, canvas_height = fig.canvas.get_width_height()
    for item in texts:
        bounds = item.get_window_extent(fig.canvas.get_renderer())
        if bounds.x0 < 14 or bounds.x1 > width - 14:
            raise RuntimeError(f"公式{index}过宽，请换行：{item.get_text()}")
        if bounds.y0 < 0 or bounds.y1 > canvas_height:
            raise RuntimeError(f"公式{index}纵向裁切")
    fig.savefig(OUT / f"公式{index:02d}.png", dpi=200)
    formula_bounds.append({"equation": index, "width": width, "height": canvas_height})
    plt.close(fig)

# 图1：同一组区间端点，分别在两个尺度上定位。
rr = 5 / 3
z = norm.ppf(.975)
se = np.sqrt(1 / 50 - 1 / 100 + 1 / 30 - 1 / 100)
log_center = np.log(rr)
log_bounds = log_center + np.array([-1, 1]) * z * se
ratio_bounds = np.exp(log_bounds)
fig, axes = plt.subplots(2, 1, figsize=(7.2, 6.5), dpi=200)
for ax, center, bounds, null, label in [
    (axes[0], log_center, log_bounds, 0, "log RR"),
    (axes[1], rr, ratio_bounds, 1, "RR (original scale)")
]:
    ax.hlines(0, bounds[0], bounds[1], color=BLUE, lw=4)
    ax.plot([bounds[0], bounds[1]], [0, 0], "|", color=BLUE, ms=17, mew=2)
    ax.plot(center, 0, "o", color=BLUE, ms=10)
    ax.axvline(null, color="#88949d", ls="--", lw=1.5)
    for x, y in [(bounds[0], .18), (center, -.3), (bounds[1], .18)]:
        ax.text(x, y, f"{x:.3f}", ha="center", fontsize=15, color=INK)
    ax.set_ylim(-.6, .6)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_position(("data", -.48))
    ax.set_xlabel(label, labelpad=12, fontsize=17)
    ax.set_xlim((-.13, 1.08) if null == 0 else (.8, 2.65))
    ax.grid(axis="x", alpha=.12)
fig.tight_layout(h_pad=2.8)
fig.savefig(OUT / "图01_区间的两个尺度.png")
plt.close(fig)

# 图2：两处log曲线及切线，无模拟；由确定性函数得到。
fig, axes = plt.subplots(2, 1, figsize=(7.2, 7.3), dpi=200)
curve_rows = []
for ax, center in zip(axes, [.50, .30]):
    p = np.linspace(center - .04, center + .04, 161)
    actual = np.log(p)
    tangent = np.log(center) + (p - center) / center
    ax.plot(p, actual, color=BLUE, lw=2.5, label="log p")
    ax.plot(p, tangent, color=ORANGE, ls="--", lw=2.2, label="Tangent line")
    ax.scatter([center], [np.log(center)], color=INK, s=28, zorder=5)
    ax.set_title(f"p = {center:.2f}: slope = {1/center:.3f}", fontsize=17,
                 loc="left", pad=12, color=INK)
    ax.set_xlabel("Probability p", fontsize=16)
    ax.set_ylabel("log p", fontsize=16)
    ax.tick_params(labelsize=13)
    ax.legend(frameon=False, loc="upper left", fontsize=13)
    ax.grid(alpha=.13)
    for x, y, t in zip(p, actual, tangent):
        curve_rows.append([center, x, y, t])
fig.tight_layout(h_pad=2.2)
fig.savefig(OUT / "图02_局部斜率.png")
plt.close(fig)

np.savetxt(ROOT / "构建源码" / "局部斜率_绘图坐标.csv", curve_rows,
           delimiter=",", header="center,p,log_p,tangent", comments="")
(ROOT / "构建源码" / "公式画布检查.json").write_text(
    json.dumps(formula_bounds, ensure_ascii=False, indent=2), encoding="utf-8")
print("已生成18幅公式及2幅示意图；公式画布无越界。")
