"""统计课堂29｜独立样本率差 A-B，UTF-8。
问题：假设A组9/20阳性，B组13/20阳性；估计总体率差的双侧区间。
依赖：numpy、pandas、scipy、statsmodels、matplotlib，见requirements.txt。
Spyder：工作目录设为项目根或代码目录，先逐节运行；F5可运行整文件。
默认控制台显示结果、图窗口显示图。SAVE_RESULTS=True才写运行结果/Python。
Newcombe与AC均调用Python包现成函数；不调用R。
"""
# %% 准备与参数：先运行本节，后面的小节沿用这里的变量
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.proportion import confint_proportions_2indep
import matplotlib.pyplot as plt

confidence = 0.95
SAVE_RESULTS = False
# Spyder分节时可能没有__file__，先从Current working directory找共享CSV。
root = Path.cwd()
if not (root / "数据/independent_counts.csv").is_file():
    root = root.parent
if not (root / "数据/independent_counts.csv").is_file() and "__file__" in globals():
    root = Path(__file__).resolve().parents[1]
data_path = root / "数据/independent_counts.csv"

# %% 问题1：两批独立的人，A率减B率是多少？
# CSV是两行假设汇总计数，没有个体记录。先看字段和行，再提取A/B。
data = pd.read_csv(data_path)
print("字段：", list(data.columns))
print(data.to_string(index=False))
if not {"group", "positive", "total"}.issubset(data.columns):
    raise ValueError("CSV需要group、positive、total字段")
if len(data) != 2 or set(data["group"]) != {"A", "B"}:
    raise ValueError("CSV必须恰有A、B各一行")
counts = data[["positive", "total"]].to_numpy(dtype=float)
if not np.isfinite(counts).all() or not np.equal(counts, np.floor(counts)).all():
    raise ValueError("计数必须是有限整数，不能缺失")
if (data["total"] <= 0).any() or (data["positive"] < 0).any():
    raise ValueError("总人数必须为正，阳性数不能为负")
if (data["positive"] > data["total"]).any() or not 0 < confidence < 1:
    raise ValueError("阳性数不能超过总人数，置信水平必须在0与1之间")
# loc选出组别相符的行；iloc[0]取这一组唯一的一行。
A = data.loc[data["group"] == "A"].iloc[0]
B = data.loc[data["group"] == "B"].iloc[0]
xA = int(A["positive"])
nA = int(A["total"])
xB = int(B["positive"])
nB = int(B["total"])
pA = xA / nA
pB = xB / nB
d = pA - pB
print(f"A阳性/阴性={xA}/{nA-xA}；B阳性/阴性={xB}/{nB-xB}")
print(f"样本率A={pA:.2%}，B={pB:.2%}；A-B={100*d:.2f}个百分点。")

# %% 问题2：Wald怎样从标准误得到区间？
# 思路：两组独立，差值方差是各自方差相加；用样本率估计未知总体率。
alpha = 1 - confidence
# ppf输入左侧累计概率；双侧95%时为.975，输出约1.96。
z = stats.norm.ppf(1 - alpha / 2)
se = np.sqrt(pA * (1-pA) / nA + pB * (1-pB) / nB)
wald_by_hand = np.array([d-z*se, d+z*se])
print(f"Wald标准误={se:.8f}，正态临界值={z:.6f}。")
print("Wald手算端点（百分点）：", 100 * wald_by_hand)

# %% 问题3：用三种现成区间函数回答同一个问题
# 四个计数参数的顺序固定A先B后；compare='diff'估计第一组减第二组。
# 返回对象是(lower, upper)二元组；alpha是1-置信水平。
wald = confint_proportions_2indep(
    xA, nA, xB, nB, method="wald", compare="diff", alpha=alpha
)
newcombe = confint_proportions_2indep(
    xA, nA, xB, nB, method="newcomb", compare="diff", alpha=alpha
)
ac = confint_proportions_2indep(
    xA, nA, xB, nB, method="agresti-caffo", compare="diff", alpha=alpha
)
# 三行各为一种方法，两列各为下限/上限。原始点估计单独保留。
intervals = np.array([wald, newcombe, ac])
method_names = ["Wald", "Newcombe", "Agresti-Caffo"]
results = pd.DataFrame(intervals, columns=["lower", "upper"])
results.insert(0, "method", method_names)
results["estimate"] = d
results["lower_pp"] = 100 * results["lower"]
results["upper_pp"] = 100 * results["upper"]
results["width_pp"] = results["upper_pp"] - results["lower_pp"]
print(results.to_string(index=False, float_format=lambda value: f"{value:.8f}"))
# AC每组各加1个阳性和1个阴性：这里只为区间运算，不改CSV的原始计数。
pA_adjusted = (xA+1) / (nA+2)
pB_adjusted = (xB+1) / (nB+2)
ac_center = pA_adjusted - pB_adjusted
print(f"AC计算中心={100*ac_center:.4f}个百分点；原始率差仍为{100*d:.2f}。")
print(f"Newcombe的{confidence:.0%}区间为{100*newcombe[0]:.2f}至{100*newcombe[1]:.2f}个百分点。")
if newcombe[0] <= 0 <= newcombe[1]:
    print("该区间包含0，尚不能确定总体率差方向；不等于证明两总体率相等。")
else:
    print("该区间不含0，方向由区间正负判断；差异大小仍需结合研究背景。")
print("这些是近似区间，不保证所有真实率下都恰好达到所标置信水平。")

# %% 问题4：把区间和原始差值放在一张图中
# 使用英文图签，便于跨平台显示；横轴已经乘100，单位为百分点。
fig, ax = plt.subplots(figsize=(8, 4.8))
y = np.arange(3)
ax.hlines(y, results["lower_pp"], results["upper_pp"], color="#0000F5", linewidth=3)
ax.scatter(np.repeat(100*d, 3), y, color="black", label="Raw A-B", zorder=3)
ax.scatter(100*ac_center, 2, marker="D", color="#E89623", label="AC center", zorder=4)
ax.axvline(0, color="gray", linestyle="--")
ax.set_yticks(y, method_names)
ax.set_ylim(2.5, -0.8)
ax.set_xlabel("Population proportion difference A-B (percentage points)")
ax.set_title(f"{confidence:.0%} two-sided CIs | A {xA}/{nA}, B {xB}/{nB}")
ax.legend(loc="best")
fig.tight_layout()
if SAVE_RESULTS:
    output = root / "运行结果/Python"
    output.mkdir(parents=True, exist_ok=True)
    results.to_csv(output / "intervals.csv", index=False)
    fig.savefig(output / "intervals.png", dpi=160)
plt.show()

# 练习：①CSV改45/100、65/100；②恢复原数据，confidence改.90；③两组均0/20。
# 每次只改一项，从准备/读取节开始按顺序运行；观察点估计、端点、图标题。
