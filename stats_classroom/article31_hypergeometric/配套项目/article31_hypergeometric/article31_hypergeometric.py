# %% 1. 题目、依赖与数据（UTF-8）
# 33名HBsAg阳性孕妇随机分组，比较新生儿HBV感染率。
# 预防组感染4未感染18；非预防组感染5未感染6。
# 来源：颜艳、王彤《医学统计学》第5版例7-4，p111。
# Spyder将项目根设为工作目录，从本节依次运行；F5可运行全文。
from pathlib import Path
from math import comb
import numpy as np
import pandas as pd
from scipy.stats import hypergeom, fisher_exact
import matplotlib.pyplot as plt
import seaborn as sns

SAVE_OUTPUT = False
project_root = Path(__file__).resolve().parent if "__file__" in globals() else Path.cwd()
hbv = pd.read_csv(project_root / "data" / "hbv_table.csv")
print(hbv)
prevention = hbv.loc[hbv["group"] == "prevention"].iloc[0]
control = hbv.loc[hbv["group"] == "control"].iloc[0]
N = int(hbv["total"].sum())
K = int(hbv["infected"].sum())
n = int(prevention["total"])
x_observed = int(prevention["infected"])
# 给定行列合计后，X表示预防组感染人数；观察值为4。

# %% 2. 手算观察表的组合数
positive_choices = comb(K, x_observed)
negative_choices = comb(N - K, n - x_observed)
all_choices = comb(N, n)
observed_probability = positive_choices * negative_choices / all_choices
print(pd.Series(dict(positive_choices=positive_choices, negative_choices=negative_choices,
                     all_choices=all_choices, observed_probability=observed_probability)))
# 分子126×134596，分母193536720，概率约.08762728。

# %% 3. 手算完整分布、矩与尾部
x = np.arange(max(0, n - (N - K)), min(n, K) + 1)
# 循环逐项数选法，每次只改变阳性计数，其余由边际确定。
probability = np.empty(len(x))
for i, count in enumerate(x):
    positive = comb(K, int(count))
    negative = comb(N - K, n - int(count))
    probability[i] = positive * negative / all_choices
mean_x = np.sum(x * probability)
variance_x = np.sum((x - mean_x)**2 * probability)
mean_formula = n * K / N
variance_formula = n * (K / N) * (1 - K / N) * (N - n) / (N - 1)
moments = pd.DataFrame([dict(probability_sum=probability.sum(), mean_x=mean_x,
                            mean_formula=mean_formula, variance_x=variance_x,
                            variance_formula=variance_formula)])
print(moments)
# 总概率1，均值6、方差1.5；方差含有限总体修正因子11/32。
left_tail = probability[x <= x_observed].sum()
right_tail = probability[x >= x_observed].sum()

# %% 4. 调用分布函数核对
# SciPy参数：M为总体总数，n为总体阳性数，N为抽取人数。
probability_function = hypergeom.pmf(x, M=N, n=K, N=n)
probability_table = pd.DataFrame(dict(x=x, probability=probability,
                                      probability_function=probability_function))
print(probability_table)
left_function = hypergeom.cdf(x_observed, M=N, n=K, N=n)
# sf(q)给P(X>q)，要包括4，所以输入4-1。
right_function = hypergeom.sf(x_observed - 1, M=N, n=K, N=n)
print(pd.DataFrame([dict(left_tail=left_tail, left_function=left_function,
                         right_tail=right_tail, right_function=right_function)]))
# 两尾均包含观察值；左尾约.10814723，右尾约.97948005。

# %% 5. 手算Fisher双侧P：概率排序
# 1e-7为浮点相等容差，纳入表概率不大于观察表的结果。
included = probability <= observed_probability * (1 + 1e-7)
probability_table["included"] = included
included_x = x[included]
fisher_manual = probability[included].sum()
print(probability_table)
print(included_x)
# 纳入0、1、2、3、4、9，双侧P约.12104475。

# %% 6. 调用完整检验函数
tab = [[int(prevention["infected"]), int(prevention["uninfected"])],
       [int(control["infected"]), int(control["uninfected"])]]
fisher_fit = fisher_exact(tab, alternative="two-sided")
fisher_check = pd.DataFrame([dict(manual=fisher_manual, function_P=fisher_fit.pvalue)])
print(fisher_check)
# alpha=.05不拒绝两总体感染率相同的H0；观察率差为-27.27个百分点。

# %% 7. 绘图：已经算出的完整离散分布
sns.set_theme(style="whitegrid")
fig, ax = plt.subplots(figsize=(7, 5))
ax.bar(x, probability, width=.72,
       color=np.where(included, "#D76B26", "#6192C3"))
ax.axvline(mean_x, color="#233449", linestyle="--")
ax.set(xlabel="X: infected in prevention group", ylabel="Probability P(X = x)", xticks=x)
fig.tight_layout()
plt.show()
# 橙柱计入双侧P，虚线为均值6；图形默认显示，保存另设开关。

# %% 8. 练习：抽11人，与留下的22人互补
probability_exercise = hypergeom.pmf(x, M=N, n=K, N=11)
remaining_positive = K - x
exercise_mean = np.sum(x * probability_exercise)
exercise_variance = np.sum((x - exercise_mean)**2 * probability_exercise)
remaining_mean = np.sum(remaining_positive * probability_exercise)
remaining_variance = np.sum((remaining_positive - remaining_mean)**2 * probability_exercise)
exercise = pd.DataFrame([dict(exercise_mean=exercise_mean, exercise_variance=exercise_variance,
                             remaining_mean=remaining_mean, remaining_variance=remaining_variance)])
print(exercise)
# 两个阳性数之和恒为9，期望3和6，方差均为1.5。

# %% 9. 可选保存：把首节开关设为True
if SAVE_OUTPUT:
    out = project_root / "results" / "python"
    out.mkdir(parents=True, exist_ok=True)
    probability_table.to_csv(out / "probability_table.csv", index=False)
    moments.to_csv(out / "moments.csv", index=False)
    fisher_check.to_csv(out / "fisher_check.csv", index=False)
    fig.savefig(out / "hypergeometric_distribution.png", dpi=180)
