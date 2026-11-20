# %% 问题与输入
# 假设8名不同患者各有一段30秒心电；同一行A/B为该记录的心率估计（次/分）。
# 目标总体平均差B−A；不同对子独立，按差值近似正态演算；本例ID唯一且无缺失。
# 在项目根目录用Spyder分节或整文件运行，不依赖__file__。
import numpy as np
import pandas as pd
from scipy import stats
rates = pd.read_csv("data/ecg_rates.csv")
print(rates)
A = rates["A_rate"].to_numpy(float)
B = rates["B_rate"].to_numpy(float)
d = B - A

# %% 手动计算：差值、标准误和配对t
n = len(d)
dbar = d.mean()
sd_d = d.std(ddof=1)                 # ddof=1表示用n−1作样本SD分母
se = sd_d / np.sqrt(n)              # 平均差标准误，单位次/分
df = n - 1
t_value = dbar / se                 # H0：总体平均差为0
p_two = 2 * stats.t.sf(abs(t_value), df)  # 对称t分布的两尾概率
ci = dbar + stats.t.ppf([.025,.975], df) * se
result = pd.DataFrame([dict(n=n, mean_diff=dbar, sd_diff=sd_d, se=se,
                           t=t_value, df=df, p=p_two, low=ci[0], high=ci[1])])
print(result.to_string(index=False))

# %% 调用函数：scipy.stats配对t
fit = stats.ttest_rel(B, A, alternative="two-sided")
print(fit)                          # statistic、pvalue、df对应手算结果
# 平均差−2、t约−6.1101、df=7、P约.0004862；均值差CI单位次/分。
# 没有金标准，这些数值不能说明算法准确性。

# %% 练习：只在练习副本中改变B的排列
print(dict(mean_B=B.mean(), mean_B_reordered=B[::-1].mean(),
           sd_diff=d.std(ddof=1), sd_diff_reordered=(B[::-1]-A).std(ddof=1)))
# B列均值不变而差值SD改变；实测分析按真实ID连接。
result.to_csv("results/Python_paired.csv", index=False)
