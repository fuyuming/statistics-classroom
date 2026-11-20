# %% 问题与输入
# 假设100名不同患者各有一段心电，两算法是否标记至少一次室早构成一对结果。
# 目标B−A总体阳性判定率；患者间独立、无缺失；无金标准，不比较准确性。
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.contingency_tables import mcnemar
tab = pd.read_csv("data/ecg_2x2.csv", index_col=0).to_numpy()
print(tab)
# 行是A阳性/阴性，列是B阳性/阴性；单元格为联合类别的人数，并非每人一行。
n = int(tab.sum())                  # 全部配对数，是率差分母
b = int(tab[0,1])                   # Python索引从0开始；A阳/B阴，贡献−1
c_count = int(tab[1,0])             # A阴/B阳，贡献＋1
m = b + c_count                    # 不一致对数，供条件二项法使用

# %% 手动计算：普通配对Wald区间
z = np.concatenate([-np.ones(b), np.zeros(n-m), np.ones(c_count)])
RD = (c_count-b)/n
se_RD = z.std(ddof=1)/np.sqrt(n)     # 样本SD分母n−1，SE用全部配对数
RD_CI = RD + stats.norm.ppf([.025,.975])*se_RD
# 标准正态临界值、无连续性校正；并非下述校正检验的区间反演。

# %% 手动计算：连续性校正McNemar与对称二项精确P
Q = max(abs(c_count-b)-1,0)**2/m
p_chi = stats.chi2.sf(Q, 1)         # 自由度1的χ²右尾概率
p_exact = min(1., 2*stats.binom.cdf(min(b,c_count), m, .5))
result = pd.DataFrame([dict(n=n,b=b,c=c_count,m=m,Q=Q,p=p_chi,exact=p_exact,
                           RD=RD,se=se_RD,low=RD_CI[0],high=RD_CI[1])])
print(result.to_string(index=False))

# %% 调用函数：statsmodels与scipy的对应检验
print(mcnemar(tab, exact=False, correction=True))
fit_exact = stats.binomtest(min(b,c_count), m, .5, alternative="two-sided")
print(fit_exact.pvalue)             # 仅核对精确P，不把条件比例CI当RD CI
# RD=.20，校正Q约12.0333，P约.0005226，精确P约.0003249。
# RD区间约[.09956,.30044]，乘100后为百分点；不一致m不是率差分母。
result.to_csv("results/Python_binary.csv", index=False)
