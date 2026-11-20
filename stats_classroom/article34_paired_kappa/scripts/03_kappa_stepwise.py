# 统计课堂34：Cohen kappa的逐步W表推导
# 2020年颜艳、王彤《医学统计学》第5版例7-3、表7-3，pp.110–111：抗核抗体检测。
# A为免疫荧光法，B为乳胶凝集法；四格计数依次为11、12、2、33。
# 目标：从三个0/1均数计算κ的方差、SE与1.96 Wald置信区间。
# 题目共有58名疑似SLE患者；X=A阳性，Y=B阳性，C=两法一致。
# 同一患者的X/Y/C有关联，不同患者相互独立；四格联合表因此用于W方差。
# 工作目录为项目根目录；默认不保存，ARTICLE34_SAVE_RESULTS=1时才写结果。

# %% 03.1 输入与完整病例
import os
import numpy as np
import pandas as pd
from statsmodels.stats.inter_rater import cohens_kappa

raw = pd.read_csv("data/paired_tests.csv")
tab_frame = raw.pivot(index="A", columns="B", values="count")
tab_frame = tab_frame.reindex(index=["positive", "negative"], columns=["positive", "negative"])
tab = tab_frame.to_numpy(dtype=float)
tab
n = tab.sum()                         # 全部58位患者，每人落入一个联合格

# %% 03.2 按四格建立X/Y/C均数
# 理论p=EX、q=EY、r=EC未知；这里以同一批配对病例的样本均数作帽号代入。
p_hat = (tab[0, 0] + tab[0, 1]) / n
q_hat = (tab[0, 0] + tab[1, 0]) / n
r_hat = (tab[0, 0] + tab[1, 1]) / n
means = pd.DataFrame([dict(p_hat=p_hat, q_hat=q_hat, r_hat=r_hat)])
means

# %% 03.3 计算三个系数
# 先用帽号值计算e_hat=p_hat*q_hat+(1-p_hat)(1-q_hat)，再按κ=(r_hat-e_hat)/(1-e_hat)归一化。
e_hat = p_hat * q_hat + (1 - p_hat) * (1 - q_hat)
kappa = (r_hat - e_hat) / (1 - e_hat)
# cr/cp/cq是C/X/Y的一阶系数；它们来自联合四格的同人记录，而非独立边际。
cr = 1 / (1 - e_hat)
cp = -(1 - kappa) * (2 * q_hat - 1) / (1 - e_hat)
cq = -(1 - kappa) * (2 * p_hat - 1) / (1 - e_hat)
coefficients = pd.DataFrame([dict(e_hat=e_hat, kappa=kappa, cr=cr, cp=cp, cq=cq)])
coefficients

# %% 03.4 联合W表和两种矩
# 行是A阳性/阴性，列是B阳性/阴性；W=cr*C+cp*X+cq*Y。
W = np.array([[cr + cp + cq, cp], [cq, cr]])
W
counts = tab
m1 = np.sum(counts * W) / n                 # E(W)：按全部患者n的经验概率矩
m2 = np.sum(counts * W ** 2) / n             # E(W²)：平方在加权前完成
moments = pd.DataFrame([dict(m1=m1, m2=m2)])
moments
# 四格计数让W同时保留X/Y/C的关联；不把三个边际当成独立样本。

# %% 03.5 方差、SE与1.96 Wald区间
# 经验矩方差为E(W²)-[E(W)]²，再除以独立患者数n；不乘n/(n-1)。
variance = (m2 - m1 ** 2) / n
se = np.sqrt(variance)
ci = kappa + np.array([-1, 1]) * 1.96 * se
stepwise_result = pd.DataFrame([dict(p_hat=p_hat, q_hat=q_hat, r_hat=r_hat, e_hat=e_hat, kappa=kappa,
    cr=cr, cp=cp, cq=cq, Wa=W[0, 0], Wb=W[0, 1], Wc=W[1, 0], Wd=W[1, 1],
    m1=m1, m2=m2, variance=variance, se=se, ci_low=ci[0], ci_high=ci[1])])
stepwise_result
print(stepwise_result.to_string(index=False))
# variance≈0.0133022442、SE≈0.1153353554；这是κ的Wald一阶近似。

# %% 03.6 完整函数核对与可选保存
# statsmodels的cohens_kappa接收完整计数表；std_kappa是一般渐近SE。
fit = cohens_kappa(tab.astype(int))
function_check = pd.DataFrame([dict(kappa_manual=kappa, kappa_statsmodels=fit.kappa,
    se_manual=se, se_statsmodels=fit.std_kappa)])
function_check
print(function_check.to_string(index=False))
# 函数核对不改变手算的W表、矩和方差分母；一阶代入不宣称κhat严格无偏。
if os.environ.get("ARTICLE34_SAVE_RESULTS") == "1":
    os.makedirs("results", exist_ok=True)
    stepwise_result.to_csv("results/Python_stepwise.csv", index=False)
