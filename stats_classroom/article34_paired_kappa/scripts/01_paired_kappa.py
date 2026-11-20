# 统计课堂34：配对率差与Cohen κ
# 58名可疑系统性红斑狼疮患者，均接受免疫荧光A和乳胶B的抗核抗体检测。
# 数据为颜艳、王彤《医学统计学》第5版例7-3的联合汇总计数。
# 问题：B−A率差、Po/Pe/κ、总体κ的95%区间、κ>0检验、固定边际可达范围。
# 条件：不同患者联合结果独立；每人两法结果保留配对。工作目录设为项目根目录。
# 第11课PDF P16、P19–25为课堂入口；本文另展开一般区间与固定边际界限。
# 使用Spyder按# %%顺序执行；SciPy/statsmodels均为第三方库，依赖见README。
# %% 问题与输入
# A=免疫荧光，B=乳胶；当前工作目录应为项目根目录，不依赖__file__。
import os
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.inter_rater import cohens_kappa
from statsmodels.stats.contingency_tables import mcnemar

raw = pd.read_csv("data/paired_tests.csv")
# CSV每行是一个联合格及其count，不能把四行误当作四位对象；按标签放回计数矩阵。
tab_frame = raw.pivot(index="A", columns="B", values="count")
tab_frame = tab_frame.reindex(index=["positive", "negative"], columns=["positive", "negative"])
tab = tab_frame.to_numpy(dtype=float)
print(tab)
n = tab.sum()
a, b, c, d = tab.ravel()

# %% 手动计算：率差、观察一致率、机会一致率与κ
# p是四格的经验联合概率；r、col是两个方法的观察边际概率。
p = tab / n
r = p.sum(axis=1)
col = p.sum(axis=0)
pA = r[0]
pB = col[0]
rd = pB - pA
po = np.trace(p)
pe = np.sum(r * col)
# κ以机会一致为基准作归一化，无量纲；rd乘100后为百分点。
kappa = (po - pe) / (1 - pe)
point_estimates = pd.DataFrame([dict(n=n, pA=pA, pB=pB, rd=rd, po=po, pe=pe, kappa=kappa)])
print(point_estimates)
# rd≈−.1724、Po≈.7586、κ≈.4550。

# %% 手动计算：κ=0对预先指定κ>0的右尾近似检验
# SE0只对应κ=0的Rosner零假设；Wald区间使用后面普通多项delta SE。
se0 = np.sqrt((pe + pe**2 - np.sum(r * col * (r + col))) / (n * (1 - pe)**2))
z = kappa / se0
# 事先方向为H1: kappa>0，因此报告标准正态右尾概率。
p_greater = stats.norm.sf(z)
null_test = pd.DataFrame([dict(se0=se0, z=z, p_greater=p_greater)])
print(null_test)
# SE0≈.120950、Z≈3.7622、右尾P≈.00008422。

# %% 手动计算：一般多项delta标准误与95%Wald区间
# g表示每个单元概率轻微改变时κ的敏感度；减去加权均值平方是多项约束的协方差项。
g = np.zeros_like(p)
for i in range(2):
    for j in range(2):
        g[i, j] = ((i == j) * (1 - pe) - (1 - po) * (r[j] + col[i])) / (1 - pe)**2
g_mean = np.sum(p * g)
g_second = np.sum(p * g**2)
variance = (g_second - g_mean**2) / n
se = np.sqrt(variance)
# 双侧95%两端各留2.5%；ppf由累计概率返回标准正态分位点。
ci = kappa + stats.norm.ppf([.025, .975]) * se
interval_result = pd.DataFrame([dict(g_mean=g_mean, g_second=g_second, variance=variance, se=se, ci_low=ci[0], ci_high=ci[1])])
print(g)
print(interval_result)
# 一般SE≈.115335；95%Wald区间≈[.228980,.681087]。

# %% 手动计算：McNemar条件二项双侧精确P
m = b + c  # 14个不一致对子；率差的分母仍是全部n=58
mcnemar_exact_p = min(1.0, 2 * stats.binom.cdf(min(b, c), int(m), .5))

# %% 调用函数：statsmodels的完整精确McNemar接口
exact_fit = mcnemar(tab.astype(int), exact=True, correction=False)
print(exact_fit)
# exact=True使用条件二项精确分布；本设定不使用卡方连续性校正。

# %% 结果核对：手算与完整函数的双侧精确P
exact_comparison = pd.DataFrame([dict(manual=mcnemar_exact_p, function_P=exact_fit.pvalue)])
print(exact_comparison)
# 两个P应均≈.012939；该接口不返回全体患者率差的置信区间。

# %% 手动枚举：固定观察边际，改变共同阳性人数
a_pos = tab[0, :].sum()
b_pos = tab[:, 0].sum()
a_min = int(max(0, a_pos + b_pos - n))
a_max = int(min(a_pos, b_pos))
rows = []
for aa in range(a_min, a_max + 1):
    # 每一轮构造满足同一行列合计的一张表；不是新增临床观察。
    bb = a_pos - aa
    cc = b_pos - aa
    dd = n - aa - bb - cc
    tt = np.array([[aa, bb], [cc, dd]], dtype=float)
    pp = tt / tt.sum()
    rr = pp.sum(axis=1)
    col_enum = pp.sum(axis=0)
    poo = np.trace(pp)
    pee = np.sum(rr * col_enum)
    rows.append(dict(a=aa, b=bb, c=cc, d=dd, rd=col_enum[0]-rr[0], po=poo, kappa=(poo-pee)/(1-pee)))
bounds = pd.DataFrame(rows)
print(bounds)
# κ范围约[−.4013,.6107]；这一枚举固定当前边际，和总体多项Wald区间的条件不同。
summary = pd.DataFrame([dict(n=n, pA=pA, pB=pB, rd=rd, po=po, pe=pe, kappa=kappa,
    se0=se0, z=z, p_greater=p_greater, se=se, ci_low=ci[0], ci_high=ci[1],
    mcnemar_exact_p=mcnemar_exact_p, kappa_min=bounds.kappa.min(),
    kappa_max=bounds.kappa.max(), po_min=bounds.po.min(), po_max=bounds.po.max())])
print(summary.to_string(index=False, float_format=lambda x: f"{x:.12g}"))

# %% 调用函数：statsmodels的完整Cohen κ接口
fit = cohens_kappa(tab.astype(int))
print(fit)
# std_kappa对应一般SE，std_kappa0对应零假设SE0；默认还列双侧P，按预设方向取单侧字段。

# %% 结果核对：κ、两类SE、区间与单侧P
function_comparison = pd.DataFrame({
    "quantity": ["kappa", "general_SE", "null_SE", "CI_low", "CI_high", "P_greater"],
    "manual": [kappa, se, se0, ci[0], ci[1], p_greater],
    "statsmodels": [fit.kappa, fit.std_kappa, fit.std_kappa0, fit.kappa_low, fit.kappa_upp, fit.pvalue_one_sided],
})
print(function_comparison)

# %% 思考：固定边际，把共同阳性a从11改为13
# 先预测B−A率差、Po、κ各怎样变化，再从bounds读取对应行。
# κ≈.455按六级参考为moderate，按Rosner三级为good；两个经验名称不改变估计值与区间。

# %% 可选保存：首次练习默认不写文件
# 若需保存，先在Spyder Console运行os.environ["ARTICLE34_SAVE_RESULTS"]="1"，再运行本节。
if os.environ.get("ARTICLE34_SAVE_RESULTS") == "1":
    os.makedirs("results", exist_ok=True)
    summary.to_csv("results/Python_summary.csv", index=False)
    bounds.to_csv("results/Python_bounds.csv", index=False)
