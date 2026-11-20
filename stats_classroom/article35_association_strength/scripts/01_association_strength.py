# 统计课堂35｜关联强度与区间。输入为第11课的假设工程试验汇总计数。
# 300台独立记录器随机分3组各100台：A硅胶、B环氧、C热缩封装。
# 弯折100次后接模拟器记录10min，无持续≥1s数据缺口判合格（课堂阈值）。
# 主表按A/B/C行、合格/不合格列排列为[[90,10],[80,20],[65,35]]。
# 在项目根目录按节顺序运行；默认不写文件，ARTICLE35_SAVE_RESULTS=1才保存。

# %% 01.1 准备与输入
import os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import contingency

raw = pd.read_csv("data/device_pass.csv")
print(raw)
tab = raw[["pass", "fail"]].to_numpy(dtype=int)  # association接口要求整数频数
assert tab.shape == (3, 2)
assert np.all(tab.sum(axis=0) > 0) and np.all(tab.sum(axis=1) > 0)
n = tab.sum()                         # 总样本量是人数/台数之和，非CSV行数
r, s = tab.shape

# %% 01.2 手算预期频数、Q、P、C与V
# axis=1给每组总数；outer把行、列合计的所有组合相乘。
expected = np.outer(tab.sum(axis=1), tab.sum(axis=0)) / n
Q = np.sum((tab - expected) ** 2 / expected)
df = (r - 1) * (s - 1)                # 自由度用于χ²参考分布
m = min(r - 1, s - 1)                 # m用于关联强度的维数归一化
p_value = stats.chi2.sf(Q, df)         # sf为右尾概率
C = np.sqrt(Q / (Q + n))
Cmax = np.sqrt(m / (m + 1))
V = np.sqrt(Q / (n * m))
point = pd.DataFrame([dict(n=n, r=r, s=s, df=df, m=m, Q=Q,
                           p=p_value, C=C, Cmax=Cmax, V=V)])
print(expected)
print(point.to_string(index=False))

# %% 01.3 A/C子表：有符号phi与0/1相关
# Python从0计数，行0是A，行2是C；此子表总台数变为200。
ac = tab[[0, 2], :]
a, b, c, d = ac.ravel()
n_sub = ac.sum()
phi = (a*d - b*c) / np.sqrt((a+b)*(c+d)*(a+c)*(b+d))
# 按已有格子计数展开类别标签，只用于核对0/1 Pearson相关。
repeats = ac.ravel().astype(int)
x = np.repeat([1, 1, 0, 0], repeats)  # A=1，C=0
y = np.repeat([1, 0, 1, 0], repeats)  # 合格=1，不合格=0
phi_pearson = np.corrcoef(x, y)[0, 1]
flipped = ac[::-1, :]                 # 实际交换两行并重新计算符号
fa, fb, fc, fd = flipped.ravel()
phi_flipped = (fa*fd-fb*fc) / np.sqrt(np.prod(flipped.sum(0))*np.prod(flipped.sum(1)))
Q_sub = stats.chi2_contingency(ac, correction=False).statistic
V_sub = contingency.association(ac, method="cramer", correction=False)
phi_result = pd.DataFrame([dict(a=a,b=b,c=c,d=d,n_sub=n_sub,Q_sub=Q_sub,
                                phi=phi,phi_pearson=phi_pearson,
                                phi_row_flipped=phi_flipped,V_sub=V_sub)])
assert np.isclose(phi, phi_pearson) and np.isclose(phi_flipped, -phi)
assert np.isclose(V_sub, abs(phi)) and np.isclose(phi**2, Q_sub/n_sub)
print(phi_result.to_string(index=False))

# %% 01.4 相同比例与不同样本量的假设试验
# 分别重算三张人数表；这些是独立设备试验的假设规模比较。
scaling_rows = []
for k in [1, 2, 4]:
    enlarged = k * tab
    nn = enlarged.sum()
    ee = np.outer(enlarged.sum(1), enlarged.sum(0)) / nn
    qq = np.sum((enlarged-ee)**2/ee)
    scaling_rows.append(dict(multiplier=k,n=nn,Q=qq,p=stats.chi2.sf(qq,df),
                             Q_over_n=qq/nn,C=np.sqrt(qq/(qq+nn)),V=np.sqrt(qq/(nn*m))))
scaling = pd.DataFrame(scaling_rows)
assert np.allclose(scaling["V"], V)
print(scaling.to_string(index=False))

# %% 01.5 完整软件核对与可选保存
# 这里真实读取SciPy输出，独立核对前面的手算量。
fit = stats.chi2_contingency(tab, correction=False)
fit_C = contingency.association(tab, method="pearson", correction=False)
fit_V = contingency.association(tab, method="cramer", correction=False)
function_comparison = pd.DataFrame(dict(quantity=["Q","P","C","V"],
    manual=[Q,p_value,C,V], software=[fit.statistic,fit.pvalue,fit_C,fit_V]))
assert np.allclose(function_comparison.manual, function_comparison.software, atol=1e-12)
print(function_comparison.to_string(index=False))
if os.environ.get("ARTICLE35_SAVE_RESULTS") == "1":
    Path("results").mkdir(exist_ok=True)
    point.to_csv("results/Python_point.csv", index=False)
    phi_result.to_csv("results/Python_phi.csv", index=False)
    scaling.to_csv("results/Python_scaling.csv", index=False)
