# %% 1 读题与读取共享数据（UTF-8；Spyder从本节按顺序运行）
# 四种孔结构各5个独立试件，每件只测一次有效压缩模量MPa。
# 20数字来自李春喜等《生物统计学》6版例7-1，仅改情境单位，非真实实验。
# 目标：样本η²、总体比例的ω²估计、偏η²分母、R²及A1−A3同时区间。
# 安装/路径见README；默认看图和打印。环境变量SAVE_RESULTS=1时另存CSV/PNG。
from pathlib import Path
import os
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
import statsmodels.api as sm

root = Path.cwd()
if not (root / 'data/scaffold.csv').is_file() and '__file__' in globals():
    root = Path(__file__).resolve().parents[1]
if not (root / 'data/scaffold.csv').is_file():
    raise FileNotFoundError('请将工作目录设为练习项目根目录（含data）。')
save_results = os.getenv('SAVE_RESULTS', '0') == '1'
out = root / '运行结果/Python'
# read_csv输入共享CSV路径，返回DataFrame：group为孔结构，modulus为MPa。
# 随后检查模量缺失并按组拆分，不自动删除试件。
scaffold = pd.read_csv(root / 'data/scaffold.csv')
print('字段：', scaffold.columns.tolist())
print(scaffold.head())
# 检查输入合法性；不固定结果，保留改数据练习的空间。
assert scaffold['modulus'].notna().all()
names = sorted(scaffold['group'].unique())
groups = []
for name in names:
    # loc挑选当前组的行，to_numpy把模量这一列变成数值数组。
    groups.append(scaffold.loc[scaffold.group == name, 'modulus'].to_numpy())
summary = scaffold.groupby('group').modulus.agg(['size', 'mean', 'std'])
print('n、均数、样本SD（n−1分母）：\n', summary)

# %% 2 先手算：平方和、误差均方与效应量
ns = np.array([len(x) for x in groups])
means = np.array([x.mean() for x in groups])
N = len(scaffold)
k = len(groups)
grand = scaffold.modulus.mean()
SSB = np.sum(ns * (means - grand)**2)
SSW = 0.0
for x in groups:
    SSW += np.sum((x - x.mean())**2)
SST = SSB + SSW
MSE = SSW / (N - k)
f_stat = (SSB / (k - 1)) / MSE
# sf计算F分布右尾；参数为观察F、组间df、误差df。
p_value = stats.f.sf(f_stat, k - 1, N - k)
eta2 = SSB / SST
omega2_raw = (SSB - (k - 1)*MSE) / (SST + MSE)
# 报告值采用effectsize相同的非负截断约定，另保留公式原值。
omega2_report = max(0.0, omega2_raw)
partial_eta2 = SSB / (SSB + SSW)
main = dict(N=N, k=k, SSB=SSB, SSW=SSW, SST=SST, MSE=MSE,
            F=f_stat, P=p_value, eta2=eta2, omega2_raw=omega2_raw,
            omega2_report=omega2_report, partial_eta2=partial_eta2)
print(pd.Series(main).to_string(float_format=lambda v: f'{v:.12g}'))
print(f'样本组间平方和占比为{eta2:.2%}；ω²是总体比例的修正点估计。')
print('整体均数检验P<0.05：', p_value < 0.05)

# %% 3 调用模型：同一带截距分类模型的R²
# C(group)明确分类变量，ols默认保留截距；fit返回拟合对象。
fit = smf.ols('modulus ~ C(group)', data=scaffold).fit()
anova_table = sm.stats.anova_lm(fit, typ=1)
print(anova_table)
main['R2'] = fit.rsquared
main['adjusted_R2'] = fit.rsquared_adj
adjusted_manual = 1 - (SSW/(N-k)) / (SST/(N-1))
assert np.isclose(main['R2'], eta2)
assert np.isclose(main['adjusted_R2'], adjusted_manual)
assert np.isclose(anova_table.loc['C(group)', 'sum_sq'], SSB)
print('R²、调整R²：', main['R2'], main['adjusted_R2'])
print('η²等于本模型未调整R²；调整R²的分母与ω²不同。')

# %% 4 A1−A3：带MPa的均数差与全部6对Tukey同时95%区间
# tukey_hsd接收4组原始数，默认共同方差；返回全部组对的差值和调整P矩阵。
# confidence_interval(.95)返回同时95%上下界；[i,j]方向是第i组减第j组。
tukey = stats.tukey_hsd(*groups)
ci = tukey.confidence_interval(confidence_level=0.95)
i = names.index('A1')
j = names.index('A3')
comparison = dict(difference=tukey.statistic[i,j], lower95=ci.low[i,j],
                  upper95=ci.high[i,j], p_adj=tukey.pvalue[i,j])
print('A1−A3（MPa），全部6对的Tukey同时95%区间：\n', pd.Series(comparison))
assert np.isclose(comparison['difference'], means[i] - means[j])

# %% 5 另一个完整小题：只理解2×2分母
# D1/D2孔结构×low/high温度，各格2个独立试件，每件只测一次模量MPa。
# 完整人工数据D1low[9,11] D2low[11,13] D1high[13,15] D2high[15,17]。
# 含孔结构、温度、交互作用的模型；仅求孔结构η²与偏η²，不扩展检验教程。
fac = pd.read_csv(root / 'data/factorial_demo.csv')
print('全部辅助例数据：\n', fac)
fac_fit = smf.ols('modulus ~ C(design) * C(temperature)', data=fac).fit()
# 平衡设计下这里Type I分解正交，模型中的各来源相加为总平方和。
fac_tab = sm.stats.anova_lm(fac_fit, typ=1)
ssa = fac_tab.loc['C(design)', 'sum_sq']
ssb = fac_tab.loc['C(temperature)', 'sum_sq']
ssab = fac_tab.loc['C(design):C(temperature)', 'sum_sq']
sse = fac_tab.loc['Residual', 'sum_sq']
fac_sst = np.sum((fac.modulus - fac.modulus.mean())**2)
factorial = dict(SS_A=ssa, SS_B=ssb, SS_AB=ssab, SSE=sse, SST=fac_sst,
                 A_eta2=ssa/fac_sst, A_partial_eta2=ssa/(ssa+sse))
assert np.isclose(ssa+ssb+ssab+sse, fac_sst)
print('孔结构效应的两种分母：\n', pd.Series(factorial))

# %% 6 看图，再把MPa改为kPa，重新拟合并核对
fig, ax = plt.subplots(figsize=(5.2,4.5))
for index, x in enumerate(groups):
    ax.scatter(index + np.linspace(-.075,.075,len(x)), x, color='#3366cc')
    ax.plot([index-.23,index+.23], [x.mean(),x.mean()], color='black', lw=2)
ax.set(xticks=range(k), xticklabels=names, xlabel='Pore design',
       ylabel='Compressive modulus (MPa)')
fig.tight_layout()
exercise = scaffold.copy()
exercise['modulus'] = exercise.modulus * 1000
ex_fit = smf.ols('modulus ~ C(group)', data=exercise).fit()
ex_tab = sm.stats.anova_lm(ex_fit, typ=1)
ex_b = ex_tab.loc['C(group)', 'sum_sq']
ex_w = ex_tab.loc['Residual', 'sum_sq']
ex_t = ex_b + ex_w
ex_mse = ex_w / (N-k)
exercise_result = dict(eta2=ex_b/ex_t,
    omega2_raw=(ex_b-(k-1)*ex_mse)/(ex_t+ex_mse),
    R2=ex_fit.rsquared, adjusted_R2=ex_fit.rsquared_adj)
ex_groups = [x * 1000 for x in groups]
ex_tukey = stats.tukey_hsd(*ex_groups)
ex_ci = ex_tukey.confidence_interval(confidence_level=.95)
ex_comparison = dict(difference=ex_tukey.statistic[i,j],
    lower95=ex_ci.low[i,j], upper95=ex_ci.high[i,j], p_adj=ex_tukey.pvalue[i,j])
for name in exercise_result:
    assert np.isclose(exercise_result[name], main[name])
for name in ['difference','lower95','upper95']:
    assert np.isclose(ex_comparison[name], 1000*comparison[name])
print('kPa重新计算的效应量及R²：\n', pd.Series(exercise_result))
print('kPa A1−A3差值和同时区间：\n', pd.Series(ex_comparison))
print('单位改变使差值和区间乘1000，无单位的效应量与R²保持不变。')
if save_results:
    out.mkdir(parents=True, exist_ok=True)
    for filename, result in [('main_results',main), ('tukey_A1_A3',comparison),
                             ('factorial_results',factorial),
                             ('exercise_results',exercise_result),
                             ('exercise_tukey_A1_A3',ex_comparison)]:
        pd.DataFrame([result]).to_csv(out / (filename+'.csv'), index=False)
    fig.savefig(out / 'scaffold_groups.png', dpi=240)
plt.show()
