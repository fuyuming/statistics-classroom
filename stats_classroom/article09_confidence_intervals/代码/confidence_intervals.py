# -*- coding: utf-8 -*-
"""统计课堂09：五个置信区间与覆盖率模拟。
Spyder：打开此文件，按F5完整运行。不要仅选择几行运行。
依赖 numpy、pandas、scipy、matplotlib；数据来自项目的“数据”文件夹。
Console显示步骤，Plots显示图，运行结果/Python保存CSV与PNG。
图中用英文以避免不同计算机的中文字体缺失；中文解释见输出和README。
"""
# %% 准备：读取共享数据，先认识字段，再进行计算
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

root = Path(__file__).resolve().parents[1]
out = root / '运行结果' / 'Python'
out.mkdir(parents=True, exist_ok=True)
inputs = {}
for file in sorted((root / '数据').glob('*.csv')):
    inputs[file.stem] = pd.read_csv(file)
    print('\n数据：', file.name)
    print(inputs[file.stem].to_string(index=False))

confidence = 0.95
alpha = 1 - confidence
# 双侧95%区间，两侧各留0.025；ppf(0.975)返回累计概率0.975对应的分位点。
upper_prob = 1 - alpha / 2
zcrit = stats.norm.ppf(upper_prob)
results = []

def report(example, estimate, lower, upper, unit, explanation):
    """保存本题数表和区间图，主程序仍按题目顺序计算。"""
    record = dict(example=example, estimate=estimate, lower=lower, upper=upper, unit=unit)
    results.append(record)
    pd.DataFrame([record]).to_csv(out / f'{example}_interval.csv', index=False)
    print(f'估计值={estimate:.8f}；95%区间=[{lower:.8f}, {upper:.8f}] {unit}')
    print(explanation)
    fig, ax = plt.subplots(figsize=(8, 2.8))
    ax.plot([lower, upper], [1, 1], lw=3, color='#1769AA')
    ax.plot(estimate, 1, 'o', color='#D65A4A')
    ax.set(yticks=[], xlabel=unit, title=f'Example {example}: 95% confidence interval', ylim=(.6, 1.4))
    ax.text(lower, .87, f'{lower:.4f}', ha='center')
    ax.text(upper, .87, f'{upper:.4f}', ha='center')
    fig.tight_layout()
    fig.savefig(out / f'{example}_interval.png', dpi=160)

# %% 例题1：A、B两批幼崽的总体平均出生重相差多少？
# 条件：两组独立，分别来自正态总体，两个总体方差相同且未知。
# 估计目标：mu_A-mu_B；使用合并方差t区间，df=n1+n2-2。
data = inputs['birth_weight_g']
a = data.loc[data['group'] == 'A', 'weight_g'].to_numpy()
b = data.loc[data['group'] == 'B', 'weight_g'].to_numpy()
n1 = len(a)
n2 = len(b)
# ddof=1使方差分母为n-1，与R的var一致。
s1_sq = np.var(a, ddof=1)
s2_sq = np.var(b, ddof=1)
difference = np.mean(a) - np.mean(b)
df = n1 + n2 - 2
pooled_variance = ((n1 - 1)*s1_sq + (n2 - 1)*s2_sq) / df
se = np.sqrt(pooled_variance*(1/n1 + 1/n2))
tcrit = stats.t.ppf(upper_prob, df)
print('\n例题1：均数差。n1,n2,方差1,方差2,合并方差,SE,df,t分位点：')
print(n1, n2, s1_sq, s2_sq, pooled_variance, se, df, tcrit)
report('01', difference, difference-tcrit*se, difference+tcrit*se, 'g',
       '区间约[-2.00,6.10]g，包含0，不能据此确定总体均数差的方向。')

# %% 例题2：总体客观缓解率是多少？
# 独立二项模型。理论上检查np和n(1-p)；p未知，实践中用观察次数作经验检查。
data = inputs['response_counts'].iloc[0]
x = int(data['events'])
n = int(data['n'])
p_hat = x / n
print('\n例题2：成功、失败次数：', x, n-x, '；是否均>=5：', min(x,n-x)>=5)
se = np.sqrt(p_hat*(1-p_hat)/n)
lower = p_hat - zcrit*se
upper = p_hat + zcrit*se
report('02', p_hat, lower, upper, 'proportion',
       'Wald近似区间约35.25%—54.75%；这里的5是经验阈值。')
# 额外比较Wilson与Clopper-Pearson区间。beta分位点用于精确二项区间反演。
denominator = 1 + zcrit**2/n
center = (p_hat + zcrit**2/(2*n))/denominator
half = zcrit*np.sqrt(p_hat*(1-p_hat)/n + zcrit**2/(4*n*n))/denominator
exact_lower = stats.beta.ppf(alpha/2, x, n-x+1) if x else 0
exact_upper = stats.beta.ppf(upper_prob, x+1, n-x) if x<n else 1
pd.DataFrame({'method':['Wald','Wilson','Clopper-Pearson'],
              'lower':[lower,center-half,exact_lower],
              'upper':[upper,center+half,exact_upper]}).to_csv(out/'02_methods.csv',index=False)

# %% 例题3：低洼地与高坡地的总体麦田发病率相差多少？
# 两组独立二项样本；各组成功、失败次数均需检查。区间用未合并标准误。
data = inputs['wheat_rust_counts']
x1, n1 = data.loc[data['group']=='lowland',['events','n']].iloc[0]
x2, n2 = data.loc[data['group']=='upland',['events','n']].iloc[0]
p1 = x1/n1
p2 = x2/n2
difference = p1-p2
se = np.sqrt(p1*(1-p1)/n1 + p2*(1-p2)/n2)
print('\n例题3：四个次数：', x1,n1-x1,x2,n2-x2)
print('四项是否均>=5：', min(x1,n1-x1,x2,n2-x2)>=5, 'SE=',se)
report('03', difference, difference-zcrit*se, difference+zcrit*se, 'proportion difference',
       '低洼地比高坡地高约11.44个百分点；区间约6.45—16.42个百分点。')

# %% 例题4：土壤铅浓度的总体方差是多少？先检查模型
# 独立性由采样设计判断；正态性结合原始点图、Q-Q图和Shapiro-Wilk检验。
x = inputs['lead_ug_g']['lead_ug_g'].to_numpy()
normality = stats.shapiro(x)
print('\n例题4：Shapiro-Wilk W、P：', normality.statistic, normality.pvalue)
print('P约0.746，未发现明显偏离正态的证据；只有8个观测，不能据此证明总体正态。')
pd.DataFrame([{'W':normality.statistic,'p_value':normality.pvalue}]).to_csv(out/'04_normality.csv',index=False)
fig, axes = plt.subplots(1,2,figsize=(9,3.5))
axes[0].plot(x,np.zeros(len(x)),'o')
axes[0].set(xlabel='Lead (ug/g)',yticks=[],title='Observed values')
# 概率绘图位置；直线穿过样本四分位数，便于和R的qqline比较。
positions = (np.arange(1,len(x)+1)-.375)/(len(x)+.25)
normal_q = stats.norm.ppf(positions)
slope = np.diff(np.quantile(x,[.25,.75]))[0]/np.diff(stats.norm.ppf([.25,.75]))[0]
intercept = np.quantile(x,.25)-slope*stats.norm.ppf(.25)
axes[1].plot(normal_q,np.sort(x),'o')
axes[1].plot(normal_q,intercept+slope*normal_q)
axes[1].set(xlabel='Theoretical normal quantile',ylabel='Lead (ug/g)',title='Normal Q-Q plot')
fig.tight_layout()
fig.savefig(out/'04_normality.png',dpi=160)
n = len(x)
s_sq = np.var(x,ddof=1)
df = n-1
# W=(n-1)S²/sigma² ~ chi-square(df)。倒数反演后，大分位点放下限分母。
q_low = stats.chi2.ppf(alpha/2,df)
q_high = stats.chi2.ppf(upper_prob,df)
lower = df*s_sq/q_high
upper = df*s_sq/q_low
print('n、样本方差、df、卡方分位点：',n,s_sq,df,q_low,q_high)
report('04',s_sq,lower,upper,'(ug/g)^2',
       '正态总体模型下的精确方差区间；开平方得到总体标准差区间。')
print('标准差区间：',np.sqrt(lower),np.sqrt(upper),'ug/g')

# %% 例题5：两个总体方差相差多少倍？
# 两组独立正态样本；不预先要求方差相等。原教材仅提供汇总量，不虚构原始数据。
data = inputs['wheat_variance_summary']
n1 = int(data.loc[0,'n'])
n2 = int(data.loc[1,'n'])
s1_sq = data.loc[0,'sample_variance']
s2_sq = data.loc[1,'sample_variance']
ratio = s1_sq/s2_sq
df1 = n1-1
df2 = n2-1
# (S1²/S2²)/(sigma1²/sigma2²) ~ F(df1,df2)。ratio是估计值，theta是总体参数。
f_low = stats.f.ppf(alpha/2,df1,df2)
f_high = stats.f.ppf(upper_prob,df1,df2)
print('\n例题5：样本方差比、自由度、F分位点：',ratio,df1,df2,f_low,f_high)
report('05',ratio,ratio/f_high,ratio/f_low,'variance ratio',
       '总体方差比区间约[1.94,31.48]，不含1；对应等尾双侧F检验在0.05水平拒绝方差相等。')
pd.DataFrame(results).to_csv(out/'all_intervals.csv',index=False)

# %% 重复抽样：95%置信水平描述的是构造区间的方法
rng = np.random.default_rng(20261012)
rounds = 1000
n1 = n2 = 10
mu1, mu2, sigma = 142,140,2.5
truth = mu1-mu2
tcrit = stats.t.ppf(upper_prob,n1+n2-2)
records = []
for i in range(rounds):
    a = rng.normal(mu1,sigma,n1)
    b = rng.normal(mu2,sigma,n2)
    estimate = a.mean()-b.mean()
    pooled = ((n1-1)*a.var(ddof=1)+(n2-1)*b.var(ddof=1))/(n1+n2-2)
    se = np.sqrt(pooled*(1/n1+1/n2))
    lower = estimate-tcrit*se
    upper = estimate+tcrit*se
    records.append([i+1,estimate,lower,upper,lower<=truth<=upper])
simulation = pd.DataFrame(records,columns=['round','estimate','lower','upper','covered'])
simulation.to_csv(out/'06_coverage.csv',index=False)
print('\n1000轮覆盖比例：',simulation['covered'].mean(),'；固定真值：',truth)
print('相同种子不保证跨语言抽到同一批随机数，因此三语言覆盖比例不必完全一致。')
fig,ax = plt.subplots(figsize=(7,8))
for row in simulation.head(40).itertuples():
    ax.plot([row.lower,row.upper],[row.round,row.round],color='#1769AA' if row.covered else '#D65A4A')
ax.axvline(truth,color='#D65A4A',ls='--')
ax.set(xlabel='Population mean difference (g)',ylabel='Round',title='95% CIs from 40 repeated samples')
ax.invert_yaxis()
fig.tight_layout()
fig.savefig(out/'06_coverage.png',dpi=160)
print('练习：把confidence改为0.99，完整重跑，观察区间宽度与覆盖比例。')
print('结果已保存到项目的 运行结果/Python。')
plt.show()
