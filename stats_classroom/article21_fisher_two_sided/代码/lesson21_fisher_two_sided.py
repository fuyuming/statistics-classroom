"""统计课堂21：Fisher双侧P值。课堂例7-4，33名孕妇的新生儿感染结局。
Spyder打开本文件后F5全文运行，或依次运行六个单元。
依赖numpy、pandas、scipy；末尾自动绘图另用matplotlib。
"""
from pathlib import Path
import runpy
import numpy as np
import pandas as pd
from scipy import stats
from scipy.special import comb

# %% 一、读题和计数检查，先用标准函数
ROOT = Path(__file__).resolve().parents[1]
data = pd.read_csv(ROOT / '数据/hbv_prevention.csv')
print(data)
counts = data[['positive', 'negative']].to_numpy(dtype=float)
if counts.shape != (2, 2) or not np.isfinite(counts).all():
    raise ValueError('需要没有缺失的2×2计数表。')
if (counts < 0).any() or not (counts == np.floor(counts)).all():
    raise ValueError('计数须为非负整数。')
O = counts.astype(int)
rows = O.sum(axis=1)  # axis=1逐行相加。
cols = O.sum(axis=0)  # axis=0逐列相加。
if (rows == 0).any() or (cols == 0).any():
    raise ValueError('每个行列合计须大于0。')
N = int(O.sum())
a = int(O[0, 0])  # Python第一个元素用索引0。
fit = stats.fisher_exact(O, alternative='two-sided')
print('标准函数的双侧P：', fit.pvalue)

# %% 二、逐表计算概率；先标准分布函数，再组合数核对
x = np.arange(max(0, rows[0]-cols[1]), min(rows[0], cols[0])+1)
# SciPy依次传入x、总人数、总阳性数、第一组人数。
q = stats.hypergeom.pmf(x, N, cols[0], rows[0])
q_obs = q[x == a][0]
q_formula = comb(cols[0], x) * comb(cols[1], rows[0]-x)
q_formula = q_formula / comb(N, rows[0])
assert np.allclose(q, q_formula, rtol=1e-12, atol=1e-14)
assert abs(q.sum()-1) < 1e-12
print(pd.DataFrame({'X':x, 'probability':q}))

# %% 三、单侧求和与两种双侧定义
included = q <= q_obs*(1+1e-7)
p_manual = float(q[included].sum())
p_left = float(q[x <= a].sum())
p_right = float(q[x >= a].sum())
p_double = min(1.0, 2*min(p_left, p_right))
assert abs(p_manual-fit.pvalue) < 1e-12
assert abs(p_left-stats.fisher_exact(O, alternative='less').pvalue) < 1e-12
assert abs(p_right-stats.fisher_exact(O, alternative='greater').pvalue) < 1e-12
print('计入双侧的X：', x[included])
print('左单侧、右单侧、双倍较小单侧、概率排序：')
print(p_left, p_right, p_double, p_manual)

# %% 四、OR与等尾95%区间；显式选条件估计
# fisher_exact.statistic给样本OR；conditional给与R相同定义的条件MLE。
conditional = stats.contingency.odds_ratio(O, kind='conditional')
interval = conditional.confidence_interval(confidence_level=0.95)
rates = O[:, 0]/rows
print('两组观察率：', rates)
print('样本OR、条件OR、区间：', fit.statistic, conditional.statistic, interval)
# 该区间是等尾区间，对应central规则；P值采用概率排序规则。

# %% 五、改变比较方向，观察双侧P与效应估计的变化

def analyze_variant(current, label):
    result = stats.fisher_exact(current, alternative='two-sided')
    left = stats.fisher_exact(current, alternative='less').pvalue
    right = stats.fisher_exact(current, alternative='greater').pvalue
    odds = stats.contingency.odds_ratio(current, kind='conditional')
    limits = odds.confidence_interval(confidence_level=0.95)
    return {
        'variant':label, 'fisher_P':result.pvalue,
        'left_P':left, 'right_P':right,
        'OR_sample':result.statistic, 'OR_conditional':odds.statistic,
        'OR_low':limits.low, 'OR_high':limits.high
    }

original = analyze_variant(O, 'original')
swap_rows = analyze_variant(O[::-1, :], 'swap_rows')
swap_columns = analyze_variant(O[:, ::-1], 'swap_columns')
exercise = pd.DataFrame([original, swap_rows, swap_columns])
assert np.allclose(exercise['fisher_P'], p_manual, rtol=1e-12)
print(exercise)
print('交换行或列后双侧P不变，左右单侧互换；OR的比较方向反转。')

# %% 六、保存三张结果表，并复现科学图
OUT = ROOT/'运行结果/Python'
OUT.mkdir(parents=True, exist_ok=True)
enumeration = pd.DataFrame({
    'X':x, 'prevention_positive':x, 'prevention_negative':rows[0]-x,
    'control_positive':cols[0]-x, 'control_negative':rows[1]-cols[0]+x,
    'probability':q, 'included':included
})
summary = pd.DataFrame([{
    'N':N, 'row1':rows[0], 'row2':rows[1], 'col1':cols[0], 'col2':cols[1],
    'observed_X':a, 'q_observed':q_obs,
    'left_P':p_left, 'right_P':p_right, 'double_min_P':p_double,
    'fisher_P':p_manual, 'OR_sample':fit.statistic,
    'OR_conditional':conditional.statistic, 'OR_low':interval.low,
    'OR_high':interval.high, 'rate1':rates[0], 'rate2':rates[1]
}])
enumeration.to_csv(OUT/'enumeration.csv', index=False)
summary.to_csv(OUT/'summary.csv', index=False)
exercise.to_csv(OUT/'swap_exercise.csv', index=False)
runpy.run_path(str(ROOT/'代码/plot_fisher_two_sided.py'), run_name='__main__')
