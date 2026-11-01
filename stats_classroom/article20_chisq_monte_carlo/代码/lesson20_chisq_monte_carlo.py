"""统计课堂20：固定边际Monte Carlo。输入为课堂249人汇总频数。
Spyder以配套根目录打开现有文件夹项目，打开本文件并F5完整运行。
依赖numpy、pandas、scipy、matplotlib；结果在运行结果/Python。
"""
# %% 一、读数据，计算六格贡献
from pathlib import Path
import runpy
import numpy as np
import pandas as pd
from scipy.special import gammaln
from scipy.stats import chi2, ncx2
from scipy.optimize import brentq
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / '运行结果/Python'
OUT.mkdir(parents=True, exist_ok=True)
data = pd.read_csv(ROOT / '数据/ace_genotype.csv')
print(data)
O = data[['DD', 'ID', 'II']].to_numpy(dtype=int)
rows = O.sum(axis=1)
cols = O.sum(axis=0)
N = int(O.sum())
E = np.outer(rows, cols) / N
contribution = (O - E)**2 / E
Q = float(contribution.sum())
asymptoticP = float(chi2.sf(Q, df=2))  # sf是右尾概率，不再乘2。
threshold = Q * (1 - 64 * np.finfo(float).eps)
print('期望表：\n', E)
print('六格贡献：\n', contribution)
pd.DataFrame({'group': np.repeat(data.group, 3),
              'genotype': ['DD', 'ID', 'II'] * 2,
              'observed': O.ravel(), 'expected': E.ravel(),
              'pearson_contribution': contribution.ravel()}).to_csv(OUT/'cells.csv', index=False)

# %% 二、把249个标签无放回抽111个，先看一张表
labels = np.repeat(np.arange(3), cols)
rng = np.random.default_rng(123)
# choice抽标签所在位置；replace=False保证一轮内每张卡片最多抽一次。
picked = rng.choice(N, size=rows[0], replace=False)
first_row = np.bincount(labels[picked], minlength=3)
first = np.vstack([first_row, cols-first_row])
print('单张模拟表：\n', first)
print('单张χ²：', ((first-E)**2/E).sum())
one = pd.DataFrame(first, columns=['DD', 'ID', 'II'])
one.insert(0, 'group', data.group)
one.to_csv(OUT/'one_table.csv', index=False)

# %% 三、多元超几何抽样批量生成第一行，第二行由列合计得到
B = 100000
rng = np.random.default_rng(123)
# colors=三种基因型总数；nsample=第一组111人；size=B次重复。
first_rows = rng.multivariate_hypergeometric(cols, rows[0], size=B)
second_rows = cols - first_rows
assert np.all(first_rows.sum(axis=1) == rows[0])
assert np.all(second_rows.sum(axis=1) == rows[1])
assert np.all(first_rows + second_rows == cols)
assert np.all(first_rows >= 0) and np.all(second_rows >= 0)
assert np.all(first_rows == np.floor(first_rows))
first_part = ((first_rows-E[0])**2/E[0]).sum(axis=1)
second_part = ((second_rows-E[1])**2/E[1]).sum(axis=1)
qs = first_part + second_part
b = int(np.count_nonzero(qs >= threshold))
P = (b+1)/(B+1)
SE = float(np.sqrt(P*(1-P)/B))
pd.DataFrame({'Q': qs}).to_csv(OUT/'simulation.csv', index=False)
print('极端次数、模拟P、模拟标准误：', b, P, SE)

# %% 四、穷举所有可能表，按条件概率累加Pearson右尾
# gammaln(x+1)=log(x!)，避免直接计算巨大的阶乘。
log_denominator = gammaln(N+1)-gammaln(rows[0]+1)-gammaln(rows[1]+1)
records = []
for dd in range(int(cols[0])+1):
    for id_count in range(int(cols[1])+1):
        ii = int(rows[0])-dd-id_count
        if 0 <= ii <= cols[2]:
            a = np.array([dd, id_count, ii])
            current = np.vstack([a, cols-a])
            current_q = float(((current-E)**2/E).sum())
            log_numerator = np.sum(gammaln(cols+1)-gammaln(a+1)-gammaln(cols-a+1))
            probability = float(np.exp(log_numerator-log_denominator))
            records.append([dd, id_count, ii, current_q, probability])
en = pd.DataFrame(records, columns=['DD', 'ID', 'II', 'Q', 'probability'])
en['pearson_tail'] = en.Q >= threshold
exact = float(en.loc[en.pearson_tail, 'probability'].sum())
assert len(en) == 4063 and abs(en.probability.sum()-1) < 1e-10
assert abs(P-exact) < 5*SE
assert np.all(qs >= 0)
en.to_csv(OUT/'enumeration.csv', index=False)
print('Pearson精确条件尾概率：', exact)

# %% 五、效应量及近似95%区间（非中心χ²反演）
# 分别求使观察Q处于累计概率.975、.025的非中心参数，再转到V尺度。
effect_v = float(np.sqrt(Q/N))
low_ncp = brentq(lambda value: ncx2.cdf(Q, 2, value)-.975, 0, 100, xtol=1e-12)
high_ncp = brentq(lambda value: ncx2.cdf(Q, 2, value)-.025, 0, 100, xtol=1e-12)
ci_low = float(np.sqrt(low_ncp/N))
ci_high = float(np.sqrt(high_ncp/N))
summary = pd.DataFrame([{'N': N, 'Qobs': Q, 'df': 2, 'asymptoticP': asymptoticP,
                         'B': B, 'b': b, 'P': P, 'SE': SE, 'exact_pearson_P': exact,
                         'effect_v': effect_v, 'ci_low': ci_low, 'ci_high': ci_high,
                         'builtinP': np.nan}])
summary.to_csv(OUT/'summary.csv', index=False)
print(summary)

# %% 六、同一序列的不同前缀：计算误差随B怎样变化？
exercise = []
sizes = sorted(set(size for size in [2000, 10000, B] if size <= B))
for size in sizes:
    count = int(np.count_nonzero(qs[:size] >= threshold))
    estimate = (count+1)/(size+1)
    error = float(np.sqrt(estimate*(1-estimate)/size))
    exercise.append([size, count, estimate, error])
pd.DataFrame(exercise, columns=['B','b','P','SE']).to_csv(OUT/'precision_exercise.csv', index=False)
print('改变种子可观察计算波动；不同语言同种子不保证相同随机序列。')
print('PYTHON_CHECKS_PASSED：边际、整数、枚举总和与Monte Carlo误差。')

# %% 七、绘图：读取本语言已保存的模拟序列，保留坐标CSV
runpy.run_path(str(ROOT/'代码/plot_chisq_monte_carlo.py'))
