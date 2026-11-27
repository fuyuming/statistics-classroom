# %% 1. 设置与共享数据（UTF-8）
# 统计课堂46：独立、正态、共同方差的四组均数比较。
# 数字来自李春喜等《生物统计学》6版例7-1；支架情境为教学改编。
# Spyder先将项目根设为工作目录，按节运行；F5运行全文件。
import os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt

ALPHA = .05
SAVE_OUTPUTS = False  # 改True才保存CSV与图片；默认在IDE中显示。
UNIT_SCALE = 1        # 改1000练习MPa转kPa。
SAVE_OUTPUTS = SAVE_OUTPUTS or os.environ.get('SAVE_OUTPUTS') == '1'
UNIT_SCALE = float(os.environ.get('UNIT_SCALE', UNIT_SCALE))
unit = 'kPa' if UNIT_SCALE == 1000 else 'MPa'
root = Path.cwd()
if not (root / 'data').exists() and '__file__' in globals():
    root = Path(__file__).resolve().parents[1]
if not (root / 'data').exists():
    raise FileNotFoundError('请把Spyder工作目录设为含data的项目根。')
dat = pd.read_csv(root / 'data' / 'scaffold.csv')
dat['modulus'] = dat['modulus'] * UNIT_SCALE
print(dat.head())
print(dat.dtypes)  # group为组名，modulus为每个独立试件的一次测量。

# %% 2. 整体模型与全部四组的共同MSE
group_names = ['A1', 'A2', 'A3', 'A4']
groups = []
for name in group_names:
    values = dat.loc[dat['group'] == name, 'modulus'].to_numpy()
    groups.append(values)
means = np.array([values.mean() for values in groups])
n = np.array([len(values) for values in groups])
k = len(groups)
df_error = len(dat) - k
ss_within = 0.0
for values in groups:
    ss_within += np.sum((values - values.mean()) ** 2)
mse = ss_within / df_error
# f_oneway是经典共同方差ANOVA，返回整体F和右尾P。
f_result = stats.f_oneway(*groups)
summary = pd.DataFrame([{'MSE': mse, 'df': df_error,
                        'F': f_result.statistic, 'F_p': f_result.pvalue}])
print(pd.DataFrame({'group': group_names, 'n': n, 'mean': means}))
print(summary.to_string(index=False))
print('整体F通过保护性LSD前置关：', f_result.pvalue < ALPHA)
print('下面LSD仍列普通P，不能作为全部组对的强FWER调整结果。')

# %% 3. LSD形式：普通双侧t、标准误与普通区间
lsd_rows = []
for i in range(k - 1):
    for j in range(i + 1, k):
        difference = means[i] - means[j]
        se = np.sqrt(mse * (1 / n[i] + 1 / n[j]))
        t_value = difference / se
        # sf计算右尾，乘2为双侧；ppf(.975,df)返回双侧5%正临界值。
        p = 2 * stats.t.sf(abs(t_value), df_error)
        half_width = stats.t.ppf(1 - ALPHA / 2, df_error) * se
        lsd_rows.append([group_names[i], group_names[j], difference, se,
                         t_value, df_error, p, difference - half_width,
                         difference + half_width, p < ALPHA])
lsd = pd.DataFrame(lsd_rows, columns=['group1', 'group2', 'difference', 'SE',
    't', 'df', 'p', 'CI_low', 'CI_high', 'reject'])
print(lsd.to_string(index=False))
print('普通P拒绝数：', lsd['reject'].sum(), '；差值和区间单位：', unit)

# %% 4. Tukey：全部六对同时推断
# SciPy返回完整矩阵，第[i,j]项方向为第i组减第j组。
tukey_result = stats.tukey_hsd(*groups)
tukey_ci = tukey_result.confidence_interval(1 - ALPHA)
tukey_rows = []
for i in range(k - 1):
    for j in range(i + 1, k):
        p = tukey_result.pvalue[i, j]
        tukey_rows.append([f'{group_names[i]}-{group_names[j]}',
            tukey_result.statistic[i, j], tukey_ci.low[i, j],
            tukey_ci.high[i, j], p, p < ALPHA])
tukey = pd.DataFrame(tukey_rows, columns=['comparison', 'difference',
    'CI_low', 'CI_high', 'adjusted_p', 'reject'])
q_critical = stats.studentized_range.ppf(1 - ALPHA, k, df_error)
print('Tukey q临界值：', q_critical)
print(tukey.to_string(index=False))
print('六比较族拒绝数：', tukey['reject'].sum())

# %% 5. Dunnett：三个候选结构各减A1
# *groups[1:]依次传入A2/A3/A4；control=groups[0]指定共同A1。
# rng固定多元t数值积分的随机源，不能保证跨版本末位完全相同。
dun_result = stats.dunnett(*groups[1:], control=groups[0],
    alternative='two-sided', rng=np.random.default_rng(20261130))
dun_ci = dun_result.confidence_interval(1 - ALPHA)
dunnett = pd.DataFrame({'comparison': ['A2-A1', 'A3-A1', 'A4-A1'],
    'difference': means[1:] - means[0], 'adjusted_p': dun_result.pvalue,
    'CI_low': dun_ci.low, 'CI_high': dun_ci.high,
    'reject': dun_result.pvalue < ALPHA})
print(dunnett.to_string(index=False))
print('共同对照族拒绝数：', dunnett['reject'].sum(), '；只保护预定三项。')
# 共同A1参照使三个差值相关；这里用共同MSE和样本量打印可核对的尺度。
v_control = mse / n[0]
v_difference = mse * (1 / n[1:] + 1 / n[0])
cov_shared = mse / n[0]
cor_shared = cov_shared / np.sqrt(v_difference[0] * v_difference[1])
dunnett_covariance = pd.DataFrame([{'v_MSE_over_n_A1': v_control,
    'var_each_treatment_minus_A1': v_difference[0],
    'cov_A2A1_A3A1': cov_shared, 'cor_A2A1_A3A1': cor_shared}])
print(dunnett_covariance.to_string(index=False))
print('说明：共同参照联合t推断使用该相关结构；MSE/n仅演示共同尺度。')

# %% 6. SNK：等组量的教学自编实现
if len(np.unique(n)) != 1:
    raise ValueError('本篇SNK演示仅适用于等组量。')
order = np.argsort(-means)  # 从大到小排序，order存原组号，i/j存排序位置。
snk_rows = []
for i in range(k - 1):
    for j in range(i + 1, k):
        span = j - i + 1   # 包含两个端点，而非只数中间均数。
        q = (means[order[i]] - means[order[j]]) / np.sqrt(mse / n[0])
        qcrit = stats.studentized_range.ppf(1 - ALPHA, span, df_error)
        raw_p = stats.studentized_range.sf(q, span, df_error)
        snk_rows.append([group_names[order[i]], group_names[order[j]],
                         i + 1, j + 1, span, q, qcrit, raw_p])
snk = pd.DataFrame(snk_rows, columns=['high', 'low', 'i', 'j', 'r', 'q',
                                    'qcrit', 'raw_range_P'])
# 某区间未拒绝，其内部比较不能拒绝。每一对须让所有包含区间都通过。
for row in snk.index:
    contains_pair = (snk['i'] <= snk.loc[row, 'i']) & (snk['j'] >= snk.loc[row, 'j'])
    snk.loc[row, 'step_P'] = snk.loc[contains_pair, 'raw_range_P'].max()
snk['reject'] = snk['step_P'] < ALPHA
# A1-A3的所有包含区间：按排序位置动态筛选，便于逐步规则核对。
a1_pos = int(np.where(order == group_names.index('A1'))[0][0]) + 1
a3_pos = int(np.where(order == group_names.index('A3'))[0][0]) + 1
contains_a1_a3 = (snk['i'] <= min(a1_pos, a3_pos)) & (snk['j'] >= max(a1_pos, a3_pos))
snk_a1_a3_ranges = snk.loc[contains_a1_a3, ['high', 'low', 'r', 'raw_range_P']]
print(snk_a1_a3_ranges.to_string(index=False))
print('A1-A3包含区间raw_range_P最大值step_P：', snk_a1_a3_ranges['raw_range_P'].max())
snk = snk.sort_values(['r', 'i'], ascending=[False, True])
print(snk.to_string(index=False))
print('SNK逐步拒绝数：', snk['reject'].sum(), '；step_P不是强FWER调整P。')

# %% 7. 显示实际观测点与均数；按需保存
fig, ax = plt.subplots(figsize=(7, 4.7))
for i, values in enumerate(groups):
    x = i + np.linspace(-.1, .1, len(values))
    ax.scatter(x, values, color='#4d8ec5')
ax.scatter(range(k), means, color='#0000F5', marker='D', label='Group mean')
ax.set(xticks=range(k), xticklabels=group_names, xlabel='Scaffold group',
       ylabel=f'Compressive modulus ({unit})', title='Individual specimens and group means')
ax.legend()
fig.tight_layout()
if SAVE_OUTPUTS:
    out = root / '运行结果' / 'Python'
    out.mkdir(parents=True, exist_ok=True)
    tables = {'summary': summary, 'LSD': lsd, 'Tukey': tukey,
              'Dunnett': dunnett, 'SNK': snk}
    for name, table in tables.items():
        table.to_csv(out / f'{name}.csv', index=False)
    fig.savefig(out / 'means_points.png', dpi=200)
plt.show()
