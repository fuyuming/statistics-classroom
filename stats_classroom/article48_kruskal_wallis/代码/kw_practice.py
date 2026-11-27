# -*- coding: utf-8 -*-
"""统计课堂48：H的构造、六瓶精确枚举、兔例Monte Carlo。
24只独立兔，每组6只；每兔两眼形成一个改善差。数据来源见README。
六瓶为教学假设数据，三条件各2瓶。所有计算读取共享CSV。
Spyder：先将工作目录设为项目根，按 # %% 顺序运行；也可F5整文件。
依赖：numpy、pandas、scipy、matplotlib。图默认显示；保存开关默认关。
"""
# %% 1 准备：路径、参数、共享数据
from pathlib import Path
import itertools
import os
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt

root = Path.cwd()
if not (root / '数据').is_dir() and '__file__' in globals():
    root = Path(__file__).resolve().parents[1]
if not (root / '数据').is_dir():
    raise FileNotFoundError('请把工作目录设为含“数据”文件夹的项目根目录。')
save_figures = os.getenv('SAVE_FIGURES', '0') == '1'
show_figures = os.getenv('KW_HEADLESS', '0') != '1'
alpha = 0.05
mc_draws = 199999
seed = 20261130
out = root / '运行结果' / 'Python'
out.mkdir(parents=True, exist_ok=True)
rabbit = pd.read_csv(root / '数据' / 'rabbit.csv')
print('每行一只兔；drug为药物，score为改善差（分）：')
print(rabbit.head())
if rabbit.isna().any().any():
    raise ValueError('本项目要求数据完整，请先处理缺失。')

# %% 2 合并编秩、B、秩方差、C与H
# rankdata用全部兔共同排序；average使相同评分取得同一平均秩。
rabbit['rank'] = stats.rankdata(rabbit['score'], method='average')
order = list(rabbit['drug'].drop_duplicates())
N = len(rabbit)
k = len(order)
center = (N + 1) / 2
sizes = []
rank_sums = []
mean_ranks = []
groups = []
for drug in order:
    selected = rabbit['drug'] == drug
    sizes.append(int(selected.sum()))
    rank_sums.append(rabbit.loc[selected, 'rank'].sum())
    mean_ranks.append(rabbit.loc[selected, 'rank'].mean())
    groups.append(rabbit.loc[selected, 'score'].to_numpy())
sizes = np.asarray(sizes)
mean_ranks = np.asarray(mean_ranks)
summary = pd.DataFrame({'drug': order, 'n': sizes,
                        'rank_sum': rank_sums, 'mean_rank': mean_ranks})
print(summary.to_string(index=False))
B = np.sum(sizes * (mean_ranks - center)**2)
rank_total_ss = np.sum((rabbit['rank'] - center)**2)
rank_variance = rank_total_ss / (N - 1)
counts = rabbit['score'].value_counts().sort_index().to_numpy()
C = 1 - np.sum(counts**3 - counts) / (N**3 - N)
if rank_variance <= 0:
    raise ValueError('全部评分相同，缺少可检验的秩变异。')
h_raw = 12 * B / (N * (N + 1))
H = B / rank_variance
# sf直接返回χ²右尾；df=k-1。kruskal返回校正后H和χ²近似P。
p_chisq = stats.chi2.sf(H, df=k-1)
native = stats.kruskal(*groups)
print(f'B={B:.9f}, SS_T={rank_total_ss:.9f}, s_r2={rank_variance:.9f}')
print(f'C={C:.12f}, H_raw={h_raw:.9f}, H={H:.9f}, chi2_P={p_chisq:.12f}')
vbase = N * (N + 1) / 12
tbase = (N - 1) * vbase
tie_values = counts
tie_losses = (tie_values**3 - tie_values) / 12
tactual = tbase - tie_losses.sum()
vactual = rank_variance
print(f'Vbase={vbase:.12f}, Tbase={tbase:.12f}')
print('tie frequencies=', tie_values.tolist(), 'losses=', tie_losses.tolist())
print(f'Tactual={tactual:.12f}, Vactual={vactual:.12f}, C=Vactual/Vbase={vactual/vbase:.12f}')
group_var = vactual * (1 / sizes - 1 / N)
group_cov = -vactual / N
print('group mean-rank variances=', group_var.tolist())
print(f'covariance between any two distinct groups={group_cov:.12f}')
if len(set(sizes)) == 1 and k == 4:
    n = sizes[0]
    m1, m2, m3, m4 = mean_ranks
    diffs = np.array([m1-m2, (m1+m2)/2-m3, (m1+m2+m3)/3-m4])
    ses = np.sqrt(np.array([2*vactual/n, 3*vactual/(2*n), 4*vactual/(3*n)]))
    zvals = diffs / ses
    contributions = zvals**2
    contrast_components = pd.DataFrame({
        'contrast': ['m1-m2', '(m1+m2)/2-m3', '(m1+m2+m3)/3-m4'],
        'difference': diffs, 'SE': ses, 'Z': zvals,
        'Z_squared_contribution': contributions})
    contrast_components.to_csv(out / 'contrast_components.csv', index=False)
    print(contrast_components.to_string(index=False))
    print(f'contrast Z^2 sum={contributions.sum():.12f}; H={H:.12f}')
    contrast_sum_z2 = contributions.sum()
else:
    contrast_sum_z2 = np.nan
    print('contrast components skipped: requires four equal-sized groups.')
print('scipy.stats.kruskal:', native)
print('χ²近似判断：', '拒绝同分布H0。' if p_chisq <= alpha else '不足以拒绝同分布H0。')
summary.to_csv(out / 'rabbit_ranks.csv', index=False)

# %% 3 六瓶教学例：观察H来自当前分组，再枚举全部合法划分
bottles = pd.read_csv(root / '数据' / 'six_bottles.csv')
print('每行一独立培养瓶，浓度单位mg/L：\n', bottles)
small_ranks = stats.rankdata(bottles['concentration_mg_L'], method='average')
bottles['rank'] = small_ranks
small_N = len(bottles)
small_center = (small_N + 1) / 2
small_var = np.var(small_ranks, ddof=1)  # ddof=1使用N-1分母。
small_summary = bottles.groupby('condition')['rank'].agg(['size', 'mean'])
if small_N != 6 or list(small_summary['size']) != [2, 2, 2] or small_var <= 0:
    raise ValueError('枚举入口限定6瓶、三条件各2瓶，且合并秩有变异。')
small_B = np.sum(small_summary['size'] * (small_summary['mean'] - small_center)**2)
small_H = small_B / small_var
all_H = []
indices = list(range(6))  # Python下标从0开始，仍对应6个不同身份的瓶。
for a in itertools.combinations(indices, 2):
    remaining = [j for j in indices if j not in a]
    for b in itertools.combinations(remaining, 2):
        cc = [j for j in remaining if j not in b]
        allocation_B = 0.0
        for group_indices in (a, b, cc):
            group_mean = np.mean(small_ranks[list(group_indices)])
            allocation_B += 2 * (group_mean - small_center)**2
        all_H.append(allocation_B / small_var)
all_H = np.asarray(all_H)
# 容差使本应相等的H计入尾部；枚举本身没有Monte Carlo误差。
tail_count = np.sum(all_H >= small_H - 1e-10)
exact_P = tail_count / len(all_H)
small_chisq_P = stats.chi2.sf(small_H, df=2)
values, frequencies = np.unique(np.round(all_H, 12), return_counts=True)
distribution = pd.DataFrame({'H': values, 'frequency': frequencies,
                             'probability': frequencies / len(all_H)})
print(distribution.to_string(index=False))
print(f'辅助 H={small_H:.9f}, total={len(all_H)}, tail={tail_count}')
print(f'exact_P={exact_P:.12f}, chi2_P={small_chisq_P:.12f}')
print('精确枚举判断：', '拒绝同分布H0。' if exact_P <= alpha else '不足以拒绝同分布H0。')
distribution.to_csv(out / 'aux_exact_distribution.csv', index=False)

# %% 4 兔例条件置换：固定每组大小，打乱有身份的单位
# default_rng创建随机数生成器；permutation均匀重排24个秩。
rng = np.random.default_rng(seed)
pooled_ranks = rabbit['rank'].to_numpy()
ends = np.cumsum(sizes)
starts = np.r_[0, ends[:-1]]
exceedances = 0
for draw in range(mc_draws):
    permuted = rng.permutation(pooled_ranks)
    perm_B = 0.0
    for i in range(k):
        perm_mean = np.mean(permuted[starts[i]:ends[i]])
        perm_B += sizes[i] * (perm_mean - center)**2
    perm_H = perm_B / rank_variance
    exceedances += perm_H >= H - 1e-10
mc_P = (exceedances + 1) / (mc_draws + 1)
mc_se = np.sqrt(mc_P * (1 - mc_P) / mc_draws)
print(f'MC B={mc_draws}, seed={seed}, exceedances={exceedances}')
print(f'P_plus1={mc_P:.9f}, MCSE={mc_se:.9f}')
print('条件置换判断：', '拒绝同分布H0。' if mc_P <= alpha else '不足以拒绝同分布H0。')
print('这是整体分布比较；具体药物组对需要另作多重比较。')
pd.DataFrame({'H':[H], 'C':[C], 'chi2_P':[p_chisq],
              'mc_draws':[mc_draws], 'tail_count':[exceedances],
              'mc_P':[mc_P], 'mc_se':[mc_se], 'small_H':[small_H],
              'small_exact_P':[exact_P], 'Vbase':[vbase], 'Tbase':[tbase],
              'Tactual':[tactual], 'Vactual':[vactual],
              'contrast_Z2_sum':[contrast_sum_z2]}).to_csv(out / 'summary.csv', index=False)

# %% 5 绘图：直接显示；SAVE_FIGURES=1才另存图片
fig1, ax1 = plt.subplots(figsize=(7, 4.5))
ax1.bar(order, mean_ranks, color='#3268c8')
ax1.axhline(center, color='#c75c24', linestyle='--', label=f'Pooled mean = {center:g}')
ax1.set_ylabel('Mean rank')
ax1.tick_params(axis='x', rotation=15)
ax1.legend()
fig1.tight_layout()
fig2, ax2 = plt.subplots(figsize=(7, 4.5))
x_step = np.r_[-0.1, values, 5.5]
y_step = np.r_[0, np.cumsum(frequencies) / len(all_H), 1]
ax2.step(x_step, y_step, where='post', label='Exact discrete CDF', linewidth=2)
grid = np.linspace(0, 5.5, 400)
ax2.plot(grid, stats.chi2.cdf(grid, df=2), label='Chi-square(2) CDF')
ax2.axvline(small_H, color='#c75c24', linestyle='--', label='Observed H')
ax2.set(xlabel='H', ylabel='Cumulative probability', ylim=(-0.02,1.05))
ax2.legend()
fig2.tight_layout()
if save_figures:
    image_dir = root / '文章配图'
    image_dir.mkdir(exist_ok=True)
    fig1.savefig(image_dir / 'python_mean_ranks.png', dpi=220)
    fig2.savefig(image_dir / 'python_aux_cdf.png', dpi=220)
if show_figures:
    plt.show()
else:
    plt.close('all')

# %% 6 自己改分组再观察
# 不改变观测和组大小：按浓度1到6，把condition改成A,B,C,C,B,A。
# 三组平均秩都为3.5，重跑第3节；预期H=0、精确P=1。
# 再运行第5节，观察竖线位置变化。不要为追求显著性挑选分组。
