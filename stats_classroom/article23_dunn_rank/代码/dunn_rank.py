"""统计课堂23：Dunn与Holm。Spyder打开本项目，在本文件按F5。
数据：第12课33只小鼠的完整存活天数；每只只属于一个组。
需要NumPy、pandas、SciPy、Matplotlib；结果写入运行结果/Python。
"""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import scipy
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / '运行结果/Python'
OUT.mkdir(parents=True, exist_ok=True)

def save_table(frame, name):
    pd.DataFrame(frame).to_csv(OUT / name, index=False, encoding='utf-8-sig')

def holm(p):
    """先排序，乘递减的比较数，累计最大值，再还原组对顺序。"""
    order = np.argsort(p)
    multipliers = np.arange(len(p), 0, -1)
    products = p[order] * multipliers
    cumulative = np.maximum.accumulate(products)
    sorted_adjusted = np.minimum(1, cumulative)
    adjusted = np.empty(len(p))
    adjusted[order] = sorted_adjusted
    return adjusted, order, products, cumulative

# %% 1. 先查看数据，再调用完整KW函数
frame = pd.read_csv(ROOT / '数据/小鼠存活天数.csv')
assert len(frame) == 33 and not frame.isna().any().any()
assert frame.id.is_unique
print(frame)
print(frame.groupby('group').size())
labels = ['A', 'B', 'C']
values = frame.days.to_numpy()
grouping = frame.group.to_numpy()
groups = [values[grouping == label] for label in labels]
standard_kw = stats.kruskal(*groups)
print('标准KW结果：', standard_kw)

# %% 2. 全体编秩，手算Dunn和Holm
# 定义可重复调用的分析函数，后面用乘10后的数据重新分析。
def calculate(x, g):
    ranks = stats.rankdata(x, method='average')
    n = np.array([np.sum(g == label) for label in labels])
    means = np.array([np.mean(ranks[g == label]) for label in labels])
    sr2 = np.var(ranks, ddof=1)  # 样本方差，除以N−1
    N = len(x)
    _, counts = np.unique(x, return_counts=True)
    tie_total = np.sum(counts**3 - counts)
    between = np.sum(n * (means - np.mean(ranks))**2)
    H = between / sr2  # 使用平均并列秩，已包含并列修正
    first = np.array([0, 0, 1])
    second = np.array([1, 2, 2])
    difference = means[second] - means[first]
    se = np.sqrt(sr2 * (1/n[first] + 1/n[second]))
    z = difference / se
    p = 2 * stats.norm.sf(np.abs(z))  # sf为右尾，乘2计入双侧
    adjusted, order, products, cumulative = holm(p)
    pair = np.array(['B-A', 'C-A', 'C-B'])
    return {
        'ranks': ranks,
        'groups': pd.DataFrame({'group': labels, 'n': n, 'rank_sum': n*means, 'mean_rank': means}),
        'kw': pd.DataFrame([{'N':N, 'H':H, 'P':stats.chi2.sf(H,2), 'sr2':sr2,
                             'tie_total':tie_total, 'epsilon2':H/(N-1)}]),
        'dunn': pd.DataFrame({'pair':pair, 'rank_diff':difference, 'SE':se,
                              'Z':z, 'P_raw':p, 'P_holm':adjusted}),
        'holm': pd.DataFrame({'pair':pair[order], 'P_raw':p[order], 'multiplier':[3,2,1],
                              'product':products, 'cumulative_max':cumulative, 'P_holm':adjusted[order]})
    }
result = calculate(values, grouping)
assert abs(result['kw'].H.iloc[0] - standard_kw.statistic) < 1e-10
assert abs(result['kw'].P.iloc[0] - standard_kw.pvalue) < 1e-12
frame['rank'] = result['ranks']
save_table(frame, 'ranks.csv')
for key, name in [('groups','groups.csv'), ('kw','kw.csv'), ('dunn','dunn.csv'), ('holm','holm_steps.csv')]:
    save_table(result[key], name)
unique, counts = np.unique(values, return_counts=True)
save_table({'days':unique, 'count':counts}, 'ties.csv')
print(result['dunn'])
print('均秩差用秩位表示，P_holm已对三对比较调整。')

# %% 3. 两两Wilcoxon/Mann–Whitney：单独编秩，近似且无连续性校正
pairs = [(0,1), (0,2), (1,2)]
w_rows = []
for index, (i,j) in enumerate(pairs):
    test = stats.mannwhitneyu(groups[i], groups[j], alternative='two-sided',
                              method='asymptotic', use_continuity=False)
    w_rows.append({'pair':result['dunn'].pair.iloc[index], 'P_no_continuity':test.pvalue})
w_table = pd.DataFrame(w_rows)
w_table['P_holm'] = holm(w_table.P_no_continuity.to_numpy())[0]
save_table(w_table, 'pairwise_wilcoxon.csv')
ab = frame[frame.group.isin(['A','B'])].copy()
ab['rank_ab'] = stats.rankdata(ab.days)
save_table(ab, 'ab_rerank.csv')
assert (ab.loc[ab.days == 4, 'rank_ab'] == 7).all()
print(w_table)

# %% 4. 置换：固定全体秩，分为A11/B10/C12
rng = np.random.default_rng(20261103)
B_permutation = 10000
difference_sim = np.empty(B_permutation)
for b in range(B_permutation):
    shuffled = rng.permutation(result['ranks'])
    mean_A = shuffled[:11].mean()
    mean_B = shuffled[11:21].mean()
    difference_sim[b] = mean_B - mean_A
save_table({'difference':difference_sim}, 'permutation_differences.csv')
permutation_summary = pd.DataFrame([{'B':B_permutation, 'seed':20261103,
    'mean':difference_sim.mean(), 'SD':difference_sim.std(ddof=1),
    'theoretical_SD':result['dunn'].SE.iloc[0]}])
save_table(permutation_summary, 'permutation_summary.csv')
print(permutation_summary)

# %% 5. 组内有放回bootstrap，给整体秩效应百分位区间
rng = np.random.default_rng(20261102)
B_bootstrap = 20000
epsilon_sim = np.empty(B_bootstrap)
for b in range(B_bootstrap):
    draw_groups = []
    for original in groups:
        draw = rng.choice(original, size=len(original), replace=True)
        draw_groups.append(draw)
    # KW会对本轮三组合并编秩并处理并列。
    H_draw = stats.kruskal(*draw_groups).statistic
    epsilon_sim[b] = H_draw / 32
interval = np.quantile(epsilon_sim, [0.025,0.975], method='linear')
bootstrap_summary = pd.DataFrame([{'B':B_bootstrap, 'seed':20261102,
    'epsilon2_observed':result['kw'].epsilon2.iloc[0],
    'percentile_2_5':interval[0], 'percentile_97_5':interval[1]}])
save_table({'epsilon2':epsilon_sim}, 'bootstrap_epsilon.csv')
save_table(bootstrap_summary, 'bootstrap_summary.csv')
print(bootstrap_summary)
print('这个区间描述整体秩效应的不确定性，不用是否跨0替代KW检验。')

# %% 6. 乘10后重新计算，而非复制原来的结果
changed = calculate(values * 10, grouping)
assert np.array_equal(result['ranks'], changed['ranks'])
assert np.allclose(result['dunn'].P_holm, changed['dunn'].P_holm, atol=1e-12)
assert abs(result['kw'].epsilon2.iloc[0] - changed['kw'].epsilon2.iloc[0]) < 1e-12
save_table({'pair':changed['dunn'].pair, 'P_raw_x10':changed['dunn'].P_raw,
    'P_holm_x10':changed['dunn'].P_holm, 'H_x10':changed['kw'].H.iloc[0],
    'epsilon2_x10':changed['kw'].epsilon2.iloc[0]}, 'unit_check.csv')
print('乘10后秩、P和秩效应不变，原始天数的中位数则乘10。')
(OUT/'versions.txt').write_text(f'Python {sys.version.split()[0]}\nNumPy {np.__version__}\npandas {pd.__version__}\nSciPy {scipy.__version__}\n')

# %% 7. 绘图；图2明确读取R置换结果，便于与公众号一致
# R输出随本项目提供；重新生成R图的数据时，先运行R主脚本。
sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_dunn import main as draw_figures
draw_figures()
