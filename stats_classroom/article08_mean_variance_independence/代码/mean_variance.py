# 明哥的微生物世界 · 均值与方差的独立性 · Python / Spyder
# 目标：核对人工例子、协方差展开、非线性依赖和正态样本的t构造。
# 数据为指定模型模拟，不是真实实验。模拟不能证明独立性。
# 下载整个项目后打开本文件，F5完整运行。UTF-8；依赖见requirements.txt。

# %% 0. 准备目录和工具
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import scipy
from scipy import stats
from plot_mean_variance import make_figures

root = Path(__file__).resolve().parents[1]
out = root / '运行结果' / 'Python'
out.mkdir(parents=True, exist_ok=True)
messages = []

def explain(message):
    print(message)
    messages.append(str(message))

def save(frame, filename):
    frame.to_csv(out / filename, index=False, float_format='%.15g')

# %% 1. 整体平移、围绕中心展开：均值和方差怎样改变？
examples = pd.read_csv(root / '数据' / 'position_examples.csv')
explain('问题1：先看共享教学数据的字段和前三行。')
explain(examples.to_string(index=False))
values = examples[['x1', 'x2', 'x3']].to_numpy()
assert np.isfinite(values).all(), '输入含缺失或非有限值，请检查数据。'
# axis=1表示对每一行计算；ddof=1让分母为n-1，与R/MATLAB一致。
examples['mean'] = values.mean(axis=1)
examples['variance'] = values.var(axis=1, ddof=1)
save(examples, 'position_results.csv')
explain(examples.to_string(index=False))
explain('答案：均值为10、20、10，样本方差为4、4、16。不能据此证明独立。')
assert np.allclose(examples['variance'], [4, 4, 16])

# %% 2. 用相同输入核对Cov(M,D)展开式，观察n=10的样本均值与方差
summary = []
conditional = []
identity_rows = []
n = 10
mu = 10
sigma = 2
# 0.975是双侧5%检验的上侧界值所对应的左侧累计概率。
t_limit = stats.t.ppf(0.975, df=n-1)
for model in ['normal', 'shifted_exponential']:
    data = pd.read_csv(root / '数据' / f'{model}_samples.csv')
    explain(f'\n问题2：{model}输入示例（每行一份n=10样本）')
    explain(data.head(3).to_string(index=False))
    x = data[[f'x{i}' for i in range(1, n+1)]].to_numpy()
    assert np.isfinite(x).all()
    mean = x.mean(axis=1)
    variance = x.var(axis=1, ddof=1)
    x1 = x[:, 0]
    x2 = x[:, 1]
    midpoint = (x1 + x2) / 2
    difference = x1 - x2
    # np.cov返回协方差矩阵，[0,1]取两列之间的协方差。
    cov_md = np.cov(midpoint, difference, ddof=1)[0, 1]
    rhs = (np.var(x1, ddof=1) - np.var(x2, ddof=1)) / 2
    pair_variance = np.var(x[:, :2], axis=1, ddof=1)
    max_pair_error = np.max(np.abs(pair_variance - difference**2 / 2))
    identity_rows.append({'model': model, 'cov_MD': cov_md, 'half_var_difference': rhs,
                          'identity_error': abs(cov_md-rhs), 'max_pair_variance_error': max_pair_error})
    assert np.isclose(cov_md, rhs, atol=1e-12)
    assert max_pair_error < 1e-12
    # 正态模型的Z和Q独立；指数模型可同样算出这些数，但没有同样的分布保证。
    z = (mean - mu) / (sigma / np.sqrt(n))
    q = (n-1) * variance / sigma**2
    t_constructed = z / np.sqrt(q / (n-1))
    t_direct = (mean - mu) / np.sqrt(variance / n)
    assert np.allclose(t_constructed, t_direct, atol=1e-12)
    save(pd.DataFrame({'sample_id': data.sample_id, 'mean': mean, 'variance': variance,
                       'M_n2': midpoint, 'D_n2': difference, 'S2_n2': pair_variance,
                       'Z': z, 'Q': q, 'T': t_direct}), f'{model}_statistics.csv')
    # 理论Cov(Xbar,S²)=总体三阶中心矩/n；正态为0，平移指数为16/10。
    theoretical_cov = 0 if model == 'normal' else 1.6
    rejection = np.mean(np.abs(t_direct) > t_limit)
    summary.append({'model': model, 'repetitions': len(x), 'n': n,
                    'mean_of_means': mean.mean(), 'mean_of_variances': variance.mean(),
                    'cov_mean_variance': np.cov(mean, variance, ddof=1)[0, 1],
                    'theoretical_cov': theoretical_cov, 'corr_mean_variance': np.corrcoef(mean, variance)[0, 1],
                    't_rejection_fraction': rejection,
                    'rejection_mcse': np.sqrt(rejection*(1-rejection)/len(x)), 't_975': t_limit})
    # 固定门槛筛选均值较低/较高的样本；门槛不由本次模拟结果决定。
    for label, selected in [('all', np.ones(len(x), dtype=bool)),
                            ('mean_lt_9_5', mean < 9.5), ('mean_gt_10_5', mean > 10.5)]:
        group_variance = variance[selected]
        conditional.append({'model': model, 'selection': label, 'count': int(selected.sum()),
                            'mean_variance': group_variance.mean(),
                            'q25_variance': np.quantile(group_variance, 0.25),
                            'median_variance': np.quantile(group_variance, 0.5),
                            'q75_variance': np.quantile(group_variance, 0.75)})

save(pd.DataFrame(identity_rows), 'covariance_identity.csv')
save(pd.DataFrame(summary), 'sampling_summary.csv')
save(pd.DataFrame(conditional), 'conditional_summary.csv')
explain(pd.DataFrame(summary).to_string(index=False))
explain(pd.DataFrame(conditional).to_string(index=False))
explain('解释：有限模拟的协方差不必恰好为0。几个分位数相近也不是独立性的证明。')
explain('正态模型理论拒绝率为5%；指数模型使用同一t界值仅作反例比较。')

# %% 3. 不相关为什么不足以证明独立？
pairs = pd.read_csv(root / '数据' / 'joint_pairs.csv')
explain('\n问题3：图2输入字段为pair_id、x、independent_y、sign。')
explain(pairs.head(3).to_string(index=False))
assert pairs['sign'].isin([-1, 1]).all()
assert np.isfinite(pairs.to_numpy()).all()
pair_rows = []
for label in ['independent', 'square', 'random_sign']:
    if label == 'independent':
        y = pairs['independent_y'].to_numpy()
    elif label == 'square':
        y = pairs['x'].to_numpy()**2
    else:
        y = pairs['sign'].to_numpy() * pairs['x'].to_numpy()
    pair_rows.append({'model': label, 'cov_XY': np.cov(pairs.x, y, ddof=1)[0, 1],
                      'corr_XY': np.corrcoef(pairs.x, y)[0, 1], 'theoretical_cov': 0})
    pairs[f'{label}_y_calculated'] = y
save(pairs, 'figure2_coordinates.csv')
save(pd.DataFrame(pair_rows), 'pair_summary.csv')
explain(pd.DataFrame(pair_rows).to_string(index=False))
explain('三种模型理论协方差都为0，后两种仍有确定的非线性约束。')

# %% 4. 自动绘图：复杂排版在配套文件中
make_figures(root, examples, pairs)
explain(f'Python {sys.version.split()[0]}；NumPy {np.__version__}；pandas {pd.__version__}；SciPy {scipy.__version__}')
explain('练习：把筛选门槛9.5、10.5改成9、11，先猜样本数及估计波动会怎样变化。')
(out / '运行记录.txt').write_text('\n'.join(messages), encoding='utf-8')
