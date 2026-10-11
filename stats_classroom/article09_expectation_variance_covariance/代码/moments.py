# %% 统计课堂09：期望、方差、协方差
# 教学模型：两独立二项样本、四类等概率多项样本；不含真实实验数据。
# 在Spyder打开本文件，F5完整运行。需要numpy、pandas、scipy、matplotlib。
# UTF-8保存；本文件自动定位项目根目录。每个问题均打印结果与解释。
from pathlib import Path
import numpy as np
import pandas as pd
from plot_moments import make_plots

root = Path(__file__).resolve().parents[1]
out = root / '运行结果' / 'Python'
out.mkdir(parents=True, exist_ok=True)
lines = []
def report(message):
    print(message)
    lines.append(str(message))

groups = pd.read_csv(root / '数据' / 'groups.csv')
categories = pd.read_csv(root / '数据' / 'categories.csv')
report('共享输入表：n为样本量，p为给定总体概率，observed为教学观察频数。')
report(groups.to_string(index=False))
report(categories.to_string(index=False))
assert (groups.n > 0).all() and ((groups.p > 0) & (groups.p < 1)).all()
assert np.isclose(categories.p.sum(), 1)
assert len(categories) == 4 and np.allclose(categories.p, 0.25)
# 以下三个比较方向专用于四类等概率。改概率时，需重新构造比较。

# %% 1. 相减的方差：分别计算两组，再相加
n1, n2 = groups.n.to_numpy()
p1, p2 = groups.p.to_numpy()
mean_d = p1 - p2
var1 = p1 * (1 - p1) / n1
var2 = p2 * (1 - p2) / n2
var_d = var1 + var2
se_d = np.sqrt(var_d)
report(f'问题1：率差中心={mean_d:.6f}，方差={var_d:.6f}，标准误={se_d:.6f}。')
report('解释：两组独立，协方差为0；标准误为9.19个百分点，绝非单份结果的标准差。')

# %% 2. 分类计数的协方差；构造整个协方差表
observed = categories.observed.to_numpy()
p = categories.p.to_numpy()
n = observed.sum()
expected = n * p
# 对角线是方差np(1-p)，非对角线是协方差-n*p_j*p_l。
cov_counts = n * (np.diag(p) - np.outer(p, p))
report('问题2：分类计数协方差矩阵（本例对角15，非对角-5）：')
report(cov_counts)
report(f'计数之和的方差=所有矩阵元素相加={cov_counts.sum():.6f}。')
contribution = (observed - expected)**2 / expected
q = contribution.sum()
# 分子分别是两两差与前两类合计减后两类合计。
z1 = (observed[0] - observed[1]) / np.sqrt(n/2)
z2 = (observed[2] - observed[3]) / np.sqrt(n/2)
z3 = (observed[0] + observed[1] - observed[2] - observed[3]) / np.sqrt(n)
z = np.array([z1, z2, z3])
report(f'问题3：Q={q:.6f}；三个方向的平方={z**2}；平方和={sum(z**2):.6f}。')
report('解释：Q=4是这一次数据的值；零假设下Q的期望是3；两者不同。')
# 此矩阵仅用于核对三个比较的方差与协方差，不要求初学者先学矩阵计算。
a = np.array([[1,-1,0,0],[0,0,1,-1],[1,1,-1,-1]], dtype=float)
a = a / np.sqrt([n/2,n/2,n])[:,None]
cov_z = a @ cov_counts @ a.T
assert np.allclose(cov_z, np.eye(3))
assert np.isclose(q, sum(z**2))
assert np.isclose(cov_counts.sum(), 0)
metrics = pd.DataFrame({'name':['mean_d','var_d','se_d','var_count','cov_count','Q','z1_sq','z2_sq','z3_sq','EQ'],
                       'value':[mean_d,var_d,se_d,cov_counts[0,0],cov_counts[0,1],q,z1*z1,z2*z2,z3*z3,len(p)-1]})
metrics.to_csv(out/'theory.csv',index=False)
pd.DataFrame(cov_counts).to_csv(out/'count_covariance.csv',index=False)

# %% 3. 重复抽样：分布模拟，不用模拟代替证明
rng = np.random.default_rng(20261011)
repeats = 20000
# 每次生成一个计数，再除以各组样本量得到样本率。
c1 = rng.binomial(n1, p1, repeats)
c2 = rng.binomial(n2, p2, repeats)
differences = c1/n1 - c2/n2
counts = rng.multinomial(n, p, repeats)  # 每行是一轮，四列为四类计数。
q_sim = ((counts - expected)**2 / expected).sum(axis=1)
pd.DataFrame({'difference':differences}).to_csv(out/'rate_simulation.csv',index=False)
sim = pd.DataFrame(counts,columns=list(categories.category))
sim['Q'] = q_sim
sim.to_csv(out/'category_simulation.csv',index=False)
report(f'模拟核对：率差均值={differences.mean():.6f}，方差={differences.var(ddof=1):.6f}，Q均值={q_sim.mean():.6f}。')
report('ddof=1采用样本方差的n-1分母，估计模拟分布的方差；有限模拟允许随机误差。')
report('练习：把四类计数改成10、20、25、25，Q应为7.5；概率仍为各0.25。')
(out/'逐题结果.txt').write_text('\n'.join(lines),encoding='utf-8')
make_plots(root, differences, q_sim, z, mean_d, se_d)
