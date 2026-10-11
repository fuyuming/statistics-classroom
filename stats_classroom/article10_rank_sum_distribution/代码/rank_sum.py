# %% 统计课堂10：秩和的精确分布。数据为人工教学评分。
# Spyder打开本文件，F5完整运行；UTF-8保存。依赖见requirements。
# 主脚本自动调用plot_rank_sum.py。无并列、并列两列来自共享输入。
from pathlib import Path
from itertools import combinations
import numpy as np
import pandas as pd
from scipy import stats
from plot_rank_sum import make_plots

root = Path(__file__).resolve().parents[1]
out = root/'运行结果'/'Python'
out.mkdir(parents=True,exist_ok=True)
messages = []
def report(value):
    print(value)
    messages.append(str(value))
data = pd.read_csv(root/'数据'/'scores.csv')
report('输入：group为组别，value为无并列评分，value_tied为并列评分，均为人工教学数据。')
report(data.to_string(index=False))
is_a = data.group.eq('A').to_numpy()
m = int(is_a.sum())
n = len(data)-m
N = m+n
assert m == 3 and n == 3
assert data.value.nunique() == N
# combinations选的是记录位置，不是去掉重复后的数值；并列时也保留所有记录。
allocations = np.array(list(combinations(range(N),m)))
mean_w = m*(N+1)/2
results = {}

# %% 1. 每种场景都逐步求秩、枚举全部20种分配
for column in ['value','value_tied']:
    ranks = stats.rankdata(data[column],method='average')
    w_obs = ranks[is_a].sum()
    u_obs = w_obs-m*(m+1)/2
    all_w = ranks[allocations].sum(axis=1)
    # 全部20种分配构成精确分布，不是20个样本；方差使用概率平均，即ddof=0。
    var_exact = ((all_w-all_w.mean())**2).mean()
    p_two = np.mean(np.abs(all_w-mean_w) >= abs(w_obs-mean_w)-1e-12)
    tied_sizes = data[column].value_counts().to_numpy()
    correction = np.sum(tied_sizes**3-tied_sizes)
    var_formula = m*n/12*(N+1-correction/(N*(N-1)))
    report(f'\n问题：{column}；秩={ranks}；W={w_obs}；U={u_obs}。')
    report(f'精确均值={all_w.mean():.6f}；方差={var_exact:.6f}；双侧P={p_two:.6f}。')
    report('P按离零假设中心至少同样远定义；每种记录分配的权重为1/20。')
    table = pd.DataFrame(allocations+1,columns=['record1','record2','record3'])
    table['W'] = all_w
    table.to_csv(out/f'{column}_allocations.csv',index=False)
    values,counts = np.unique(all_w,return_counts=True)
    pmf = pd.DataFrame({'W':values,'count':counts,'probability':counts/len(all_w)})
    pmf.to_csv(out/f'{column}_pmf.csv',index=False)
    report(pmf.to_string(index=False))
    assert np.isclose(var_exact,var_formula)
    assert np.isclose(all_w.mean(),mean_w)
    results[column] = (w_obs,u_obs,var_exact,p_two,pmf)

# %% 2. 核对软件统计量名称：SciPy返回U而非原始秩和
x = data.loc[is_a,'value']
y = data.loc[~is_a,'value']
test = stats.mannwhitneyu(x,y,alternative='two-sided',method='exact')
assert np.isclose(test.statistic,results['value'][1])
assert np.isclose(test.pvalue,results['value'][3])
report(f'\nSciPy无并列精确检验：U={test.statistic}，P={test.pvalue}。')
report('有并列值时这里自行枚举条件分布，不套用SciPy的无并列exact算法。')
# 无并列W为整数；左尾在w+0.5处用连续曲线近似。
w = results['value'][0]
var_w = results['value'][2]
p_normal = 2*stats.norm.cdf((w+.5-mean_w)/np.sqrt(var_w))
report(f'本例连续性校正的双侧正态近似={p_normal:.6f}，精确值=0.2。')
metrics = pd.DataFrame({'name':['W','U','mean_W','var_W','p_exact','var_tied','p_tied','p_normal_cc'],
    'value':[w,results['value'][1],mean_w,var_w,results['value'][3],results['value_tied'][2],results['value_tied'][3],p_normal]})
metrics.to_csv(out/'theory.csv',index=False)

# %% 3. 大样本模拟：均匀选取20个名次，非原始实验数据
rng = np.random.default_rng(20261011)
m_large = 20
n_large = 20
N_large = m_large+n_large
repeats = 20000
# 每行独立均匀随机数的排序索引构成一个随机排列；取前m个是均匀不放回选择。
permutations = np.argsort(rng.random((repeats,N_large)),axis=1)+1
w_sim = permutations[:,:m_large].sum(axis=1)
center = m_large*(N_large+1)/2
sd = np.sqrt(m_large*n_large*(N_large+1)/12)
z_sim = (w_sim-center)/sd
pd.DataFrame({'W':w_sim,'Z':z_sim}).to_csv(out/'large_simulation.csv',index=False)
report(f'模拟m=n=20：理论W中心={center:.6f}，SD={sd:.6f}；模拟Z均值={z_sim.mean():.6f}，方差={z_sim.var(ddof=1):.6f}。')
report('模拟仅检查分布表现，不是正态极限的证明。把m_large和n_large改小，可观察离散性。')
(out/'逐题结果.txt').write_text('\n'.join(messages),encoding='utf-8')
make_plots(root,data,results['value'][4],z_sim,m_large,n_large)
