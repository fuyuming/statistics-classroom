# %% 统计课堂08：从精确二项分布到样本率、率差与Q
# Spyder中打开本文件，F5完整运行。UTF-8保存，依赖见requirements.txt。
# 输入全是教学模型。图2和图3是理论计算，不是模拟实验。
from pathlib import Path
from math import gcd
import numpy as np
import pandas as pd
from scipy.stats import binom, norm
from plot_expectation_variance_covariance import make_plots, make_reference_plots
root = Path(__file__).resolve().parents[1]
out = root/'运行结果'/'v2'/'Python'
out.mkdir(parents=True, exist_ok=True)
lines, results = [], {}
def report(text):
    print(text)
    lines.append(str(text))

def save(frame, name):
    frame.to_csv(out/name, index=False, encoding='utf-8-sig')

# %% 1. 单样本率：先精确计算，再比较正态近似
# 每一行是一组给定的n、p；k是所问事件“计数至少k”的边界。
cases = pd.read_csv(root/'数据'/'binomial_cases_v2.csv')
curves = []
for row in cases.itertuples(index=False):
    n, p = int(row.n), float(row.p)
    x = np.arange(n+1)
    se = np.sqrt(p*(1-p)/n)
    exact = binom.pmf(x,n,p)
    # 用小区间概率与离散概率比较，不把正态密度误当作概率。
    normal_mass = norm.cdf((x+.5)/n,p,se)-norm.cdf((x-.5)/n,p,se)
    curves.append(pd.DataFrame(dict(case=row.case,n=n,p=p,x=x,proportion=x/n,
                                    exact=exact,normal_mass=normal_mass)))
curve = pd.concat(curves, ignore_index=True)
save(curve,'single_distribution.csv')
n, p, k = 10, .1, 3
se = np.sqrt(p*(1-p)/n)
# sf(k-1) = P(X>k-1) = P(X>=k)，不能写成sf(k)。
results.update(single_mean=p,single_var=se**2,single_se=se,
               single_exact=binom.sf(k-1,n,p),
               single_normal=norm.sf(k/n,p,se),
               single_cc=norm.sf((k-.5)/n,p,se))
report('问题1：n=10、p=0.1，样本率至少0.3有多大概率？')
report(f"精确={results['single_exact']:.9f}；未校正={results['single_normal']:.9f}；校正={results['single_cc']:.9f}。")
report('解释：边界3移动到2.5，除以10得到0.25；np=1，连续性校正不能消除偏斜造成的误差。')

# %% 2. 独立两组：枚举所有计数组合，合并相同的率差
G = pd.read_csv(root/'数据'/'groups.csv')
n1,n2 = G.n.to_numpy(dtype=int)
p1,p2 = G.p.to_numpy(dtype=float)
assert min(n1,n2)>0 and 0<p1<1 and 0<p2<1
mean_d = p1-p2
var_d = p1*(1-p1)/n1+p2*(1-p2)/n2
se_d = np.sqrt(var_d)
# 整数键避免浮点数将同一个差值拆成两个格点。
# n1=40,n2=50时，公分母L=200，key=5*x1-4*x2。
L = int(n1*n2//gcd(int(n1),int(n2)))
x1,x2 = np.meshgrid(np.arange(n1+1),np.arange(n2+1),indexing='ij')
key = x1*(L//n1)-x2*(L//n2)
prob = binom.pmf(x1,n1,p1)*binom.pmf(x2,n2,p2)
joint = pd.DataFrame(dict(x1=x1.ravel(),x2=x2.ravel(),key=key.ravel(),
                          difference=key.ravel()/L,probability=prob.ravel()))
save(joint,'joint_pairs.csv')
pmf = joint.groupby('key',as_index=False).probability.sum()
pmf['difference'] = pmf.key/L
save(pmf,'difference_pmf.csv')
# 文章固定问题D>=0.2。使用整数键比较，避开0.2的二进制表示误差。
threshold = .2
cut_key = int(np.ceil(threshold*L-1e-12))
exact_tail = pmf.loc[pmf.key>=cut_key,'probability'].sum()
# 本例阈值附近格点连续；这是应用半格校正所检查的局部条件。
assert all(v in set(pmf.key) for v in [cut_key-1,cut_key,cut_key+1])
cc_boundary = (cut_key-.5)/L
results.update(mean_d=mean_d,var_d=var_d,se_d=se_d,
               pair_count=len(joint),difference_exact=exact_tail,
               difference_normal=norm.sf(threshold,mean_d,se_d),
               difference_cc=norm.sf(cc_boundary,mean_d,se_d),cc_boundary=cc_boundary)
report('\n问题2：率差的中心与波动怎样计算？')
report(f'期望={mean_d:.9f}；方差={var_d:.9f}；标准误={se_d:.9f}。')
report('解释：期望线性相减；方差先保留-2Cov，因两组独立才消去。')
report(f'问题3：{len(joint)}种联合组合；D>=0.2的精确概率={exact_tail:.9f}。')
report(f"未校正={results['difference_normal']:.9f}；边界{cc_boundary:.6f}校正后={results['difference_cc']:.9f}。")
report('解释：总体率差设为0.1，这是抽样尾概率，不是检验两总体率相等的P值。')
# 精确分布不只检查概率和，还要核对解析均值、方差。
assert np.isclose(pmf.probability.sum(),1,atol=1e-12)
assert np.isclose((pmf.difference*pmf.probability).sum(),mean_d,atol=1e-12)
assert np.isclose(((pmf.difference-mean_d)**2*pmf.probability).sum(),var_d,atol=1e-12)

# %% 3. 分类计数：协方差与Q；这部分不需要正态假设
C = pd.read_csv(root/'数据'/'categories.csv')
observed = C.observed.to_numpy()
p = C.p.to_numpy()
assert len(p)==4 and np.allclose(p,.25)
n = observed.sum()
expected = n*p
cov_counts = n*(np.diag(p)-np.outer(p,p))
q = np.sum((observed-expected)**2/expected)
# 四类等概率时的三个比较方向；不能原样用于任意概率。
a = np.array([[1,-1,0,0],[0,0,1,-1],[1,1,-1,-1]])/np.sqrt([n/2,n/2,n])[:,None]
z = a@observed
assert np.allclose(a@cov_counts@a.T,np.eye(3))
assert np.isclose(q,np.sum(z*z))
results.update(var_count=cov_counts[0,0],cov_count=cov_counts[0,1],Q=q,
               z1_sq=z[0]**2,z2_sq=z[1]**2,z3_sq=z[2]**2,EQ=len(p)-1)
save(pd.DataFrame(cov_counts),'count_covariance.csv')
report('\n问题4：四分类协方差矩阵：\n'+str(cov_counts))
report(f'计数总和方差={cov_counts.sum():.9f}；Q={q:.9f}；三个平方={z*z}；E(Q)=3。')
report('解释：总数固定，分类计数负相关；Q=三个比较的平方和是代数恒等式，卡方参照仍是大样本近似。')
# 保存理论结果，便于三语言独立实现后交叉核对。
save(pd.DataFrame(dict(name=list(results),value=list(results.values()))),'theory.csv')
(out/'逐题结果.txt').write_text('\n'.join(lines),encoding='utf-8')
make_plots(root,curve,pmf,results)

# %% 4. 复现协方差曲线与Q模拟图；模拟用于观察近似，不代替证明。
rng = np.random.default_rng(20261011)
counts = rng.multinomial(int(n), p, 20000)
q_sim = ((counts-expected)**2/expected).sum(axis=1)
save(pd.DataFrame({"Q":q_sim}),"category_simulation.csv")
make_reference_plots(root,q_sim,z)
