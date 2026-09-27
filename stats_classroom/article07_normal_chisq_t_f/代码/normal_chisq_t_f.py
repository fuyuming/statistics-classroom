"""统计课堂07。安装 numpy scipy matplotlib；运行本文件即可。
UTF-8；数据与输出均根据脚本位置定位。数据为人工教学数据。
"""
from pathlib import Path
import csv, gzip, sys
import numpy as np
import scipy
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / '运行结果' / 'Python'
FIG = ROOT / '文章配图' / 'Python'
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

def save(name, names, columns):
    with (OUT / name).open('w', encoding='utf-8', newline='') as f:
        w = csv.writer(f); w.writerow(names); w.writerows(zip(*columns))

def finish(fig, name):
    fig.savefig(FIG / name, dpi=180); plt.close(fig)

# 1. 理论临界值与错误率。
dfs = np.array([1, 5, 10, 30, 100])
q = stats.t.ppf(.975, dfs); fq = stats.f.ppf(.95, 1, dfs)
critical = np.column_stack([dfs, dfs+1, q, q*q, fq, 2*stats.t.sf(1.96, dfs)])
assert np.allclose(q*q, fq, rtol=1e-10)
save('critical_values.csv', ['df','n_one_sample','t_975','t_975_squared','f_95','rejection_at_1_96'], critical.T)

# 2. 三张分布图及全部坐标，英文标注避免依赖本机中文字体。
colors = ['#178C88','#4144A5','#E36854','#333333']
x = np.linspace(-5,5,1001)
y = [stats.norm.pdf(x)] + [stats.t.pdf(x,k) for k in [1,5,30]]
save('t_normal_curves.csv',['x','normal','t_df1','t_df5','t_df30'],[x]+y)
fig,ax=plt.subplots(figsize=(10,6),layout='constrained')
for vals,label,c,ls in zip(y,['N(0,1)','t: df=1','t: df=5','t: df=30'],['#222222']+colors[:3],['--','-','-','-']):
    ax.plot(x,vals,label=label,color=c,ls=ls)
ax.set(xlabel='Statistic value',ylabel='Density',title='Student t and standard normal distributions');ax.legend()
finish(fig,'02-t_normal.png')
x=np.linspace(.001,55,3001); ks=[1,2,3,5,10,30]
save('chisq_curves.csv',['x']+[f'df_{k}' for k in ks],[x]+[stats.chi2.pdf(x,k) for k in ks])
fig,axes=plt.subplots(2,1,figsize=(10,9),layout='constrained')
for ax, subset in zip(axes,[[1,2],[3,5,10,30]]):
    for k,c,ls in zip(subset,colors,['-','--',':','-.']): ax.plot(x,stats.chi2.pdf(x,k),label=f'df={k}',color=c,ls=ls)
    ax.set(xlabel='Chi-square value',ylabel='Density'); ax.legend()
axes[0].set(xlim=(0,8),ylim=(0,1.5),title='df=1, 2: df=1 density diverges at zero; upper part clipped')
axes[1].set_title('Chi-square: increasing degrees of freedom')
finish(fig,'03_chisq.png')
x=np.linspace(.001,6,2001); pairs=[(3,10),(5,10),(10,10),(5,5),(5,30)]
save('f_curves.csv',['x']+[f'F_{a}_{b}' for a,b in pairs],[x]+[stats.f.pdf(x,a,b) for a,b in pairs])
fig,axes=plt.subplots(2,1,figsize=(10,9),layout='constrained')
for ax,ps,title in zip(axes,[pairs[:3],[pairs[3],pairs[1],pairs[4]]],['Fixed denominator df=10','Fixed numerator df=5']):
    for (a,b),c,ls in zip(ps,colors,['-','--',':']): ax.plot(x,stats.f.pdf(x,a,b),label=f'F({a},{b})',color=c,ls=ls)
    ax.axvline(1,color='gray',ls=':',lw=1);ax.set(xlabel='F value',ylabel='Density',title=title);ax.legend()
finish(fig,'04_f.png')

# 3. 同一教学数据：手写 OLS、现成 t 检验、ANOVA 相互核对。
with (ROOT/'数据'/'teaching_data.csv').open(encoding='utf-8') as f: rows=list(csv.DictReader(f))
y=np.array([float(r['response']) for r in rows]); x=np.array([float(r['x']) for r in rows]); group=np.array([r['group']=='B' for r in rows])
def ols_test(predictor):
    X=np.column_stack([np.ones(len(y)),predictor]); beta=np.linalg.lstsq(X,y,rcond=None)[0]
    residual=y-X@beta; df=len(y)-2; sse=residual@residual
    cov=(sse/df)*np.linalg.inv(X.T@X)
    t=beta[1]/np.sqrt(cov[1,1]); F=((np.sum((y-y.mean())**2)-sse))/(sse/df)
    return t,F,df
res=stats.ttest_ind(y[~group],y[group],equal_var=True); a=stats.f_oneway(y[~group],y[group]); gt,gf,gdf=ols_test(group.astype(float))
values=[res.statistic**2,a.statistic,gt**2]; ps=[res.pvalue,a.pvalue,2*stats.t.sf(abs(gt),gdf)]
assert np.allclose(values,values[0],rtol=1e-10) and np.allclose(ps,ps[0],rtol=1e-10)
save('t_anova_regression_equivalence.csv',['method','statistic','p_value'],[['pooled_t_squared','one_way_ANOVA_F','regression_group_t_squared'],values,ps])
t,F,df=ols_test(x); assert np.isclose(t*t,F,rtol=1e-10)
save('simple_regression_equivalence.csv',['slope_t','slope_t_squared','overall_F','residual_df'],[[t],[t*t],[F],[df]])

# 4. 独立正态平方和构造卡方，再构造 t、F。跨语言随机流不同。
rng=np.random.default_rng(20260927); B=50000; nu=10
z=rng.normal(size=B); u=np.sum(rng.normal(size=(B,nu))**2,axis=1)
tsim=z/np.sqrt(u/nu); fsim=z*z/(u/nu)
assert np.allclose(tsim**2,fsim)
probs=np.array([.025,.5,.975])
save('construction_simulation_quantiles.csv',['probability','empirical_t','theoretical_t'],[probs,np.quantile(tsim,probs),stats.t.ppf(probs,nu)])
with gzip.open(OUT/'construction_draws.csv.gz','wt',encoding='utf-8',newline='') as f:
    w=csv.writer(f);w.writerow(['z','u_chisq10','t10','f1_10']);w.writerows(zip(z,u,tsim,fsim))
report=f'All checks passed.\ncritical values:\n{critical}\nGroup t^2 / ANOVA F: {values}\nRegression slope t^2 / overall F: {t*t}, {F}\nPython {sys.version}\nNumPy {np.__version__}; SciPy {scipy.__version__}; Matplotlib {matplotlib.__version__}\n'
(OUT/'运行记录.txt').write_text(report,encoding='utf-8');print(report)
