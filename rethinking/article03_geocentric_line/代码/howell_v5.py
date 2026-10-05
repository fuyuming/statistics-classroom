"""精读03 v5：作者的中心化模型；三维网格与二次近似。
Spyder 打开后 F5，或 python 代码/howell_v5.py。依赖 numpy/scipy/matplotlib。
Python 的 quadratic 是同原理实现，不冒称调用了 R 的 quap。
先运行 R 可额外核对真实 quap；没有 R 也能独立生成全部图表。
"""
# %% 共享数据：每一行是一位成年人，身高厘米、体重千克。
from pathlib import Path
import json
import numpy as np
from scipy.special import logsumexp
from scipy.stats import norm, lognorm, qmc
from scipy.optimize import minimize
ROOT=Path(__file__).resolve().parent.parent
RES=ROOT/'运行结果/v5'
RES.mkdir(exist_ok=True,parents=True)
d=np.genfromtxt(ROOT/'运行结果/00_教材数据_Howell1成年人.csv',names=True,delimiter=',')
h=d['height']; y=d['weight']; hc=h.mean(); x=h-hc; n=len(y)
print('前五位的身高、体重：',np.column_stack([h,y])[:5])
# 充分统计量只用来加速；与逐个计算残差平方后相加相同。
sy=y.sum(); sx=x.sum(); sxx=x@x; sxy=x@y; syy=y@y

def sse(a,b):
    return syy-2*a*sy-2*b*sxy+n*a*a+2*a*b*sx+b*b*sxx

# %% 三维等距网格：计算窗口不是先验边界。
def grid_fit(k=101,wide=False):
    ranges=[(42,48),(.4,.85),(3,6)] if wide else [(43,47),(.45,.8),(3.3,5.5)]
    axes=[np.linspace(lo,hi,k) for lo,hi in ranges]
    a,b,s=np.meshgrid(*axes,indexing='ij')
    lp=-n*np.log(s)-sse(a,b)/(2*s*s)+norm.logpdf(a,60,10)+lognorm.logpdf(b,1)
    weights=np.exp(lp-logsumexp(lp))
    means=np.array([np.sum(weights*v) for v in [a,b,s]])
    edge=sum(weights.take([0,-1],axis=i).sum() for i in range(3))
    return axes,weights,means,float(edge)
axes,weights,means,edge=grid_fit()
# 加密及扩大范围检查：不是把后验集中区误当成新的均匀先验。
_,_,dense_means,dense_edge=grid_fit(151)
_,_,wide_means,wide_edge=grid_fit(151,True)
assert np.max(abs(means-dense_means))<1e-5
assert np.max(abs(means-wide_means))<1e-5
assert edge<1e-8
marginals=[weights.sum(axis=tuple(j for j in range(3) if j!=i)) for i in range(3)]
np.savetxt(RES/'Python_grid_means.csv',means[None,:],delimiter=',',header='a,b,sigma',comments='')
np.savez_compressed(RES/'grid_posterior.npz',a=axes[0],b=axes[1],sigma=axes[2],weights=weights)
# 离散网格的等尾分位数；区间端点受网格步长限制。
def wquantile(v,w,p):
    order=np.argsort(v)
    return np.interp(p,np.cumsum(w[order])-.5*w[order],v[order])

# %% 二次近似：求后验峰值，再以峰附近曲率构造协方差。
def objective(t):
    a,b,s=t
    if b<=0 or s<=0 or s>=10: return 1e100
    return n*np.log(s)+sse(a,b)/(2*s*s)+(a-60)**2/200+np.log(b)+np.log(b)**2/2

def gradient(t):
    a,b,s=t
    if b<=0 or s<=0 or s>=10: return np.zeros(3)
    da=n*a+b*sx-sy; db=a*sx+b*sxx-sxy
    return np.array([da/s**2+(a-60)/100,db/s**2+(1+np.log(b))/b,n/s-sse(a,b)/s**3])
fit=minimize(objective,[45,.63,4.2],jac=gradient,method='BFGS',options={'gtol':1e-7})
mode=fit.x
a,b,s=mode
hes=np.array([[n/s**2+.01,sx/s**2,-2*(n*a+b*sx-sy)/s**3],
 [sx/s**2,sxx/s**2-np.log(b)/b**2,-2*(a*sx+b*sxx-sxy)/s**3],
 [-2*(n*a+b*sx-sy)/s**3,-2*(a*sx+b*sxx-sxy)/s**3,-n/s**2+3*sse(a,b)/s**4]])
cov=np.linalg.inv(hes)
assert np.max(abs(gradient(mode)))<1e-4
rfile=RES/'R_quap_parameters.csv'
if rfile.exists():
    r=np.genfromtxt(rfile,delimiter=',',names=True,dtype=None,encoding='utf8')
    assert np.max(abs(r['mean']-mode))<1e-4
    rcov=np.genfromtxt(RES/'R_quap_covariance.csv',delimiter=',',skip_header=1)
    assert np.max(abs(rcov-cov))<1e-5
    # 优先使用已实跑 R quap 的参数和协方差；缺 R 时由上述同原理计算独立生成。
    mode=np.asarray(r['mean'],dtype=float)
    cov=rcov
np.savetxt(RES/'Python_quadratic_parameters.csv',np.column_stack([mode,np.sqrt(np.diag(cov))]),delimiter=',',header='mean,sd',comments='')
# %% 用低差异点减少数值波动。它是积分辅助，不是新的观测数据。
# 网格采用条件逆CDF抽取，保留三个参数间的依赖；不能分别独立抽边际。
u=qmc.Sobol(4,scramble=True,seed=20261005).random_base2(18)
ps=weights.sum(axis=(0,1))
isig=np.searchsorted(np.cumsum(ps),u[:,0])
samples=np.empty((len(u),3))
for k in np.unique(isig):
    rows=np.where(isig==k)[0]
    joint=weights[:,:,k]/ps[k]
    pb=joint.sum(axis=0)
    ib=np.searchsorted(np.cumsum(pb),u[rows,1])
    for j in np.unique(ib):
        rr=rows[ib==j]
        pa=joint[:,j]/pb[j]
        ia=np.searchsorted(np.cumsum(pa),u[rr,2])
        samples[rr,0]=axes[0][np.minimum(ia,len(axes[0])-1)]
        samples[rr,1]=axes[1][j]
        samples[rr,2]=axes[2][k]
qsamples=mode+norm.ppf(u[:,:3])@np.linalg.cholesky(cov).T
assert np.all(qsamples[:,1:]>0)
z=norm.ppf(u[:,3])

def predictions(draws,heights):
    rows=[]
    for H in heights:
        mu=draws[:,0]+draws[:,1]*(H-hc)
        new=mu+draws[:,2]*z
        rows.append([H,mu.mean(),*np.quantile(mu,[.055,.945]),*np.quantile(new,[.055,.945])])
    return np.array(rows)
heights=np.sort(np.unique(np.r_[np.linspace(h.min(),h.max(),81),160]))
grid=predictions(samples,heights)
qgrid=predictions(qsamples,heights)
for name,values in [('grid_predictions',grid),('quap_predictions',qgrid)]:
    np.savetxt(RES/(name+'.csv'),values,delimiter=',',header='height,mean,mean_lower89,mean_upper89,new_lower89,new_upper89',comments='')
summary=[]
for i,name in enumerate(['alpha_centered','beta','sigma']):
    summary.append({'parameter':name,'grid_mean':means[i],'grid_lower89':wquantile(axes[i],marginals[i],.055),'grid_upper89':wquantile(axes[i],marginals[i],.945),'quap_mean':mode[i],'quap_lower89':mode[i]+norm.ppf(.055)*np.sqrt(cov[i,i]),'quap_upper89':mode[i]+norm.ppf(.945)*np.sqrt(cov[i,i])})
import csv
with (RES/'comparison.csv').open('w') as f:
    wr=csv.DictWriter(f,fieldnames=summary[0],lineterminator="\n"); wr.writeheader(); wr.writerows(summary)
half_mu=samples[:len(u)//2,0]+samples[:len(u)//2,1]*(160-hc)
half_new=half_mu+samples[:len(u)//2,2]*z[:len(u)//2]
full160=grid[heights==160][0]
integration_diff=float(np.max(abs(np.r_[np.quantile(half_mu,[.055,.945]),np.quantile(half_new,[.055,.945])]-full160[2:])))
assert integration_diff<.02
checks={'integration_half_full_maxdiff_kg':integration_diff,'n':n,'center':hc,'edge_mass':edge,'dense_mean_maxdiff':float(max(abs(means-dense_means))),'wide_mean_maxdiff':float(max(abs(means-wide_means))),'grid_at160':grid[heights==160][0].tolist(),'quap_at160':qgrid[heights==160][0].tolist()}
(RES/'checks.json').write_text(json.dumps(checks,indent=2))
print(summary); print(checks)
# 图形单独保存：正文与绘图共用以上结果，避免手填数字。
import runpy
runpy.run_path(str(ROOT/'代码/howell_plots_v5.py'),init_globals=globals())
