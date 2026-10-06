"""精读04：作者分类模型和年龄—身高样条，Python二次近似复算。
Spyder F5；依赖numpy、scipy、matplotlib。自动画六张中文图，保留GUI并保存PNG。
输入数据与R生成的样条基函数表已随项目提供，无需先运行R。
"""
from pathlib import Path
import numpy as np
from scipy.optimize import minimize, minimize_scalar
from scipy.stats import norm, qmc
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'运行结果/v3'; OUT.mkdir(exist_ok=True,parents=True)
d=np.genfromtxt(ROOT/'数据/v3/Howell1_adults.csv',names=True,delimiter=',')
h=d['height']; y=d['weight']; sex=d['male'].astype(int); hc=h.mean()
print('人数：女性',sum(sex==0),'男性',sum(sex==1),'平均身高',hc)
G=np.column_stack([sex==0,sex==1]).astype(float)
X2=np.column_stack([G,G*(h-hc)[:,None]])
# 正态似然、先验与作者完全相同；sigma是未知参数。
def fit_category(X):
    k=X.shape[1]; n=len(y)
    def obj(t):
        b=t[:k]; s=t[-1]
        if s<=0 or s>=10 or (k==4 and np.any((b[2:]<=0)|(b[2:]>=1))): return 1e100
        residual=X@b-y
        val=n*np.log(s)+(residual@residual)/(2*s*s)+np.sum((b[:2]-60)**2)/200
        # PPT第59页斜率在(0,1)内先验密度为常数。
        return val
    def grad(t):
        b=t[:k]; s=t[-1]
        if s<=0 or s>=10 or (k==4 and np.any((b[2:]<=0)|(b[2:]>=1))): return np.zeros(k+1)
        residual=X@b-y
        gb=X.T@residual/s**2
        gb[:2]+=(b[:2]-60)/100

        return np.r_[gb,n/s-(residual@residual)/s**3]
    start=[42,49,5.5] if k==2 else [45,45,.65,.61,4.2]
    f=minimize(obj,start,jac=grad,method='BFGS',options={'gtol':1e-6})
    b=f.x[:k]; s=f.x[-1]; residual=X@b-y
    prior=np.r_[.01,.01] if k==2 else np.r_[.01,.01,0.,0.]
    H=np.zeros((k+1,k+1)); H[:k,:k]=X.T@X/s**2+np.diag(prior)
    H[:k,k]=-2*(X.T@residual)/s**3; H[k,:k]=H[:k,k]
    H[k,k]=-n/s**2+3*(residual@residual)/s**4
    assert max(abs(grad(f.x)))<.002
    return f.x,np.linalg.inv(H)
m1,v1=fit_category(G); m2,v3=fit_category(X2)
# 样条：共享basis是作者R splines::bs生成；不同语言拟合同一设计矩阵。
ch=np.genfromtxt(ROOT/'数据/v3/age_height.csv',names=True,delimiter=',')
B=np.loadtxt(ROOT/'数据/v3/age_basis.csv',delimiter=',',skiprows=1)
bpred=np.loadtxt(ROOT/'数据/v3/age_basis_prediction.csv',delimiter=',',skiprows=1)
year=bpred[:,0]; BP=bpred[:,1:]
X=np.column_stack([np.ones(len(B)),B]); Y=ch['height']; k=X.shape[1]
pmean=np.r_[120,np.zeros(k-1)]; precision=np.r_[1,np.full(k-1,1/625)]
# 给定log_sigma，正态先验下系数峰值可直接求解；外层搜索log_sigma。
def beta_at(l):
    v=np.exp(-2*l)
    return np.linalg.solve(v*(X.T@X)+np.diag(precision),v*(X.T@Y)+precision*pmean)
def profile(l):
    b=beta_at(l); res=X@b-Y
    return len(Y)*l+(res@res)*np.exp(-2*l)/2+np.sum(precision*(b-pmean)**2)/2+2*l*l
f=minimize_scalar(profile,bounds=(-1,4),method='bounded',options={'xatol':1e-11})
l=f.x; b=beta_at(l); res=X@b-Y; v=np.exp(-2*l)
H=np.zeros((k+1,k+1)); H[:k,:k]=v*X.T@X+np.diag(precision)
H[:k,k]=-2*v*(X.T@res); H[k,:k]=H[:k,k]; H[k,k]=2*v*(res@res)+4
ms=np.r_[b,l]; vs=np.linalg.inv(H)
for name,m,V in [('m1',m1,v1),('m2',m2,v3),('ms',ms,vs)]:
    np.savetxt(OUT/f'Python_{name}.csv',np.column_stack([m,np.sqrt(np.diag(V))]),delimiter=',',header='mean,sd',comments='')
    ref=np.genfromtxt(OUT/f'R_{name}.csv',delimiter=',',skip_header=1,usecols=(1,2))
    # R quap的数值Hessian与本脚本解析Hessian允许小幅数值差异。
    assert np.max(abs(m-ref[:,0]))<.02,(name,np.max(abs(m-ref[:,0])))
# 两个不同问题：均值之差与独立抽取两个新个体之差。
z=norm.ppf(.945)
c=np.array([-1,1,0]); delta=c@m1; sd=np.sqrt(c@v1@c)
u=qmc.Sobol(5,scramble=True,seed=20261006).random_base2(17)
q=norm.ppf(u)
draw=m1+q[:,:3]@np.linalg.cholesky(v1).T
individual=draw[:,1]-draw[:,0]+draw[:,2]*(q[:,4]-q[:,3])
c160=np.array([1,-1,160-hc,-(160-hc),0])
d160=c160@m2; sd160=np.sqrt(c160@v3@c160)
summary=np.array([[delta,delta-z*sd,delta+z*sd,norm.cdf(delta/sd)],
 [individual.mean(),*np.quantile(individual,[.055,.945]),np.mean(individual>0)],
 [d160,d160-z*sd160,d160+z*sd160,norm.cdf(d160/sd160)]])
np.savetxt(OUT/'contrast_summary.csv',summary,delimiter=',',header='mean,lower89,upper89,prob_positive',comments='')
XP=np.column_stack([np.ones(len(BP)),BP]); curve=XP@ms[:-1]
curve_sd=np.sqrt(np.einsum('ij,jk,ik->i',XP,vs[:-1,:-1],XP))
np.savetxt(OUT/'age_curve.csv',np.column_stack([year,curve,curve-z*curve_sd,curve+z*curve_sd]),delimiter=',',header='age,mean,lower89,upper89',comments='')
print('三行分别为均值差、新个体差、160厘米处女性减男性均值差：\n',summary)
print('年龄—身高记录数',len(Y),'年龄范围',min(ch['age']),max(ch['age']),'basis列数',B.shape[1])
import runpy
runpy.run_path(str(ROOT/'代码/lesson04_plots_v3.py'),init_globals=globals())
