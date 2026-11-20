# 统计课堂34：配对κ Bootstrap、百分位区间与BCa区间
# 背景：第11课例4的58名疑似系统性红斑狼疮患者，各接受免疫荧光法A与乳胶凝集法B检测抗核抗体。
# 题源：颜艳、王彤《医学统计学》第5版（2020）例7-3、表7-3，教材pp110–111；四格为[11,12;2,33]。
# 目标是给定样本条件下的κ波动，不把Monte Carlo区间称为严格精确区间。
# 共享bootstrap_counts.csv每行是58个联合结果类别的重抽样计数，不是虚构病历。
# %% 02.1 导入与输入（不依赖__file__，请在项目根目录按分节执行）
# n是58名患者，B是50000次共享联合计数重抽样，seed记录参考生成条件。
import os, sys
import numpy as np
import pandas as pd
from scipy.stats import norm
sys.path.insert(0, "helpers")
from kappa_bootstrap_python import kappa_from_counts, norm_inter_r, bca_levels
raw = pd.read_csv("data/paired_tests.csv")
counts = raw["count"].to_numpy(dtype=np.int64); n = int(counts.sum()); B=50000; seed=20261120
assert np.array_equal(counts, [11,12,2,33]) and n == 58
shared = pd.read_csv("data/bootstrap_counts.csv").to_numpy(dtype=np.int64)
assert shared.shape == (B,4) and np.all(shared.sum(axis=1)==n)

# %% 02.2 手动计算：原样本κ、bootstrap κ*及基本统计量
# shared每行是四种联合结果的计数[a,b,c,d]，不是一行病历；以下显式展开每行κ公式。
a,b,c,d = counts
kappa_hat = 2*(a*d-b*c)/((a+b)*(b+d)+(c+d)*(a+c))
sa,sb,sc,sd = (shared[:,0],shared[:,1],shared[:,2],shared[:,3])
den = (sa+sb)*(sb+sd)+(sc+sd)*(sa+sc)
kappa_star = 2*(sa*sd-sb*sc)/den
invalid = int(np.sum(~np.isfinite(kappa_star)))
if invalid != 0: raise ValueError("本例出现不可定义样本，不能静默删除")
bootstrap_mean=float(kappa_star.mean()); bootstrap_bias=bootstrap_mean-kappa_hat
bootstrap_variance=float(kappa_star.var(ddof=1)); bootstrap_se=float(kappa_star.std(ddof=1))
print(pd.DataFrame([dict(a=counts[0],b=counts[1],c=counts[2],d=counts[3],n=n,kappa=kappa_hat)]))

# %% 02.3 手动计算：58行delete-one jackknife、z0、acceleration与R norm.inter端点
# cell有58个类别标签；np.delete删去第i名患者，np.bincount恢复四格计数。
cell=np.repeat(np.arange(4), counts)
jack=np.array([kappa_from_counts(np.bincount(np.delete(cell,i), minlength=4)) for i in range(n)])
if np.any(~np.isfinite(jack)): raise ValueError("jackknife出现不可定义样本")
L=(n-1)*(jack.mean()-jack)
bias_probability=float(np.mean(kappa_star < kappa_hat)) # 严格小于
z0=float(norm.ppf(bias_probability)); acc=float(np.sum(L**3)/(6*np.sum(L**2)**1.5))
alpha_bca=bca_levels(z0,acc)
print(dict(bias_probability=bias_probability,z0=z0,acceleration=acc,alpha_bca=alpha_bca))

# %% 02.4 调用函数：自编κ helper
# helper的kappa_from_counts接收[a,b,c,d]，返回单个κ；这里逐行调用核对κ*。
assert np.allclose(kappa_star, np.array([kappa_from_counts(row) for row in shared]))

# %% 02.5 调用函数：自编norm.inter与BCa helper
# norm_inter_r按R boot::norm.inter的正态分位尺度插值，不使用NumPy默认线性分位插值。
percentile=norm_inter_r(kappa_star,[.025,.975])
bca=norm_inter_r(kappa_star,alpha_bca)
assert np.allclose(norm_inter_r(kappa_star,[.025,.975]), percentile, rtol=0, atol=0)

# %% 02.6 结果核对：宽列与可选保存
result=pd.DataFrame([dict(B=B,seed=seed,invalid=invalid,kappa=kappa_hat,
  bootstrap_mean=bootstrap_mean,bootstrap_bias=bootstrap_bias,bootstrap_variance=bootstrap_variance,
  bootstrap_se=bootstrap_se,bias_probability=bias_probability,z0=z0,acceleration=acc,
  alpha_low=alpha_bca[0],alpha_high=alpha_bca[1],percentile_low=percentile[0],
  percentile_high=percentile[1],bca_low=bca[0],bca_high=bca[1])])
print(result.to_string(index=False, float_format=lambda x:f"{x:.14g}"))
if os.environ.get("ARTICLE34_SAVE_RESULTS")=="1":
    os.makedirs("results", exist_ok=True); result.to_csv("results/Python_bootstrap.csv",index=False)
# bootstrap均值/偏差是当前样本条件重抽样的描述，不是总体无偏性证明；BCa调整分位点而非加伪计数。
