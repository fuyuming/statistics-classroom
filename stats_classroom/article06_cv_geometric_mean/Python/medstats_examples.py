"""统计课堂 06（附篇）：变异系数与几何均数的 Python 实现。

安装：python -m pip install -r requirements.txt
本文件保存共用函数；请运行同目录中同名的入口文件（2-9-几何均数.py、2-10-变异系数.py）。
数据：data/05.txt（10 名患者来源分离株的环丙沙星 MIC）、data/newborn_2025.csv（10 名新生儿身高体重）。
口径：标准差一律为样本标准差（ddof=1）；四分位数用位置线性插值，与 R 默认 type=7 一致。
"""
from pathlib import Path
import json
try:
    import numpy as np
    import pandas as pd
    from scipy import stats
    import matplotlib.pyplot as plt
except ImportError as exc:
    raise SystemExit('请先安装依赖：python -m pip install -r requirements.txt') from exc
ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data'
OUT = ROOT / 'output' / 'Python'
OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams['font.sans-serif'] = ['Arial Unicode MS', 'PingFang SC', 'Microsoft YaHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

def savefig(name):
    """图形与控制台表格互相核对；保存后关闭，批处理不会堆积窗口。"""
    plt.tight_layout(); plt.savefig(OUT / (name+'.png'), dpi=160); plt.close()

def shape(x):
    """返回DescTools三种算法；m_k分母n，样本方差s²分母n-1。
    所有峰度均为超额峰度（正态基线0），不是普通峰度（基线3）。
    """
    x=np.asarray(x,float); n=len(x); dx=x-x.mean(); m2=np.mean(dx**2)
    if n<4 or m2==0: raise ValueError('形状比较要求n>=4且样本方差>0')
    g1=np.mean(dx**3)/m2**1.5; g2=np.mean(dx**4)/m2**2-3
    return {'g1':g1,'G1':g1*np.sqrt(n*(n-1))/(n-2),'b1':g1*((n-1)/n)**1.5,
            'g2':g2,'G2':((n+1)*g2+6)*(n-1)/((n-2)*(n-3)),'b2':(g2+3)*((n-1)/n)**2-3}

def describe(x):
    """均数=加法中心；SD=个体离散；中位数/IQR=排序摘要。不要自动以P值选指标。"""
    x=np.asarray(x,float)
    if not np.all(np.isfinite(x)): raise ValueError('先说明并处理缺失/非有限值，不能默默改变分母')
    q1,med,q3=np.quantile(x,[.25,.5,.75],method='linear')
    return dict(n=len(x),mean=x.mean(),median=med,sd=x.std(ddof=1),q1=q1,q3=q3,iqr=q3-q1,
                cv_percent=100*x.std(ddof=1)/x.mean(),**shape(x))

def run(example):
    """每个入口只执行对应例题。输出位于 output/Python。"""
    result={}
    if example=='2-9':
        x=pd.read_csv(DATA/'05.txt',sep=r'\s+')['Concentration'].to_numpy();assert np.all(x>0)
        result=describe(x);logx=np.log(x);gm=np.exp(logx.mean());gsd=np.exp(logx.std(ddof=1));result.update(GM=gm,GSD=gsd,GM_over_GSD=gm/gsd,GM_times_GSD=gm*gsd,log_shape=shape(logx))
        # 备查：SciPy 的 D'Agostino 变换与 moments::agostino.test 对应双侧 P 值。
        test=stats.skewtest(x);result['skew_test_backup']={'z':float(test.statistic),'p':float(test.pvalue)}
        fig,axes=plt.subplots(1,2,figsize=(10,4));axes[0].hist(x,bins=10);axes[0].axvline(x.mean(),c='red');axes[0].set_xlabel('MIC (mg/L)');axes[1].hist(logx,bins=8);axes[1].set_xlabel('ln(MIC)');savefig(example)
        # GM/GSD到GM*GSD是倍数范围，不能当成均值置信区间。
        assert np.allclose([x.mean(),np.median(x),gm,gsd],[.225,.0625,.10153155,3.70265996])
    elif example=='2-10':
        d=pd.read_csv(DATA/'newborn_2025.csv');result={key:describe(d[key]) for key in ['height_cm','weight_kg']}
        # CV=SD/mean，只对零点有实质意义且均数为正远离0的变量解释。
        print(pd.DataFrame(result).loc[['n','mean','median','sd','cv_percent']]);plt.bar(['身高','体重'],[result[k]['cv_percent'] for k in result]);plt.ylabel('CV (%)');savefig(example)
        assert np.allclose([result[k]['cv_percent'] for k in result],[8.09682,11.28283])
    else:raise ValueError(example)
    (OUT/(example+'_results.json')).write_text(json.dumps(result,ensure_ascii=False,indent=2,default=lambda v:v.item()),encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,default=lambda v:v.item()))
    return result
