"""统计课堂09 v2：理论概率配图；计算入口为moments_v2.py。"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import font_manager
from scipy.stats import norm

def make_plots(root, curve, pmf, results):
    available = {f.name for f in font_manager.fontManager.ttflist}
    font = next((f for f in ['Microsoft YaHei','PingFang SC','Noto Sans CJK SC','SimHei','Arial Unicode MS'] if f in available),None)
    if font is None:
        raise RuntimeError('计算结果已保存。绘图需要中文字体，请安装思源黑体后重新运行。')
    plt.rcParams.update({'font.family':font,'axes.unicode_minus':False,'font.size':12,
                         'axes.spines.top':False,'axes.spines.right':False,
                         'savefig.facecolor':'white'})
    blue,orange = '#277E91','#D67D38'
    dest=root/'文章配图'
    fig,axes=plt.subplots(2,2,figsize=(11,8.5),layout='constrained')
    for ax, name in zip(axes.ravel(), ['small','larger','rare','balanced']):
        part=curve[curve.case==name]
        n,p=int(part.n.iloc[0]),float(part.p.iloc[0])
        # 稀有事件用局部横轴，保留所有有视觉意义的概率质量。
        se=np.sqrt(p*(1-p)/n)
        ax.bar(part.proportion,part.exact,width=.72/n,color=blue,alpha=.78,label='二项精确概率')
        ax.plot(part.proportion,part.normal_mass,'o-',markersize=3,color=orange,label='正态小区间概率')
        ax.set(xlim=(max(-.6/n,p-4*se),min(1+.6/n,p+4*se)),
               title=f'n = {n}，p = {p:g}；np = {n*p:g}',xlabel='样本率',ylabel='概率')
        ax.grid(axis='y',alpha=.15)
    axes[0,0].legend(frameon=False,fontsize=10)
    fig.suptitle('人数多，不代表稀有事件就适合正态近似',fontsize=18)
    fig.savefig(dest/'v2_02_样本率与正态近似.png',dpi=210)
    curve.to_csv(dest/'v2_02_曲线坐标.csv',index=False,encoding='utf-8-sig')
    fig,axes=plt.subplots(1,2,figsize=(11,4.7),layout='constrained',gridspec_kw={'width_ratios':[1.5,1]})
    d=pmf.difference.to_numpy()
    mass=norm.cdf(d+.0025,results['mean_d'],results['se_d'])-norm.cdf(d-.0025,results['mean_d'],results['se_d'])
    axes[0].vlines(d,0,pmf.probability,color=blue,alpha=.8,label='精确概率')
    axes[0].plot(d,mass,color=orange,lw=2,label='正态小区间概率')
    axes[0].axvline(.2,color='#5b576f',ls='--',label='率差0.2')
    axes[0].set(xlim=(-.25,.45),xlabel='样本率之差',ylabel='概率',title='2091种计数组合，合并成率差分布')
    axes[0].legend(frameon=False,fontsize=10)
    vals=[results[v] for v in ['difference_exact','difference_normal','difference_cc']]
    bars=axes[1].bar(['精确','正态\n未校正','正态\n半格校正'],vals,color=[blue,'#acb3b8',orange],width=.65)
    axes[1].bar_label(bars,labels=[f'{x:.6f}' for x in vals],padding=5,fontsize=11)
    axes[1].set(ylim=(0,.18),ylabel='尾概率',title='率差至少0.2的概率')
    for ax in axes: ax.grid(axis='y',alpha=.15)
    fig.savefig(dest/'v2_03_率差精确分布.png',dpi=210)
    pd.DataFrame(dict(difference=d,exact=pmf.probability,normal_mass=mass)).to_csv(dest/'v2_03_曲线坐标.csv',index=False,encoding='utf-8-sig')
    if plt.get_backend().lower()!='agg': plt.show()
    else: plt.close('all')
