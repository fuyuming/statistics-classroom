"""统计课堂08：理论概率配图；计算入口为expectation_variance_covariance.py。"""
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
    blue,orange = '#56b4e9','#D67D38'
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


from scipy.stats import chi2

def make_reference_plots(root, q_sim, z):
    fonts = {f.name for f in font_manager.fontManager.ttflist}
    available = [x for x in ['Microsoft YaHei','PingFang SC','Noto Sans CJK SC','SimHei','Arial Unicode MS'] if x in fonts]
    if not available:
        raise RuntimeError('计算结果已保存。绘图请先安装中文字体（如Noto Sans CJK SC）并重启Python。')
    plt.rcParams.update({'font.family':available[0], 'axes.unicode_minus':False,
                         'font.size':13,'axes.spines.top':False,'axes.spines.right':False,
                         'figure.facecolor':'#fbfaf6','axes.facecolor':'#fbfaf6'})
    folder = root/'文章配图'
    folder.mkdir(exist_ok=True)
    rho = np.linspace(-1,1,201)
    variance = 2-2*rho
    fig,ax=plt.subplots(figsize=(9,5.3),layout='constrained')
    ax.plot(rho,variance,color='#56b4e9',lw=2.5)
    ax.scatter([0,.8],[2,.4],color='#bb6a38',s=45,zorder=3)
    ax.annotate('不相关：差的方差为2',(0,2),xytext=(-.4,3.0),bbox={'facecolor':'#fbfaf6','edgecolor':'none','pad':3},arrowprops={'arrowstyle':'->'})
    ax.annotate('正相关：共同波动在差中抵消',(.8,.4),xytext=(-.2,.9),arrowprops={'arrowstyle':'->'})
    ax.set(xlabel='相关系数 ρ',ylabel='Var(X − Y)',title='两个量相减，波动有多大？\n固定 Var(X) = Var(Y) = 1')
    fig.savefig(folder/'01_协方差与差值.png',dpi=180)
    pd.DataFrame({'rho':rho,'variance_of_difference':variance}).to_csv(folder/'01_curve.csv',index=False)
    fig,axes=plt.subplots(2,1,figsize=(9,8),layout='constrained')
    axes[0].bar([r'$Z_1^2$',r'$Z_2^2$',r'$Z_3^2$'],z**2,color=['#56b4e9','#5264a7','#bb6a38'])
    for i,value in enumerate(z**2):axes[0].text(i,value+.07,f'{value:.1f}',ha='center')
    axes[0].set(ylabel='该方向的平方',ylim=(0,3.8),title='三个比较方向：0.4 + 0.4 + 3.2 = Q = 4')
    x=np.linspace(.001,20,501)
    axes[1].hist(q_sim,bins=np.arange(0,np.ceil(q_sim.max())+1,.5),density=True,color='#aedcf5',label='多项模型模拟（20,000次）')
    axes[1].plot(x,chi2.pdf(x,3),color='#5264a7',lw=2,label='χ²(3) 理论密度')
    axes[1].set(xlabel='Pearson Q',ylabel='密度',xlim=(0,20),title='零假设：n = 80，四类概率各0.25')
    axes[1].legend(fontsize=11)
    fig.savefig(folder/'03_Q与三个方向.png',dpi=180)
    pd.DataFrame({'Q':x,'chisq3_density':chi2.pdf(x,3)}).to_csv(folder/'03_curve.csv',index=False)
    # 在Spyder中显示；命令行验证可在外部设置MPLBACKEND=Agg。
    if plt.get_backend().lower() != 'agg':
        plt.show()
