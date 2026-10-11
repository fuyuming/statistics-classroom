"""配套排版与作图。主脚本自动调用；只用已计算结果和理论曲线。"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import font_manager
from scipy.stats import chi2

def make_plots(root, differences, q_sim, z, center, se):
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
    ax.plot(rho,variance,color='#19837b',lw=2.5)
    ax.scatter([0,.8],[2,.4],color='#bb6a38',s=45,zorder=3)
    ax.annotate('不相关：差的方差为2',(0,2),xytext=(-.4,3.0),bbox={'facecolor':'#fbfaf6','edgecolor':'none','pad':3},arrowprops={'arrowstyle':'->'})
    ax.annotate('正相关：共同波动在差中抵消',(.8,.4),xytext=(-.2,.9),arrowprops={'arrowstyle':'->'})
    ax.set(xlabel='相关系数 ρ',ylabel='Var(X − Y)',title='两个量相减，波动有多大？\n固定 Var(X) = Var(Y) = 1')
    fig.savefig(folder/'01_协方差与差值.png',dpi=180)
    pd.DataFrame({'rho':rho,'variance_of_difference':variance}).to_csv(folder/'01_curve.csv',index=False)
    fig,ax=plt.subplots(figsize=(9,5.3),layout='constrained')
    ax.hist(differences,bins=np.linspace(-.5,.6,45),density=True,color='#75b5ad',edgecolor='white')
    ax.axvline(center,color='#b66c3e',ls='--',lw=2,label=f'理论中心 {center:.2f}')
    ax.axvline(center-se,color='#5264a7',ls=':',label=f'中心 ± 1 SE；SE = {se:.5f}')
    ax.axvline(center+se,color='#5264a7',ls=':')
    ax.set(xlabel='样本率之差（比例单位）',ylabel='模拟密度',title='率差是统计量，每次抽样都会变化\nn1 = 40，p1 = 0.3；n2 = 50，p2 = 0.2；两组独立')
    ax.legend(fontsize=11)
    fig.savefig(folder/'02_率差的抽样波动.png',dpi=180)
    fig,axes=plt.subplots(2,1,figsize=(9,8),layout='constrained')
    axes[0].bar([r'$Z_1^2$',r'$Z_2^2$',r'$Z_3^2$'],z**2,color=['#19837b','#5264a7','#bb6a38'])
    for i,value in enumerate(z**2):axes[0].text(i,value+.07,f'{value:.1f}',ha='center')
    axes[0].set(ylabel='该方向的平方',ylim=(0,3.8),title='三个比较方向：0.4 + 0.4 + 3.2 = Q = 4')
    x=np.linspace(.001,20,501)
    axes[1].hist(q_sim,bins=np.arange(0,np.ceil(q_sim.max())+1,.5),density=True,color='#a7cdc7',label='多项模型模拟（20,000次）')
    axes[1].plot(x,chi2.pdf(x,3),color='#5264a7',lw=2,label='χ²(3) 理论密度')
    axes[1].set(xlabel='Pearson Q',ylabel='密度',xlim=(0,20),title='零假设：n = 80，四类概率各0.25')
    axes[1].legend(fontsize=11)
    fig.savefig(folder/'03_Q与三个方向.png',dpi=180)
    pd.DataFrame({'Q':x,'chisq3_density':chi2.pdf(x,3)}).to_csv(folder/'03_curve.csv',index=False)
    # 在Spyder中显示；命令行验证可在外部设置MPLBACKEND=Agg。
    if plt.get_backend().lower() != 'agg':
        plt.show()
