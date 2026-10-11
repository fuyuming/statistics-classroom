"""中文配图由主文件调用；图1为名次示意、图2为精确概率、图3为模拟密度。"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import font_manager
from scipy.stats import norm

def make_plots(root,data,pmf,z_sim,m,n):
    fonts={f.name for f in font_manager.fontManager.ttflist}
    found=[x for x in ['Microsoft YaHei','PingFang SC','Noto Sans CJK SC','SimHei','Arial Unicode MS'] if x in fonts]
    if not found:
        raise RuntimeError('计算已保存。请安装中文字体Noto Sans CJK SC等，重启Python再绘图。')
    plt.rcParams.update({'font.family':found[0],'font.size':13,'axes.unicode_minus':False,
       'figure.facecolor':'#fbfaf6','axes.facecolor':'#fbfaf6','axes.spines.top':False,'axes.spines.right':False})
    folder=root/'文章配图'
    folder.mkdir(exist_ok=True)
    ordered=data.sort_values('value').reset_index(drop=True)
    fig,ax=plt.subplots(figsize=(9,4.8),layout='constrained')
    for i,row in ordered.iterrows():
        color='#19837b' if row.group=='A' else '#ba723f'
        ax.scatter(i+1,1,s=800,color=color)
        ax.text(i+1,1,row.group,color='white',ha='center',va='center',fontsize=17)
        ax.text(i+1,1.37,f'数值 {row.value}',ha='center')
        ax.text(i+1,.64,f'秩 {i+1}',ha='center')
    ax.text(3.5,.05,'A组秩和 W = 1 + 2 + 4 = 7',ha='center',fontsize=20,color='#19837b')
    ax.set(xlim=(.4,6.6),ylim=(-.3,1.8),title='合并排序，才有可比较的名次')
    ax.axis('off')
    fig.savefig(folder/'01_排序与秩和.png',dpi=180)
    fig,ax=plt.subplots(figsize=(9,5.3),layout='constrained')
    colors=['#ba723f' if x<=7 or x>=14 else '#19837b' for x in pmf.W]
    ax.bar(pmf.W,pmf.probability,color=colors,width=.7)
    for row in pmf.itertuples():
        ax.text(row.W,row.probability+.004,f'{row.count}/20',ha='center',fontsize=11)
    ax.axvline(10.5,ls='--',color='#5363aa',label='零假设中心10.5')
    ax.set(xticks=np.arange(6,16),ylim=(0,.195),xlabel='A组原始秩和 W',ylabel='精确概率',title='20种等可能分配，10个可能秩和\n橙色：与W = 7至少同样偏离中心，双侧P = 4/20')
    ax.legend(fontsize=11)
    fig.savefig(folder/'02_精确分布.png',dpi=180)
    fig,ax=plt.subplots(figsize=(9,5.3),layout='constrained')
    ax.hist(z_sim,bins=np.arange(-6,6.2,.2),density=True,color='#a6cdc6',label='20,000轮秩分配模拟')
    x=np.linspace(-4.5,4.5,501)
    ax.plot(x,norm.pdf(x),color='#5363aa',lw=2,label='标准正态理论密度')
    ax.set(xlim=(-4.5,4.5),xlabel='Z = [W − E(W)] / SD(W)',ylabel='密度',title=f'样本量增大后，秩和怎样接近正态？\nm = {m}，n = {n}；无并列值的零假设')
    ax.legend(fontsize=11)
    fig.savefig(folder/'03_正态近似.png',dpi=180)
    pd.DataFrame({'Z':x,'normal_density':norm.pdf(x)}).to_csv(folder/'03_curve.csv',index=False)
    if plt.get_backend().lower()!='agg':
        plt.show()
