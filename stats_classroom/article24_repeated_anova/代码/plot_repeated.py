"""绘制第24篇数据图。直接运行本文件，或在主分析脚本运行后调用。
数据来自课堂8人×4条件例题。QQ图采用正态理论分位数；直线通过第一、第三四分位点。
"""
from pathlib import Path
import itertools
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties, fontManager
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'文章配图'
OUT.mkdir(exist_ok=True)
font = next((f.name for f in fontManager.ttflist if f.name in ['PingFang SC', 'Noto Sans CJK SC']), 'Arial Unicode MS')
plt.rcParams.update({'font.family': font, 'font.size': 13, 'axes.titlesize': 15,
 'axes.labelsize': 13, 'xtick.labelsize': 12, 'ytick.labelsize': 12,
 'axes.spines.top': False, 'axes.spines.right': False, 'axes.unicode_minus': False,
 'savefig.facecolor':'white'})
x = pd.read_csv(ROOT/'数据/四种声音条件反应率.csv').set_index('subject')[list('ABCD')].to_numpy(float)
colors = ['#4776B4','#D25F46','#789241','#906BB0','#C89832','#2A9092','#B5658C','#647485']

def save(fig, name):
    fig.savefig(OUT/name, dpi=220, bbox_inches='tight', pad_inches=.18)
    plt.close(fig)

fig, axes=plt.subplots(2,1,figsize=(5.5,8.3),layout='constrained')
for ax,y in zip(axes,[x,x-x.mean(1,keepdims=True)]):
    for i in range(8): ax.plot(range(4),y[i],'-o',color=colors[i],lw=1.2,ms=4,alpha=.85,label=f'受试者{i+1}')
    ax.plot(range(4),y.mean(0),'-s',color='#131F2A',lw=2.8,ms=6,label='条件均值')
    ax.set_xticks(range(4),list('ABCD')); ax.set_xlabel('声音条件'); ax.grid(axis='y',alpha=.18)
axes[0].set_title('原始记录：每条细线连接同一个人',loc='left',pad=12)
axes[0].set_ylabel('反应率（%）')
axes[1].set_title('每人减去自己的四次平均值',loc='left',pad=12)
axes[1].set_ylabel('偏离个人平均值（百分点）')
axes[1].axhline(0,color='#aaa',ls='--',lw=1)
axes[0].legend(ncol=3,fontsize=10,loc='upper center',bbox_to_anchor=(.5,1.02),frameon=False)
axes[0].set_ylim(7.3,14.5)
save(fig,'01-个体轨迹与个人中心化.png')

fig,ax=plt.subplots(figsize=(5.5,4.3),layout='constrained')
vals=[32.38,14.3875,10.7875]; labels=['人际平均水平差异','声音条件均值差异','剩余变异']; cs=['#8DAACB','#39A6A3','#DCA960']
ax.barh([2,1,0],vals,color=cs,height=.58)
for yy,v,label in zip([2,1,0],vals,labels):
    ax.text(.3, yy+.32,label,ha='left',va='bottom',fontsize=12)
    ax.text(v+.4,yy,f'{v:g}',va='center',fontsize=12)
ax.set_xlim(0,40);ax.set_yticks([]);ax.set_xlabel('平方和（百分点²）')
ax.set_title('总平方和 57.555 分成三部分',loc='left',pad=14)
ax.spines['left'].set_visible(False)
save(fig,'02-三部分平方和.png')

fig,axes=plt.subplots(6,1,figsize=(5.4,18),layout='constrained')
for ax,(i,j) in zip(axes,itertools.combinations(range(4),2)):
    d=x[:,j]-x[:,i]
    q=stats.norm.ppf((np.arange(1,9)-.375)/(8+.25))
    ax.scatter(q,np.sort(d),s=35,color='#237BA3',zorder=3)
    q1,q3=stats.norm.ppf([.25,.75]); y1,y3=np.quantile(d,[.25,.75]); slope=(y3-y1)/(q3-q1)
    ax.plot(q,y1+slope*(q-q1),color='#C7664D',lw=1.4)
    ax.set_title(f'{"ABCD"[j]}−{"ABCD"[i]}：Shapiro–Wilk P = {stats.shapiro(d).pvalue:.3f}',loc='left',fontsize=14)
    ax.set_xlabel('标准正态理论分位数');ax.set_ylabel('差值（百分点）');ax.grid(alpha=.15)
save(fig,'03-六对差值QQ图.png')

fig,ax=plt.subplots(figsize=(5.5,5.1),layout='constrained')
pairs=list(itertools.combinations(range(4),2))
for yy,(i,j) in enumerate(pairs[::-1]):
    d=x[:,j]-x[:,i]; m=d.mean(); half=stats.t.ppf(1-.05/12,7)*stats.sem(d)
    ax.errorbar(m,yy,xerr=half,fmt='o',color='#207FA5',capsize=4,lw=2,ms=6)
ax.axvline(0,color='#999',ls='--',lw=1)
ax.set_yticks(range(6),[f'{"ABCD"[j]}−{"ABCD"[i]}' for i,j in pairs[::-1]])
ax.set_xlabel('平均反应率之差（百分点）')
ax.set_title('六对均值差与同时置信区间',loc='left',pad=30,fontsize=15)
ax.text(0,1.02,'Bonferroni；整组置信水平≥95%',transform=ax.transAxes,fontsize=12)
ax.grid(axis='x',alpha=.18)
save(fig,'04-配对均值差与同时区间.png')
print('已生成4张教学数据图。')
