"""重建本文默认题目的科学配图。改参练习请使用配套三语言脚本。"""
from pathlib import Path
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import numpy as np
from statsmodels.stats.proportion import confint_proportions_2indep

root=Path(__file__).resolve().parents[1]
out=root/'文章配图'
out.mkdir(exist_ok=True)
from matplotlib import font_manager
candidates=['Hiragino Sans GB','PingFang SC','Microsoft YaHei','SimHei','Noto Sans CJK SC','WenQuanYi Zen Hei']
available={f.name for f in font_manager.fontManager.ttflist}
chosen=next((f for f in candidates if f in available),None)
if chosen is None:
    raise RuntimeError('请先安装一种中文字体：Noto Sans CJK SC / 思源黑体等，再运行。')
font=FontProperties(family=chosen)
plt.rcParams.update({'font.family':font.get_name(),'axes.unicode_minus':False,'mathtext.fontset':'dejavusans'})
p=root/'配套项目/数据/independent_counts.csv'
if p.exists():
    with p.open(encoding='utf-8-sig') as f: rows=list(csv.DictReader(f))
    xA,nA=int(rows[0]['positive']),int(rows[0]['total'])
    xB,nB=int(rows[1]['positive']),int(rows[1]['total'])
else:
    xA,nA,xB,nB=9,20,13,20
assert (xA,nA,xB,nB)==(9,20,13,20), '此入口重建默认文章配图；改参请用配套脚本'
confidence=.95
blue='#0000F5'
# 图1：只表示两批独立的人；不使用连线，排成各自独立的计数框。
fig,axs=plt.subplots(2,1,figsize=(7.2,5.6),dpi=150)
for ax,label,x,n in zip(axs,['A','B'],[xA,xB],[nA,nB]):
    ax.set_xlim(-.7,10.1);ax.set_ylim(-1,2.1);ax.axis('off')
    ax.text(-.5,1.8,f'{label}组  {x}/{n} 阳性（{x/n:.0%}）',fontproperties=font,fontsize=19,color=blue,weight='bold')
    for i in range(n):
        ax.scatter(i%10, .75-(i//10)*.8,s=220,facecolor=blue if i<x else '#DDE3EE',edgecolor='white',linewidth=1)
fig.text(.06,.04,'独立抽样 · 两组共40人 · 每人检测一次',fontproperties=font,fontsize=17,color='#243449')
fig.subplots_adjust(left=.05,right=.96,bottom=.11,top=.98,hspace=.18)
fig.savefig(out/'01_两批独立样本.png',facecolor='white')
plt.close(fig)
# 核心公式：固定1080px宽、相同主要数学字号；不用tight裁切导致手机缩放不一。
fig=plt.figure(figsize=(9,2.6),dpi=120)
fig.patch.set_facecolor('white')
fig.text(.5,.69,r'$\widehat{SE}=\sqrt{\frac{\hat p_A(1-\hat p_A)}{n_A}+\frac{\hat p_B(1-\hat p_B)}{n_B}}$',ha='center',va='center',fontsize=32,color='#172A44')
fig.text(.5,.22,r'$d\;\pm\;1.96\,\widehat{SE}$',ha='center',va='center',fontsize=32,color='#172A44')
fig.savefig(out/'公式01_Wald.png',facecolor='white')
plt.close(fig)
# 图2：端点全部由官方Python函数计算，原始点估计与AC区间中心分开。
methods=['wald','newcomb','agresti-caffo']
labels=['Wald','Newcombe','Agresti–Caffo']
ci=[confint_proportions_2indep(xA,nA,xB,nB,method=m,compare='diff',alpha=1-confidence) for m in methods]
d=xA/nA-xB/nB
adjusted=(xA+1)/(nA+2)-(xB+1)/(nB+2)
fig,ax=plt.subplots(figsize=(7.2,5.9),dpi=150)
for y,interval,label in zip([2,1,0],ci,labels):
    lo,hi=np.array(interval)*100
    ax.plot([lo,hi],[y,y],color=blue,lw=4,solid_capstyle='round')
    ax.scatter([d*100],[y],color=blue,s=58,zorder=3)
    ax.text((lo+hi)/2,y+.28,f'{label}    [{lo:.2f}, {hi:.2f}]',ha='center',fontsize=17,color='#172A44')
ax.scatter([adjusted*100],[0],marker='D',s=58,color='#F29F39',zorder=4)
ax.axvline(0,color='#6D7D90',linestyle='--',lw=1.2)
ax.set_xlim(-57,18);ax.set_ylim(-.65,2.85)
ax.set_xticks([-50,-40,-30,-20,-10,0,10])
ax.tick_params(axis='x',labelsize=14)
ax.set_yticks([])
ax.set_xlabel('A减B的总体率差（百分点）',fontproperties=font,fontsize=18,labelpad=12)
ax.set_title(f'三种{confidence:.0%}区间 · A {xA}/{nA}，B {xB}/{nB}',fontproperties=font,fontsize=21,pad=16,color=blue)
for spine in ['left','right','top']:ax.spines[spine].set_visible(False)
ax.grid(axis='x',color='#EDF0F7',zorder=0)
fig.text(.10,.045,'● 原始样本率差    ◆ 调整后的区间中心',fontproperties=font,fontsize=16,color='#243449')
fig.subplots_adjust(left=.08,right=.96,bottom=.20,top=.86)
fig.savefig(out/'02_三种区间比较.png',facecolor='white')
plt.close(fig)
with (out/'区间图坐标.csv').open('w',encoding='utf-8',newline='') as f:
    w=csv.writer(f);w.writerow(['method','raw_difference','lower','upper','interval_center','confidence'])
    for m,v in zip(methods,ci):w.writerow([m,d,*v,adjusted if m=='agresti-caffo' else (sum(v)/2),confidence])
print('图1、图2、式1及图2坐标已重建')
