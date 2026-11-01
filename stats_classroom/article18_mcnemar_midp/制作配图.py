"""依据三语言配套R的真实枚举输出绘制条件二项图；无随机模拟。
依赖matplotlib、numpy、pandas；先运行配套R，再 python 制作配图.py。
"""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Rectangle,Patch
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'文章配图'
OUT.mkdir(exist_ok=True)
data=pd.read_csv(ROOT/'运行结果/R/all_outcomes.csv')
assert (data.m==14).all()
fonts={f.name for f in font_manager.fontManager.ttflist}
for name in ['PingFang SC','Heiti TC','Arial Unicode MS','Noto Sans CJK SC']:
    if name in fonts:
        plt.rcParams['font.sans-serif']=[name]
        break
else:raise RuntimeError('需要中文字体')
plt.rcParams['axes.unicode_minus']=False
plt.rcParams['font.size']=12
plt.rcParams['text.color']='#27343a'
BLUE='#1769aa';ORANGE='#c95516';GOLD='#daa327';PALE='#c5d6e5'
k=data.k.to_numpy();prob=data.prob.to_numpy()
tail=(k<=2)|(k>=12)
boundary=(k==2)|(k==12)
colors=np.where(boundary,GOLD,np.where(tail,ORANGE,PALE))
def style(ax):
    ax.set_axisbelow(True);ax.grid(axis='y',color='#e3e9ee',lw=.7)
    for side in ['top','right']:ax.spines[side].set_visible(False)
    for side in ['left','bottom']:ax.spines[side].set_color('#9cabb8')
    ax.tick_params(labelsize=12)
# 图1: 全分布与尾部局部，用不同刻度但明确标注。
fig,axes=plt.subplots(2,1,figsize=(5.2,7.4),dpi=240)
axes[0].bar(k,prob,color=colors,width=.78)
axes[0].set_xticks(range(0,15,2));axes[0].set_xlim(-.65,14.65)
axes[0].set_ylabel('概率');axes[0].set_ylim(0,.24)
axes[0].set_xlabel('B：仅荧光阳性的不一致人数')
axes[0].set_title('全部15种结果：中心在7',loc='left',fontsize=14,pad=12)
axes[0].legend(handles=[Patch(facecolor=ORANGE,label='更极端'),Patch(facecolor=GOLD,label='同样极端')],frameon=False,ncol=2,fontsize=11,loc='upper right')
kt=k[tail];pt=prob[tail];positions=np.arange(len(kt))
axes[1].bar(positions,pt,color=colors[tail],width=.65)
axes[1].set_xticks(positions,[str(x) for x in kt]);axes[1].set_ylabel('概率（放大刻度）')
axes[1].set_ylim(0,.0067);axes[1].set_yticks([0,.002,.004,.006])
axes[1].set_xlabel('仅显示双尾的六个取值')
axes[1].set_title('双尾放大：2与12是边界',loc='left',fontsize=14,pad=12)
for ax in axes:style(ax)
fig.text(.14,.975,'给定 m = 14，H$_0$：方向概率各为1/2',fontsize=14,va='top')
fig.subplots_adjust(left=.17,right=.96,top=.87,bottom=.10,hspace=.65)
fig.savefig(OUT/'01-二项概率与双尾.png',facecolor='white');plt.close(fig)
data[['m','k','prob']].assign(in_tail=tail,is_boundary=boundary).to_csv(OUT/'01-理论概率.csv',index=False)
# 图2: 两种权重并列，虚框只画midP未计入的半根边界。
weighted=pt.copy();weighted[boundary[tail]]*=.5
fig,ax=plt.subplots(figsize=(5.2,4.8),dpi=240)
w=.34
ax.bar(positions-w/2,pt,width=w,color=BLUE,label='精确：边界全计')
ax.bar(positions+w/2,weighted,width=w,color=ORANGE,label='mid-P：边界取半')
for i in np.where(boundary[tail])[0]:
    ax.add_patch(Rectangle((positions[i],weighted[i]),w,pt[i]-weighted[i],facecolor='#fff9eb',edgecolor=GOLD,linestyle='--',linewidth=1.5))
ax.set_xticks(positions,[str(x) for x in kt]);ax.set_xlabel('B的尾部取值（省略3至11）',labelpad=10)
ax.set_ylabel('计入P值的概率');ax.set_ylim(0,.0065)
ax.set_yticks([0,.002,.004,.006]);style(ax)
ax.set_title('边界各取一半，其余尾部不变',loc='left',fontsize=15,pad=55)
ax.legend(frameon=False,loc='lower left',bbox_to_anchor=(-.02,1.01),fontsize=12)
fig.subplots_adjust(left=.19,right=.96,top=.71,bottom=.17)
fig.savefig(OUT/'02-边界概率的一半.png',facecolor='white');plt.close(fig)
pd.DataFrame({'B':kt,'exact_weighted_probability':pt,'midP_weighted_probability':weighted}).to_csv(OUT/'02-边界权重.csv',index=False)
# 图3: 离散取值点，连线只作阅读辅助，不补造非整数数据。
plot=data[data.k>=10]
fig,ax=plt.subplots(figsize=(5.2,4.8),dpi=240)
ax.plot(plot.k,plot.exact_P,'o-',color=BLUE,label='条件精确P',lw=1.7,markersize=6)
ax.plot(plot.k,plot.midP,'s--',color=ORANGE,label='mid-P',lw=1.5,markersize=5)
ax.axhline(.05,color='#67727d',ls=':',lw=1.3)
ax.text(13.35,.055,'α = 0.05',fontsize=11,color='#67727d')
ax.annotate('0.05737',xy=(11,.057373046875),xytext=(11.5,.10),fontsize=12,color=BLUE,arrowprops={'arrowstyle':'->','color':BLUE})
ax.annotate('0.03516',xy=(11,.03515625),xytext=(11.65,.037),fontsize=12,color=ORANGE,arrowprops={'arrowstyle':'->','color':ORANGE})
ax.text(12.05,.016,'0.01294',color=BLUE,fontsize=11)
ax.set_xticks(range(10,15));ax.set_xlim(9.8,14.35);ax.set_ylim(0,.205)
ax.set_yticks([0,.05,.10,.15,.20]);ax.set_ylabel('双侧P值');ax.set_xlabel('b：仅荧光阳性的不一致人数',labelpad=10)
ax.set_title('固定 m = 14，P值跨过0.05',loc='left',fontsize=15,pad=45)
ax.legend(loc='lower left',bbox_to_anchor=(-.02,1.01),ncol=2,frameon=False,fontsize=12)
style(ax);fig.subplots_adjust(left=.16,right=.96,top=.75,bottom=.17)
fig.savefig(OUT/'03-P值跳跃与阈值.png',facecolor='white');plt.close(fig)
plot[['m','k','exact_P','midP']].to_csv(OUT/'03-P值坐标.csv',index=False)
print('3张配图与坐标已生成，均为给定m=14的理论量。')
