"""由lesson04_v3.py调用。保留六张中文Figure，并输出公众号配图。"""
import matplotlib.pyplot as plt
from matplotlib import font_manager
FONTS={f.name for f in font_manager.fontManager.ttflist}
font=next((f for f in ['PingFang SC','Microsoft YaHei','SimHei','Noto Sans CJK SC','Heiti SC'] if f in FONTS),None)
if font is None: raise RuntimeError('请安装 Noto Sans CJK SC 中文字体并重启内核。')
plt.rcParams.update({'font.family':font,'axes.unicode_minus':False,'font.size':13})
FIG=ROOT/'文章配图_v3';FIG.mkdir(exist_ok=True)
TEAL='#167e83'; ORANGE='#ca6841'; INK='#193c46'; BG='#faf8f3'
def setup(num,title,sub):
    fig=plt.figure(figsize=(8,8),facecolor=BG)
    fig.text(.10,.95,'STATISTICAL RETHINKING / 04',color=TEAL,fontsize=11)
    fig.text(.10,.895,f'{num}  {title}',fontsize=23,color=INK,weight='bold')
    fig.text(.10,.85,sub,fontsize=13,color=INK)
    return fig
def axes_at(fig,rect):
    ax=fig.add_axes(rect,facecolor=BG)
    ax.spines[['top','right']].set_visible(False)
    ax.grid(axis='y',alpha=.15);return ax
def save(fig,name,foot):
    fig.text(.10,.075,foot,fontsize=12,color=INK,linespacing=1.6)
    fig.text(.10,.025,'明哥的微生物世界 · 原课数据 · 配套代码可复算',fontsize=10,color=TEAL)
    fig.savefig(FIG/name,dpi=180,facecolor=BG)
# 01 点是个体，粗线是组均值的不确定性。
fig=setup('01','先分清：人，与组平均','每个小点是一个人；大点和短线描述组平均。')
ax=axes_at(fig,[.13,.23,.8,.55]);rng=np.random.default_rng(20261006)
for g,col,label in [(0,TEAL,'女性'),(1,ORANGE,'男性')]:
    ax.scatter(g+rng.uniform(-.13,.13,sum(sex==g)),y[sex==g],s=15,alpha=.45,color=col)
    ax.errorbar(g+.23,m1[g],yerr=z*np.sqrt(v1[g,g]),fmt='o',capsize=6,color=INK,lw=3)
    ax.text(g+.28,m1[g],f'{m1[g]:.2f}',color=INK)
ax.set(xticks=[0,1],xticklabels=[f'女性（{sum(sex==0)}人）',f'男性（{sum(sex==1)}人）'],ylabel='体重（千克）',xlim=(-.4,1.65))
save(fig,'01-个体与组平均.png','组平均有差异，个体体重仍有大量重叠。\n短线为组均值的89%可信区间，不是个体体重范围。')
# 02 同尺度看均值差与新个体差。
fig=setup('02','平均差很多，人人都如此吗？','两张图横轴相同；都计算“男性 − 女性”。')
xx=np.linspace(-20,35,600)
for rect,values,title in [([.13,.52,.8,.25],norm.pdf(xx,delta,sd),'组平均体重之差'),([.13,.2,.8,.25],None,'各独立抽一位新个体的体重之差')]:
    ax=axes_at(fig,rect)
    if values is None:
        ax.hist(individual,bins=np.linspace(-20,35,100),density=True,color=TEAL,alpha=.65)
    else:ax.plot(xx,values,color=TEAL,lw=2);ax.fill_between(xx,values,color=TEAL,alpha=.2)
    ax.axvline(0,color=ORANGE,ls='--');ax.set(xlim=(-20,35),ylabel='概率密度',title=title)
ax.set_xlabel('体重差（千克）；小于0表示女性更重')
save(fig,'02-均值差与个体差.png',f'均值差约 {delta:.2f} 千克；但随机一对新个体中，\n男性更重的模型预测概率约 {summary[1,3]:.0%}。上下纵轴尺度不同。')
# 03 分类后每组一条线。
fig=setup('03','如果身高一样，还差多少？','把问题改成：比较相同身高处的组平均体重。')
ax=axes_at(fig,[.13,.23,.8,.55])
for g,col,label in [(0,TEAL,'女性'),(1,ORANGE,'男性')]:
    hh=np.linspace(min(h[sex==g]),max(h[sex==g]),100)
    ax.scatter(h[sex==g],y[sex==g],s=15,color=col,alpha=.35)
    ax.plot(hh,m2[g]+m2[g+2]*(hh-hc),color=col,lw=2,label=label)
ax.axvline(160,color=INK,ls='--',lw=1);ax.set(xlabel='身高（厘米）',ylabel='体重（千克）');ax.legend(frameon=False)
save(fig,'03-每组一条线.png','这里比较的是同身高处的平均体重。\n它与“不限制身高的组平均差”回答不同问题。')
# 04 同一身高的后验对比曲线，限制共同观测身高范围。
fig=setup('04','组间差，也可能随身高变化','差值 = 女性平均体重 − 男性平均体重。')
ax=axes_at(fig,[.13,.23,.8,.55])
hh=np.linspace(max(min(h[sex==0]),min(h[sex==1])),min(max(h[sex==0]),max(h[sex==1])),100)
C=np.column_stack([np.ones(len(hh)),-np.ones(len(hh)),hh-hc,-(hh-hc),np.zeros(len(hh))])
mm=C@m2;ss=np.sqrt(np.einsum('ij,jk,ik->i',C,v3,C))
ax.fill_between(hh,mm-z*ss,mm+z*ss,color=TEAL,alpha=.2);ax.plot(hh,mm,color=TEAL,lw=2)
ax.axhline(0,color=ORANGE,ls='--');ax.set(xlabel='共同观测范围内的身高（厘米）',ylabel='同身高处的平均体重差（千克）')
save(fig,'04-同身高的平均差.png',f'160厘米处：均值差约 {d160:.2f} 千克，\n89%可信区间约 {summary[2,1]:.2f} 至 {summary[2,2]:.2f} 千克。')
# 05 三个实际三次B样条基函数组成的教学示意，不伪称拟合结果。
fig=setup('05','曲线，可以用局部积木拼起来','取年龄模型中相邻的三条基函数作示意。')
ix=[8,9,10]; weights_demo=np.array([6.,-3.,5.])
ax=axes_at(fig,[.13,.53,.8,.24])
for j,c in zip(ix,[TEAL,ORANGE,INK]):ax.plot(year,BP[:,j],color=c,label=f'积木 {j-7}')
ax.set(xlim=(20,55),ylabel='基函数的值',title='每块积木只在一段年龄内起作用');ax.legend(fontsize=10,frameon=False)
ax=axes_at(fig,[.13,.21,.8,.23]);ax.plot(year,BP[:,ix]@weights_demo,color=TEAL,lw=3)
ax.set(xlim=(20,55),xlabel='年龄（岁）',ylabel='加权和',title='示意：6 × 积木1 − 3 × 积木2 + 5 × 积木3')
save(fig,'05-样条积木示意.png','这里的6、−3、5是讲解用权重，不是拟合结果。\n真正拟合时，所有积木的权重一起由先验与数据确定。')
# 06 作者樱花模型，均值的后验区间。
fig=setup('06','年龄与身高，需要一条曲线','点是不同年龄的人；并非同一个人的生长记录。')
ax=axes_at(fig,[.13,.23,.8,.55]);ax.scatter(ch['age'],ch['height'],s=9,color='#91a5a6',alpha=.45)
ax.fill_between(year,curve-z*curve_sd,curve+z*curve_sd,color=TEAL,alpha=.25)
ax.plot(year,curve,color=TEAL,lw=2);ax.set(xlabel='年龄（岁）',ylabel='身高（厘米）')
save(fig,'06-年龄身高曲线.png','深线为平均身高，阴影为均值的89%可信区间。\n这不是某一个人身高的预测区间。')
print('六张中文图已保存，交互环境中显示并保留。')
plt.show()
