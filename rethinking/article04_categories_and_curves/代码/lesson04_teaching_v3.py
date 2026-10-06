"""PPT主线的补充教学图：因果图、已知真值的模拟、全样本曲线。
这些图为本文重绘，非作者原图截图。Run后原生显示并保存中文PNG。
多项式仅用最小二乘演示函数形状，不充当PPT贝叶斯拟合的精确复刻。
"""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import font_manager
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'文章配图_v3'
font=next((x for x in ['PingFang SC','Microsoft YaHei','Noto Sans CJK SC','Heiti SC'] if x in {f.name for f in font_manager.fontManager.ttflist}),None)
if not font:raise RuntimeError('请安装 Noto Sans CJK SC 中文字体')
plt.rcParams.update({'font.family':font,'axes.unicode_minus':False,'font.size':13})
bg='#faf8f3';ink='#193c46';teal='#167e83';orange='#ca6841'
def page(title,sub):
 f=plt.figure(figsize=(9,8),facecolor=bg);f.text(.08,.94,title,fontsize=23,color=ink,weight='bold');f.text(.08,.89,sub,fontsize=13,color=teal);return f
def save(f,name,note):
 f.text(.08,.065,note,fontsize=12,color=ink,linespacing=1.5);f.savefig(OUT/name,dpi=180,facecolor=bg)
f=page('先决定：想知道哪一种影响？','对应PPT第9—28页；箭头是建模假设，需要科学知识支持。')
a=f.add_axes([.08,.4,.84,.43]);a.axis('off');a.set(xlim=(0,1),ylim=(0,1))
pos={'S':(.18,.18),'H':(.18,.8),'W':(.79,.8)}
for name,label in [('S','S：性别'),('H','H：身高'),('W','W：体重')]:
 x,y=pos[name];a.text(x,y,label,ha='center',va='center',fontsize=20,color=ink,bbox=dict(boxstyle='round,pad=.55',fc='white',ec=teal,lw=2))
for u,v in [('S','H'),('H','W'),('S','W')]:
 a.annotate('',xy=pos[v],xytext=pos[u],arrowprops=dict(arrowstyle='-|>',lw=3,color=orange,shrinkA=40,shrinkB=45))
f.text(.12,.32,'总效应：同时保留 S → H → W 与 S → W。',fontsize=17,color=ink)
f.text(.12,.24,'直接效应：在指定身高处比较 S → W。',fontsize=17,color=ink)
f.text(.12,.16,'身高的效应：在同一性别中比较身高变化。',fontsize=17,color=ink)
save(f,'A-因果问题.png','散点图本身不能决定箭头方向。\n以下因果解释都以这套简化结构及相应识别假设成立为条件。')
f=page('先造一个知道答案的世界','对应PPT第21—38、50页；两组模拟用于回答不同问题。')
rng=np.random.default_rng(20261005);S=np.repeat([0,1],50);H=np.where(S==0,150.,160.)+rng.normal(0,5,100)
for j,(aa,bb,title) in enumerate([(np.array([0.,0.]),np.array([.5,.6]),'例1：总体均值差为21千克'),(np.array([0.,10.]),np.array([.5,.5]),'例2：总差15，直接差10')]):
 W=aa[S]+bb[S]*H+rng.normal(0,5,100);ax=f.add_axes([.10+j*.46,.26,.37,.53],facecolor=bg)
 for g,c,label in [(0,teal,'女性'),(1,orange,'男性')]:
  ax.scatter(H[S==g],W[S==g],s=15,color=c,alpha=.4);hh=np.linspace(135,175,60);ax.plot(hh,aa[g]+bb[g]*hh,color=c,label=label)
 ax.set(xlabel='身高（厘米）',ylabel='体重（千克）',title=title,ylim=(60,120));ax.legend(frameon=False,fontsize=10);ax.spines[['right','top']].set_visible(False)
 np.savetxt(ROOT/f'数据/v3/synthetic_case{j+1}.csv',np.column_stack([S+1,H,W]),delimiter=',',header='S,H,W',comments='')
save(f,'B-模拟检验.png','图中的体重由作者的教学规则生成，不是真实人体数据。\n例2：平均身高差10厘米带来5千克，再加直接差10千克。')
f=page('把儿童放回来，一条直线不够了','对应PPT第65—71页；544人的身高与体重。')
d=np.genfromtxt(ROOT/'数据/v3/Howell1_all.csv',delimiter=',',names=True);h=d['height'];w=d['weight'];xx=np.linspace(50,200,400)
ax=f.add_axes([.12,.24,.81,.56],facecolor=bg);ax.scatter(h,w,s=13,color='#8ea9ac',alpha=.4)
for degree,c in [(1,ink),(2,teal),(4,orange)]:
 co=np.polynomial.polynomial.polyfit((h-h.mean())/h.std(),w,degree)
 yy=np.polynomial.polynomial.polyval((xx-h.mean())/h.std(),co)
 ax.plot(xx,yy,color=c,lw=2,label={1:'直线',2:'二次多项式',4:'四次多项式'}[degree])
ax.axvspan(50,h.min(),color='grey',alpha=.1);ax.axvspan(h.max(),200,color='grey',alpha=.1)
ax.set(xlabel='身高（厘米）',ylabel='体重（千克）',xlim=(50,200),ylim=(0,100));ax.legend(frameon=False);ax.spines[['right','top']].set_visible(False)
save(f,'C-多项式形状.png','三条线用最小二乘作形状演示，未画后验区间。\n灰色为样本身高范围之外；贴近现有点不等于解释了生长机制。')
print('三张PPT主线教学图已保存；GUI环境显示。')
plt.show()
