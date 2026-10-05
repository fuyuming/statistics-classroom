"""由 howell_v5.py 自动调用；所有图都使用 v5 同一三参数模型。"""
import matplotlib.pyplot as plt
from matplotlib import font_manager
OUT=ROOT/'文章配图_v5'
OUT.mkdir(exist_ok=True)
ma,mb,sigma=means
at160=grid[grid[:,0]==160][0]
a=samples[:,0]; b=samples[:,1]
post=np.full(len(a),1/len(a))
font=next((v for v in ['PingFang SC','Microsoft YaHei','Noto Sans CJK SC'] if v in {q.name for q in font_manager.fontManager.ttflist}), 'DejaVu Sans')
plt.rcParams.update({'font.family':font,'font.size':15,'axes.unicode_minus':False,'axes.labelsize':15,'xtick.labelsize':13,'ytick.labelsize':13,'savefig.facecolor':'#faf8f3'})
BG='#faf8f3'; INK='#193c46'; TEAL='#167e83'; ORANGE='#ca6841'; PALE='#cce3df'; GRAY='#91a5a6'
def base(num,title,sub):
    fig=plt.figure(figsize=(8,8),facecolor=BG)
    fig.text(.10,.95,'STATISTICAL RETHINKING  /  03',fontsize=11,color=TEAL,weight='bold')
    fig.text(.10,.895,f'{num}  {title}',fontsize=23,color=INK,weight='bold')
    fig.text(.10,.85,sub,fontsize=14,color=INK)
    return fig
def axis(fig,rect):
    ax=fig.add_axes(rect,facecolor=BG)
    for side in ['top','right']: ax.spines[side].set_visible(False)
    for side in ['bottom','left']: ax.spines[side].set_color('#c5cecc')
    ax.tick_params(colors=INK); ax.grid(axis='y',color='#dde3df',lw=.6); ax.set_axisbelow(True)
    return ax
def footer(fig,text):
    fig.text(.10,.075,text,fontsize=13,color=INK,linespacing=1.6)
    fig.text(.10,.025,'明哥的微生物世界  ·  Howell1 成年人  ·  三参数网格与 quap 对照',fontsize=10,color=TEAL)
def save(fig,name):
    fig.savefig(OUT/name,dpi=200); plt.close(fig)
# 01 全部真实点 + 斜率的几何含义
fig=base('01','一条线，说的是平均关系','一个点是一位成年人；点不必落在线上。')
ax=axis(fig,[.12,.23,.82,.55]); ax.scatter(h,y,s=24,color=GRAY,alpha=.65,edgecolor='none')
ax.plot(grid[:,0],grid[:,1],color=TEAL,lw=3)
u=ma+mb*(150-hc); v=ma+mb*(160-hc)
ax.plot([150,160,160],[u,u,v],color=ORANGE,lw=2,ls='--'); ax.text(155,u-2.3,'+10 厘米',ha='center',color=ORANGE,fontsize=13)
ax.text(161,(u+v)/2,'约 +6.29 千克',color=ORANGE,fontsize=13)
ax.set(xlabel='身高（厘米）',ylabel='体重（千克）',xlim=(134,181),ylim=(29,65))
footer(fig,'这批数据中：身高相差 1 厘米，平均体重相差约 0.63 千克。\n这是群体关联，不能当作同一个人的增重规律。'); save(fig,'01-身高体重平均关系.png')
# 02 一位真实个体的残差
idx=int(np.argmin(np.abs(h-160)+.03*np.abs(y-(ma+mb*x)-6)))
H=float(h[idx]); W=float(y[idx]); M=ma+mb*(H-hc)
fig=base('02','点离直线多远？','残差 = 实际体重 − 模型给出的平均体重')
ax=axis(fig,[.12,.26,.82,.52]); ax.scatter(h,y,s=19,color=GRAY,alpha=.3,edgecolor='none'); ax.plot(grid[:,0],grid[:,1],color=TEAL,lw=2.5)
ax.plot([H,H],[M,W],color=ORANGE,lw=3); ax.scatter([H,H],[M,W],c=[TEAL,ORANGE],s=90,zorder=5)
ax.annotate(f'实际读数 {W:.2f}',(H,W),xytext=(H+4,W+3),color=ORANGE,fontsize=14,arrowprops={'arrowstyle':'-','color':ORANGE})
ax.annotate(f'线上均值 {M:.2f}',(H,M),xytext=(H+5,M-4),color=TEAL,fontsize=14,arrowprops={'arrowstyle':'-','color':TEAL})
ax.set(xlabel='身高（厘米）',ylabel='体重（千克）',xlim=(134,183),ylim=(29,67))
footer(fig,f'示例取自真实记录：身高 {H:.2f} 厘米。\n残差 = {W:.2f} − {M:.2f} ≈ {W-M:+.2f} 千克。'); save(fig,'02-一个人的残差.png')
# 03 条件分布
fig=base('03','同样身高，体重仍会不同','把身高固定在 160 厘米，看看模型允许哪些体重。')
ax=axis(fig,[.12,.25,.82,.52]); z=np.linspace(32,73,401); mu=at160[1]; d=norm.pdf(z,mu,sigma)
ax.fill_between(z,d,color=PALE); ax.plot(z,d,color=TEAL,lw=3); ax.axvline(mu,color=ORANGE,lw=2,ls='--')
ax.text(mu+.6,d.max()*.9,f'中心约 {mu:.2f} 千克',color=INK,fontsize=14)
ax.annotate('',(mu-sigma,.025),(mu+sigma,.025),arrowprops={'arrowstyle':'<->','color':INK,'lw':1.5}); ax.text(mu,.017,'中心左右各一个 σ',ha='center',fontsize=13,color=INK)
ax.set(xlabel='可能的体重（千克）',ylabel='概率密度',ylim=(0,.11),xlim=(32,73))
footer(fig,f'σ ≈ {sigma:.2f} 千克：描述个体读数的散布，不是最大误差。\n这是固定平均参数后的模型曲线；概率要看曲线下的面积。'); save(fig,'03-同身高的体重波动.png')
np.savetxt(RES/'160cm_条件正态密度.csv',np.column_stack([z,d]),delimiter=',',header='weight_kg,density',comments='')
# 04 分别展示本文分析先验与后验，不冒充原课完整先验
fig=base('04','数据怎样筛选候选直线？','上下图尺度相同；每一条细线都代表一组参数。')
rng=np.random.default_rng(20261005)
pa=rng.normal(60,10,40); pb=rng.lognormal(0,1,40)
ix=rng.choice(len(post),40,p=post)
for rect,A,B,label in [([.12,.51,.82,.26],pa,pb,'看数据前：作者模型的先验'),([.12,.18,.82,.26],a[ix],b[ix],'看数据后：按后验权重抽出的直线')]:
    ax=axis(fig,rect)
    for aa,bb in zip(A,B): ax.plot(grid[:,0],aa+bb*(grid[:,0]-hc),color=TEAL,alpha=.23,lw=1.2)
    ax.scatter(h,y,s=10,color=ORANGE,alpha=.3,edgecolor='none'); ax.set(xlim=(134,181),ylim=(-150,300),ylabel='体重（千克）'); ax.set_title(label,fontsize=14,loc='left',color=INK,pad=8)
ax.set_xlabel('身高（厘米）')
fig.text(.10,.067,'先验：α ~ Normal(60,10)，β ~ LogNormal(0,1)。',fontsize=12,color=INK)
fig.text(.10,.025,'均值线示意；未加入个体波动，超出画幅的线段被截断。',fontsize=11,color=TEAL)
save(fig,'04-先验到后验直线.png')
np.savetxt(RES/'先验后验代表线.csv',np.column_stack([pa,pb,a[ix],b[ix]]),delimiter=',',header='prior_alpha,prior_beta,posterior_alpha,posterior_beta',comments='')
# 05 区间和固定160切片对照
fig=base('06','平均趋势确定，个人仍有差异','同一个身高，要分清两个不同的问题。')
ax=axis(fig,[.12,.42,.82,.36]); ax.fill_between(grid[:,0],grid[:,4],grid[:,5],color=PALE,label='新个体：89% 预测区间'); ax.fill_between(grid[:,0],grid[:,2],grid[:,3],color=TEAL,alpha=.75,label='平均体重：89% 可信区间'); ax.plot(grid[:,0],grid[:,1],color=INK,lw=1.4)
ax.axvline(160,color=ORANGE,lw=1,ls='--'); ax.set(xlabel='身高（厘米）',ylabel='体重（千克）',xlim=(134,181)); ax.legend(loc='upper left',fontsize=11,frameon=False)
ax=axis(fig,[.17,.18,.73,.13]); ax.grid(False); ax.set_xlim(40,60); ax.set_ylim(-.6,1.7)
for yy,lo,hi,col,label in [(1,at160[2],at160[3],TEAL,'平均体重'),(0,at160[4],at160[5],GRAY,'新个体')]:
    ax.plot([lo,hi],[yy,yy],color=col,lw=9,solid_capstyle='round'); ax.text(57,yy+.27,f'{lo:.2f}—{hi:.2f}',ha='right',va='bottom',fontsize=11,color=INK)
ax.set_yticks([0,1],['新个体','平均体重']); ax.set_xlabel('身高 160 厘米时的体重（千克）',fontsize=12); ax.spines['left'].set_visible(False)
fig.text(.10,.075,'宽区间包含两层不确定性：平均线在哪里 ＋ 个人偏离多少。',fontsize=12,color=INK)
fig.text(.10,.025,'同一 Howell1 模型  ·  α、β、σ 均未知  ·  低差异数值积分',fontsize=11,color=TEAL)
save(fig,'05-平均与个人预测区间.png')
print('generated five figures; beta=',mb,'160cm=',at160)

# 07 两种计算的边际后验，密度面积约为1；quap 保留联合协方差。
fig=base('05','同一个模型，两种计算','实线：三维网格；虚线：quap 的正态近似。')
for i,(label,unit) in enumerate([('α：平均身高处的平均体重','千克'),('β：身高相差1厘米的平均体重差','千克/厘米'),('σ：同身高体重的分散程度','千克')]):
    ax=axis(fig,[.14,.60-i*.22,.78,.14])
    dx=axes[i][1]-axes[i][0]
    ax.plot(axes[i],marginals[i]/dx,color=TEAL,lw=2)
    zz=np.linspace(mode[i]-4*np.sqrt(cov[i,i]),mode[i]+4*np.sqrt(cov[i,i]),301)
    ax.plot(zz,norm.pdf(zz,mode[i],np.sqrt(cov[i,i])),color=ORANGE,lw=2,ls='--')
    ax.set_title(label,loc='left',fontsize=12,color=INK)
    ax.set_ylabel('密度',fontsize=11)
    ax.set_xlabel(unit,fontsize=11)
fig.text(.10,.055,'前两项几乎重合；σ 的网格后验略向右偏。',fontsize=13,color=INK)
save(fig,'07-网格与quap后验对照.png')
