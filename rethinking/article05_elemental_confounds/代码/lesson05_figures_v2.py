"""重绘作者因果图与公众号公式图片；不依赖图片生成模型绘制科学内容。"""
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib import font_manager
ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/'文章配图_v2'
fonts={f.name for f in font_manager.fontManager.ttflist}
font=next(f for f in ['PingFang SC','Microsoft YaHei','Noto Sans CJK SC','Heiti SC'] if f in fonts)
plt.rcParams.update({'font.family':font,'axes.unicode_minus':False})
bg='#faf8f3';ink='#173f50';teal='#167e83';ochre='#c4853f'
def dag(name,title,nodes,edges,note,selected=None):
 fig,ax=plt.subplots(figsize=(8,4.7),facecolor=bg);ax.set(xlim=(0,1),ylim=(0,1));ax.axis('off')
 fig.text(.06,.91,title,size=20,color=ink,weight='bold')
 boxes={}
 for k,(x,y,label) in nodes.items():
  text=ax.text(x,y,label,ha='center',va='center',size=16,color=ink,bbox=dict(boxstyle='round,pad=.5',fc='white',ec=ochre if k==selected else teal,lw=2 if k==selected else 1.4))
  boxes[k]=text.get_bbox_patch()
 fig.canvas.draw()
 for u,v in edges:
  ax.annotate('',xy=nodes[v][:2],xytext=nodes[u][:2],arrowprops=dict(arrowstyle='-|>',lw=2.5,color=teal,patchA=boxes[u],patchB=boxes[v],shrinkA=5,shrinkB=5))
 fig.text(.06,.06,note,size=12,color=ink)
 fig.savefig(OUT/name,dpi=180,facecolor=bg);plt.close(fig)
dag('01-叉.png','叉：同一个原因，影响两边',{'X':(.18,.3,'X'),'Z':(.5,.7,'Z：共同原因'),'Y':(.82,.3,'Y')},[('Z','X'),('Z','Y')],'本图只有这两条箭头；给定Z后，X与Y独立生成。')
dag('02-婚姻.png','先画假设：结婚年龄可能影响两个比例',{'M':(.18,.3,'M：结婚率'),'A':(.5,.7,'A：结婚年龄'),'D':(.82,.3,'D：离婚率')},[('A','M'),('A','D'),('M','D')],'研究M对D的影响：需要处理经过A的共同原因路径。\n箭头是科学假设；散点和回归不能单独证明它。')
dag('03-植物.png','管：处理可以通过减少真菌，让植物长高',{'T':(.15,.65,'T：处理'),'F':(.38,.22,'F：真菌'),'H1':(.66,.65,'H1：末次株高'),'H0':(.88,.22,'H0：初始株高')},[('T','F'),('F','H1'),('T','H1'),('H0','H1')],'按PPT第48页保留可能的直接通路。原书的教学生成规则中，\n处理只通过真菌起作用；总体效果包括这条经过真菌的通路。')
dag('04-对撞.png','只看获资助者，两项评价的搭配变了',{'N':(.15,.7,'新颖性'),'T':(.85,.7,'可信度'),'A':(.5,.25,'获得资助')},[('N','A'),('T','A')],'青色框：两项评价；赭色框：用来筛选的共同结果。\n新颖性低仍能入选者，更多依靠较高的可信度。',selected='A')
dag('05-幸福.png','按婚姻状态分组，组里的人发生了变化',{'A':(.15,.7,'年龄'),'H':(.85,.7,'幸福程度'),'M':(.5,.25,'婚姻状态')},[('A','M'),('H','M')],'青色框：两个原因；赭色框：用来分组的共同结果。\n幸福程度较高的人更早进入已婚组，逐渐离开未婚组。',selected='M')
dag('06-后代结构.png','检测结果，是实际真菌状态的一个线索',{'X':(.12,.7,'X\n（处理）'),'Z':(.5,.7,'Z\n（真菌状态）'),'Y':(.88,.7,'Y\n（生长）'),'A':(.5,.17,'A\n（检测结果）')},[('X','Z'),('Z','Y'),('Z','A')],'沿用作者的箭头结构；括号内为植物例子的补充解释。\n按检测结果分组，也会部分地按实际真菌状态分组。')
dag('07-祖辈.png','同一个变量，在不同路径上角色不同',{'G':(.12,.7,'G：祖辈'),'P':(.5,.7,'P：父母'),'U':(.88,.45,'U：未测原因'),'C':(.5,.16,'C：孙辈')},[('G','P'),('P','C'),('G','C'),('U','P'),('U','C')],'调整P会阻断G→P→C，同时可能打开G→P←U→C。')
formulas={
 '公式01-标准化.png':[r'$x_{\mathrm{std}}=\frac{x-\bar{x}}{s_x}$'],
 '公式02-婚姻模型.png':[r'$D_i\sim\mathrm{Normal}(\mu_i,\sigma)$',r'$\mu_i=\alpha+\beta_M M_i+\beta_A A_i$'],
 '公式03-先验.png':[r'$\alpha\sim\mathrm{Normal}(0,0.2)$',r'$\beta_M,\beta_A\sim\mathrm{Normal}(0,0.5)$',r'$\sigma\sim\mathrm{Exponential}(1)$'],
 '公式04-指数先验.png':[r'$p(\sigma\mid\lambda)=\lambda e^{-\lambda\sigma},\quad \sigma>0$',r'$E[\sigma]=\frac{1}{\lambda},\quad \lambda=1\Rightarrow E[\sigma]=1$'],
 '公式05-均值差.png':[r'$\mu(M=1,A=a)-\mu(M=0,A=a)=\beta_M$']}
for name,lines in formulas.items():
 fig=plt.figure(figsize=(8,1.0+len(lines)*.65),facecolor=bg)
 for i,line in enumerate(lines):fig.text(.5,1-(i+1)/(len(lines)+1),line,ha='center',va='center',fontsize=24,color=ink)
 fig.savefig(OUT/name,dpi=200,facecolor=bg,bbox_inches='tight',pad_inches=.25);plt.close(fig)
print('7张因果图、5张公式图片已导出。')
