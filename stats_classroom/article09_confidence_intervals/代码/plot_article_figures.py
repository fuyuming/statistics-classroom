# 文章配图复现：先运行R总入口生成覆盖率CSV，再运行此文件。
# 数学公式图直接按公式绘制；封面来自课件，不由本脚本生成。
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
import pandas as pd
from scipy.stats import t
base=Path(__file__).resolve().parents[1]
out=base/'文章配图'
# 自动选择可用中文字体；没有字体时请安装Noto Sans CJK SC后重跑。
font_candidates = ['Microsoft YaHei', 'PingFang SC', 'Noto Sans CJK SC', 'SimHei']
font_names = {f.name for f in font_manager.fontManager.ttflist}
chosen = next((f for f in font_candidates if f in font_names), 'DejaVu Sans')
plt.rcParams.update({'font.family':chosen,'axes.unicode_minus':False,'font.size':13,'savefig.facecolor':'white','mathtext.fontset':'stix'})
blue='#1769aa';red='#cf5549'
formulas={
'02_合并方差':[r's_p^2=\frac{(n_1-1)s_1^2+(n_2-1)s_2^2}{n_1+n_2-2}',r'=\frac{3\times8.25+4\times5.20}{7}=6.507143',r'\mathrm{SE}=s_p\sqrt{\frac{1}{n_1}+\frac{1}{n_2}}=1.711203'],
'03_t反解':[r'P\!\left(-t\leq\frac{D-\delta}{\mathrm{SE}}\leq t\right)=0.95',r'P(-t\,\mathrm{SE}\leq D-\delta\leq t\,\mathrm{SE})=0.95',r'P(D-t\,\mathrm{SE}\leq\delta\leq D+t\,\mathrm{SE})=0.95',r't=t_{0.975,7}=2.364624'],
'04_率区间':[r'\widehat{\mathrm{SE}}=\sqrt{\frac{\hat p(1-\hat p)}{n}}',r'P\!\left(-1.96\leq\frac{\hat p-p}{\widehat{\mathrm{SE}}}\leq1.96\right)\approx0.95',r'P(\hat p-1.96\widehat{\mathrm{SE}}\leq p\leq\hat p+1.96\widehat{\mathrm{SE}})\approx0.95'],
'05_率差区间':[r'\widehat{\mathrm{SE}}=\sqrt{\frac{\hat p_1(1-\hat p_1)}{n_1}+\frac{\hat p_2(1-\hat p_2)}{n_2}}',r'(\hat p_1-\hat p_2)\ \pm\ 1.96\widehat{\mathrm{SE}}'],
'06_方差反解':[r'P\!\left(q_{0.025}\leq\frac{(n-1)S^2}{\sigma^2}\leq q_{0.975}\right)=0.95',r'P\!\left(\frac{(n-1)S^2}{q_{0.975}}\leq\sigma^2\leq\frac{(n-1)S^2}{q_{0.025}}\right)=0.95'],
'07_方差比反解':[r'\frac{R}{\theta}=\frac{S_1^2/\sigma_1^2}{S_2^2/\sigma_2^2}\sim F(n_1-1,n_2-1)',r'P\!\left(f_{0.025}\leq\frac{R}{\theta}\leq f_{0.975}\right)=0.95',r'P\!\left(\frac{R}{f_{0.975}}\leq\theta\leq\frac{R}{f_{0.025}}\right)=0.95']}
for name,lines in formulas.items():
 fig=plt.figure(figsize=(7.2,len(lines)*.8+.2))
 for i,line in enumerate(lines):fig.text(.5,1-(i+.6)/(len(lines)+.2),'$'+line+'$',ha='center',va='center',fontsize=22,color=blue)
 fig.savefig(out/(name+'.png'),dpi=200,bbox_inches='tight',pad_inches=.18);plt.close(fig)
fig,ax=plt.subplots(figsize=(7,4.4));x=np.linspace(-5,5,1600);a=t.ppf(.025,7);b=-a;y=t.pdf(x,7)
ax.plot(x,y,color=blue,lw=2);ax.fill_between(x,y,where=(x>=a)&(x<=b),color='#b5dff7');ax.fill_between(x,y,where=(x<a)|(x>b),color='#f4cdc7')
ax.axvline(a,color=red,ls='--');ax.axvline(b,color=red,ls='--');ax.text(0,.16,'中间面积 95%',ha='center',color=blue,fontsize=15)
ax.text(-3.6,.04,'2.5%',ha='center');ax.text(3.6,.04,'2.5%',ha='center')
ax.set(xticks=[a,0,b],xticklabels=['−2.3646','0','2.3646'],xlabel='t 统计量',ylabel='概率密度',title='先取抽样分布中间的95%：t(7)')
ax.spines[['top','right']].set_visible(False);fig.tight_layout();fig.savefig(out/'01_中间95.png',dpi=200);plt.close(fig)
df=pd.read_csv(base/'运行结果/R/06_重复抽样.csv').iloc[:40]
fig,ax=plt.subplots(figsize=(7,8));color=np.where(df.covered,blue,red)
ax.hlines(df['round'],df.lower,df.upper,colors=color,lw=2);ax.scatter(df.estimate,df['round'],c=color,s=12)
ax.axvline(2,color='#252d3e',ls='--',label='固定真值：2 g');ax.set(ylim=(41,0),xlabel='总体均数差的95%置信区间（g）',ylabel='抽样轮次',title='同一总体、同一方法，换一轮样本就换一个区间')
ax.legend(loc='lower right');ax.spines[['top','right']].set_visible(False);fig.tight_layout();fig.savefig(out/'08_重复抽样.png',dpi=200);plt.close(fig)
# 保存一幅阅读导航图：不同参数不混在同一数轴。
rows=[('均数差','−2.00 至 6.10 g'),('总体率','35.25% 至 54.75%'),('率差','6.45 至 16.42 个百分点'),('总体方差','0.06549 至 0.62061 (μg/g)²'),('方差比','1.94 至 31.48')]
fig=plt.figure(figsize=(7.2,4.7));fig.text(.06,.92,'五个问题，五个95%置信区间',fontsize=22,color=blue,weight='bold')
for i,(name,res) in enumerate(rows):
 y=.77-i*.15;fig.text(.06,y,name,fontsize=17,color=blue);fig.text(.30,y,res,fontsize=16)
fig.savefig(out/'09_结果速览.png',dpi=200,bbox_inches='tight');plt.close(fig)

items={'10_样本率抽样分布':[r'X\sim\mathrm{Bin}(n,p),\quad\hat p=X/n',r'\hat p\ \overset{\mathrm{approx}}{\sim}\ N\!\left(p,\frac{p(1-p)}{n}\right)',r'\frac{\hat p-p}{\sqrt{p(1-p)/n}}\ \overset{\mathrm{approx}}{\sim}\ N(0,1)'], '11_率差抽样分布':[r'D=\hat p_1-\hat p_2',r'D\ \overset{\mathrm{approx}}{\sim}\ N\!\left(p_1-p_2,\frac{p_1(1-p_1)}{n_1}+\frac{p_2(1-p_2)}{n_2}\right)',r'\frac{D-(p_1-p_2)}{\sqrt{p_1(1-p_1)/n_1+p_2(1-p_2)/n_2}}',r'\overset{\mathrm{approx}}{\sim}\ N(0,1)'], '12_均数差抽样分布':[r'S_p^2=\frac{(n_1-1)S_1^2+(n_2-1)S_2^2}{n_1+n_2-2}',r'T=\frac{(\bar X_1-\bar X_2)-(\mu_1-\mu_2)}{S_p\sqrt{1/n_1+1/n_2}}',r'\sim t(n_1+n_2-2)'], '13_方差抽样分布':[r'S^2=\frac{1}{n-1}\sum_{i=1}^{n}(X_i-\bar X)^2',r'W=\frac{(n-1)S^2}{\sigma^2}\sim\chi^2(n-1)']}
for name,lines in items.items():
 fig=plt.figure(figsize=(7.5,len(lines)*.9))
 for i,line in enumerate(lines):fig.text(.5,1-(i+.5)/len(lines),'$'+line+'$',ha='center',va='center',fontsize=21,color='#1769aa')
 fig.savefig(base/'文章配图'/(name+'.png'),dpi=200,bbox_inches='tight',pad_inches=.2);plt.close(fig)

from scipy import stats
x=np.array([4.2,4.5,3.6,4.7,4.0,3.8,3.7,4.2])
fig,axs=plt.subplots(1,2,figsize=(10,4.6),layout='constrained')
axs[0].scatter(x,[0,0,0,0,0,0,0,.04],s=70,color='#1769AA')
axs[0].set(xlabel='铅浓度（μg/g）',yticks=[],ylim=(-.1,.15),title='先看8个原始观测')
q=stats.norm.ppf((np.arange(1,9)-.375)/(8+.25))
slope=np.diff(np.quantile(x,[.25,.75]))[0]/np.diff(stats.norm.ppf([.25,.75]))[0]
intercept=np.quantile(x,.25)-slope*stats.norm.ppf(.25)
axs[1].scatter(q,np.sort(x),s=60,color='#1769AA')
axs[1].plot(q,intercept+slope*q,color='#D65A4A')
axs[1].set(xlabel='标准正态理论分位数',ylabel='样本铅浓度（μg/g）',title='正态Q–Q图：看点是否大致沿直线')
fig.suptitle('土壤铅浓度：先检查正态模型是否合适',fontsize=17,color='#253858')
fig.supxlabel('Shapiro–Wilk：W = 0.95342，P = 0.74561；未发现明显偏离，不能据此证明正态。',fontsize=11)
fig.savefig(base/'文章配图/14_铅浓度正态性.png',dpi=180)
