from pathlib import Path
import json, shutil
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
base=Path(__file__).resolve().parent.parent if Path(__file__).resolve().parent.name=='配图源' else Path('outputs/统计课堂46_多组均数方法选择')
assets=base/'文章配图'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':13,'axes.spines.top':False,'axes.spines.right':False})
formulas=[
('01',[r'SE_{ij}=\sqrt{MSE(1/n_i+1/n_j)}',r't_{ij}=\frac{\bar y_i-\bar y_j}{SE_{ij}},\quad df=N-k'], '均数差标准误与共同误差自由度'),
('02',[r'LSD=t_{0.975,16}\,SE=0.8777\ \mathrm{MPa}',r'1.10\pm0.8777=[0.2223,\ 1.9777]'], 'LSD临界差与A1减A3的普通95%置信区间，单位MPa'),
('03',[r'q_{ij}=\frac{|\bar y_i-\bar y_j|}{\sqrt{MSE/n}}=\sqrt{2}|t_{ij}|',r'HSD=q_{0.95;4,16}\sqrt{MSE/5}',r'=4.0461\times0.29275=1.1845\ \mathrm{MPa}'], '等组量Tukey统计量与临界差'),
('04',[r'P\!\left(\max_{i=2,3,4}|T_{i1}|\leq d\right)=0.95',r'(\bar y_i-\bar y_1)\pm d\,SE_{i1}'], 'Dunnett三项共同对照比较的联合临界值与同时区间'),
('05',[r'r=j-i+1,\quad q_{ij}=\frac{\bar y_{(i)}-\bar y_{(j)}}{\sqrt{MSE/n}}',r'CR_r=q_{0.95;r,16}\sqrt{MSE/5}'], 'SNK排序跨度含两个端点，按跨度使用临界极差')]
manifest=[]
for num,lines,alt in formulas:
 fig=plt.figure(figsize=(7,0.6+0.70*len(lines)),dpi=160,facecolor='white')
 for z,line in enumerate(lines): fig.text(.5,1-(0.65+0.70*z)/(0.6+0.70*len(lines)),f'${line}$',ha='center',va='center',fontsize=26)
 fig.savefig(assets/f'公式{num}.png',dpi=320,facecolor='white')
 plt.close(fig)
 manifest.append({'number':num,'latex':lines,'alt':alt,'font_pt':26,'canvas_width_inches':7})
(base/'配图源/公式源.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2))
project=base/'article46_posthoc_selection'
if not project.is_dir():
 project=base
p=project/'data/scaffold.csv'
data=pd.read_csv(p)
fig,ax=plt.subplots(figsize=(6,4.2),layout='constrained')
for pos,group in enumerate(['A1','A2','A3','A4']):
 y=data.loc[data.group==group,'modulus'].values
 ax.scatter(pos+np.linspace(-.12,.12,5),y,s=38,color='#4d8ec5',alpha=.8)
 ax.plot([pos-.21,pos+.21],[y.mean()]*2,color='#0000F5',lw=3)
 ax.text(pos,19.2,f'{y.mean():.2f}',ha='center',color='#0000F5')
ax.set(xticks=range(4),xticklabels=['A1','A2','A3','A4'],ylabel='Compressive modulus (MPa)',ylim=(14.8,19.8))
ax.set_title('Individual specimens and group means',loc='left',fontsize=14)
fig.savefig(assets/'图1_数据.png',dpi=240);plt.close(fig)
fig,ax=plt.subplots(figsize=(6,3.5),layout='constrained')
from scipy import stats
groups=[data.loc[data.group==g,'modulus'].to_numpy() for g in ['A1','A2','A3','A4']]
mse=sum(np.sum((x-x.mean())**2) for x in groups)/(len(data)-len(groups))
se=np.sqrt(mse*(1/len(groups[0])+1/len(groups[2])))
delta=groups[0].mean()-groups[2].mean()
dun=stats.dunnett(*groups[1:],control=groups[0],rng=np.random.default_rng(20261130))
dun_ci=dun.confidence_interval(.95)
half=[stats.t.ppf(.975,16)*se,stats.studentized_range.ppf(.95,4,16)*np.sqrt(mse/5),((dun_ci.high-dun_ci.low)/2)[1]]
for i,(name,h) in enumerate(zip(['LSD: individual 95%','Tukey: 6 pairs, simultaneous 95%','Dunnett: 3 vs A1, simultaneous 95%'],half)):
 ax.errorbar(delta,2-i,xerr=h,fmt='o',color=['#64748b','#0000F5','#14806c'][i],capsize=5,lw=2)
ax.axvline(0,color='#555',ls='--',lw=1)
ax.set(yticks=[2,1,0],yticklabels=['LSD','Tukey','Dunnett'],xlabel='A1 minus A3 (MPa)',ylim=(-.6,2.6),xlim=(-.3,2.5))
ax.set_title('Same difference, different inference families',loc='left',fontsize=13)
fig.savefig(assets/'图2_区间.png',dpi=240);plt.close(fig)
target=base/'配图源/build_article_assets.py'
if Path(__file__).resolve()!=target.resolve():
 text=Path(__file__).read_text()
 start=text.index('if not p.exists():')
 end=text.index('data=pd.read_csv(p)',start)
 text=text[:start]+text[end:]
 target.write_text(text)
