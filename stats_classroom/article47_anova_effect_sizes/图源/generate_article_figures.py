"""从共享教学数据重建正文科学图及公式。运行方式见文章目录README。"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
BASE = Path(__file__).resolve().parents[1]
IMG = BASE / '文章配图'
F = IMG / '公式'
F.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','mathtext.fontset':'stix', 'axes.spines.top':False,'axes.spines.right':False})
# 所有公式同一画布宽度和数学字号，不逐图缩小；分行避免长式挤压。
formulas = {
'01': [r'SS_B=\sum_{i=1}^{k}n_i(\bar y_i-\bar y)^2',r'SS_W=\sum_{i=1}^{k}\sum_{j=1}^{n_i}(y_{ij}-\bar y_i)^2',r'SS_T=\sum_{i,j}(y_{ij}-\bar y)^2=SS_B+SS_W'],
'02': [r'\eta^2=\frac{SS_B}{SS_T}',r'=\frac{15.4495}{22.3055}\approx0.692632'],
'03': [r'\mu=\frac{1}{4}\sum_{i=1}^{4}\mu_i,\quad V_\mu=\frac{1}{4}\sum_{i=1}^{4}(\mu_i-\mu)^2',r'\eta^2_{\mathrm{pop}}=\frac{V_\mu}{V_\mu+\sigma^2}'],
'04': [r'MSE=\frac{SS_W}{N-k}=\frac{6.8560}{16}=0.4285',r'\mathrm{E}(SS_B)=NV_\mu+(k-1)\sigma^2',r'\mathrm{E}(SS_T)=NV_\mu+(N-1)\sigma^2'],
'04a': [r'SS_B=\sum_i n_i(\bar y_i-\mu)^2-N(\bar y-\mu)^2',r'\mathrm{E}\!\left[\sum_i n_i(\bar Y_i-\mu)^2\right]=NV_\mu+k\sigma^2',r'\mathrm{E}\!\left[N(\bar Y-\mu)^2\right]=\sigma^2'],
'05': [r'\omega^2=\frac{SS_B-(k-1)MSE}{SS_T+MSE}',r'=\frac{15.4495-3\times0.4285}{22.3055+0.4285}',r'=\frac{14.1640}{22.7340}\approx0.623032'],
'06': [r'\eta_p^2=\frac{SS_{\mathrm{effect}}}{SS_{\mathrm{effect}}+SS_{\mathrm{error}}}'],
'07a': [r'SS_A=4[(12-13)^2+(14-13)^2]=8',r'SS_B=4[(11-13)^2+(15-13)^2]=32',r'SS_E=4(1^2+1^2)=8,\quad SS_T=48',r'SS_{AB}=48-8-32-8=0'],
'07b': [r'\eta_A^2=\frac{SS_A}{SS_T}=\frac{8}{48}\approx0.167',r'\eta_{p,A}^2=\frac{SS_A}{SS_A+SS_E}=\frac{8}{8+8}=0.500'],
'08': [r'R^2=1-\frac{SS_W}{SS_T}',r'=\frac{SS_B}{SS_T}=\eta^2\approx0.692632'],
'09': [r'R^2_{\mathrm{adj}}=1-\frac{SS_W/(N-k)}{SS_T/(N-1)}',r'=\frac{SS_B-(k-1)MSE}{SS_T}\approx0.635000'],
'10': [r'\eta^2_{\mathrm{pop}}=\frac{f^2}{1+f^2}',r'f=0.10\quad\Rightarrow\quad\eta^2_{\mathrm{pop}}\approx0.009901',r'f=0.25\quad\Rightarrow\quad\eta^2_{\mathrm{pop}}\approx0.058824',r'f=0.40\quad\Rightarrow\quad\eta^2_{\mathrm{pop}}\approx0.137931']}
for num, lines in formulas.items():
    height = 110 + 155*len(lines)
    fig=plt.figure(figsize=(4,height/300),dpi=300,facecolor='white')
    texts=[]
    for j,line in enumerate(lines):
        y = 1-(65+155*j+50)/height
        texts.append(fig.text(.5,y,'$'+line+'$',ha='center',va='center',fontsize=18,color='#202329'))
    fig.canvas.draw()
    for tx in texts:
        bb=tx.get_window_extent()
        assert bb.width < 1160, (num, bb.width)
    fig.savefig(F/f'式{num}.png',dpi=300,facecolor='white')
    plt.close(fig)
(F/'formulas.json').write_text(json.dumps(formulas,ensure_ascii=False,indent=2),encoding='utf-8')
(F/'formulas.tex').write_text('\n\n'.join('% 式（'+n+'）\n\\begin{aligned}\n'+' \\\\\n'.join(lines)+'\n\\end{aligned}' for n,lines in formulas.items()),encoding='utf-8')
d=pd.read_csv(BASE/'练习项目/data/scaffold.csv')
groups = [d.loc[d.group==g,'modulus'].to_numpy() for g in ['A1','A2','A3','A4']]
means=np.array([g.mean() for g in groups]); grand=d.modulus.mean()
b=sum(len(g)*(g.mean()-grand)**2 for g in groups); w=sum(sum((g-g.mean())**2) for g in groups)
fig,ax=plt.subplots(figsize=(4,3.4),dpi=300)
for i,g in enumerate(groups):
    ax.scatter(i+np.linspace(-.075,.075,len(g)),g,s=40,c='#3366cc',zorder=3,edgecolors='white',linewidths=.5)
    ax.plot([i-.23,i+.23],[g.mean(),g.mean()],lw=2.5,color='#111827',zorder=4)
ax.set(xticks=range(4),xticklabels=['A1','A2','A3','A4'],ylabel='Compressive modulus (MPa)',xlabel='Pore design',ylim=(14.8,19.5))
ax.tick_params(labelsize=13);ax.xaxis.label.set_size(14);ax.yaxis.label.set_size(14)
ax.grid(axis='y',color='#eeeeee',zorder=0)
fig.tight_layout();fig.savefig(IMG/'图1_试件与均数.png',dpi=240);plt.close(fig)
# 比例示意使用真实平方和；标签均用英语，中文解释位于正文和图注。
fig,ax=plt.subplots(figsize=(4,2.6),dpi=300)
ax.barh([0],[b],color='#0000f5',height=.5)
ax.barh([0],[w],left=b,color='#e1e5ed',height=.5)
ax.text(b/2,0,f'SSB\n{b:.4f}\n{b/(b+w):.2%}',ha='center',va='center',fontsize=15,color='white',linespacing=1.4)
ax.text(b+w/2,0,f'SSW\n{w:.4f}\n{w/(b+w):.2%}',ha='center',va='center',fontsize=14,color='#202329',linespacing=1.4)
ax.set(yticks=[],ylim=(-.45,.48),xlim=(0,b+w),xlabel=r'Sum of squares (MPa$^2$)')
ax.set_title(f'SST = {b+w:.4f}',fontsize=16,pad=12)
ax.set_xticks([0,5,10,15,20]);ax.tick_params(labelsize=12);ax.xaxis.label.set_size(13)
ax.spines['left'].set_visible(False)
fig.tight_layout();fig.savefig(IMG/'图2_平方和比例.png',dpi=240);plt.close(fig)
print('Rebuilt 12 formula PNGs, LaTeX sources and 2 scientific figures.')
