"""从数学表达式和精确组合数重建文章图式；封面沿用课堂PDF。"""
from pathlib import Path
import os
import tempfile
os.environ.setdefault('MPLCONFIGDIR', tempfile.mkdtemp(prefix='hypergeom_figures_'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd
from math import comb
import fitz

OUT = Path(__file__).resolve().parent.parent
FIG = OUT / '文章配图'
FIG.mkdir(parents=True, exist_ok=True)
BLUE = '#0000F5'
INK = '#233449'
ORANGE = '#D76B26'
plt.rcParams.update({'font.family':'DejaVu Sans', 'font.size':15,
                     'mathtext.fontset':'dejavusans', 'axes.spines.top':False,
                     'axes.spines.right':False, 'savefig.facecolor':'white'})

equations = [
 [r'$\binom{33}{22}=\dfrac{33!}{22!\,11!}=193536720$'],
 [r'$\Pr(X=4)=\dfrac{\binom{9}{4}\binom{24}{18}}{\binom{33}{22}}$',
  r'$=\dfrac{16959096}{193536720}\approx0.08762728$'],
 [r'$\Pr(X=x)=\dfrac{\binom{K}{x}\binom{N-K}{n-x}}{\binom{N}{n}}$'],
 [r'$\max\{0,\ n-(N-K)\}\leq x$', r'$x\leq\min\{n,K\},\quad x\in\mathbb{Z}$'],
 [r'$\mathrm{E}(X)=n\dfrac{K}{N}$',
  r'$\mathrm{Var}(X)=n\dfrac{K}{N}\left(1-\dfrac{K}{N}\right)\dfrac{N-n}{N-1}$'],
 [r'$\Pr(X_1=x\mid S=9)$',
  r'$=\dfrac{\Pr(X_1=x,\ X_0=9-x)}{\Pr(S=9)}$',
  r'$=\dfrac{\binom{22}{x}\binom{11}{9-x}\,p^9(1-p)^{24}}{\binom{33}{9}\,p^9(1-p)^{24}}$',
  r'$=\dfrac{\binom{22}{x}\binom{11}{9-x}}{\binom{33}{9}}$'],
 [r'$\dfrac{\binom{22}{x}\binom{11}{9-x}}{\binom{33}{9}}$',
  r'$=\dfrac{\binom{9}{x}\binom{24}{22-x}}{\binom{33}{22}}$'],
 [r'$D=\dfrac{X}{22}-\dfrac{9-X}{11}=\dfrac{3X}{22}-\dfrac{9}{11}$',
  r'$X=4:\quad D=-\dfrac{3}{11}\approx-0.2727$',
  r'$\mathrm{E}(D\mid S=9)=0$',
  r'$\mathrm{Var}(D\mid S=9)=\dfrac{27}{968}\approx0.027893$'],

 [r'$\hat{p}=\dfrac{50+30}{100+100}=0.40$',
  r'$\widehat{\mathrm{SE}}_0=\sqrt{\hat{p}(1-\hat{p})\left(\dfrac{1}{n_1}+\dfrac{1}{n_0}\right)}$',
  r'$=\sqrt{0.40\times0.60\left(\dfrac{1}{100}+\dfrac{1}{100}\right)}$',
  r'$\approx0.069282$',
  r'$Z=\dfrac{(\hat{p}_1-\hat{p}_0)-0}{\widehat{\mathrm{SE}}_0}$',
  r'$=\dfrac{(0.50-0.30)-0}{0.0692820323}\approx2.886751$',
  r'$P_{\mathrm{normal}}=2\Phi(-|Z|)$',
  r'$\approx0.00389242$'],
 [r'$\Pr(A=a)=\dfrac{\binom{80}{a}\binom{120}{100-a}}{\binom{200}{100}}$',
  r'$a=0,1,\ldots,80$',
  r'$\Pr(A=50)\approx0.001798976$']

]
for number, lines in enumerate(equations, 1):
    # 所有公式固定画布宽及32pt数学字体；逐行换行避免长式缩字号。
    height = 1.65 * len(lines) + .65
    fig = plt.figure(figsize=(9,height))
    for i,line in enumerate(lines):
        fig.text(.5, 1-(i+.62)/(len(lines)+.50), line,
                 ha='center', va='center', fontsize=32, color=INK)
    fig.savefig(FIG / f'公式{number:02}.png', dpi=180)
    plt.close(fig)

# 条件概率定义与原式6共用式号，按6a/6b分步阅读。
lines = [r'$\Pr(E\mid B)=\dfrac{\Pr(E\cap B)}{\Pr(B)}$', r'$\Pr(B)>0$']
fig = plt.figure(figsize=(9,3.95))
for i,line in enumerate(lines):
    fig.text(.5,1-(i+.62)/2.50,line,ha='center',va='center',fontsize=32,color=INK)
fig.savefig(FIG/'公式06a_条件概率.png',dpi=180)
plt.close(fig)

# 图1按实际卡片数量画出观察结果；每个点代表不同对象。
fig, ax = plt.subplots(figsize=(9,9))
ax.set(xlim=(0,10),ylim=(0,10))
ax.axis('off')
for y, h, title in [(5.5,3.9,'Selected: 22 cards'),(.3,3.9,'Remaining: 11 cards')]:
    ax.add_patch(patches.FancyBboxPatch((.35,y),9.3,h,boxstyle='round,pad=0.02',
                                      facecolor='#F6F8FC',edgecolor='#CED8E5',linewidth=1.5))
    ax.text(.75,y+h-.5,title,fontsize=24,color=INK,weight='bold')
def dots(count,start_x,y,positive):
    for i in range(count):
        col = i % 9
        row = i // 9
        ax.add_patch(patches.Circle((start_x+.82*col,y-.60*row),.17,
                  facecolor=ORANGE if positive else '#E4EAF2',
                  edgecolor=ORANGE if positive else '#72859C',linewidth=1.5))
dots(4,1.0,7.8,True)
ax.text(5.1,7.7,'4 positive',fontsize=21,color=INK)
dots(18,1.0,6.9,False)
ax.text(1.0,5.85,'18 negative',fontsize=21,color=INK)
dots(5,1.0,2.6,True)
ax.text(5.7,2.5,'5 positive',fontsize=21,color=INK)
dots(6,1.0,1.6,False)
ax.text(1.0,.65,'6 negative',fontsize=21,color=INK)
ax.text(5,4.7,'33 cards = 9 positive + 24 negative',ha='center',fontsize=21,color=INK)
fig.subplots_adjust(left=.01,right=.99,bottom=.01,top=.99)
fig.savefig(FIG/'01_不放回抽取.png',dpi=180)
plt.close(fig)

x = np.arange(10)
q = np.array([comb(9,int(a))*comb(24,22-int(a))/comb(33,22) for a in x])
included = q <= q[4]*(1+1e-7)
pd.DataFrame({'x':x,'probability':q,'included_two_sided':included}).to_csv(
    OUT/'图式源码'/'分布坐标.csv',index=False)
fig, ax = plt.subplots(figsize=(9,7))
ax.bar(x,q,color=[ORANGE if a else '#6192C3' for a in included],width=.73)
ax.bar([4],[q[4]],color=ORANGE,edgecolor=INK,linewidth=2.5,width=.73)
ax.axvline(6,color=INK,ls='--',lw=1.5)
ax.text(6.15,.346,'Mean = 6',fontsize=18,color=INK)
ax.annotate('Observed: x = 4',xy=(4,q[4]),xytext=(.1,.21),fontsize=19,color=INK,
            arrowprops={'arrowstyle':'->','color':INK})
ax.set(xlabel='X: infected in prevention group',ylabel='Probability P(X = x)',
       xticks=x,ylim=(0,.38))
ax.legend(handles=[patches.Patch(color=ORANGE,label='Included in two-sided P'),
                   patches.Patch(color='#6192C3',label='Other tables')],
          loc='upper left',frameon=False,fontsize=16)
# 参数已在中文图注交代，保留数据区空白。
ax.tick_params(labelsize=17)
fig.tight_layout()
fig.savefig(FIG/'02_超几何分布与P值.png',dpi=180)
plt.close(fig)

# 封面保留为独立无姓名PDF，可直接渲染重建PNG。
cover_pdf=OUT/'图式源码'/'课堂封面_无姓名.pdf'
if cover_pdf.is_file():
    cover=fitz.open(cover_pdf)
    cover[0].get_pixmap(matrix=fitz.Matrix(1.8,1.8)).save(FIG/'00_课堂封面.png')
else:
    assert (FIG/'00_课堂封面.png').is_file(), '请保留配套课堂封面PNG。'
print('Built cover, 11 equation images, 2 figures')
