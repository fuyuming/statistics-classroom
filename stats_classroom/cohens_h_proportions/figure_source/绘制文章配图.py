"""重建文章中的确定性理论图和公式；假设教学资料，无随机模拟。"""
from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / '文章配图'
OUT.mkdir(exist_ok=True)
available = {f.name for f in font_manager.fontManager.ttflist}
font = next((s for s in ['Hiragino Sans GB', 'PingFang SC', 'Noto Sans CJK SC', 'Microsoft YaHei', 'SimHei'] if s in available), None)
if font is None:
    raise RuntimeError('需要中文字体：思源黑体或其他中文字体。')
plt.rcParams.update({'font.family':font, 'axes.unicode_minus':False, 'mathtext.fontset':'dejavusans', 'font.size':16})
BLUE = '#0000F5'
INK = '#172A44'
ORANGE = '#D48120'

def formula(name, lines):
    fig = plt.figure(figsize=(7.2, 0.88 * len(lines) + 0.35), dpi=150)
    for i, line in enumerate(lines):
        fig.text(.5, 1 - (i + .65) / (len(lines) + .35), line, ha='center', va='center', fontsize=26, color=INK)
    fig.savefig(OUT / name, facecolor='white')
    plt.close(fig)

formula('公式01_h.png', [r'$g(p)=2\arcsin\sqrt{p}$', r'$\hat h=g(\hat p_1)-g(\hat p_0)$', r'$=1.57080-1.15928=0.41152$'])
formula('公式02_方差与标准误代入.png', [r'$p=0.5:\quad \mathrm{Var}(\hat p)=\frac{0.5(1-0.5)}{100}=0.0025$', r'$SE(\hat p)=\sqrt{0.0025}=0.05$', r'$p=0.1:\quad \mathrm{Var}(\hat p)=\frac{0.1(1-0.1)}{100}=0.0009$', r'$SE(\hat p)=\sqrt{0.0009}=0.03$'])
formula('公式03_局部直线.png', [r'$g(\hat p)\approx g(p)+g\prime(p)(\hat p-p)$', r'$g(\hat p)-g(p)\approx g\prime(p)(\hat p-p)$'])
formula('公式04_方差稳定.png', [r'$g\prime(p)=2\times\frac{1}{\sqrt{1-p}}\times\frac{1}{2\sqrt{p}}$', r'$=\frac{1}{\sqrt{p(1-p)}}$', r'$\mathrm{Var}[g(\hat p)]\approx[g\prime(p)]^2\mathrm{Var}(\hat p)$', r'$=\frac{1}{p(1-p)}\times\frac{p(1-p)}{n}=\frac{1}{n}$'])
formula('公式05_h区间.png', [r'$\hat h\ \overset{\mathrm{approx}}{\sim}\ N\!\left(h,\frac{1}{n_1}+\frac{1}{n_0}\right)$', r'$SE(\hat h)\approx\sqrt{\frac{1}{n_1}+\frac{1}{n_0}}$', r'$Z=\frac{\hat h-h}{SE(\hat h)}\ \overset{\mathrm{approx}}{\sim}\ N(0,1)$', r'$95\%\ \mathrm{CI}:\quad\hat h\ \pm\ 1.96\,SE(\hat h)$'])

# 图1：两个率分别做角变换，再相减得到 h；RD 留在原始率层。
g = lambda x: 2*np.arcsin(np.sqrt(x))
p0, p1 = .30, .50
g0, g1 = g(p0), g(p1)
rd, hval = p1 - p0, g1 - g0
fig = plt.figure(figsize=(7.2, 8.0), dpi=150)
ax = fig.add_axes([0, 0, 1, 1]); ax.set_axis_off()

def card(x, y, w, h, title, value, accent=BLUE, fill='#F5F7FF'):
    patch = plt.matplotlib.patches.FancyBboxPatch(
        (x, y), w, h, boxstyle='round,pad=0.012,rounding_size=0.018',
        transform=ax.transAxes, facecolor=fill, edgecolor=accent, linewidth=2.2)
    ax.add_patch(patch)
    ax.text(x + w/2, y + h*.69, title, transform=ax.transAxes,
            ha='center', va='center', fontsize=18, color=INK, fontweight='bold')
    ax.text(x + w/2, y + h*.31, value, transform=ax.transAxes,
            ha='center', va='center', fontsize=23, color=accent, fontweight='bold')

ax.text(.5, .965, '两个率分别变换，再相减得到 h', transform=ax.transAxes,
        ha='center', va='top', fontsize=25, color=BLUE, fontweight='bold')
ax.text(.5, .925, r'$g(p)=2\arcsin\sqrt{p}$', transform=ax.transAxes,
        ha='center', va='top', fontsize=21, color=INK)
ax.text(.5, .835, '原始率相减：  RD = 0.50 - 0.30 = 0.20（20个百分点）',
        transform=ax.transAxes, ha='center', va='center', fontsize=18, color=ORANGE,
        bbox=dict(boxstyle='round,pad=.35', facecolor='#FFF9F0', edgecolor=ORANGE, linewidth=1.5))
card(.08, .66, .34, .12, r'新方案  $\hat p_1$', '0.50（50%）', BLUE)
card(.58, .66, .34, .12, r'常规对照  $\hat p_0$', '0.30（30%）', ORANGE, '#FFF9F0')

# 两条独立的变换箭头，避免把 RD 误读为 h 的直接输入。
for x in (.25, .75):
    ax.annotate('', xy=(x, .505), xytext=(x, .65), xycoords=ax.transAxes,
                arrowprops=dict(arrowstyle='-|>', color=BLUE, lw=2.4))
ax.text(.29, .575, '各自做 g 变换', transform=ax.transAxes, ha='left', va='center',
        fontsize=16, color=INK)
ax.text(.79, .575, '各自做 g 变换', transform=ax.transAxes, ha='left', va='center',
        fontsize=16, color=INK)
card(.08, .38, .34, .12, r'$g(\hat p_1)$', '1.57080 弧度', BLUE)
card(.58, .38, .34, .12, r'$g(\hat p_0)$', '1.15928 弧度', ORANGE, '#FFF9F0')

ax.annotate('', xy=(.5, .225), xytext=(.5, .355), xycoords=ax.transAxes,
            arrowprops=dict(arrowstyle='-|>', color=BLUE, lw=2.6))
ax.text(.5, .295, r'新方案减对照：$g(\hat p_1)-g(\hat p_0)$', transform=ax.transAxes,
        ha='center', va='center', fontsize=19, color=INK,
        bbox=dict(boxstyle='round,pad=.28', facecolor='white', edgecolor=BLUE, linewidth=1.4))
patch = plt.matplotlib.patches.FancyBboxPatch(
    (.10, .07), .80, .14, boxstyle='round,pad=0.014,rounding_size=0.025',
    transform=ax.transAxes, facecolor='#EEF0FF', edgecolor=BLUE, linewidth=2.8)
ax.add_patch(patch)
ax.text(.5, .16, r'$\hat h=1.57080-1.15928=0.41152$', transform=ax.transAxes,
        ha='center', va='center', fontsize=21, color=BLUE, fontweight='bold')
ax.text(.5, .105, 'Cohen h（弧度尺度的效应大小）', transform=ax.transAxes,
        ha='center', va='center', fontsize=17, color=INK)
fig.savefig(OUT/'01_两率到h_图解.png', facecolor='white')
plt.close(fig)
with (OUT/'图1坐标.csv').open('w', encoding='utf-8', newline='') as f:
    w = csv.writer(f, lineterminator='\n')
    w.writerow(['group','p','g_p','metric','value'])
    w.writerow(['new_plan',p1,g1,'RD',rd])
    w.writerow(['control',p0,g0,'RD',rd])
    w.writerow(['difference','', '', 'h', hval])

# 图2：上下两个相同坐标范围的面板；RD是横向差，h是纵向差。
p = np.linspace(.001,.999,600)
g = lambda x: 2*np.arcsin(np.sqrt(x))
examples = [(.30,.50,BLUE,'课件原例：30% → 50%'),(.10,.30,ORANGE,'变式：10% → 30%')]
fig, axes = plt.subplots(2,1,figsize=(7.2,9.8),dpi=150,sharex=True,sharey=True)
coord_rows=[]
for ax,(p0,p1,color,title) in zip(axes,examples):
    y0,y1 = g(p0),g(p1); hval=y1-y0; rd=p1-p0
    ax.plot(p,g(p),color=INK,lw=3,zorder=1)
    ax.scatter([p0,p1],[y0,y1],color=color,s=65,zorder=4)
    # 投影虚线：端点分别投影到横向 RD 箭头和纵向 h 箭头。
    y_rd=max(.10,y0-.22); x_h=min(.88,p1+.12)
    ax.vlines([p0,p1],y_rd,[y0,y1],colors=color,linestyles=(0,(3,3)),lw=1.5,alpha=.8,zorder=2)
    ax.hlines([y0,y1],[p0,p1],x_h,colors=color,linestyles=(0,(3,3)),lw=1.5,alpha=.8,zorder=2)
    ax.annotate('',xy=(p1,y_rd),xytext=(p0,y_rd),
                arrowprops=dict(arrowstyle='<->',color=color,lw=2.4),zorder=3)
    ax.annotate('',xy=(x_h,y1),xytext=(x_h,y0),
                arrowprops=dict(arrowstyle='<->',color=color,lw=2.4),zorder=3)
    ax.text(max(.23,(p0+p1)/2),y_rd-.08,f'RD = {rd:.2f}（20个百分点）',ha='center',va='top',color=color,fontsize=15)
    ax.text(x_h+.025,(y0+y1)/2,f'h = {hval:.3f}',ha='left',va='center',color=color,fontsize=16,fontweight='bold')
    ax.text(.03,.93,title,transform=ax.transAxes,color=color,fontsize=18,fontweight='bold',va='top')
    ax.set_xlim(0,1); ax.set_ylim(0,np.pi)
    ax.set_xticks([0,.2,.4,.6,.8,1],['0%','20%','40%','60%','80%','100%'])
    ax.grid(alpha=.14); ax.spines[['top','right']].set_visible(False)
    coord_rows.append([p0,p1,round(rd,2),y0,y1,hval])
axes[0].set_ylabel('g(p)（弧度）',fontsize=16)
axes[1].set_ylabel('g(p)（弧度）',fontsize=16)
axes[1].set_xlabel('概率 p（原始尺度）',fontsize=16)
fig.suptitle('RD 是横向差；Cohen h 是纵向角度差',color=BLUE,fontsize=21,y=.985)
fig.subplots_adjust(left=.15,right=.91,bottom=.075,top=.95,hspace=.24)
fig.savefig(OUT/'02_角变换.png',facecolor='white')
plt.close(fig)
with (OUT/'图2坐标.csv').open('w',encoding='utf-8',newline='') as f:
    w=csv.writer(f, lineterminator='\n'); w.writerow(['p0','p1','RD','g0','g1','h']); w.writerows(coord_rows)

# 图3：h的正态近似区间，与Pearson P分开计算。
h=2*np.arcsin(np.sqrt(.5))-2*np.arcsin(np.sqrt(.3))
fig,ax=plt.subplots(figsize=(7.2,5.7),dpi=150)
rows=[]
for y,n in zip([2,1,0],[20,100,1000]):
    se=np.sqrt(2/n);lo=h-stats.norm.ppf(.975)*se;hi=h+stats.norm.ppf(.975)*se
    pv=stats.chi2_contingency([[n*.5,n*.5],[n*.3,n*.7]],correction=False).pvalue
    ax.plot([lo,hi],[y,y],color=BLUE,lw=4);ax.scatter([h],[y],color=BLUE,s=60,zorder=3)
    ax.text((lo+hi)/2,y+.20,f'每组{n}人  [{lo:.3f}, {hi:.3f}]',ha='center',fontsize=18,color=INK)
    rows.append([n,h,se,lo,hi,pv])
ax.axvline(0,color='#75869A',ls='--');ax.set_xlim(-.3,1.1);ax.set_ylim(-.5,2.55);ax.set_yticks([])
ax.set_xlabel('h（角变换后的差值）',fontsize=19,labelpad=12)
ax.set_title('样本率相同，人数改变估计精度',fontsize=22,color=BLUE,pad=18)
ax.spines[['top','right','left']].set_visible(False);ax.grid(axis='x',alpha=.14)
fig.subplots_adjust(left=.08,right=.97,bottom=.19,top=.84)
fig.savefig(OUT/'03_幅度与精度.png',facecolor='white')
plt.close(fig)
with (OUT/'图3坐标.csv').open('w',encoding='utf-8',newline='') as f:
    w=csv.writer(f, lineterminator='\n');w.writerow(['n_per_group','h','SE','lower','upper','Pearson_P']);w.writerows(rows)
print('5幅公式和3幅图解/数据图已生成；字体：',font)
