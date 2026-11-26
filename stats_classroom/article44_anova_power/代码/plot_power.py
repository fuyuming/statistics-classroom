"""从共享规划参数重建文章数值图。UTF-8。
直接运行显示图形；ARTICLE44_SAVE_RESULTS=1时保存PNG和坐标CSV。
图2（两类错误示意图）沿用文章图片，不由本数值绘图脚本生成。
"""
# %% 01 输入、字体与概率函数
from pathlib import Path
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import font_manager
from scipy import stats

try:
    root = Path(__file__).resolve().parents[1]
except NameError:
    root = Path.cwd()
if root.name == '代码':
    root = root.parent
raw = pd.read_csv(root/'数据/planning_scenarios.csv')
settings = pd.read_csv(root/'数据/settings.csv').set_index('setting')['value']
alpha = float(settings['alpha'])
target = float(settings['target_power'])
n_check = int(settings['n_check'])
cases = {}
for name in ['large','small']:
    part = raw[raw.scenario == name].sort_values('group')
    cases[name] = (part.mean_MPa.to_numpy(), float(part.sd_MPa.iloc[0]))
fonts = {item.name for item in font_manager.fontManager.ttflist}
preferred = [name for name in ['Arial Unicode MS','PingFang SC','Noto Sans CJK SC','SimHei'] if name in fonts]
if not preferred:
    raise RuntimeError('请先安装Noto Sans CJK SC等中文字体。')
plt.rcParams.update({'font.sans-serif':preferred,'axes.unicode_minus':False,
                     'font.size':17,'axes.spines.top':False,'axes.spines.right':False})
blue, orange, green, purple = '#366FAD','#E77B28','#479458','#8865A8'
save_results = os.getenv('ARTICLE44_SAVE_RESULTS','0') == '1'
figure_dir = root/'文章配图'
coords_dir = root/'运行结果/绘图坐标'
if save_results:
    figure_dir.mkdir(exist_ok=True)
    coords_dir.mkdir(parents=True,exist_ok=True)


def effect(means,sd):
    return np.sqrt(np.mean((means-means.mean())**2))/sd


def probability(n,f,k=4):
    critical = stats.f.ppf(1-alpha,k-1,k*(n-1))
    if f == 0:
        return stats.f.sf(critical,k-1,k*(n-1))
    return stats.ncf.sf(critical,k-1,k*(n-1),k*n*f*f)


def required_n(f):
    n = 2
    while probability(n,f) < target:
        n += 1
    return n


def finish(fig,name,coordinates,rect=None):
    fig.tight_layout(rect=rect)
    if save_results:
        fig.savefig(figure_dir/name,dpi=220,facecolor='white')
        coordinates.to_csv(coords_dir/name.replace('.png','.csv'),index=False)
    plt.show()


# %% 02 图1：先看原始量纲下的差异和波动，两个面板共用同一坐标
all_means = np.concatenate([cases[name][0] for name in cases])
max_sd = max(sd for means,sd in cases.values())
x = np.linspace(all_means.min()-3*max_sd,all_means.max()+3*max_sd,1000)
fig,axes = plt.subplots(2,1,figsize=(7,7.1),sharex=True,sharey=True)
rows = []
for ax,(name,(means,sd)) in zip(axes,cases.items()):
    for label,mean,color in zip(['A','B','C','D'],means,[blue,orange,green,purple]):
        density = stats.norm.pdf(x,mean,sd)
        ax.plot(x,density,lw=2.4,color=color,label=label)
        for point,height in zip(x,density):
            rows.append({'scenario':name,'group':label,'x_MPa':point,'density':height,
                         'mean_MPa':mean,'sd_MPa':sd})
    mu_text = '、'.join(f'{value:g}' for value in means)
    ax.set_title(f'均数：{mu_text}；SD＝{sd:g}',fontsize=15)
    ax.set_ylabel('密度')
    ax.set_xlim(x.min(),x.max())
fig.suptitle('图1  比较两种均数差异',fontsize=20,y=.98)
handles,labels = axes[0].get_legend_handles_labels()
fig.legend(handles,labels,ncol=4,loc='upper center',bbox_to_anchor=(.5,.95),fontsize=14,frameon=False)
axes[1].set_xlabel('压缩模量 / MPa')
finish(fig,'01_两种总体分布.png',pd.DataFrame(rows),rect=(0,0,1,.96))

# %% 03 图3：在相同组数、alpha和目标功效下，f越小需要的n越多
case_f = {name:effect(means,sd) for name,(means,sd) in cases.items()}
f_min = min(.2,min(case_f.values())*.9)
f_max = max(1.2,max(case_f.values())*1.05)
f_values = np.linspace(f_min,f_max,121)
counts = np.array([required_n(value) for value in f_values])
fig,ax = plt.subplots(figsize=(7,4.9))
ax.plot(f_values,counts,lw=2.5,color=blue)
for name,color,y_fraction in [('large',orange,.30),('small',green,.78)]:
    f_value = case_f[name]
    n = required_n(f_value)
    ax.scatter([f_value],[n],s=70,color=color,zorder=3)
    label = '较大差异' if name == 'large' else '较小差异'
    ax.annotate(f'{label}\nf＝{f_value:.3f}，每组{n}个',xy=(f_value,n),
                xytext=(.48,y_fraction),textcoords='axes fraction',fontsize=15,
                arrowprops={'arrowstyle':'->','color':color},color=color)
ax.set(xlabel='ANOVA效应量 f',ylabel='所需的每组独立试件数',
       title=f'图3  达到{target:.0%}功效需要多少试件？',ylim=(0,counts.max()*1.08))
finish(fig,'03_效应量与样本量.png',pd.DataFrame({'f':f_values,'n_per_group':counts,
       'target_power':target,'alpha':alpha,'groups':4}))

# %% 04 图4：保持小差异总体不变，增加n会提高检出机会
small_mu,small_sd = cases['small']
f_small = case_f['small']
n_plan = required_n(f_small)
sizes = np.arange(2,max(100,int(n_plan*1.5))+1)
powers = np.array([probability(n,f_small) for n in sizes])
fig,ax = plt.subplots(figsize=(7,4.8))
ax.plot(sizes,powers,lw=2.5,color=blue)
ax.axhline(target,color='#888888',ls='--',label=f'目标{target:.0%}')
for n,color,position in [(n_check,orange,(.07,.61)),(n_plan,green,(.53,.65))]:
    value = probability(n,f_small)
    ax.scatter([n],[value],color=color,s=65,zorder=3)
    ax.annotate(f'每组{n}个\n功效{value:.1%}',xy=(n,value),xytext=position,
                textcoords='axes fraction',fontsize=15,
                arrowprops={'arrowstyle':'->','color':color},color=color)
ax.set(xlabel='每组独立试件数 n',ylabel='功效',title='图4  同一效应，增加独立样本量',ylim=(0,1.03))
ax.legend(loc='lower right',fontsize=14,frameon=False)
finish(fig,'04_小差异功效曲线.png',pd.DataFrame({'n_per_group':sizes,'power':powers,
       'f':f_small,'sigma_MPa':small_sd,'alpha':alpha}))

