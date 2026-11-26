"""从共享CSV与明确的总体模型重建图1和图2。需要numpy/pandas/scipy/seaborn/matplotlib。"""
from pathlib import Path
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import font_manager
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / '文章配图'
DATA = ROOT / '配套项目/数据/residual_cases.csv'
FIG.mkdir(exist_ok=True)
sns.set_theme(style='whitegrid')
available = {f.name for f in font_manager.fontManager.ttflist}
family = next((f for f in ['PingFang SC', 'Heiti TC', 'Arial Unicode MS', 'Noto Sans CJK SC'] if f in available), 'DejaVu Sans')
plt.rcParams.update({'font.family': family, 'font.size': 18, 'axes.labelsize': 18,
                     'xtick.labelsize': 16, 'ytick.labelsize': 16,
                     'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.unicode_minus': False})
colors = ['#3155ad', '#008e91', '#ae6373']
names = ['A1', 'A2', 'A3']
mu = np.array([10., 20., 30.])
data = pd.read_csv(DATA)
data['fitted'] = data.groupby(['scenario', 'group'])['value'].transform('mean')
data['residual'] = data['value'] - data['fitted']

def save(fig, name):
    fig.savefig(FIG/name, dpi=220, facecolor='white')
    plt.close(fig)

# 图1为已知总体模型的密度，而不是从15个样本估计密度。
fig, axes = plt.subplots(2, 1, figsize=(6.4, 7.1), constrained_layout=True)
x = np.linspace(2, 38, 1000)
components = np.array([stats.norm.pdf(x, m, 2) for m in mu])
for i, name in enumerate(names):
    axes[0].plot(x, components[i], color=colors[i], label=name, lw=1.8)
axes[0].plot(x, components.mean(axis=0), color='#202b3d', ls='--', lw=2.4, label='等权混合')
axes[0].set_title('A  各组正态，混合后可以多峰', loc='left', fontsize=17, pad=12)
axes[0].set(xlabel='压缩模量（MPa）', ylabel='密度', ylim=(0,.26))
axes[0].legend(ncol=4, fontsize=11, loc='upper center', frameon=False)
err = np.linspace(-8, 8, 700)
axes[1].plot(err, stats.norm.pdf(err, 0, 2), color=colors[0], lw=2.4)
axes[1].set_title('B  减去各自总体均数后的真实误差', loc='left', fontsize=17, pad=12)
axes[1].set(xlabel='真实误差（MPa）', ylabel='密度')
save(fig, '01_组内正态与混合分布.png')

def residual_panel(ax, scenario):
    d = data.loc[data['scenario'] == scenario]
    for i, name in enumerate(names):
        g = d.loc[d['group'] == name]
        ax.scatter(g['fitted'], g['residual'], color=colors[i], s=57, label=name, alpha=.9)
    ax.axhline(0, color='#66717d', ls='--', lw=1.2)
    ax.set(xlabel='拟合组均数（MPa）', ylabel='残差（MPa）', xlim=(5, 34), ylim=(-9.2, 9.2))
    ax.legend(ncol=3, fontsize=12, frameon=False, loc='upper center')

fig, axes = plt.subplots(2, 1, figsize=(6.4, 7.4), constrained_layout=True)
residual_panel(axes[0], 'A')
axes[0].set_title('A  场景A：残差与拟合值', loc='left', fontsize=17, pad=12)
a = data.loc[data['scenario'] == 'A'].sort_values('residual').copy()
N = len(a)
a['theoretical'] = stats.norm.ppf((np.arange(1, N+1)-.5)/N)
for i, name in enumerate(names):
    g = a.loc[a['group'] == name]
    axes[1].scatter(g['theoretical'], g['residual'], color=colors[i], s=57, label=name)
# 与R qqline默认相同：样本分位数type=7，直线经过25%和75%分位对应点。
yq = np.quantile(a['residual'], [.25, .75], method='linear')
zq = stats.norm.ppf([.25, .75])
slope = (yq[1]-yq[0])/(zq[1]-zq[0])
intercept = yq[0]-slope*zq[0]
xx = np.array([a['theoretical'].min(), a['theoretical'].max()])
axes[1].plot(xx, intercept+slope*xx, color='#66717d', ls='--', lw=1.4)
axes[1].set_title('B  合并15个残差的正态Q–Q图', loc='left', fontsize=17, pad=12)
axes[1].set(xlabel='标准正态理论分位数', ylabel='残差（MPa）')
save(fig, '02_残差图与合并QQ图.png')

print('已按共享CSV与给定总体参数生成图1、图2。')
