"""图：只改变一个观测值时，均值M和差值D的变化方向。
人工教学例子，不是随机样本；箭头用于解释系数符号，不是协方差估计。
可单独运行，也由主脚本自动调用。输入为数据/direction_examples.csv。
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import font_manager

def make_direction_figure(root):
    data = pd.read_csv(root / '数据' / 'direction_examples.csv')
    # 每次只改变一个观测值，另一个保持10不变。
    data['M'] = (data.x1 + data.x2) / 2
    data['D'] = data.x1 - data.x2
    assert np.allclose(data.M, [10, 11, 10, 11])
    assert np.allclose(data.D, [0, 2, 0, -2])
    out = root / '运行结果' / 'Python'
    out.mkdir(parents=True, exist_ok=True)
    data.to_csv(out / 'direction_coordinates.csv', index=False)
    print('变化方向示例（人工数据，不是协方差估计）：')
    print(data.to_string(index=False))
    fonts = {f.name for f in font_manager.fontManager.ttflist}
    font = next((f for f in ['PingFang SC','Arial Unicode MS','Noto Sans CJK SC',
                             'Microsoft YaHei','SimHei'] if f in fonts), None)
    zh = font is not None
    plt.rcParams.update({'font.family': font or 'DejaVu Sans', 'font.size': 15,
                         'axes.unicode_minus': False})
    teal, orange, navy = '#167d80', '#cc773e', '#18394b'
    fig, axes = plt.subplots(2, 1, figsize=(8.4, 11.7))
    fig.subplots_adjust(top=.82,bottom=.24,left=.18,right=.89,hspace=.72)
    fig.suptitle('都推高均值，对差值的影响却相反' if zh else
                 'Both raise the mean; their effects on D oppose',
                 fontsize=21,color=navy,y=.975)
    fig.text(.5,.932,'M = (X₁ + X₂) / 2       D = X₁ − X₂',ha='center',fontsize=18,color=navy)
    fig.text(.5,.895,'从同一起点出发，每次只增加一个观测值' if zh else
             'Same starting point; change only one observation at a time',ha='center',fontsize=14)
    for i, (scenario, color) in enumerate([('x1_increases',teal),('x2_increases',orange)]):
        ax = axes[i]
        frame = data[data.scenario==scenario]
        start, end = frame.iloc[0], frame.iloc[1]
        if zh:
            title = ['① 只增加 X₁：10 → 12，X₂ 固定为 10',
                     '② 只增加 X₂：10 → 12，X₁ 固定为 10'][i]
        else:
            title = ['Only X1: 10 to 12; X2 stays at 10',
                     'Only X2: 10 to 12; X1 stays at 10'][i]
        ax.set_title(title,fontsize=15,loc='left',pad=17,color=color)
        ax.set_xlim(9.7,11.6)
        ax.set_ylim(-2.8,2.8)
        ax.set_xticks([10,10.5,11,11.5])
        ax.set_yticks([-2,0,2])
        ax.set_xlabel('均值 M（两幅图共用刻度）' if zh else 'Mean M (same scales)',fontsize=14)
        ax.set_ylabel('差值 D' if zh else 'Difference D',fontsize=14)
        ax.grid(alpha=.18)
        ax.spines[['top','right']].set_visible(False)
        ax.annotate('',xy=(end.M,end.D),xytext=(start.M,start.D),
                    arrowprops={'arrowstyle':'-|>','color':color,'lw':3,'mutation_scale':23})
        ax.scatter([start.M],[start.D],s=70,color=navy,zorder=5)
        ax.scatter([end.M],[end.D],s=90,color=color,zorder=5)
        ax.text(9.78,.35,'起点 (10, 0)' if zh else 'Start (10, 0)',fontsize=13,color=navy)
        ax.text(11.06,end.D,f'(11, {int(end.D)})',va='center',fontsize=14,color=color)
        ax.text(.98,.14 if i==0 else .86,
                'M +1，D +2' if i==0 else 'M +1，D −2',transform=ax.transAxes,
                ha='right',va='center',fontsize=16,color=color,weight='bold')
    fig.text(.5,.161,'独立时，两项对 Cov(M,D) 的贡献为：' if zh else
             'Under independence, the covariance contributions are:',ha='center',fontsize=15,color=navy)
    fig.text(.5,.119,r'$\frac{1}{2}\mathrm{Var}(X_1)\quad-\quad\frac{1}{2}\mathrm{Var}(X_2)$',
             ha='center',fontsize=23,color=navy)
    fig.text(.5,.075,'同方差 → 两项大小相等、符号相反 → 总和为 0' if zh else
             'Equal variances: equal magnitudes, opposite signs, sum = 0',ha='center',fontsize=15,color=teal)
    fig.text(.5,.031,'箭头展示变化方向；协方差的结论来自公式与条件。' if zh else
             'Arrows illustrate directions; the covariance claim follows from the formula.',
             ha='center',fontsize=12,color='#657a7d')
    (root / '文章配图').mkdir(exist_ok=True)
    fig.savefig(root / '文章配图' / '03_均值与差值的变化方向.png',dpi=180)

if __name__ == '__main__':
    make_direction_figure(Path(__file__).resolve().parents[1])
