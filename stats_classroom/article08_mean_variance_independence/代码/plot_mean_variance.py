"""本课正文图。坐标由共享CSV产生；不强制后端，Spyder可查看图窗。"""
import matplotlib.pyplot as plt
from matplotlib import font_manager

def make_figures(root, examples, pairs):
    available = {font.name for font in font_manager.fontManager.ttflist}
    chinese_font = next((f for f in ['PingFang SC', 'Arial Unicode MS', 'Noto Sans CJK SC',
                                    'Microsoft YaHei', 'SimHei'] if f in available), None)
    chinese = chinese_font is not None
    plt.rcParams.update({'font.family': chinese_font or 'DejaVu Sans', 'font.size': 13,
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'axes.unicode_minus': False, 'figure.facecolor': 'white'})
    folder = root / '文章配图'
    folder.mkdir(exist_ok=True)
    teal, orange, navy = '#167d80', '#cc773e', '#18394b'
    fig, axes = plt.subplots(3, 1, figsize=(8.4, 7.2), sharex=True)
    fig.subplots_adjust(top=.84, bottom=.12, hspace=.62, left=.12, right=.94)
    fig.suptitle('整体位置与分散程度，可以分别变化' if chinese else 'Position and spread can change separately',
                 fontsize=19, color=navy, y=.96)
    for i, ax in enumerate(axes):
        row = examples.iloc[i]
        xs = row[['x1','x2','x3']].to_numpy(dtype=float)
        ax.axhline(0, color='#ccd9dd', lw=1)
        ax.axvline(row['mean'], color=orange, ls='--', lw=1.5)
        ax.scatter(xs, [0,0,0], s=130, color=teal, zorder=3)
        for value in xs:
            ax.text(value, .21, str(int(value)), ha='center', color=navy)
        title = f"{row['group']}   "
        title += (f"均值 = {row['mean']:g}，样本方差 = {row['variance']:g}" if chinese
                  else f"Mean = {row['mean']:g}, sample variance = {row['variance']:g}")
        ax.set_title(title, loc='left', fontsize=14, pad=12)
        ax.set_ylim(-.45,.65)
        ax.set_xlim(4,24)
        ax.set_yticks([])
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
    axes[-1].set_xlabel('观测值（无量纲；三组共用刻度）' if chinese else 'Value (dimensionless; shared scale)')
    fig.savefig(folder / '01_位置与分散.png', dpi=180)

    fig, axes = plt.subplots(3, 1, figsize=(8.4, 13))
    fig.subplots_adjust(top=.90,bottom=.07,left=.15,right=.93,hspace=.60)
    fig.suptitle('协方差为零，依赖关系仍可能存在' if chinese else 'Zero covariance does not ensure independence',
                 fontsize=19, color=navy, y=.97)
    fig.text(.5,.932,'三个模型的理论协方差均为0；每图1200个模拟点' if chinese else
             'All theoretical covariances = 0; 1,200 simulated pairs per panel',ha='center',fontsize=12)
    titles_zh=['独立正态：联合正态，且独立', '平方关系：不相关，但不独立', '随机变号：各自正态，但不联合正态']
    titles_en=['Independent normals: jointly normal and independent',
               'Y = X²: uncorrelated but dependent',
               'Y = ±X: normal marginals, not jointly normal']
    keys=['independent','square','random_sign']
    for i, ax in enumerate(axes):
        ax.scatter(pairs.x, pairs[f'{keys[i]}_y_calculated'], s=9, alpha=.35, color=teal)
        ax.set_title((titles_zh if chinese else titles_en)[i],loc='left',fontsize=14,pad=14)
        ax.set_xlim(-3.4,3.4)
        ax.set_ylim((-.3,10) if i == 1 else (-3.4,3.4))
        ax.set_xlabel('X')
        ax.set_ylabel(['Y','Y = X²','Y = ±X'][i])
        ax.grid(alpha=.18)
    fig.savefig(folder / '02_不相关与独立.png',dpi=180)
