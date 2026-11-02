"""从共享数据与计算表生成文章图；由tukey_q.py自动调用。
独立运行时先运行主脚本。科学曲线与对应CSV一并保存。
图形保留供Spyder查看；批处理核对时可通过MPLBACKEND=Agg选择无界面后端。
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib import font_manager

def make_figures(root, wide, pairs, nu):
    output = root / '文章配图'
    output.mkdir(exist_ok=True)
    available = {f.name for f in font_manager.fontManager.ttflist}
    preferred = ['PingFang SC', 'Microsoft YaHei', 'Noto Sans CJK SC', 'SimHei', 'Arial Unicode MS']
    font = next((f for f in preferred if f in available), 'DejaVu Sans')
    plt.rcParams.update({'font.sans-serif': [font], 'axes.unicode_minus': False,
                         'font.size': 16, 'axes.spines.top': False, 'axes.spines.right': False})
    values = wide.to_numpy(dtype=float)
    means = values.mean(axis=0)
    colors = ['#2759a5', '#df8540', '#15998f']
    fig, ax = plt.subplots(figsize=(6.5, 5.5), layout='constrained')
    for i in range(3):
        offsets = np.linspace(-0.12, 0.12, 10)
        ax.scatter(i + offsets, values[:, i], s=44, color=colors[i], alpha=.65)
        ax.hlines(means[i], i-.28, i+.28, color=colors[i], linewidth=3)
        ax.text(i, 849, f'均值 {means[i]:.0f}', ha='center', fontsize=15)
    ax.annotate('', xy=(2.65,756), xytext=(2.65,666),
                arrowprops={'arrowstyle':'<->', 'color':'#26364a', 'lw':1.8})
    ax.text(2.72,711,'极差\n90',va='center',fontsize=15)
    ax.axhline(756, xmax=.87, color='#aaa', ls=':', lw=1)
    ax.axhline(666, xmax=.87, color='#aaa', ls=':', lw=1)
    ax.set(xticks=[0,1,2], xticklabels=['品种1','品种2','品种3'], ylabel='地块产量（原题数值）',
           xlim=(-.5,3.22), ylim=(575,875), title='三个均值中，最大差异为90')
    fig.savefig(output/'01-产量与均值极差.png', dpi=240)

    x = np.linspace(0,7,281)
    coordinate = pd.DataFrame({'q':x})
    fig, ax = plt.subplots(figsize=(6.5,5.3),layout='constrained')
    for j,k in enumerate([2,3,5]):
        tail = stats.studentized_range.sf(x,k,nu)
        cutoff = stats.studentized_range.ppf(.95,k,nu)
        coordinate[f'tail_k{k}'] = tail
        ax.plot(x,tail,color=colors[j],lw=2.5,label=f'k＝{k}；q*＝{cutoff:.4f}')
        ax.scatter([cutoff],[.05],color=colors[j],s=38,zorder=4)
        ax.vlines(cutoff,0,.05,color=colors[j],lw=1.5)
    ax.axhline(.05,color='#555',ls='--',lw=1)
    ax.text(5.55,.082,'α＝0.05',fontsize=14)
    ax.set(xlabel='标准化极差 q',ylabel='右尾概率 Pr(Q > q)',xlim=(0,7),ylim=(0,1.02),
           title='组数增加，最大差异的门槛提高\n残差自由度均为27')
    ax.legend(loc='upper right',frameon=False,fontsize=14)
    fig.savefig(output/'02-q分布与组数.png',dpi=240)
    coordinate.to_csv(output/'02-q分布曲线.csv',index=False,encoding='utf-8-sig')

    fig, ax = plt.subplots(figsize=(6.5,4.9),layout='constrained')
    for r,row in pairs.iterrows():
        ax.plot([row['lower'],row['upper']],[2-r,2-r],lw=3,color=colors[r])
        ax.scatter([row['difference']],[2-r],s=65,color=colors[r],zorder=3)
        label=f"{row['difference']:.0f}  [{row['lower']:.2f}, {row['upper']:.2f}]"
        ax.text(row['difference'],2-r+.23,label,ha='center',fontsize=14,color=colors[r])
    ax.axvline(0,color='#777',ls='--',lw=1.3)
    ax.set(yticks=[2,1,0],yticklabels=['品种2−1','品种3−1','品种3−2'],
           xlabel='均值差（原题数值）',xlim=(-165,145),ylim=(-.55,2.65),
           title='均值差与Tukey\n同时95%置信区间')
    fig.savefig(output/'03-均值差与同时区间.png',dpi=240)
    pairs.to_csv(output/'03-均值差与同时区间.csv',index=False,encoding='utf-8-sig')
    wide.to_csv(output/'01-地块产量.csv',index=False,encoding='utf-8-sig')

    # 与R qqnorm/qqline相同：n<=10用a=.375，参考线通过第1/3四分位数。
    fig,axes = plt.subplots(3,1,figsize=(6.5,10),layout='constrained')
    qq_rows = []
    n = len(values)
    a = 0.375 if n <= 10 else 0.5
    theoretical = stats.norm.ppf((np.arange(1,n+1)-a)/(n+1-2*a))
    normal_quartiles = stats.norm.ppf([.25,.75])
    for i,ax in enumerate(axes):
        observed = np.sort(values[:,i])
        observed_quartiles = np.quantile(observed,[.25,.75])
        slope = np.diff(observed_quartiles)[0]/np.diff(normal_quartiles)[0]
        intercept = observed_quartiles[0]-slope*normal_quartiles[0]
        reference = intercept+slope*theoretical
        ax.scatter(theoretical,observed,s=36,color=colors[i])
        ax.plot(theoretical,reference,color='#2759a5',lw=1.6)
        ax.set(title=f'品种{i+1}',xlabel='标准正态理论分位数',ylabel='排序后的产量')
        for x,y,line in zip(theoretical,observed,reference):
            qq_rows.append({'group':i+1,'theoretical':x,'observed':y,'reference':line})
    fig.savefig(output/'04-正态QQ诊断.png',dpi=240)
    fig.savefig(root/'运行结果/Python/qq_diagnostics.png',dpi=150)
    pd.DataFrame(qq_rows).to_csv(output/'04-QQ坐标.csv',index=False,encoding='utf-8-sig')
    # 在Spyder保留所有窗口；脚本批处理环境由调用方自行设置后端。
    if plt.get_backend().lower() != 'agg':
        plt.show(block=False)

def make_null_family_figure(root, simulation):
    """绘制全局H0下两种门槛的整组误报教学模拟图，并保存坐标CSV。"""
    output = root / '文章配图'
    output.mkdir(exist_ok=True)
    available = {f.name for f in font_manager.fontManager.ttflist}
    preferred = ['PingFang SC', 'Microsoft YaHei', 'Noto Sans CJK SC', 'SimHei', 'Arial Unicode MS']
    font = next((f for f in preferred if f in available), 'DejaVu Sans')
    plt.rcParams.update({'font.sans-serif': [font], 'axes.unicode_minus': False,
                         'font.size': 16})
    rows = []
    for _, row in simulation.iterrows():
        false_count = int(row['false_positive_trials'])
        no_count = int(row['trials']) - false_count
        for status, count in [('有误报', false_count), ('无误报', no_count)]:
            rows.append({'method': row['method'], 'status': status, 'count': count,
                         'proportion': count / int(row['trials'])})
    coordinates = pd.DataFrame(rows)
    coordinates.to_csv(output / '05-整组误报模拟-Python.csv', index=False, encoding='utf-8-sig')
    # 固定画布区域，文字放在柱条上方，避免短红柱装不下计数。
    fig = plt.figure(figsize=(6.5, 6.4), layout=None)
    fig.set_layout_engine(None)
    ax = fig.add_axes([.08, .25, .86, .49])
    false_counts = simulation['false_positive_trials'].astype(int).to_numpy()
    totals = simulation['trials'].astype(int).to_numpy()
    no_counts = totals - false_counts
    y_positions = [1.1, .1]
    ax.barh(y_positions, false_counts, height=.28, color='#d95f59', label='至少一对误报')
    ax.barh(y_positions, no_counts, left=false_counts, height=.28,
            color='#a8d5ba', label='三对均未误报')
    for label, y, false_count, no_count, total in zip(
            ['逐对t检验','Tukey检验'], y_positions, false_counts, no_counts, totals):
        ax.text(0, y+.23, label, fontsize=16, fontweight='bold', va='bottom')
        ax.text(10000, y+.23, f'误报{false_count:,}轮（{false_count/total:.2%}）',
                ha='right', va='bottom', fontsize=14, color='#b34b46')
        ax.text(false_count+no_count/2, y, f'{no_count:,}轮未误报',
                ha='center', va='center', fontsize=15, color='#244b3a')
    ax.set(xlim=(0,10000), ylim=(-.25,1.65), yticks=[])
    ax.set_xticks([0,5000,10000])
    ax.tick_params(axis='x', labelsize=14)
    ax.spines['left'].set_visible(False)
    ax.set_xlabel('试验轮数（每轮同时比较三对）', fontsize=14, labelpad=10)
    fig.text(.08,.955,'10000轮试验，有多少轮误报？',fontsize=18,fontweight='bold',va='top')
    fig.text(.08,.89,'教学模拟：三组总体均值均为720',fontsize=14,va='top')
    fig.text(.08,.845,'共同方差1842，每组10个独立正态观测',fontsize=14,va='top')
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles,labels,loc='lower left',bbox_to_anchor=(.06,.035),
               frameon=False,ncol=1,fontsize=14)
    for destination in [output]:
        destination.mkdir(exist_ok=True)
        fig.savefig(destination/'05-整组误报模拟-Python.png',dpi=240)
    if plt.get_backend().lower() != 'agg':
        plt.show(block=False)

if __name__ == '__main__':
    root = Path(__file__).resolve().parents[1]
    wide = pd.read_csv(root/'数据/小麦产量.csv',encoding='utf-8-sig')
    pairs = pd.read_csv(root/'运行结果/Python/pairs.csv')
    make_figures(root,wide,pairs,27)
