"""绘制第23篇科学图。用共享数据及R实际置换结果生成PNG与坐标CSV。
在配套项目根目录运行代码/plot_dunn.py；也可由dunn_rank.py调用。
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import rankdata
import matplotlib.pyplot as plt
from matplotlib import font_manager

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / '文章配图'
COLORS = {'A': '#2367a5', 'B': '#d86c37', 'C': '#319483'}

def font_setup():
    available = {f.name for f in font_manager.fontManager.ttflist}
    names = ['PingFang SC', 'Noto Sans CJK SC', 'Microsoft YaHei', 'SimHei', 'DejaVu Sans']
    plt.rcParams.update({'font.family': [n for n in names if n in available],
        'font.size': 12, 'axes.unicode_minus': False, 'axes.spines.top': False,
        'axes.spines.right': False, 'axes.titleweight': 'bold', 'figure.facecolor': 'white'})

def save(fig, name):
    OUT.mkdir(exist_ok=True)
    fig.savefig(OUT / name, dpi=260, facecolor='white')
    plt.close(fig)

def offsets(values):
    # 同值点对称错开，不用随机抖动改变图形。
    result = np.zeros(len(values))
    for value in np.unique(values):
        idx = np.flatnonzero(np.asarray(values) == value)
        result[idx] = np.linspace(-.18, .18, len(idx)) if len(idx)>1 else 0
    return result

def main():
    font_setup()
    OUT.mkdir(exist_ok=True)
    dat = pd.read_csv(ROOT/'数据/小鼠存活天数.csv')
    dat['rank'] = rankdata(dat.days, method='average')
    dat.to_csv(OUT/'01-天数与统一秩-坐标.csv', index=False)
    fig, axes = plt.subplots(2, 1, figsize=(5.6, 7.4))
    fig.subplots_adjust(left=.25, right=.94, top=.91, bottom=.09, hspace=.55)
    fig.suptitle('同一批小鼠：从天数换成秩', fontsize=16, y=.98, weight='bold')
    for ax, field, label in zip(axes, ['days','rank'], ['存活天数（天）','全体33个观察的秩']):
        for i,g in enumerate('ABC'):
            vals=dat.loc[dat.group==g,field].to_numpy()
            ax.scatter(vals,i+offsets(vals),s=34,color=COLORS[g],alpha=.85,zorder=3)
            if field=='rank':
                mean=vals.mean()
                ax.plot([mean,mean],[i-.3,i+.3],color='#202c35',lw=2,zorder=4)
                ax.text(mean,i+.45,f'均秩 {mean:.2f}',ha='center',fontsize=11)
        ax.set_yticks([0,1,2], ['A（11只）','B（10只）','C（12只）'])
        ax.set_ylim(2.8,-.6)
        ax.set_xlabel(label)
        ax.grid(axis='x',color='#e3e8eb',zorder=0)
    axes[0].set_xlim(0,23); axes[0].set_xticks([0,5,10,15,20])
    axes[1].set_xlim(0,35);axes[1].set_xticks([0,5,10,15,20,25,30,35])
    save(fig,'01-天数与统一秩.png')

    perm_path=ROOT/'运行结果/R/permutation_differences.csv'
    if not perm_path.exists():
        raise FileNotFoundError('请先运行R主脚本，生成实际置换结果：'+str(perm_path))
    perm=pd.read_csv(perm_path)['difference'].to_numpy()
    ranks=dat['rank'].to_numpy(); observed=dat.loc[dat.group=='B','rank'].mean()-dat.loc[dat.group=='A','rank'].mean()
    se=np.sqrt(np.var(ranks,ddof=1)*(1/11+1/10))
    fig,axes=plt.subplots(2,1,figsize=(5.6,7.7))
    fig.subplots_adjust(left=.16,right=.95,top=.90,bottom=.09,hspace=.67)
    fig.suptitle('单个秩与平均秩差，波动不同',fontsize=16,y=.98,weight='bold')
    values,counts=np.unique(ranks,return_counts=True)
    axes[0].bar(values,counts,width=.95,color='#4784b6')
    axes[0].set(xlabel='单个秩',ylabel='小鼠只数',xlim=(0,34),ylim=(0,6.4))
    axes[0].set_title('33个秩固定不变',fontsize=13,pad=10)
    axes[0].text(.98,.92,'秩的方差 = 92.65625',transform=axes[0].transAxes,ha='right',fontsize=12)
    bins=np.linspace(-17,17,35)
    c,edges,_=axes[1].hist(perm,bins=bins,color='#41958b',alpha=.9)
    axes[1].axvline(0,color='#718087',ls=':',lw=1)
    axes[1].axvline(observed,color='#d76b36',ls='--',lw=2)
    axes[1].set(xlabel='随机分组后的平均秩差（B−A）',ylabel='模拟轮数',xlim=(-17,17))
    axes[1].set_title('每轮随机分成11、10、12个，重复10000轮',fontsize=12,pad=10)
    axes[1].text(.02,.96,f'理论标准误 {se:.3f}\nR模拟标准差 {perm.std(ddof=1):.3f}',transform=axes[1].transAxes,va='top',fontsize=11)
    axes[1].annotate(f'观察差\n{observed:.2f}',xy=(observed,c.max()*.45),xytext=(15,c.max()*.8),ha='center',fontsize=11,color='#b44f23',arrowprops={'arrowstyle':'->','color':'#b44f23'})
    axes[1].set_ylim(0,c.max()*1.32)
    save(fig,'02-秩与平均秩差的波动.png')
    pd.DataFrame({'rank':values,'count':counts}).to_csv(OUT/'02-单个秩-坐标.csv',index=False)
    pd.DataFrame({'bin_left':edges[:-1],'bin_right':edges[1:],'count':c.astype(int)}).to_csv(OUT/'02-平均秩差-直方图.csv',index=False)
    (OUT/'02-参数.json').write_text(json.dumps({'source':'R/permutation_differences.csv','rounds':len(perm),'n':[11,10,12],'rank_variance':float(np.var(ranks,ddof=1)),'theoretical_SE':float(se),'observed_difference':float(observed),'simulated_SD':float(perm.std(ddof=1))},ensure_ascii=False,indent=2))

    fig,axs=plt.subplots(2,1,figsize=(5.6,5.8))
    fig.subplots_adjust(top=.83,bottom=.13,left=.13,right=.94,hspace=.85)
    fig.suptitle('同样三个“4天”，排序集合变了',fontsize=15,y=.98,weight='bold')
    for ax,positions,title,mean in zip(axs,[[7,8,9],[6,7,8]],['全体33个观察一起排','只取A、B的21个观察重排'],[8,7]):
        ax.set_xlim(4.5,10.5);ax.set_ylim(-.55,1.2)
        ax.spines[['left','bottom']].set_visible(False)
        ax.set_yticks([]);ax.set_xticks([5,6,7,8,9,10]);ax.set_xlabel('排序位置')
        for x in range(5,11):ax.scatter(x,0,s=480,color='#e8edf1',zorder=1)
        for x in positions:
            ax.scatter(x,0,s=480,color='#3784b9',zorder=2)
            ax.text(x,0,'4',ha='center',va='center',color='white',weight='bold',fontsize=15)
        ax.set_title(title,fontsize=13,pad=12)
        ax.text(7.5,.8,'占第'+ '、'.join(map(str,positions))+'位',ha='center',fontsize=12)
        ax.text(7.5,.42,f'共同平均秩 = {mean}',ha='center',fontsize=13,weight='bold',color='#1769aa')
    save(fig,'03-全体秩与两组重编秩.png')
    pd.DataFrame({'ranking_set':['all33','AB21'],'n':[33,21],'first_position':[7,6],'last_position':[9,8],'average_rank':[8,7]}).to_csv(OUT/'03-重编秩-坐标.csv',index=False)
    print('3张科学图及坐标/参数已生成。')

if __name__=='__main__':
    main()
