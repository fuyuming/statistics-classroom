"""配套绘图模块。主脚本最后调用 draw_figures；初学者可稍后阅读。
每张图先保存坐标 CSV，再保存 PNG；文字说明和读图练习在主脚本。
"""
from pathlib import Path
import numpy as np
from scipy import stats
import pandas as pd
import matplotlib.pyplot as plt


def draw_figures(root):
    """root 是项目目录；这里只处理绘图，不改变教学原始数据。"""
    OUT = root / '运行结果' / 'Python'
    FIG = root / '文章配图' / 'Python'
    FIG.mkdir(parents=True, exist_ok=True)
    dfs = np.array([1, 5, 10, 30, 100])
    q = stats.t.ppf(0.975, dfs)

    def save(name, names, columns):
        # zip 将每列数组转为逐行记录，columns 的顺序与 names 对应。
        table = pd.DataFrame(dict(zip(names, columns)))
        table.to_csv(OUT / name, index=False, encoding='utf-8-sig')

    def finish(fig, name):
        fig.savefig(FIG / name, dpi=180)
        # 保留窗口，便于 Spyder 的 Plots 窗格查看；批处理由主脚本关闭。

    # 2. 三张分布图及全部坐标，英文标注避免依赖本机中文字体。
    colors = ['#178C88','#4144A5','#E36854','#333333']
    x = np.linspace(-5,5,1001)
    y = [stats.norm.pdf(x)] + [stats.t.pdf(x,k) for k in [1,5,30]]
    save('t_normal_curves.csv',['x','normal','t_df1','t_df5','t_df30'],[x]+y)
    fig,ax=plt.subplots(figsize=(10,6),layout='constrained')
    for vals,label,c,ls in zip(y,['N(0,1)','t: df=1','t: df=5','t: df=30'],['#222222']+colors[:3],['--','-','-','-']):
        ax.plot(x,vals,label=label,color=c,ls=ls)
    ax.set(xlabel='Statistic value',ylabel='Density',title='Student t and standard normal distributions')
    ax.legend()
    finish(fig,'02-t_normal.png')
    x=np.linspace(.001,55,3001)
    ks=[1,2,3,5,10,30]
    save('chisq_curves.csv',['x']+[f'df_{k}' for k in ks],[x]+[stats.chi2.pdf(x,k) for k in ks])
    fig,axes=plt.subplots(2,1,figsize=(10,9),layout='constrained')
    for ax, subset in zip(axes,[[1,2],[3,5,10,30]]):
        for k,c,ls in zip(subset,colors,['-','--',':','-.']): ax.plot(x,stats.chi2.pdf(x,k),label=f'df={k}',color=c,ls=ls)
        ax.set(xlabel='Chi-square value',ylabel='Density')
        ax.legend()
    axes[0].set(xlim=(0,8),ylim=(0,1.5),title='df=1, 2: df=1 density diverges at zero; upper part clipped')
    axes[1].set_title('Chi-square: increasing degrees of freedom')
    finish(fig,'03_chisq.png')
    x=np.linspace(.001,6,2001)
    pairs=[(3,10),(5,10),(10,10),(5,5),(5,30)]
    save('f_curves.csv',['x']+[f'F_{a}_{b}' for a,b in pairs],[x]+[stats.f.pdf(x,a,b) for a,b in pairs])
    fig,axes=plt.subplots(2,1,figsize=(10,9),layout='constrained')
    for ax,ps,title in zip(axes,[pairs[:3],[pairs[3],pairs[1],pairs[4]]],['Fixed denominator df=10','Fixed numerator df=5']):
        for (a,b),c,ls in zip(ps,colors,['-','--',':']): ax.plot(x,stats.f.pdf(x,a,b),label=f'F({a},{b})',color=c,ls=ls)
        ax.axvline(1,color='gray',ls=':',lw=1)
        ax.set(xlabel='F value',ylabel='Density',title=title)
        ax.legend()
    finish(fig,'04_f.png')
    
    # 构造关系：数学条件写在每一行，分布变量不能与原始资料混淆。
    fig,ax=plt.subplots(figsize=(10,8),layout='constrained')
    ax.set_axis_off()
    ax.set_title('Normal, chi-square, t and F: construction',fontsize=17,pad=20)
    rows=[
     (r'Independent $Z_1,\ldots,Z_k\sim N(0,1)$',r'$U=\sum_{i=1}^{k} Z_i^2\sim\chi^2(k)$'),
     (r'$Z\sim N(0,1),\ U\sim\chi^2(\nu)$; independent',r'$T=\frac{Z}{\sqrt{U/\nu}}\sim t(\nu)$'),
     (r'$U\sim\chi^2(d_1),\ V\sim\chi^2(d_2)$; independent',r'$F=\frac{U/d_1}{V/d_2}\sim F(d_1,d_2)$')]
    for y0,(inputs,formula) in zip([.84,.54,.24],rows):
        ax.text(.5,y0,inputs,ha='center',va='center',fontsize=13,transform=ax.transAxes)
        ax.text(.5,y0-.12,formula,ha='center',va='center',fontsize=22,color='#167d80',
                bbox=dict(boxstyle='round,pad=.5',facecolor='#f0f7f6',edgecolor='#c3ddda'),transform=ax.transAxes)
    ax.text(.5,.01,r'$T\sim t(\nu)\quad\Rightarrow\quad T^2\sim F(1,\nu)$',ha='center',fontsize=16,transform=ax.transAxes)
    finish(fig,'01_construction.png')
    correct=2*stats.t.sf(q,dfs)
    wrong=2*stats.t.sf(1.96,dfs)
    assert np.allclose(correct,.05,atol=1e-10)
    save('error_rates.csv',['df','n_one_sample','correct_t','wrong_1_96'],[dfs,dfs+1,correct,wrong])
    fig,ax=plt.subplots(figsize=(10,6),layout='constrained')
    pos=np.arange(len(dfs))
    for offset,vals,label,c in [(-.19,correct,'Correct t cutoff','#178C88'),(.19,wrong,'Wrong cutoff: 1.96','#E36854')]:
        bars=ax.bar(pos+offset,100*vals,width=.38,label=label,color=c)
        ax.bar_label(bars,fmt='%.2f',padding=3,fontsize=10)
    ax.axhline(5,color='gray',ls=':',lw=1)
    ax.set_xticks(pos,[f'df={k}' for k in dfs])
    ax.set(ylim=(0,38),xlabel='Degrees of freedom',ylabel='Type I error (%)',
           title='Exact probabilities under H0: t statistic compared with 1.96')
    ax.legend()
    finish(fig,'05_error_rates.png')
    
