# -*- coding: utf-8 -*-
"""第14篇的解释图；由small_sample_proportion.py调用。
使用理论分布作图，不生成或冒充新的受试者数据。
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.mathtext import math_to_image
from matplotlib.font_manager import FontProperties
from scipy.stats import norm, binom


def draw_explanations(figure_dir, output_dir):
    figure_dir, output_dir = Path(figure_dir), Path(output_dir)
    z = norm.ppf(0.975)
    grid = np.linspace(-4, 4, 2001)
    density = norm.pdf(grid)
    fig, ax = plt.subplots(figsize=(6.4, 4.5))
    ax.plot(grid, density, color='#194e7b', linewidth=2)
    ax.fill_between(grid, 0, density, where=np.abs(grid) <= z, color='#a9d8f1')
    ax.fill_between(grid, 0, density, where=np.abs(grid) >= z, color='#efb09b')
    ax.vlines([-z, z], 0, norm.pdf(z), color='#c24c32', linestyles='--')
    ax.text(0, .16, '中间 95%', ha='center', fontsize=19, color='#164d75')
    for x0, text_x in [(-z, -2.8), (z, 2.8)]:
        ax.annotate('2.5%', xy=(np.sign(x0)*2.3, .015), xytext=(text_x, .13),
                    ha='center', fontsize=15, arrowprops={'arrowstyle':'->','color':'#b54125'})
    ax.set_xticks([-z, 0, z], ['−1.96', '0', '1.96'])
    ax.set(xlim=(-4, 4), ylim=(0, .45), xlabel='标准化后的数值 Z', ylabel='概率密度',
           title='1.96：标准正态分布的97.5%分位点')
    ax.spines[['top','right']].set_visible(False)
    fig.tight_layout()
    fig.savefig(figure_dir/'03_标准正态与196.png', dpi=200, bbox_inches='tight')
    plt.close(fig)
    pd.DataFrame({'z':grid,'density':density}).to_csv(output_dir/'标准正态图数据.csv',index=False)

    # 用连续正态在每个整数格内的面积，与二项概率进行同量纲比较。
    # 第一组p=6/13仅为示意，不能从这一张图断言真实p就是样本率。
    settings=[(13,6/13,'本例的形状演示：n=13，p=6/13'),
              (13,.1,'同样13次，成功率较低：n=13，p=0.1'),
              (100,.1,'保持p=0.1，增加到100次')]
    fig, axes = plt.subplots(3, 1, figsize=(6.4, 11.2))
    records=[]
    for ax, (n, p, title) in zip(axes, settings):
        k=np.arange(n+1)
        mass=binom.pmf(k,n,p)
        mu=n*p; se=np.sqrt(n*p*(1-p))
        normal_cells=norm.cdf(k+.5,mu,se)-norm.cdf(k-.5,mu,se)
        # 100次示例仅显示0—22次，其余右尾极小；图注明确截断。
        max_show=n if n==13 else 22
        mask=k<=max_show
        ax.bar(k[mask],mass[mask],width=.72,color='#99cdeb',label='二项：恰好k次的概率')
        ax.plot(k[mask],normal_cells[mask],color='#ca633a',marker='o',markersize=3,
                linewidth=1.4,label='正态：每个整数格的面积')
        ax.set(title=title,xlabel='成功次数 k',ylabel='概率',xlim=(-.6,max_show+.6))
        ax.set_xticks(np.arange(0,max_show+1,2))
        ax.legend(frameon=False,fontsize=11)
        ax.spines[['top','right']].set_visible(False)
        for kk, pm, nc in zip(k,mass,normal_cells):
            records.append({'n':n,'p':p,'k':kk,'binomial_mass':pm,'normal_cell_area':nc})
    fig.tight_layout(h_pad=2.2)
    fig.savefig(figure_dir/'04_二项与正态近似.png',dpi=200,bbox_inches='tight')
    plt.close(fig)
    pd.DataFrame(records).to_csv(output_dir/'二项与正态逐格概率.csv',index=False)

    # 四张公式图在手机上分开显示，保留真正的分式与根号。
    equations={
        '公式01_样本率正态近似.png':[
            r'$\hat p=\frac{X}{n}\ \overset{\mathrm{approx}}{\sim}\ N\!\left(p,\frac{p(1-p)}{n}\right)$',
            r'$Z=\frac{\hat p-p}{\sqrt{p(1-p)/n}}\ \overset{\mathrm{approx}}{\sim}\ N(0,1)$'],
        '公式02_Wald推导.png':[
            r'$\widehat{\mathrm{SE}}=\sqrt{\frac{\hat p(1-\hat p)}{n}}$',
            r'$\Pr\!\left(-1.96\leq\frac{\hat p-p}{\widehat{\mathrm{SE}}}\leq1.96\right)\approx0.95$',
            r'$\left[\hat p-1.96\widehat{\mathrm{SE}},\ \hat p+1.96\widehat{\mathrm{SE}}\right]$'],
        '公式03_Wilson得分.png':[
            r'$Z(p_0)=\frac{\hat p-p_0}{\sqrt{p_0(1-p_0)/n}}$',
            r'$|Z(p_0)|\leq z\quad\Longleftrightarrow\quad n(\hat p-p_0)^2\leq z^2p_0(1-p_0)$'],
        '公式04_Wilson二次式.png':[
            r'$(n+z^2)p_0^2-(2n\hat p+z^2)p_0+n\hat p^2\leq0$']}
    for filename, expressions in equations.items():
        fig, ax=plt.subplots(figsize=(7.0,1.18*len(expressions)))
        ax.axis('off')
        for i, expression in enumerate(expressions):
            ax.text(.5,1-(i+.5)/len(expressions),expression,ha='center',va='center',fontsize=19,color='#174f7e',transform=ax.transAxes)
        fig.savefig(figure_dir/filename,dpi=240,bbox_inches='tight',facecolor='white')
        plt.close(fig)
