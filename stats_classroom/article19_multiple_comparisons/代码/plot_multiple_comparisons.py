"""从已计算结果绘图；坐标逐行保存，图形本身不重新计算检验。"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
root = Path(__file__).resolve().parents[1]
out = root / '文章配图'
out.mkdir(exist_ok=True)
plt.rcParams['font.sans-serif'] = ['Arial Unicode MS', 'Noto Sans CJK SC', 'DejaVu Sans']
plt.rcParams['font.size'] = 12
m = np.arange(1, 8)
prob = 1 - .95 ** m
fig, ax = plt.subplots(figsize=(5.2, 4.2))
ax.plot(m, prob, 'o-', color='#1769aa')
ax.axhline(.05, ls='--', color='gray')
ax.annotate('3次：14.26%', (3, prob[2]), xytext=(3.2, .11), fontsize=12)
ax.set(xlabel='独立检验次数', ylabel='至少一次误报的概率', ylim=(0,.34))
ax.set_title('各原假设均为真，检验相互独立\n每次实际误报率为5%', fontsize=13)
ax.set_yticks([0,.05,.1,.15,.2,.25,.3], ['0','5%','10%','15%','20%','25%','30%'])
fig.tight_layout()
fig.savefig(out/'01-多次检验的误报概率.png', dpi=220)
pd.DataFrame({'m':m,'probability':prob}).to_csv(out/'01-多次检验的误报概率.csv', index=False)

results = pd.read_csv(root/'运行结果/Python/pair_results.csv')
fig, axes = plt.subplots(3, 1, figsize=(5.4, 8.5))
pairs = ['物理—药物', '物理—膏药', '药物—膏药']
methods = ['原始P', 'Bonferroni', 'Holm']
colors = ['#777777','#d98c00','#008f6b']
coords = []
for i, ax in enumerate(axes):
    values = results.loc[i, ['raw_p_yates','bonferroni_p','holm_p']].to_numpy(dtype=float)
    ax.axvline(.05, color='#c73737', ls='--', lw=1.3)
    for j, value in enumerate(values):
        y = 2-j
        ax.scatter(value, y, c=colors[j], s=42, zorder=4)
        # 数值放在对应点上方，不与参考线或点重叠。
        ax.text(max(.004, min(value, .10)), y+.24, f'{value:.7g}', fontsize=12, va='bottom')
        coords.append({'pair':results.loc[i,'pair'],'method':methods[j],'p':value})
    ax.set(xlim=(-.004,.18),ylim=(-.35,2.75),title=pairs[i])
    ax.set_yticks([2,1,0],methods)
    ax.set_xticks([0,.05,.10,.15],['0','0.05','0.10','0.15'])
    ax.grid(axis='x',alpha=.16)
    ax.spines[['top','right']].set_visible(False)
axes[-1].set_xlabel('P值（红虚线为0.05）')
fig.tight_layout(h_pad=1.4)
fig.savefig(out/'02-三对比较的校正结果.png',dpi=220)
pd.DataFrame(coords).to_csv(out/'02-三对比较的校正结果.csv', index=False)
# Spyder中图保留在Plots窗口；自动批量运行可设置MPLBACKEND=Agg。
