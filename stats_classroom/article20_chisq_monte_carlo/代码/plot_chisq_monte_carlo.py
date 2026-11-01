"""从所选语言的模拟记录生成三图；保留对应计数、直方图和累计P的坐标。"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[1]
LANGUAGE = globals().get('SOURCE_LANGUAGE', 'Python')
source = ROOT / '运行结果' / LANGUAGE / 'simulation.csv'
out = ROOT / '文章配图'
out.mkdir(exist_ok=True)
d = pd.read_csv(source)
q = d['Q'].to_numpy()
summary = pd.read_csv(ROOT / '运行结果' / LANGUAGE / 'summary.csv').iloc[0]
q_obs = float(summary['Qobs'])
p_exact = float(summary['exact_pearson_P'])
cutoff = q_obs * (1 - 64 * np.finfo(float).eps)
extreme = q >= cutoff
B = len(q)
b = int(extreme.sum())
plt.rcParams['font.sans-serif'] = ['Arial Unicode MS', 'Noto Sans CJK SC', 'DejaVu Sans']
plt.rcParams['font.size'] = 12
blue = '#2477b5'
red = '#c94a3a'
# 对应图：原表与一张模拟表，固定的边际和变化的格内数采用不同颜色。
observed = pd.read_csv(ROOT / '数据/ace_genotype.csv')
one = pd.read_csv(ROOT / '运行结果' / LANGUAGE / 'one_table.csv')
original_counts = observed[['DD', 'ID', 'II']].to_numpy(dtype=int)
simulated_counts = one[['DD', 'ID', 'II']].to_numpy(dtype=int)
assert np.array_equal(original_counts.sum(0), simulated_counts.sum(0))
assert np.array_equal(original_counts.sum(1), simulated_counts.sum(1))
fig, ax = plt.subplots(figsize=(5.4, 6.5))
ax.set_axis_off()
fig.subplots_adjust(left=.02, right=.98, top=.98, bottom=.02)
ax.text(.5, .975, '固定边际，改变格内计数', ha='center', va='top', fontsize=16, weight='bold', color=blue)
blue_fill = '#e2effb'
orange_fill = '#fff0da'
coords = []
for counts, title, y, label in [(original_counts,'观察表',.59,'observed'),
                                (simulated_counts,'一张模拟表',.18,'simulated')]:
    row_sums = counts.sum(axis=1)
    col_sums = counts.sum(axis=0)
    display = [['DN', *counts[0], row_sums[0]],
               ['无DN', *counts[1], row_sums[1]],
               ['合计', *col_sums, int(counts.sum())]]
    ax.text(.5,y+.285,title,ha='center',fontsize=14,weight='bold')
    table = ax.table(cellText=display, colLabels=['组别','DD','ID','II','合计'],
                     cellLoc='center', bbox=[.03,y,.94,.255],
                     colWidths=[.20,.19,.19,.19,.23])
    table.auto_set_font_size(False)
    table.set_fontsize(15)
    for (i,j), cell in table.get_celld().items():
        cell.set_edgecolor('white')
        cell.set_linewidth(2)
        if i==0:
            cell.set_facecolor('#f1f4f7')
            cell.get_text().set_weight('bold')
        elif i==3 or j==4:
            cell.set_facecolor(blue_fill)
            cell.get_text().set_color(blue)
            cell.get_text().set_weight('bold')
        elif j==0:
            cell.set_facecolor('#f1f4f7')
        else:
            cell.set_facecolor(orange_fill)
            cell.get_text().set_color('#9b5400')
    for row in display:
        coords.append([label,*row])
ax.annotate('', xy=(.5,.545), xytext=(.5,.58),
            arrowprops={'arrowstyle':'-|>', 'color':'#687a88','lw':2})
ax.text(.5,.515,'按原假设随机重新分组',ha='center',fontsize=13,color='#465762')
ax.text(.5,.115,'蓝色：行列合计保持不变',ha='center',fontsize=13,color=blue,weight='bold')
ax.text(.5,.06,'橙色：格内人数随抽签结果变化',ha='center',fontsize=13,color='#9b5400')
fig.savefig(out/'00-固定边际与格内变化.png',dpi=240)
pd.DataFrame(coords,columns=['table','group','DD','ID','II','row_total']).to_csv(out/'00-固定边际与格内变化.csv',index=False)

# 把观察统计量作为箱边界，保证红柱正好覆盖右尾。
edges = np.concatenate([np.arange(0, cutoff, .5), np.arange(cutoff, q.max() + .51, .5)])
counts, edges = np.histogram(q, bins=edges)
colors = [red if x >= cutoff else blue for x in edges[:-1]]
fig, (ax, tail_ax) = plt.subplots(2, 1, figsize=(5.4, 7.2))
ax.bar(edges[:-1], counts, width=np.diff(edges), align='edge', color=colors, edgecolor='white', linewidth=.4)
ax.axvline(q_obs, color=red, lw=1.6)
ax.text(.97, .94, f'观察 χ² = {q_obs:.4f}\n全右尾 {b:,} / {B:,} 次', transform=ax.transAxes, ha='right', va='top', fontsize=12)
ax.set(xlabel='每张模拟表的 Pearson χ²', ylabel='模拟次数', xlim=(0, edges[-1]))
ax.set_title(f'{LANGUAGE}：固定边际，模拟{B:,}次', fontsize=14)
# 第二面板放大靠近观察值的右尾；第一面板保留全部模拟范围。
mask = edges[:-1] >= cutoff
tail_ax.bar(edges[:-1][mask], counts[mask], width=np.diff(edges)[mask], align='edge', color=red, edgecolor='white', linewidth=.4)
tail_ax.axvline(q_obs, color=red, lw=1.6)
tail_ax.set(xlabel='模拟 Pearson χ²（右尾局部）', ylabel='模拟次数', xlim=(7.5,16), title='右尾放大：χ²约7.91至16')
for panel in (ax, tail_ax):
    panel.spines[['top', 'right']].set_visible(False)
    panel.ticklabel_format(axis='y', style='plain')
fig.tight_layout(h_pad=1.5)
fig.savefig(out / '01-模拟卡方分布与右尾.png', dpi=240)
pd.DataFrame({'left':edges[:-1], 'right':edges[1:], 'count':counts, 'in_right_tail':edges[:-1]>=cutoff}).to_csv(out/'01-直方图坐标.csv',index=False)
assert counts.sum() == B and counts[edges[:-1]>=cutoff].sum() == b
# 保留完整逐次结果；图中取稠密的对数间隔，避免10万点挤在左端。
iterations = np.unique(np.r_[np.geomspace(min(100,B),B,700).astype(int), [v for v in [2000,10000,B] if v<=B]])
p_run = (np.cumsum(extreme)[iterations-1] + 1) / (iterations+1)
fig, ax = plt.subplots(figsize=(5.4, 4.7))
ax.plot(iterations, p_run, color=blue, lw=1.4, label='累计模拟P')
ax.axhline(p_exact,color=red,ls='--',lw=1.3,label='穷举条件尾概率')
ax.scatter([B], [p_run[-1]], color=blue, s=28, zorder=3)
ax.annotate(f'{B:,}次：{p_run[-1]:.5f}',(B,p_run[-1]),xytext=(.97,.88),textcoords='axes fraction',ha='right',arrowprops={'arrowstyle':'-','color':blue},fontsize=12)
ax.set_xscale('log')
ax.set(xlabel='已完成的模拟次数（对数刻度）',ylabel='累计模拟P值',xlim=(min(100,B),B*1.1))
ax.set_xticks([100,1000,10000,100000], ['100','1000','1万','10万'])
ax.legend(loc='lower right',frameon=False,fontsize=11)
ax.grid(axis='y',alpha=.15)
ax.spines[['top','right']].set_visible(False)
fig.tight_layout()
fig.savefig(out/'02-累计模拟P值.png',dpi=240)
pd.DataFrame({'B':iterations,'p_mc':p_run,'p_exact':p_exact}).to_csv(out/'02-累计P坐标.csv',index=False)
# 自动核验可设置MPLBACKEND=Agg；交互环境保留图供查看。
print('三图及数据/坐标已生成。')
