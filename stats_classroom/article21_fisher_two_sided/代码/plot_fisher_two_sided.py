"""读取共享计数，复现Fisher概率图和三张表；保存每张图的数值依据。"""
from pathlib import Path
from math import comb
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / '文章配图'
OUT.mkdir(exist_ok=True)
data = pd.read_csv(ROOT / '数据/hbv_prevention.csv')
counts = data[['positive', 'negative']].to_numpy()
if counts.shape != (2, 2) or not np.isfinite(counts).all():
    raise ValueError('请输入完整的2×2计数表。')
if (counts < 0).any() or not (counts == np.floor(counts)).all():
    raise ValueError('计数必须是非负整数。')
counts = counts.astype(int)
row1, row2 = counts.sum(axis=1)
col1, col2 = counts.sum(axis=0)
N = counts.sum()
a = counts[0, 0]
xs = np.arange(max(0, row1-col2), min(row1, col1)+1)
q = np.array([comb(int(col1), int(x))*comb(int(col2), int(row1-x))/comb(int(N), int(row1)) for x in xs])
q_obs = q[xs == a][0]
included = q <= q_obs*(1+1e-7)
assert np.isclose(q.sum(), 1)
plt.rcParams['font.sans-serif'] = ['Arial Unicode MS', 'Noto Sans CJK SC', 'DejaVu Sans']
plt.rcParams['font.size'] = 14
red = '#c94a3a'
blue = '#2477b5'

# 图1：保留完整分布，另把观察值和右侧最接近门槛的表放大比较。
fig, axes = plt.subplots(2, 1, figsize=(5.4, 7.5), gridspec_kw={'height_ratios':[1.25,1]})
colors = np.where(included, red, blue)
axes[0].bar(xs, q, color=colors, width=.72)
axes[0].axhline(q_obs, color='#48565f', linestyle='--', linewidth=1.4)
axes[0].text(.02, .96, '红色：计入双侧P\n蓝色：不计入', transform=axes[0].transAxes,
             va='top', fontsize=13)
axes[0].text(.98, .56, f'观察表概率\nq = {q_obs:.5f}', transform=axes[0].transAxes,
             ha='right', fontsize=12, bbox={'facecolor':'white','edgecolor':'none','alpha':.9})
axes[0].set(xticks=xs, xlabel='预防组阳性人数 x', ylabel='原假设下的概率', title='给定边际：每张表有多大概率？')
# 原例选4、8、9；改数据后改为观察值与右端两个支持值，避免超出支持集。
selected = list(dict.fromkeys([int(a), int(xs[-2]), int(xs[-1])])) if len(xs)>1 else [int(a)]
selected_q = np.array([q[xs==x][0] for x in selected])
selected_in = selected_q <= q_obs*(1+1e-7)
positions = np.arange(len(selected))
axes[1].bar(positions, selected_q, color=np.where(selected_in, red, blue), width=.55)
axes[1].axhline(q_obs, color='#48565f', linestyle='--', linewidth=1.4)
for i, v in enumerate(selected_q):
    axes[1].text(i, v+max(selected_q)*.035, f'{v:.8f}', ha='center', fontsize=12)
labels = [f'x = {x}\n'+('计入' if inc else '不计入') for x, inc in zip(selected, selected_in)]
axes[1].set(xticks=positions, xticklabels=labels, ylabel='原假设下的概率',
            title='单独比较观察表与右端两表', ylim=(0,max(selected_q)*1.28))
for ax in axes:
    ax.spines[['top','right']].set_visible(False)
fig.tight_layout(h_pad=2)
fig.savefig(OUT/'01-超几何概率与双侧范围.png', dpi=240)
pd.DataFrame({'X':xs,'probability':q,'included':included,'observed_probability':q_obs}).to_csv(OUT/'01-概率图数据.csv', index=False)

# 图2：同一边际下三张表的对应变化；合计统一浅蓝底。
fig, ax = plt.subplots(figsize=(5.4, 8.2))
fig.subplots_adjust(left=.02, right=.98, top=.98, bottom=.02)
ax.axis('off')
ax.text(.5,.98,'同一组行列合计，三张可能的表',ha='center',va='top',fontsize=16,color=blue,weight='bold')
rows = []
for idx, x in enumerate(selected):
    y = .70 - idx*.30
    qi = q[xs==x][0]
    take = qi<=q_obs*(1+1e-7)
    label = '观察表' if x == a else '右侧表'
    verdict = '计入双侧P' if take else '不计入双侧P'
    ax.text(.5,y+.215,f'{label} x＝{x}：{verdict}',ha='center',fontsize=14,color=red if take else blue,weight='bold')
    ax.text(.5,y+.179,f'概率 {qi:.8f}',ha='center',fontsize=12,color='#465762')
    table_data = [['预防组',x,row1-x,row1],['非预防组',col1-x,row2-col1+x,row2],['合计',col1,col2,N]]
    table = ax.table(cellText=table_data,colLabels=['组别','阳性','阴性','合计'],cellLoc='center',
                     bbox=[.03,y-.01,.94,.16],colWidths=[.31,.23,.23,.23])
    table.auto_set_font_size(False)
    table.set_fontsize(14)
    for (i,j), cell in table.get_celld().items():
        cell.set_edgecolor('white')
        cell.set_linewidth(1.5)
        if i==0 or j==0:
            cell.set_facecolor('#f1f4f7')
        if i==3 or j==3:
            cell.set_facecolor('#e2effb')
            cell.get_text().set_color(blue)
        if i in [1,2] and j in [1,2]:
            cell.set_facecolor('#fff0e8' if take else '#edf2f6')
    rows.append({'X':x,'prevention_positive':x,'prevention_negative':row1-x,
                 'control_positive':col1-x,'control_negative':row2-col1+x,'probability':qi,'included':take})
ax.text(.5,.045,'蓝色合计始终为：22、11、9、24、33' if (row1,row2,col1,col2)==(22,11,9,24) else '蓝色行列合计保持观察值不变',ha='center',fontsize=12,color=blue)
fig.savefig(OUT/'02-观察表与右侧两表.png',dpi=240)
pd.DataFrame(rows).to_csv(OUT/'02-三张表数据.csv',index=False)
print('Fisher两幅图及对应数据已生成。')
