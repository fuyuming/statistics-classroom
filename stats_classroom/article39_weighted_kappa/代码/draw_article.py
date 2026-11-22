"""重建正文公式与分歧位置图。项目根/代码目录均可运行；不依赖分析结果缓存。"""
from pathlib import Path
import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties, findfont
from matplotlib.patches import Rectangle

root = Path(__file__).resolve().parents[1]
out = root / '文章配图'
out.mkdir(exist_ok=True)
# STIX用于可重建的数学公式，中文图采用系统字体回退。
plt.rcParams.update({'font.family': 'STIXGeneral', 'mathtext.fontset': 'stix'})
formulae = [
('公式01_未加权.png', [r'$\widehat{\kappa}=\dfrac{P_o-P_e}{1-P_e}$', r'$=\dfrac{0.68-0.25}{1-0.25}=0.573333\ldots$']),
('公式02_权重定义.png', [r'$\ell^{L}_{ij}=\dfrac{|i-j|}{3},\quad w^{L}_{ij}=1-\ell^{L}_{ij}$', r'$\ell^{Q}_{ij}=\left(\dfrac{|i-j|}{3}\right)^2,\quad w^{Q}_{ij}=1-\ell^{Q}_{ij}$']),
('公式03_加权构造.png', [r'$P_{o,w}=\sum_i\sum_j w_{ij}p_{ij}$', r'$P_{e,w}=\sum_i\sum_j w_{ij}r_i c_j$',r'$\widehat{\kappa}_w=\dfrac{P_{o,w}-P_{e,w}}{1-P_{e,w}}$']),
('公式04_观察加权.png', [r'$P_{o,w}=\dfrac{68\times1+24\times\frac{2}{3}}{100}$',r'$\qquad+\dfrac{8\times\frac{1}{3}+0\times0}{100}$',r'$=\dfrac{13}{15}\approx0.866667$']),
('公式05_机会加权.png', [r'$P_{e,w}=\dfrac{4\times1+6\times\frac{2}{3}}{16}$',r'$\qquad+\dfrac{4\times\frac{1}{3}+2\times0}{16}$',r'$=\dfrac{7}{12}\approx0.583333$']),
('公式06_线性Kappa.png', [r'$\widehat{\kappa}_{L}=\dfrac{13/15-7/12}{1-7/12}=0.680$']),
('公式07_损失形式.png', [r'$\widehat{\kappa}_w=1-\dfrac{L_o}{L_e}$',r'$\widehat{\kappa}_{L}=1-\dfrac{2/15}{5/12}=0.680$'])]
# 以统一源字号绘制后，按实际文字包围框紧凑裁切；保留24源像素白边。
for filename, lines in formulae:
    fig = plt.figure(figsize=(8.2, 1.13 * len(lines)), dpi=180, facecolor='white')
    texts = []
    for i, line in enumerate(lines):
        texts.append(fig.text(.5, 1-(i+.5)/len(lines), line, ha='center', va='center', fontsize=32, color='#142d4e'))
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    text_bbox = plt.matplotlib.transforms.Bbox.union([text.get_window_extent(renderer) for text in texts])
    bbox_inches = text_bbox.transformed(fig.dpi_scale_trans.inverted())
    pad_inches = 24 / fig.dpi
    bbox_inches = plt.matplotlib.transforms.Bbox.from_extents(
        bbox_inches.x0 - pad_inches, bbox_inches.y0 - pad_inches,
        bbox_inches.x1 + pad_inches, bbox_inches.y1 + pad_inches)
    fig.savefig(out / filename, dpi=180, facecolor='white', bbox_inches=bbox_inches, pad_inches=0)
    # 公式是学生重建材料的一部分；正文只引用同一份生成结果。
    plt.close(fig)

chinese = FontProperties(family=['PingFang SC', 'Heiti TC', 'Noto Sans CJK SC', 'Arial Unicode MS'])
# 以首个实际可用中文字体绘图，防止逐字回退出现缺字。
font_path = findfont(chinese)
cn = FontProperties(fname=font_path)
with (root / '数据' / 'paired_tables.csv').open(encoding='utf-8-sig') as f:
    data = list(csv.DictReader(f))
tables = {'main': np.zeros((4, 4), dtype=int), 'control': np.zeros((4, 4), dtype=int)}
for row in data:
    tables[row['case']][int(row['A']), int(row['B'])] = int(row['count'])
# 数据图的标题从当前输入计算；改变练习数据后，不沿用固定的“8人”。
far_main = int(tables['main'][0, 3] + tables['main'][3, 0])
far_control = int(tables['control'][0, 3] + tables['control'][3, 0])
far_difference = far_control - far_main
title_main = '主例：分歧集中在近处' if far_main == 0 else f'主例：有{far_main}名最远分歧'
if far_difference > 0:
    title_control = f'对照：多出{far_difference}名最远分歧'
elif far_difference < 0:
    title_control = f'对照：减少{-far_difference}名最远分歧'
else:
    title_control = '对照：最远分歧人数相同'
colour_max = max(1, int(max(tables['main'].max(), tables['control'].max())))
fig, axes = plt.subplots(2, 1, figsize=(6.8, 11), dpi=200)
for ax, case, title in zip(axes, ['main', 'control'], [title_main, title_control]):
    tab = tables[case]
    ax.imshow(tab, cmap='Blues', vmin=0, vmax=colour_max)
    for i in range(4):
        for j in range(4):
            ax.text(j, i, str(tab[i,j]), ha='center', va='center', fontsize=20, color='white' if tab[i,j]>=colour_max*2/3 else '#142d4e')
            if i == j:
                ax.add_patch(Rectangle((j-.5, i-.5),1,1,fill=False,edgecolor='#0000a8',linewidth=2.6))
            if abs(i-j) == 3:
                ax.add_patch(Rectangle((j-.5, i-.5),1,1,fill=False,edgecolor='#ec7c26',linewidth=2.6))
    ax.set_xticks(range(4), ['无','轻','中','重'], fontproperties=cn, fontsize=17)
    ax.set_yticks(range(4), ['无','轻','中','重'], fontproperties=cn, fontsize=17)
    ax.set_xlabel('评定者B（列）', fontproperties=cn, fontsize=17)
    ax.set_ylabel('评定者A（行）', fontproperties=cn, fontsize=17)
    ax.set_title(title, fontproperties=cn, fontsize=21, pad=15, color='#0000a8')
    ax.set_xticks(np.arange(-.5, 4, 1), minor=True)
    ax.set_yticks(np.arange(-.5, 4, 1), minor=True)
    ax.grid(which='minor', color='white', linewidth=1.5)
    ax.tick_params(which='minor', bottom=False, left=False)
fig.subplots_adjust(left=.18, right=.90, bottom=.065, top=.95, hspace=.42)
fig.savefig(out / '图01_分歧位置.png', facecolor='white')
plt.close(fig)
print('已重建7张公式与1张数据图。')
