"""从实际三语言结果CSV重建κ区间比较图；先运行项目01和02脚本并保存结果。"""
from pathlib import Path
import csv
import matplotlib.pyplot as plt

out = Path(__file__).resolve().parent
results = out.parent / '配套项目/article34_paired_kappa/results'
with (results / 'R_summary.csv').open() as f:
    ordinary = next(csv.DictReader(f))
with (results / 'R_bootstrap.csv').open() as f:
    bootstrap = next(csv.DictReader(f))
# 原始点估计相同；读取真正计算的端点，不在绘图函数里重新抽样或构造区间。
point = float(bootstrap['kappa'])
assert abs(point - float(ordinary['kappa'])) < 1e-12
rows = [
    ('Delta-Wald', float(ordinary['ci_low']), float(ordinary['ci_high'])),
    ('Paired bootstrap: percentile', float(bootstrap['percentile_low']), float(bootstrap['percentile_high'])),
    ('Paired bootstrap: BCa', float(bootstrap['bca_low']), float(bootstrap['bca_high'])),
]
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 14,
    'axes.spines.top': False, 'axes.spines.right': False, 'axes.spines.left': False})
fig, ax = plt.subplots(figsize=(6.4, 5.3), layout='constrained')
colors = ['#1769aa', '#087f8c', '#df7357']
for y, (label, low, high), color in zip([2, 1, 0], rows, colors):
    ax.hlines(y, low, high, color=color, linewidth=4)
    ax.vlines([low, high], y-.07, y+.07, color=color, linewidth=2)
    ax.scatter(point, y, color=color, s=60, zorder=3)
    ax.text(.16, y+.23, label, color=color, fontsize=13)
    ax.text(low, y-.15, f'{low:.4f}', ha='center', va='top', fontsize=12)
    ax.text(high, y-.15, f'{high:.4f}', ha='center', va='top', fontsize=12)
ax.axvline(point, color='#b0bec5', linewidth=1, linestyle=':', zorder=0)
ax.set(xlim=(.14, .76), ylim=(-.5, 2.65), yticks=[], xlabel='Population kappa')
ax.set_title('Three 95% confidence intervals', fontsize=17, pad=15)
ax.text(.5, 1.01, f'n = 58; sample kappa = {point:.4f}', transform=ax.transAxes,
        ha='center', fontsize=12)
ax.grid(axis='x', alpha=.15)
fig.savefig(out / '01_三种区间比较.png', dpi=200, bbox_inches='tight', facecolor='white')
plt.close(fig)
print('Rebuilt CI comparison from actual R result CSVs.')
