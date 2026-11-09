# 配套科学图，从Python主入口保存的精确分布坐标绘制（UTF-8）。
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / '运行结果' / 'Python' / '图像'
FIG.mkdir(parents=True, exist_ok=True)
RESULT = ROOT / '运行结果' / 'Python'
plt.rcParams['font.sans-serif'] = ['STHeiti', 'PingFang SC', 'Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['font.size'] = 12


def draw(ax, name, title, limits):
    data = pd.read_csv(RESULT / ('distribution_' + name + '.csv'))
    shown = data.D.between(*limits)
    values = data[shown]
    width = float(np.diff(data.D).min()) * 0.82
    ax.bar(values.D, values.exact_probability, width=width, color='#4C78A8', alpha=0.8, label='精确概率')
    ax.plot(values.D, values.normal_cell_probability, 'o', ms=3.5, color='#E45756', label='正态单格概率')
    ax.set(xlim=limits, xlabel='率差 D=k/n（后测减前测）', ylabel='单个取值的概率', title=title)
    ax.set_ylim(0, max(values.exact_probability.max(), values.normal_cell_probability.max()) * 1.25)
    ax.grid(axis='y', alpha=0.2)
    ax.legend(loc='upper right', fontsize=10)


fig, axes = plt.subplots(2, 1, figsize=(6.6, 8.5), dpi=190)
draw(axes[0], 'normal_n20', '给定变化概率，n=20', (-0.65, 0.25))
draw(axes[1], 'normal_n1000', '同一变化概率，n=1000（横轴局部放大）', (-0.27, -0.13))
fig.tight_layout(h_pad=2.4)
fig.savefig(FIG / '图1_配对率差精确与近似.png')
fig.savefig(FIG / '图1_配对率差精确与近似.pdf')
plt.close(fig)

fig, ax = plt.subplots(figsize=(6.6, 4.8), dpi=190)
draw(ax, 'rare_n1000', 'n=1000，平均只有3人发生变化', (-0.009, 0.007))
fig.tight_layout()
fig.savefig(FIG / '图2_变化稀少.png')
fig.savefig(FIG / '图2_变化稀少.pdf')
plt.close(fig)

# 图3继续使用Python逐项计算141组三项概率，不读取R曲线。
import runpy
runpy.run_path(str(ROOT / "代码" / "plot_grid_search.py"), run_name="__main__")
