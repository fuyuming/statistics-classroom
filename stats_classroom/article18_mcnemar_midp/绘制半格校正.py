"""绘制半格边界与整数刻度示意图。依赖 matplotlib、numpy。"""
from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np

ROOT = Path(__file__).resolve().parent
OUT = ROOT / '文章配图'
OUT.mkdir(exist_ok=True)
available = {f.name for f in font_manager.fontManager.ttflist}
for candidate in ['PingFang SC', 'Heiti TC', 'Arial Unicode MS', 'Noto Sans CJK SC']:
    if candidate in available:
        plt.rcParams['font.sans-serif'] = [candidate]
        break
else:
    raise RuntimeError('需要中文字体')
plt.rcParams['axes.unicode_minus'] = False
blue = '#1769aa'
orange = '#ce641e'
pale = '#fff0e4'

# 同一线性变换下的整数格子：B间隔1，D=2B-14间隔2。
fig, axes = plt.subplots(2, 1, figsize=(5.2, 5.9), dpi=240)
for ax, values, edge, observed, title, edge_label, delta_label in [
    (axes[0], np.arange(10, 15), 11.5, 12, '二项计数 B：每格相隔1', '边界11.5', '半格 = 0.5'),
    (axes[1], np.arange(6, 15, 2), 9, 10, '人数差 D：每格相隔2', '边界9', '半格 = 1'),
]:
    step = values[1] - values[0]
    ax.set_xlim(values[0] - .65 * step, values[-1] + .65 * step)
    ax.set_ylim(-.85, 1.4)
    ax.axvspan(edge, values[-1] + .5 * step, color=pale, ymin=.25, ymax=.72)
    ax.hlines(0, values[0] - .5 * step, values[-1] + .5 * step, color='#a3b1bf', lw=1)
    for v in values:
        ax.vlines(v - .5 * step, -.2, .65, color='#c5ccd2', lw=.8)
        ax.scatter(v, 0, s=65, color=orange if v >= observed else blue, zorder=5)
        ax.text(v, -.36, str(v), ha='center', fontsize=13)
    ax.vlines(values[-1] + .5 * step, -.2, .65, color='#c5ccd2', lw=.8)
    ax.vlines(edge, -.25, .84, color=orange, ls='--', lw=1.5)
    ax.annotate('', xy=(edge, .9), xytext=(observed, .9),
                arrowprops={'arrowstyle': '<->', 'color': orange, 'lw': 1.3})
    ax.text((edge + observed) / 2, 1.1, delta_label, ha='center', fontsize=13, color=orange)
    ax.text(edge, -.72, edge_label, ha='center', fontsize=12, color=orange)
    ax.text(observed + 1.1 * step, .42, '右尾方向 →', ha='center', fontsize=12)
    ax.set_title(title, loc='left', fontsize=15, pad=12)
    ax.axis('off')
fig.subplots_adjust(left=.07, right=.98, top=.82, bottom=.15, hspace=.65)
fig.text(.07, .975, '同一边界，换一个刻度', fontsize=17, va='top', weight='bold')
fig.text(.07, .025, 'D = B − (14 − B) = 2B − 14\n本例：12 − 2 = 10', fontsize=12, linespacing=1.6)
fig.savefig(OUT / '04-半格到减一.png', facecolor='white')
plt.close(fig)

params = {
    'm': 14,
    'observed_B': 12,
    'B_tail_boundary': 11.5,
    'D_definition': 'B minus the other discordant count (14-B)',
    'D_transform': 'B-(14-B)=2*B-14',
    'observed_D': 10,
    'D_tail_boundary': 9,
    'diagram_type': 'lattice boundaries, not probability density',
}
(OUT / '04-刻度参数.json').write_text(json.dumps(params, ensure_ascii=False, indent=2))
print('半格校正示意图与刻度参数已生成。')
