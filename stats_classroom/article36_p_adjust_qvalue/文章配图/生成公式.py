from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--out', type=Path, default=Path(__file__).resolve().parent)
out = parser.parse_args().out
out.mkdir(parents=True, exist_ok=True)
panels = [
    ('公式01_FWER.png', [r'$\mathrm{FWER}=\Pr(V\geq 1)$'], (5.2, .8), (.5,)),
    ('公式02_FDR.png', [r'$\mathrm{FDR}=E\left[\frac{V}{\max(R,1)}\right]$'], (5.2, .8), (.5,)),
    ('公式03_BH规则.png', [
        r'$k=\max\{i:p_{(i)}\leq \frac{i}{m}\alpha\}$',
        r'$P^{\mathrm{BH}}_{(i)}=\min\left\{1,\min_{j\geq i}\frac{m\,p_{(j)}}{j}\right\}$'], (5.2, 1.55), (.74, .27)),
    ('公式03_pi0与q.png', [
        r'$\widehat{\pi}_0=\min\left\{1,\frac{\#\{p_i\geq\lambda\}}{m(1-\lambda)}\right\}$',
        r'$\widehat q_{(i)}=\widehat{\pi}_0\,P^{\mathrm{BH}}_{(i)}$'], (5.2, 1.55), (.74, .27)),
]
for name, lines, size, ys in panels:
    fig = plt.figure(figsize=size, dpi=240, facecolor='white')
    for i, line in enumerate(lines):
        fig.text(.5, ys[i], line, ha='center', va='center', fontsize=19, color='#27343a')
    fig.savefig(out/name, bbox_inches='tight', pad_inches=.1)
    plt.close(fig)
print('Generated 4 formula panels (2 current + 2 archived/optional).')
