from pathlib import Path
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('--out', type=Path, default=Path(__file__).resolve().parent)
out = parser.parse_args().out
out.mkdir(parents=True, exist_ok=True)
panels = [
    ('公式01_FWER.png', r'$\mathrm{FWER}=\Pr(V\geq 1)$'),
    ('公式02_FDR.png', r'$\mathrm{FDR}=E\left[\frac{V}{\max(R,1)}\right]$'),
]
for name, formula in panels:
    fig = plt.figure(figsize=(5.2, .8), dpi=240, facecolor='white')
    fig.text(.5, .5, formula, ha='center', va='center', fontsize=19, color='#27343a')
    fig.savefig(out/name, bbox_inches='tight', pad_inches=.1)
    plt.close(fig)
print('Generated the 2 formula panels used by this article.')
