"""Rebuild the numbered mathematical panels; optional arguments select filenames."""
from pathlib import Path
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

out = Path(__file__).resolve().parent
plt.rcParams.update({'font.family': 'STIXGeneral', 'mathtext.fontset': 'stix'})
items = [('公式01_期望数与卡方.png',
  ['$E_{ij}=\\dfrac{N_{i+}N_{+j}}{n}$',
   '$Q=\\sum_{i=1}^{r}\\sum_{j=1}^{s}\\dfrac{(N_{ij}-E_{ij})^2}{E_{ij}}$',
   '$\\nu=(r-1)(s-1)$']),
 ('公式02_CramerV.png',
  ['$m=\\min(r-1,s-1)$',
   '$\\widehat V=\\sqrt{\\dfrac{Q}{nm}},\\qquad 0\\leq\\widehat V\\leq1$',
   '$\\widehat V=\\sqrt{\\dfrac{18.657938}{300\\times1}}\\approx0.249385$'])]
selected = set(sys.argv[1:])
for name, lines in items:
    if selected and name not in selected:
        continue
    fig = plt.figure(figsize=(8.2, 1.15 * len(lines)), dpi=162, facecolor='white')
    for i, line in enumerate(lines):
        fig.text(.5, 1-(i+.5)/len(lines), line, ha='center', va='center', fontsize=24, color='#142d4e')
    fig.savefig(out / name, dpi=162, facecolor='white', bbox_inches='tight', pad_inches=.20)
    plt.close(fig)
print('Rebuilt mathematical panels:', len(selected) if selected else len(items))
