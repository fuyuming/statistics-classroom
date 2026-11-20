from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

out = Path(__file__).resolve().parent
out.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({'mathtext.fontset': 'stix', 'font.family': 'STIXGeneral'})
items = [
    ('公式01_独立.png', [
        r'$P(X\leq x,\,Y\leq y)=P(X\leq x)\,P(Y\leq y)$',
        r'$F_{X,Y}(x,y)=F_X(x)\,F_Y(y)\quad (\mathrm{all}\ x,y)$'
    ]),
    ('公式02_配对t.png', [
        r'$d_i=B_i-A_i,\qquad SE(\bar d)=\dfrac{s_d}{\sqrt{n}}$',
        r'$t=\dfrac{\bar d-0}{s_d/\sqrt{n}},\quad df=n-1$',
        r'$95\%\ \mathrm{CI}:\quad \bar d\ \pm\ t_{0.975,n-1}\dfrac{s_d}{\sqrt{n}}$'
    ]),
    ('公式03_差值方差.png', [
        r'$\mathrm{Var}(B-A)=\mathrm{Var}(B)+\mathrm{Var}(A)$',
        r'$-\ 2\,\mathrm{Cov}(A,B)$'
    ])
]
for name, lines in items:
    height = 1.20 * len(lines)
    fig = plt.figure(figsize=(9.6, height), dpi=160, facecolor='white')
    for i, line in enumerate(lines):
        fig.text(.5, 1-(i+.5)/len(lines), line, ha='center', va='center', fontsize=26, color='#142d4e')
    fig.savefig(out / name, dpi=180, facecolor='white', bbox_inches='tight', pad_inches=.24)
    plt.close(fig)
print(out)
