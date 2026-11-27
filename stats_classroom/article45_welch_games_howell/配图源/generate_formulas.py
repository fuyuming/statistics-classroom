"""将LaTeX渲染为高清公式PNG。数据图由generate_figures.R统一绘制。"""
from pathlib import Path
import json
import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / '文章配图'
OUT.mkdir(exist_ok=True)
# 固定画布宽3.66in，300dpi；手机显示366px时主公式约18px。
FORMULAS = {
 '01': [r'v_i=\frac{s_i^2}{n_i},\qquad \mathrm{SE}(\bar y_i)=\sqrt{v_i}', r'w_i=\frac{1}{v_i}=\frac{n_i}{s_i^2}'],
 '02': [r'W=\sum_{i=1}^{k}w_i,\qquad \bar y_w=\frac{\sum_{i=1}^{k}w_i\bar y_i}{W}'],
 '03': [r'A=\frac{\sum_{i=1}^{k}w_i(\bar y_i-\bar y_w)^2}{k-1}'],
 '04': [r'C=\sum_{i=1}^{k}\frac{(1-w_i/W)^2}{n_i-1}',r'F_W=\frac{A}{1+2(k-2)C/(k^2-1)}',r'\nu=\frac{k^2-1}{3C}'],
 '05': [r'F_W=\frac{12.377002}{1+0.171352/4}=11.868577',r'\nu=\frac{8}{3\times0.171352}=15.562533',r'H_0:\quad F_W\ \overset{\mathrm{approx}}{\sim}\ F(2,15.562533)'],
 '06': [r'\mathrm{SE}_{ij}=\sqrt{v_i+v_j}',r'\nu_{ij}=\frac{(v_i+v_j)^2}{v_i^2/(n_i-1)+v_j^2/(n_j-1)}'],
 '07': [r't_{ij}=\frac{d_{ij}}{\mathrm{SE}_{ij}}',r'q_{ij}=\sqrt{2}|t_{ij}|=\frac{\sqrt{2}|d_{ij}|}{\mathrm{SE}_{ij}}',r'=\frac{|d_{ij}|}{\sqrt{(v_i+v_j)/2}}'],
 '07a': [r'Z_1,\ldots,Z_k\ \sim N(0,1)',r'U\sim\chi_\nu^2,\quad U\perp(Z_1,\ldots,Z_k)',r'Q_{k,\nu}=\frac{\max_i Z_i-\min_i Z_i}{\sqrt{U/\nu}}'],
 '08': [r'P_{\mathrm{adj},ij}=\Pr(Q_{k,\nu_{ij}}\geq q_{ij})',r'\mathrm{CI}_{ij}:\quad d_{ij}\ \pm\ q_{0.95;k,\nu_{ij}}\frac{\mathrm{SE}_{ij}}{\sqrt{2}}'],
 '09': [r'8.156667\ \pm\ 3.805083\times\frac{1.884437}{\sqrt{2}}',r'=8.156667\ \pm\ 5.070266',r'=[3.086400,\ 13.226933]\ \mathrm{mmHg}']
}
plt.rcParams.update({'mathtext.fontset':'stix','font.family':'DejaVu Sans','font.size':13})
for no, lines in FORMULAS.items():
    height = .22 + .60*len(lines)
    fig = plt.figure(figsize=(3.66,height), dpi=300, facecolor='white')
    for i, line in enumerate(lines):
        fig.text(.5, 1-(.11+.60*(i+.5))/height, '$'+line+'$', ha='center', va='center',fontsize=13,color='#172d45')
    fig.savefig(OUT/f'公式{no}.png',dpi=300,facecolor='white')
    plt.close(fig)
(ROOT/'配图源'/'formulas.json').write_text(json.dumps(FORMULAS,ensure_ascii=False,indent=2),encoding='utf-8')
print('10张高清公式已生成；正文两张数据图请运行generate_figures.R。')
