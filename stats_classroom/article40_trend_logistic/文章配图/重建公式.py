"""按统一手机视觉字号生成本篇5幅公式，并保存明确显示尺寸。"""
from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from PIL import Image

root = Path(__file__).resolve().parent
d = pd.read_csv(root.parent / '配套项目/article40_trend_logistic/数据/共享教学数据.csv')
p0 = d.failed.sum() / d.n.sum()
sbar = (d.n * d.score).sum() / d.n.sum()
t = d.score - sbar
U = (t * (d.failed - d.n*p0)).sum()
variance = p0 * (1-p0) * (d.n * t**2).sum()
z = U / variance**.5
q = U**2 / variance

# 沿用第35篇的STIX数学字体与手机显示尺度。
# 24pt × 162dpi / 72 / 3 = 18 CSS像素；裁切后仍使用同一3:1显示比例。
# 不把短式拉满整行，也不单独缩小长式。超宽时应先按数学步骤分行。
FONT_PT = 24
DPI = 162
PIXELS_PER_CSS = 3
MAX_DISPLAY_WIDTH = 300
plt.rcParams.update({'font.family':'STIXGeneral', 'mathtext.fontset':'stix'})
items = [
 ('公式01_胜算与logit.png', [r'$\mathrm{odds}_i=\dfrac{p_i}{1-p_i}$',r'$\mathrm{logit}(p_i)=\log\!\left(\dfrac{p_i}{1-p_i}\right)$']),
 ('公式02_Logistic概率模型.png',[r'$\log\!\left(\dfrac{p_i}{1-p_i}\right)=\alpha+\beta s_i$',r'$p_i=\dfrac{\exp(\alpha+\beta s_i)}{1+\exp(\alpha+\beta s_i)}$']),
 ('公式03_OR与区间.png',[r'$\mathrm{OR}=\exp(\beta)$',r'$\mathrm{CI}_{95\%}=[\exp(\widehat\beta-z_{0.975}\,\mathrm{SE}),$'+'\n'+r'$\qquad\qquad\exp(\widehat\beta+z_{0.975}\,\mathrm{SE})]$']),
 ('公式04_趋势得分与方差.png',[r'$U=\sum_i(s_i-\bar s)(x_i-n_i\bar p)$',r'$V_0=\bar p(1-\bar p)\sum_i n_i(s_i-\bar s)^2$']),
 ('公式05_趋势Q.png',[rf'$Z=\dfrac{{U}}{{\sqrt{{V_0}}}}\approx {z:.4f}$',rf'$Q=Z^2=\dfrac{{U^2}}{{V_0}}\approx {q:.4f}$',r'$Q\ \overset{d}{\longrightarrow}\ \chi^2_1\quad(H_0)$'])]
layout = {'font_family':'STIXGeneral', 'mathtext_fontset':'stix', 'font_pt':FONT_PT,
          'dpi':DPI, 'pixels_per_css':PIXELS_PER_CSS,
          'nominal_math_css_px':FONT_PT*DPI/72/PIXELS_PER_CSS, 'images':{}}
for name, lines in items:
    heights = [1.05 if '\n' in line else .90 for line in lines]
    height = sum(heights) + .25
    fig = plt.figure(figsize=(7.5,height),dpi=DPI,facecolor='white')
    cursor = height-.13
    for line, line_height in zip(lines,heights):
        fig.text(.5,(cursor-line_height/2)/height,line,ha='center',va='center',
                 fontsize=FONT_PT,color='#142d4e',linespacing=1.55)
        cursor -= line_height
    # 只裁掉多余空白，所有公式保留相同物理字号和边距。
    fig.savefig(root/name,dpi=DPI,facecolor='white',bbox_inches='tight',pad_inches=.20)
    plt.close(fig)
    width, height_px = Image.open(root/name).size
    css_width = width/PIXELS_PER_CSS
    assert css_width <= MAX_DISPLAY_WIDTH, f'{name}过宽，应分行：{css_width:.1f}px'
    layout['images'][name]={'pixel_width':width, 'pixel_height':height_px,
                           'display_width_css_px':css_width,
                           'display_height_css_px':height_px/PIXELS_PER_CSS}
(root/'公式显示尺寸.json').write_text(json.dumps(layout,ensure_ascii=False,indent=2))
print('Generated 5 panels at a uniform 18 CSS px nominal math size.')
for name,info in layout['images'].items():print(name,info)
