"""重建本篇7张公式图（式1、1a及式2至6）。需latex2mathml、playwright与本机Chrome。
采用同一数学字体和字号，按内容定宽；输出为可离线显示的PNG。
"""
from pathlib import Path
import os
import json
import xml.etree.ElementTree as ET
from latex2mathml.converter import convert
from playwright.sync_api import sync_playwright

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / '文章配图'
LAYOUT = json.loads(Path(__file__).with_name('公式样式.json').read_text())
# 原数学定义及表1的具体代入。array用等号对齐；相邻求和/连乘符号加数学间距。
FORMULAS = {
    1: r'''\begin{array}{rl}
\bar{x}&=\dfrac{\sum_i(n_{i+}\times x_i)}{N}\\
\bar{y}&=\dfrac{\sum_j(n_{+j}\times y_j)}{N}\\
L_{xx}&=\sum_i\left[n_{i+}\times(x_i-\bar{x})^2\right]\\
L_{yy}&=\sum_j\left[n_{+j}\times(y_j-\bar{y})^2\right]\\
T&=\sum_i\,\sum_j\left[n_{ij}\times(x_i-\bar{x})\times(y_j-\bar{y})\right]\\
r&=\dfrac{T}{\sqrt{L_{xx}\times L_{yy}}}
\end{array}''',
    '1a': r'''\begin{array}{rl}
\bar{x}&=\dfrac{98\times1+63\times2+59\times3+58\times4}{278}\\
&=\dfrac{633}{278}\approx2.277\\
\bar{y}&=\dfrac{122\times1+89\times2+41\times3+26\times4}{278}\\
&=\dfrac{527}{278}\approx1.896\\
L_{xx}&=98\times(1-\bar{x})^2+63\times(2-\bar{x})^2\\
&\quad+59\times(3-\bar{x})^2+58\times(4-\bar{x})^2\\
&\approx367.673\\
L_{yy}&=122\times(1-\bar{y})^2+89\times(2-\bar{y})^2\\
&\quad+41\times(3-\bar{y})^2+26\times(4-\bar{y})^2\\
&\approx263.975\\
T&=70\times(1-\bar{x})\times(1-\bar{y})\\
&\quad+\cdots+14\times(4-\bar{x})\times(4-\bar{y})\\
&\approx149.032\\
r&\approx\dfrac{149.032}{\sqrt{367.673\times263.975}}\\
&\approx0.4784
\end{array}''',
    2: r'''\begin{array}{rl}
\operatorname{E}_0(T\mid\mathcal{M})&=0\\
\operatorname{Var}_0(T\mid\mathcal{M})&=\dfrac{L_{xx}\times L_{yy}}{N-1}\\
Z&=\dfrac{T}{\sqrt{\dfrac{L_{xx}\times L_{yy}}{N-1}}}\\
Q&=Z^2=(N-1)\times r^2
\end{array}''',
    3: r'''P_0(\mathbf{n}\mid\mathcal{M})=
\dfrac{\prod_i n_{i+}!\;\times\;\prod_j n_{+j}!}{N!\;\times\;\prod_i\,\prod_j n_{ij}!}''',
    4: r'''P_{\mathrm{cond}}=
\sum_{\mathbf{n}:\,|T(\mathbf{n})|\ge|T_{\mathrm{obs}}|}
P_0(\mathbf{n}\mid\mathcal{M})''',
    5: r'''\begin{array}{rl}
\hat{p}&=\dfrac{k}{B}\\
p_{\mathrm{MC}}&=\dfrac{k+1}{B+1}\\
\mathrm{MCSE}&=\sqrt{\dfrac{\hat{p}\times(1-\hat{p})}{B}}
\end{array}''',
    6: r'''\begin{array}{rl}
P(K=0\mid p)&=(1-p)^B\\
p_{\mathrm{U}}&=1-0.05^{1/B}
\end{array}''',
}
STYLE = f'''body{{margin:0;background:#fff}}
.formula{{display:inline-block;box-sizing:border-box;width:max-content;
max-width:{LAYOUT['canvas_max_css']}px;padding:{LAYOUT['padding_css']}px;
font-size:{LAYOUT['font_size_css']}px;color:#20252b}}
.formula math+math{{margin-top:24px}}
math{{font-family:{LAYOUT['font_family']};math-style:normal}}
'''

def mathml(tex):
    root = ET.fromstring(convert(tex, display='block'))
    for node in root.iter():
        if node.tag.endswith('mtable'):
            node.set('displaystyle', 'true')
            node.set('rowspacing', '0.55em')
            node.set('columnspacing', '0.15em')
        if node.tag.endswith('mtd'):
            node.set('style', 'text-align:' + node.get('columnalign', 'center'))
        if node.tag.endswith('mover') and len(node) == 2 and node[1].text == '^':
            node.set('accent', 'true')
            node[1].set('stretchy', 'true')
            node[1].text = '\u02c6'
    for row in root.iter():
        if row.tag.endswith('mtr'):
            for i, cell in enumerate(row):
                cell.set('style', 'text-align:' + ('right' if i == 0 else 'left'))
    return ET.tostring(root, encoding='unicode').replace('ns0:', '').replace(':ns0', '')

if __name__ == '__main__':
    OUT.mkdir(exist_ok=True)
    with sync_playwright() as runtime:
        executable = os.environ.get('CHROMIUM_EXECUTABLE')
        browser = runtime.chromium.launch(executable_path=executable) if executable else runtime.chromium.launch(channel='chrome')
        page = browser.new_page(viewport={'width':720,'height':1200},
                                device_scale_factor=LAYOUT['device_scale_factor'])
        for number, tex in FORMULAS.items():
            blocks = tex if isinstance(tex, list) else [tex]
            document = '<meta charset="utf-8"><style>' + STYLE + '</style><div class="formula">' + ''.join(mathml(block) for block in blocks) + '</div>'
            page.set_content(document)
            page.evaluate('document.fonts.ready')
            bounds = page.locator('.formula').evaluate('''el => {
                const box = el.getBoundingClientRect();
                const math = [...el.querySelectorAll('math')].map(m => m.getBoundingClientRect());
                return {width:box.width, scrollWidth:el.scrollWidth, mathRight:Math.max(...math.map(m=>m.right)), right:box.right};
            }''')
            if (bounds['scrollWidth'] > bounds['width'] + 1
                    or bounds['mathRight'] > bounds['right'] + 1):
                raise RuntimeError(f'公式{number}超出字号对应画布，请分行：{bounds}')
            label = f'{number:02d}' if isinstance(number, int) else str(number)
            page.locator('.formula').screenshot(path=str(OUT / f'公式{label}.png'))
        browser.close()
