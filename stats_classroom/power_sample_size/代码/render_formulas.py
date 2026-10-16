# 文章公式排版；保存原始LaTeX，供重新生成。无统计参数拟合。
from pathlib import Path
import matplotlib.pyplot as plt
root = Path(__file__).resolve().parents[1]
out = root / '文章配图'
out.mkdir(exist_ok=True)
formulas = [('公式_样本均数标准误', '\\mathrm{SE}(\\bar{X})=\\frac{\\sigma}{\\sqrt{n}}=\\frac{1.58}{\\sqrt{100}}=0.158\\ \\mathrm{cm}'), ('公式_统计功效', '1-\\beta=1-\\Phi\\!\\left(z_{1-\\alpha}-d\\sqrt{n}\\right)'), ('公式_标准化效应量', 'd=\\frac{\\delta}{\\sigma}=\\frac{0.40}{1.58}\\approx0.253'), ('公式_样本量', 'n_{\\min}=\\left\\lceil\\left[\\frac{(z_{1-\\alpha}+z_{1-\\beta})\\sigma}{\\delta}\\right]^2\\right\\rceil'), ('公式_样本量代入', '\\left[\\frac{(1.6449+0.8416)\\times1.58}{0.40}\\right]^2\\approx96.46')]
for name, formula in formulas:
    fig = plt.figure(figsize=(8,1.25),facecolor='white')
    fig.text(.5,.5,'$'+formula+'$',ha='center',va='center',fontsize=25,color='#183b56')
    fig.savefig(out/(name+'.png'),dpi=220,bbox_inches='tight',pad_inches=.15)
    plt.close(fig)
