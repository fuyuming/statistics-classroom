"""由第25篇抽样分布延伸至估计：固定1080px画布、32pt数学字号。"""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(__file__).resolve().parents[1]
output=root/'文章配图'
plt.rcParams['mathtext.fontset']='dejavusans'
def export(name,lines):
    fig=plt.figure(figsize=(9,2.9),dpi=120,facecolor='white')
    text_items=[]
    for y,line in zip([.72,.24],lines):
        text_items.append(fig.text(.5,y,line,ha='center',va='center',fontsize=32,color='#172A44'))
    fig.canvas.draw()
    widths=[t.get_window_extent().width for t in text_items]
    assert max(widths)<1040,(name,widths)
    fig.savefig(output/name,facecolor='white')
    plt.close(fig)
    print(name,1080,'px','32pt',widths)
export('公式_率差抽样分布.png',[
 r'$E(D)=\delta$',
 r'$\mathrm{Var}(D)=\frac{p_A(1-p_A)}{n_A}+\frac{p_B(1-p_B)}{n_B}$'
])
export('公式_从概率到区间.png',[
 r'$P\!\left(-1.96\leq\frac{D-\delta}{SE(D)}\leq1.96\right)\approx0.95$',
 r'$D-1.96\,SE(D)\leq\delta\leq D+1.96\,SE(D)$'
])
