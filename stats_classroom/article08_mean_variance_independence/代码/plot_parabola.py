"""Y=X²的理论关系图。X~N(0,1)，E[X]=0、E[Y]=1。
曲线坐标是函数网格，不是抽样数据；协方差的期望覆盖整个实数轴。
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import font_manager

def make_parabola_figure(root):
    x = np.linspace(-3,3,601)
    y = x**2
    centered_product = x*(y-1)
    density = np.exp(-x**2/2)/np.sqrt(2*np.pi)
    curve = pd.DataFrame({'x':x,'y':y,'centered_product':centered_product,
                          'normal_density':density,'weighted_integrand':centered_product*density})
    assert np.allclose(centered_product,-centered_product[::-1],atol=1e-12)
    out = root / '运行结果' / 'Python'
    out.mkdir(parents=True,exist_ok=True)
    curve.to_csv(out / 'parabola_coordinates.csv',index=False,float_format='%.15g')
    print('抛物线：X=-2、2时Y都为4，中心化乘积分别为-6、6；理论协方差为0。')
    fonts={f.name for f in font_manager.fontManager.ttflist}
    font=next((f for f in ['PingFang SC','Arial Unicode MS','Noto Sans CJK SC','Microsoft YaHei','SimHei'] if f in fonts),None)
    zh=font is not None
    plt.rcParams.update({'font.family':font or 'DejaVu Sans','font.size':15,'axes.unicode_minus':False})
    teal,orange,navy='#167d80','#cc773e','#18394b'
    fig,ax=plt.subplots(figsize=(8.4,8.7))
    fig.subplots_adjust(left=.12,right=.95,top=.81,bottom=.26)
    fig.suptitle('不相关，也可以有确定的曲线关系' if zh else
                 'Uncorrelated, yet connected by a curve',fontsize=21,color=navy,y=.975)
    fig.text(.5,.917,'X ∼ N(0, 1)       Y = X²',ha='center',fontsize=20,color=navy)
    fig.text(.5,.87,'X 在负方向或正方向远离 0，Y 都会增大' if zh else
             'Large negative or positive X both give large Y',ha='center',fontsize=15)
    ax.plot(x,y,color=teal,lw=3)
    ax.axvline(0,color='#9cabad',lw=1)
    ax.axhline(1,color='#839196',ls='--',lw=1.4,label='E[Y] = 1')
    ax.scatter([-2,2],[4,4],s=100,color=[orange,teal],zorder=5)
    ax.plot([-2,2],[4,4],ls=':',color='#9cabad',lw=1)
    ax.annotate('X = −2，Y = 4\nX(Y − 1) = −6',xy=(-2,4),xytext=(-1.55,7.2),
                ha='center',fontsize=15,color=orange,
                arrowprops={'arrowstyle':'->','color':orange,'lw':1.5})
    ax.annotate('X = 2，Y = 4\nX(Y − 1) = 6',xy=(2,4),xytext=(1.55,7.2),
                ha='center',fontsize=15,color=teal,
                arrowprops={'arrowstyle':'->','color':teal,'lw':1.5})
    ax.set(xlim=(-3.2,3.2),ylim=(-.3,10),xlabel='X',ylabel='Y = X²')
    ax.set_xticks([-3,-2,-1,0,1,2,3])
    ax.set_yticks([0,1,4,9])
    ax.spines[['top','right']].set_visible(False)
    ax.grid(alpha=.13)
    ax.legend(loc='upper center',frameon=False,fontsize=13)
    fig.text(.5,.182,'相反的中心化乘积 × 相同的对称概率密度' if zh else
             'Opposite centered products × equal symmetric densities',ha='center',fontsize=15,color=navy)
    fig.text(.5,.125,'Cov(X,Y) = E[X(X² − 1)] = 0',ha='center',fontsize=21,color=teal)
    fig.text(.5,.067,'每一对 x 与 −x 都如此；但知道 X，仍能完全确定 Y。' if zh else
             'Every pair x and -x cancels; knowing X still determines Y.',ha='center',fontsize=13,color=navy)
    fig.text(.5,.025,'图示范围 −3 ≤ X ≤ 3；协方差的期望覆盖整个实数轴。' if zh else
             'Plot range: -3 to 3. The expectation covers the full real line.',ha='center',fontsize=12,color='#657a7d')
    fig.savefig(root / '文章配图' / '04_抛物线与零协方差.png',dpi=180)

if __name__=='__main__':
    make_parabola_figure(Path(__file__).resolve().parents[1])
