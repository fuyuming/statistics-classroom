from pathlib import Path
import numpy as np
from scipy.stats import norm
import matplotlib
import matplotlib.pyplot as plt
r=Path(__file__).resolve().parents[1]
(r/'文章配图').mkdir(exist_ok=True)
plt.rcParams.update({'font.size':12,'axes.spines.top':False,'axes.spines.right':False})
fig,axes=plt.subplots(2,1,figsize=(8,8),sharex=True,sharey=True)
se=1.58/np.sqrt(100)
for ax,a in zip(axes,[.05,.01]):
 c=7.25+norm.ppf(1-a)*se
 beta=norm.cdf((c-7.65)/se)
 x=np.sort(np.append(np.linspace(6.65,8.3,2500),c))
 y0=norm.pdf(x,7.25,se);y1=norm.pdf(x,7.65,se)
 ax.plot(x,y0,color='#285CA0',lw=2,label='True mean = 7.25 cm')
 ax.plot(x,y1,color='#C17A21',lw=2,label='True mean = 7.65 cm')
 ax.fill_between(x,0,y0,where=x>=c,color='#72B9E3',alpha=.65,label=f'False positive: {a:.0%}')
 ax.fill_between(x,0,y1,where=x<=c,color='#D16A5A',alpha=.4,label=f'Miss: {beta:.2%}')
 ax.axvline(c,color='#333333',ls='--')
 ax.text(c+.015,2.8,f'Cutoff\n{c:.4f}',fontsize=11)
 ax.set(title=f'n = 100, alpha = {a:.2f}, power = {1-beta:.2%}',ylabel='Density',ylim=(0,3.6))
 ax.legend(loc='upper left',fontsize=10,frameon=False)
axes[-1].set_xlabel('Sample mean of 100 fish (cm)')
fig.tight_layout()
fig.savefig(r/'文章配图/04_降低alpha与漏判.png',dpi=180)
plt.close(fig)
