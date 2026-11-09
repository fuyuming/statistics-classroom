"""固定候选总体率差，比较未知相同概率的粗网格与细网格搜索。"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.special import gammaln, xlogy

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / '数据'
RESULT = ROOT / '运行结果' / 'Python'
FIG = ROOT / '运行结果' / 'Python' / '图像'
FIG.mkdir(parents=True, exist_ok=True)
RESULT.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({'font.sans-serif': ['STHeiti','PingFang SC','Arial Unicode MS','DejaVu Sans'],
                     'axes.unicode_minus': False, 'font.size': 14})

# 输入为主例在Wang排序下、计算下限所用的真实尾部计数组合。
# 它包含观察组合自身及更前面的组合，不改用简单的D排序。
counts = pd.read_csv(DATA / 'wang_lower_tail_n20.csv')[['n_plus','n_zero','n_minus']].to_numpy(dtype=int)
n = int(counts[0].sum())
assert np.all(counts.sum(axis=1) == n)
theta = -0.45
upper = 1 - abs(theta)
coeff = gammaln(n+1) - gammaln(counts+1).sum(axis=1)


def tail_probability(p_t):
    # theta固定后，每个p_T确定一组三类概率。
    values = np.atleast_1d(p_t)
    probs = np.column_stack(((1-values+theta)/2, values, (1-values-theta)/2))
    assert probs.min() >= -1e-12
    probs = np.maximum(probs, 0)
    # 逐个三类人数组合求概率，再对整个尾部相加；xlogy处理0人数与0概率。
    terms = coeff[:,None] + xlogy(counts[:,None,:], probs[None,:,:]).sum(axis=2)
    return np.exp(terms).sum(axis=0)


# 教学图先试11点，再在粗峰的左右相邻点之间试21点。
# 此处点数为便于读图而选，正文PairedCI实际调用使用60/40。
coarse = np.linspace(0, upper, 11)
y_coarse = tail_probability(coarse)
best = int(np.argmax(y_coarse))
left = coarse[max(0,best-1)]
right = coarse[min(len(coarse)-1,best+1)]
fine = np.linspace(left, right, 21)
y_fine = tail_probability(fine)
fine_best = int(np.argmax(y_fine))
reference = np.linspace(0, upper, 4001)
y_reference = tail_probability(reference)

for name, values, probs in [('coarse',coarse,y_coarse),('fine',fine,y_fine),('reference',reference,y_reference)]:
    pd.DataFrame({'theta':theta,'p_t':values,'tail_probability':probs}).to_csv(RESULT / f'wang_grid_{name}.csv',index=False)
summary = {'theta':theta, 'tail_count':len(counts), 'p_t_range':[0,upper],
           'coarse_points':len(coarse), 'fine_points':len(fine),
           'coarse_step':float(coarse[1]-coarse[0]),'fine_step':float(fine[1]-fine[0]),
           'coarse_peak_p_t':float(coarse[best]),'coarse_peak_probability':float(y_coarse[best]),
           'fine_range':[float(left),float(right)],'fine_peak_p_t':float(fine[fine_best]),
           'fine_peak_probability':float(y_fine[fine_best])}
(RESULT / 'wang_grid_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')

fig, axes = plt.subplots(2,1,figsize=(7.1,8.8),dpi=210)
ax=axes[0]
ax.plot(reference,y_reference,color='#6f7c88',lw=1.7,label='密集计算的参考曲线')
ax.scatter(coarse,y_coarse,color='#1266ab',s=43,zorder=3,label='粗网格：11个试算点')
ax.scatter([coarse[best]],[y_coarse[best]],facecolors='none',edgecolors='#df6b25',s=160,lw=2,zorder=4)
ax.axvspan(left,right,color='#df6b25',alpha=.12,label='接下来细搜的区域')
ax.set(title=f'上：固定候选 θ = {theta:.2f}，先在整个范围粗搜',xlim=(0,upper),ylim=(0,max(y_reference)*1.28))
ax.annotate(f'粗网格最高点\n$p_T$={coarse[best]:.3f}',xy=(coarse[best],y_coarse[best]),
            xytext=(.60,.87),textcoords='axes fraction',fontsize=12,
            arrowprops={'arrowstyle':'->','color':'#df6b25'},ha='center')
ax.legend(loc='lower left',fontsize=11,frameon=False)
ax=axes[1]
local_x=np.linspace(left,right,1601)
ax.plot(local_x,tail_probability(local_x),color='#6f7c88',lw=1.7)
ax.scatter(fine,y_fine,color='#df6b25',s=33,zorder=3,label='细网格：21个试算点')
ax.scatter(coarse[max(0,best-1):min(len(coarse),best+2)],y_coarse[max(0,best-1):min(len(coarse),best+2)],
           color='#1266ab',s=55,zorder=4,label='这个区域内的粗网格点')
ax.scatter([fine[fine_best]],[y_fine[fine_best]],facecolors='none',edgecolors='#a42d2d',s=170,lw=2,zorder=5)
spread=max(float(y_fine.max()-y_fine.min()),1e-6)
ax.set(title=f'下：放大 [{left:.3f}, {right:.3f}]，再细搜',xlim=(left,right),
       ylim=(max(0,y_fine.min()-.18*spread),y_fine.max()+.48*spread))
ax.annotate(f'细网格最高点\n$p_T$={fine[fine_best]:.4f}\n尾部概率={y_fine[fine_best]:.6f}',
            xy=(fine[fine_best],y_fine[fine_best]),xytext=(.77,.83),textcoords='axes fraction',fontsize=12,
            arrowprops={'arrowstyle':'->','color':'#a42d2d'},ha='center')
ax.legend(loc='lower left',fontsize=11,frameon=False)
for ax in axes:
    ax.set_xlabel('两次结果相同的总体概率 $p_T$')
    ax.set_ylabel('单侧尾部概率')
    ax.grid(alpha=.17)
    ax.spines[['top','right']].set_visible(False)
    ax.ticklabel_format(axis='y',style='plain',useOffset=False)
fig.tight_layout(h_pad=2)
fig.savefig(FIG / '图3_未知概率的两步网格搜索.png')
fig.savefig(FIG / '图3_未知概率的两步网格搜索.pdf')
plt.close(fig)
print(json.dumps(summary,ensure_ascii=False))
