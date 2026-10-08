"""精读05：按作者2023年第05讲计算。Python / Spyder 中按 F5。
输入：rethinking包的WaffleDivorce；输出：运行结果/v2/Python，四张原生中文图。
主模型先验与PPT31页一致；二次近似不等于调用R quap。二元例子精确枚举，非抽样。
"""
# %% 1. 路径和数据：一行代表一个地区，不能当作一个人
from pathlib import Path
import itertools
import numpy as np
from scipy.optimize import minimize
from scipy.stats import norm
import matplotlib.pyplot as plt
from matplotlib import font_manager
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / '运行结果/v2/Python'
OUT.mkdir(parents=True, exist_ok=True)
d = np.genfromtxt(ROOT / '数据/v2/WaffleDivorce.csv', delimiter=',', names=True, encoding='utf-8')
raw = np.column_stack([d['Divorce'], d['Marriage'], d['MedianAgeMarriage']])
print('前5行：离婚率、结婚率、结婚年龄中位数\n', raw[:5])
# ddof=1：采用样本标准差，分母是人数减1；与R的scale一致。
standard = (raw - raw.mean(axis=0)) / raw.std(axis=0, ddof=1)
y = standard[:, 0]
X = np.column_stack([np.ones(len(y)), standard[:, 1:]])
prior_sd = np.array([0.2, 0.5, 0.5])

# %% 2. 后验：似然乘先验；在峰顶附近作二次近似
# 前三个参数是截距、结婚率斜率、年龄斜率，最后一个是残差标准差。
def objective(theta):
    beta = theta[:3]
    sigma = theta[3]
    if sigma <= 0:
        return np.inf
    residual = y - X @ beta
    return len(y)*np.log(sigma) + np.sum(residual**2)/(2*sigma**2) + np.sum((beta/prior_sd)**2)/2 + sigma
fit = minimize(objective, [0, 0, -0.5, 0.8], method='BFGS', options={'gtol': 1e-6})
if not fit.success:
    raise RuntimeError(fit.message)
theta = fit.x
beta = theta[:3]
sigma = theta[3]
residual = y - X @ beta
# Hessian：峰顶曲率决定近似后验的宽度，也保留参数间的相关性。
H = np.zeros((4,4))
H[:3,:3] = X.T @ X / sigma**2 + np.diag(1/prior_sd**2)
H[:3,3] = 2*X.T @ residual / sigma**3
H[3,:3] = H[:3,3]
H[3,3] = -len(y)/sigma**2 + 3*np.sum(residual**2)/sigma**4
cov = np.linalg.inv(H)
sd = np.sqrt(np.diag(cov))
# 中间89%的区间，两端各留5.5%。
z89 = norm.ppf(0.945)
summary = np.column_stack([theta, sd, theta-z89*sd, theta+z89*sd])
np.savetxt(OUT/'posterior.csv', summary, delimiter=',', header='mean,sd,lower89,upper89', comments='')
np.savetxt(OUT/'covariance.csv', cov, delimiter=',')
print('行顺序 alpha, beta_M, beta_A, sigma；列 mean,sd,lower89,upper89\n', summary)
print('beta_M：年龄相同、结婚率增加1个标准差时，平均离婚率的变化。')

# %% 3. 枚举作者二元生成规则的所有可能组合
# 0和1只是状态标签。伯努利概率：状态是1取p，是0取1-p。
def bern(value, probability):
    return probability if value == 1 else 1-probability
rows = []
for model in range(1,5):
    for x,z,yv,a in itertools.product([0,1], repeat=4):
        if model == 1: # 叉：Z分别影响X和Y；给定Z后独立生成
            w = .5 * bern(x,.1+.8*z) * bern(yv,.1+.8*z) * .5
        elif model == 2: # 管：X影响Z，Z影响Y
            w = .5 * bern(z,.1+.8*x) * bern(yv,.1+.8*z) * .5
        elif model == 3: # 对撞：X、Y独立，任一为1都会提高Z=1的概率
            w = .25 * bern(z,.9 if x+yv>0 else .2) * .5
        else: # 中介的后代A：90%机会与Z一致
            w = .5 * bern(z,.1+.8*x) * bern(yv,.1+.8*z) * bern(a,.1+.8*z)
        rows.append([model,x,z,yv,a,w])
rows = np.array(rows)
np.savetxt(OUT/'joint.csv', rows, delimiter=',', header='model,X,Z,Y,A,probability', comments='')
# 对每种结构，比较全体、Z=0、Z=1、A=0、A=1。
answers=[]
for model in range(1,5):
    r=rows[rows[:,0]==model]
    for group in range(5):
        if group == 0:
            sub=r
        else:
            col=2 if group<3 else 4
            val=(group-1)%2
            sub=r[r[:,col]==val]
        w=sub[:,5]/sub[:,5].sum()
        x=sub[:,1]; yy=sub[:,3]
        ex=np.sum(w*x); ey=np.sum(w*yy)
        corr=(np.sum(w*x*yy)-ex*ey)/np.sqrt(ex*(1-ex)*ey*(1-ey))
        answers.append([model,group,corr])
answers=np.array(answers)
np.savetxt(OUT/'associations.csv', answers, delimiter=',', header='model,group,correlation', comments='')
print('理论相关：模型1叉、2管、3对撞、4中介的后代；组0全体、1/2按Z、3/4按A\n',answers)

# 原书171页：有真菌少长3，无真菌平均长5。
plant = np.array([[0,.5,(1-.5)*5+.5*2],[1,.1,(1-.1)*5+.1*2]])
np.savetxt(OUT/'plant_expectation.csv',plant,delimiter=',',header='treatment,fungus_probability,mean_growth',comments='')
print('植物总效应（期望）：',plant[1,2]-plant[0,2])

# %% 4. 原生中文绘图：保存并保留窗口；无界面核验可外设MPLBACKEND=Agg
fonts={f.name for f in font_manager.fontManager.ttflist}
font=next((s for s in ['PingFang SC','Microsoft YaHei','Noto Sans CJK SC','Heiti SC'] if s in fonts),None)
if font is None:
    raise RuntimeError('请安装 Noto Sans CJK SC 中文字体后重启。')
plt.rcParams.update({'font.family':font,'axes.unicode_minus':False,'font.size':12})
colors=['#167e83','#c4853f']; bg='#faf8f3'
def save(fig,name):
    fig.tight_layout()
    fig.savefig(OUT/name,dpi=180,facecolor=bg)
fig,ax=plt.subplots(figsize=(8,4.5),facecolor=bg)
ax.errorbar(theta[1:3], [1,0], xerr=z89*sd[1:3], fmt='o', color=colors[0], capsize=5)
ax.axvline(0,color='gray',linestyle='--');ax.set(yticks=[1,0],yticklabels=['结婚率 βM','结婚年龄 βA'],xlabel='标准化后的斜率（点：近似均值；线：89%区间）',title='年龄相同以后，结婚率的斜率接近零')
save(fig,'01-婚姻模型.png')
fig,axes=plt.subplots(1,2,figsize=(9,4.5),facecolor=bg)
rng=np.random.default_rng(20261008)
# 同一批标准正态数，单独展示先验尺度的影响。固定M=0，仅画条件均值线。
seeds=rng.normal(size=(20,2));grid=np.linspace(-2,2,50)
for ax,scales,title in zip(axes,[(10,10),(.2,.5)],['过宽先验','PPT采用的先验']):
    for u,v in seeds:
        ax.plot(grid,scales[0]*u+scales[1]*v*grid,color=colors[0],alpha=.3)
    ax.set(xlabel='标准化结婚年龄',ylabel='标准化离婚率均值',title=title,ylim=(-50,50))
save(fig,'02-先验均值线.png')
fig,ax=plt.subplots(figsize=(8,4.5),facecolor=bg)
xs=np.linspace(-3.5,3.5,500)
# 平均变化的后验与新增两个独立残差后的预测差别，后者作后验混合。
draws=rng.multivariate_normal(theta,cov,size=10000)
draws=draws[draws[:,3]>0]
predictive=np.mean(norm.pdf(xs[:,None],draws[:,1],np.sqrt(2)*draws[:,3]),axis=1)
ax.plot(xs,norm.pdf(xs,theta[1],sd[1]),label='平均变化的不确定性',color=colors[0])
ax.plot(xs,predictive,label='两次独立预测结果之差',color=colors[1])
ax.set(xlabel='离婚率变化（标准差单位）',ylabel='概率密度',title='预测结果的差，比平均变化更分散');ax.legend()
save(fig,'03-干预比较.png')
fig,ax=plt.subplots(figsize=(8,4.5),facecolor=bg)
selected=answers[(answers[:,0]==4)&np.isin(answers[:,1],[0,1,3])]
ax.bar(['全体','直接固定中介 Z','固定后代 A'],selected[:,2],color=[colors[0],colors[1],colors[0]])
ax.set(ylim=(0,.75),ylabel='X与Y的理论相关',title='看后代的值，会间接获得中介的信息')
for i,v in enumerate(selected[:,2]):ax.text(i,v+.025,f'{v:.3f}',ha='center')
save(fig,'04-后代.png')
plt.show()
