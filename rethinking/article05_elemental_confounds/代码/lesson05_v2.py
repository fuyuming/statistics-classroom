# 【课堂问题】结婚率较高的地区，离婚率也较高。这种关联是否混入了结婚年龄的影响？
# 【数据】作者WaffleDivorce数据，49州和哥伦比亚特区，共50行。
# 【计算思路】①统一尺度；②写出似然与先验；③近似联合后验；④解释相同年龄下的斜率。
# 【第二个问题】为什么控制第三个变量，有时消除关联，有时反而制造关联？
# 【核对办法】按四种结构的生成规则枚举0/1状态，比较全体与分组后的相关。
# 建议先运行第一部分并读后验结果，再进入第二部分；数值输出和图形要结合当前问题阅读。

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
# 列顺序决定后续模型：Divorce离婚率、Marriage结婚率（均为每千名成年人）、
# MedianAgeMarriage结婚年龄中位数（岁）；每行一个州或特区。
# 换数据先核对这三个字段、单位和缺失值。
raw = np.column_stack([d['Divorce'], d['Marriage'], d['MedianAgeMarriage']])
print('前5行：离婚率、结婚率、结婚年龄中位数\n', raw[:5])
# ddof=1：采用样本标准差，分母是观测数减1；与R的scale一致。
standard = (raw - raw.mean(axis=0)) / raw.std(axis=0, ddof=1)
y = standard[:, 0]
# 这里X是设计矩阵（不是后文因果图中的二元X）：三列是常数1、M、A。
# 常数列乘截距；另两列乘各自斜率。y是标准化离婚率。
X = np.column_stack([np.ones(len(y)), standard[:, 1:]])
# 先验怎么填：依次是alpha、beta_M、beta_A的标准差，不能填方差。
# 这里三者先验均值都为0；sigma的指数先验另写在objective最后一项中。
prior_sd = np.array([0.2, 0.5, 0.5])

# %% 2. 后验：似然乘先验；在峰顶附近作二次近似
# 前三个参数是截距、结婚率斜率、年龄斜率，最后一个是残差标准差。
def objective(theta):
    beta = theta[:3]
    sigma = theta[3]
    if sigma <= 0:
        return np.inf
    residual = y - X @ beta
    # 去掉不依赖参数的常数后，负对数后验分成四项：
    # n*log(sigma)及残差平方项来自正态似然；beta项来自正态先验；
    # 最后的sigma来自速率为1的指数先验。优化整个式子得到MAP。
    return len(y)*np.log(sigma) + np.sum(residual**2)/(2*sigma**2) + np.sum((beta/prior_sd)**2)/2 + sigma
# [0,0,-0.5,0.8]是优化初值，顺序alpha/beta_M/beta_A/sigma；sigma须为正。
# 初值不是先验，修改先验尺度应改prior_sd。若要改指数先验速率lambda，
# 负对数后验最后一项需改为lambda*sigma（与参数无关的常数可略）。
# 本例解析Hessian与零均值正态/指数先验配套，换分布需重新推导，不能只改标签。
fit = minimize(objective, [0, 0, -0.5, 0.8], method='BFGS', options={'gtol': 1e-6})
if not fit.success:
    raise RuntimeError(fit.message)
theta = fit.x  # 后验峰值：四参数的联合近似以此为中心
beta = theta[:3]  # alpha、beta_M、beta_A，Python索引从0开始
sigma = theta[3]  # 第4个参数为残差标准差
residual = y - X @ beta
# Hessian：峰顶曲率决定近似后验的宽度，也保留参数间的相关性。
H = np.zeros((4,4))
H[:3,:3] = X.T @ X / sigma**2 + np.diag(1/prior_sd**2)
H[:3,3] = 2*X.T @ residual / sigma**3
H[3,:3] = H[:3,3]
H[3,3] = -len(y)/sigma**2 + 3*np.sum(residual**2)/sigma**4
# 逆Hessian是联合正态近似的协方差；不是只输出一组最优参数。
# 此处theta既是MAP，也用作近似分布的均值；没有进行MCMC。
cov = np.linalg.inv(H)  # 方差在对角线，非对角线保留参数间协方差
sd = np.sqrt(np.diag(cov))
# 中间89%的区间，两端各留5.5%。
# 区间参数：中央概率c使用(1+c)/2分位点；89%填0.945，95%填0.975。
# 修改覆盖概率时须同步输出列名、图注及所有使用z89的绘图说明。
z89 = norm.ppf(0.945)
# 【结果核对】四行依次对应四参数，四列依次为均值、SD、下限、上限。
summary = np.column_stack([theta, sd, theta-z89*sd, theta+z89*sd])
np.savetxt(OUT/'posterior.csv', summary, delimiter=',', header='mean,sd,lower89,upper89', comments='')
np.savetxt(OUT/'covariance.csv', cov, delimiter=',')
# 读posterior.csv：mean是近似后验均值，sd是参数的后验标准差。
# sigma的mean约0.785描述数据剩余波动，sigma的sd描述对该波动大小的不确定性。
# beta_M约-0.065，89%区间约[-0.306,0.176]，方向仍不确定，不等于证明作用为0。
# beta_A约-0.614：固定M时，A增加1个标准差，平均D降低0.614个标准差。
print('行顺序 alpha, beta_M, beta_A, sigma；列 mean,sd,lower89,upper89\n', summary)
print('beta_M：年龄相同、结婚率增加1个标准差时，平均离婚率的变化。')

# %% 3. 枚举作者二元生成规则的所有可能组合
# 0和1只是状态标签。伯努利概率：状态是1取p，是0取1-p。
def bern(value, probability):
    return probability if value == 1 else 1-probability
# 此处四个字母都是二元状态，A不再指结婚年龄，X不再指前面的设计矩阵。
# 0.1+0.8*z在z=0/1时分别给出0.1/0.9的概率。
# 枚举的是可能数据状态，概率已知；这与网格近似未知参数的后验是两回事。
rows = []
for model in range(1,5):
    for x,z,yv,a in itertools.product([0,1], repeat=4):
        if model == 1: # 叉：Z分别影响X和Y；给定Z后独立生成
            w = .5 * bern(x,.1+.8*z) * bern(yv,.1+.8*z) * .5
        elif model == 2: # 管：X影响Z，Z影响Y
            w = .5 * bern(z,.1+.8*x) * bern(yv,.1+.8*z) * .5
        elif model == 3: # 对撞：X、Y独立，任一为1都会提高Z=1的概率
            w = .25 * bern(z,.9 if x+yv>0 else .2) * .5
        else: # 中间变量的后代A：90%机会与Z一致
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
        # 分组后重新归一化，得到组内条件概率；0/1变量的方差是E(X)*(1-E(X))。
        w=sub[:,5]/sub[:,5].sum()
        x=sub[:,1]; yy=sub[:,3]
        ex=np.sum(w*x)  # 组内X=1的概率
        ey=np.sum(w*yy)  # 组内Y=1的概率
        corr=(np.sum(w*x*yy)-ex*ey)/np.sqrt(ex*(1-ex)*ey*(1-ey))
        answers.append([model,group,corr])
# 练习只改后代模型最后的bern(a,.1+.8*z)为bern(a,.5)。
# 这让A不再携带Z的信息：固定A后的相关应从约0.390回到全体的0.640。
answers=np.array(answers)
np.savetxt(OUT/'associations.csv', answers, delimiter=',', header='model,group,correlation', comments='')
print('理论相关：模型1叉、2管、3对撞、4中间变量的后代；组0全体、1/2按Z、3/4按A\n',answers)

# 原书171页：有真菌少长3，无真菌平均长5。
# 每行依次为处理状态、真菌概率、平均增长。无真菌长5，有真菌长2，
# 用各状态概率加权得到3.5和4.7；差1.2是理论期望，不是后验拟合结果。
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
# 从联合正态近似抽10000组参数，保留协方差；这是近似后验抽样，不是MCMC。
# 两个独立预测残差相减，方差相加为2*sigma²，因此下面标准差用sqrt(2)*sigma。
# 青色线看平均变化的不确定性；赭色线还包含未来观测的随机波动。
draws=rng.multivariate_normal(theta,cov,size=10000)
draws=draws[draws[:,3]>0]
# 【把do(M)展开为两次计算】每行用同一组参数、同一个年龄，只改变M。
# 从观测年龄中有放回抽取，保留本数据的年龄构成；独立种子便于复算。
intervention_rng = np.random.default_rng(20261009)
age_same = intervention_rng.choice(standard[:, 2], size=len(draws), replace=True)
M_before = 0.0  # 标准化结婚率0：原始均值约20.11/千名成年人
M_after = 1.0   # 标准化结婚率1：原始均值加1个SD，约23.91/千名成年人
baseline = draws[:, 0] + draws[:, 2] * age_same  # alpha + beta_A*A，两次共用
mu_before = baseline + draws[:, 1] * M_before
mu_after = baseline + draws[:, 1] * M_after
mean_change = mu_after - mu_before
# 把M的差改成其他值时，均值差应为beta_M乘这个差；这里等于beta_M。
assert np.allclose(mean_change, draws[:, 1] * (M_after-M_before))
# PPT还分别生成两次观测值，各加一个独立正态残差。
# 正态差的密度可直接计算：均值为mean_change，SD为sqrt(2)*sigma。
# 图5用密度平均得到平滑曲线，与反复模拟两次观测再相减针对同一分布。
np.savetxt(OUT/'intervention_pairs.csv',
           np.column_stack([age_same, mu_before, mu_after, mean_change]), delimiter=',',
           header='age_standardized,mu_M0,mu_M1,mean_change', comments='')
print('同一年龄下两种情境的前5行：年龄、M=0均值、M=1均值、均值差\n',
      np.column_stack([age_same, mu_before, mu_after, mean_change])[:5])
predictive=np.mean(norm.pdf(xs[:,None],mean_change,np.sqrt(2)*draws[:,3]),axis=1)
ax.plot(xs,norm.pdf(xs,theta[1],sd[1]),label='平均变化的不确定性',color=colors[0])
ax.plot(xs,predictive,label='两次独立预测结果之差',color=colors[1])
ax.set(xlabel='离婚率变化（标准差单位）',ylabel='概率密度',title='预测结果的差，比平均变化更分散');ax.legend()
save(fig,'03-干预比较.png')
fig,ax=plt.subplots(figsize=(8,4.5),facecolor=bg)
selected=answers[(answers[:,0]==4)&np.isin(answers[:,1],[0,1,3])]
ax.bar(['全体','直接固定中间变量 Z','固定后代 A'],selected[:,2],color=[colors[0],colors[1],colors[0]])
ax.set(ylim=(0,.75),ylabel='X与Y的理论相关',title='看后代的值，会间接获得中间变量的信息')
for i,v in enumerate(selected[:,2]):ax.text(i,v+.025,f'{v:.3f}',ha='center')
save(fig,'04-后代.png')
plt.show()
