# %% 问题：真实提高0.40 cm时，100尾的功效是多少？至少需要多少尾？
# 来源：李春喜《生物统计学》第6版例4-1的数值，另拟右侧检验设计。
# Python依赖见requirements.txt。Spyder打开本文件，F5运行。
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import norm
import matplotlib.pyplot as plt
root = Path(__file__).resolve().parents[1]
out = root / '运行结果' / 'Python'
out.mkdir(parents=True, exist_ok=True)
parameters = pd.read_csv(root / '数据' / '设计参数.csv')
print('教学设计参数（cm），不是原始测量数据：')
print(parameters)
mu0 = parameters.loc[0, 'mu0']
mu1 = parameters.loc[0, 'mu1']
sigma = parameters.loc[0, 'sigma']
alpha = parameters.loc[0, 'alpha']
target = parameters.loc[0, 'target_power']
delta = mu1 - mu0
d = delta / sigma

# %% 1. 先用H0边界确定拒绝门槛，再在指定H1下计算漏判概率
n = 100
se = sigma / np.sqrt(n)
# ppf输入左侧累计概率，输出标准正态分位数。
cutoff = mu0 + norm.ppf(1 - alpha) * se
# cdf计算左侧概率：真实均数为mu1时，未超过门槛的机会。
beta = norm.cdf((cutoff - mu1) / se)
power = 1 - beta

# %% 2. 反求样本量，并检查相邻整数
n_continuous = ((norm.ppf(1-alpha) + norm.ppf(target)) / d)**2
n_min = int(np.ceil(n_continuous))
# sf计算右侧概率。功效是指定真实差异下的拒绝概率，不是H1为真的概率。
power_previous = norm.sf(norm.ppf(1-alpha) - d*np.sqrt(n_min-1))
power_minimum = norm.sf(norm.ppf(1-alpha) - d*np.sqrt(n_min))
values = [se, cutoff, beta, power, d, n_continuous, n_min, power_previous, power_minimum]
names = ['se','cutoff','beta','power','d','n_continuous','n_min','power_previous','power_minimum']
result = pd.DataFrame({'quantity': names, 'value': values})
result.to_csv(out / '结果.csv', index=False)
print(result)
print('n_min是当前设计条件下达到目标功效的最小整数样本量。')
assert power_previous < target <= power_minimum

# %% 3. 区分设计功效与观察结果：三个假设研究，观察差值均为0.40 cm
observed_n = np.array([25, 100, 400])
observed_difference = 0.40
observed_se = sigma / np.sqrt(observed_n)
z = observed_difference / observed_se
p_two_sided = 2 * norm.sf(np.abs(z))
ci_low = observed_difference - norm.ppf(.975) * observed_se
ci_high = observed_difference + norm.ppf(.975) * observed_se
observed = pd.DataFrame({'n':observed_n, 'difference':observed_difference, 'p_two_sided':p_two_sided, 'ci_low':ci_low, 'ci_high':ci_high})
observed.to_csv(out / '观察结果情境.csv', index=False)
print(observed)
print('这一节改为双侧问题；三种情境的观察差值相同，P值不同。')

# %% 4. 生成文章图。排版放在配套文件中，计算步骤保留在本入口。
exec((root / '代码' / 'plot_power.py').read_text(encoding='utf-8'))
exec((root / '代码' / 'plot_threshold.py').read_text(encoding='utf-8'))
exec((root / '代码' / 'render_formulas.py').read_text(encoding='utf-8'))
plt.show()
