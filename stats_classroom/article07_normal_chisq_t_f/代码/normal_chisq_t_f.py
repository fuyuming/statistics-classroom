# 明哥的微生物世界 · 统计课堂 07 · Python / Spyder
# 学习目标：①认识抽样分布；②理解χ²、t、F的构造；③核对t²=F。
# 本课数据全部为教学设定，不是真实实验，不用结果推断生物学结论。
# 完整下载项目，在 Spyder 打开本文件后按 F5。文件保存为 UTF-8。
# 先读本脚本，复杂的绘图排版在 plot_distributions.py，最后会自动调用。
# 需要 numpy、pandas、scipy、matplotlib，安装步骤见项目 README。

# %% 0. 准备工具与目录：第一次先完整运行，以后可按分节从上到下练习
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import scipy
from scipy import stats
import matplotlib
import matplotlib.pyplot as plt

# __file__ 是本脚本；上一级是“代码”，再上一级才是项目目录。
project_dir = Path(__file__).resolve().parents[1]
out_dir = project_dir / '运行结果' / 'Python'
out_dir.mkdir(parents=True, exist_ok=True)

# 一个小工具：同时把说明显示在控制台，并收集起来留作运行记录。
messages = []
def explain(message):
    print(message)
    messages.append(str(message))

# %% 1. 均数的“抽样波动”，和个体的“离散程度”有什么区别？
# 设一个正态总体：总体标准差 sigma=2，独立样本量 n=10。
# 标准误描述样本均数的波动，不是单个观测的波动。
sigma = 2.0
sample_size = 10
standard_error = sigma / np.sqrt(sample_size)
explain('\n问题1：总体标准差为2，取10个独立观测，均数的标准误是多少？')
explain(f'答案：2 / sqrt(10) = {standard_error:.6f}。这里总体标准差是假设已知的。')

# %% 2. 为什么自由度10的t，不能随手拿1.96作界值？
# 双侧检验水平alpha=0.05：左右各留0.025，中间覆盖0.95。
# ppf 是分位点函数：输入左侧累计概率，返回横轴上的界值。
# sf 是右尾概率函数：sf(a) = P(T > a)。
degrees_freedom = 10
alpha = 0.05
t_cutoff = stats.t.ppf(1 - alpha / 2, df=degrees_freedom)
f_cutoff = stats.f.ppf(1 - alpha, dfn=1, dfd=degrees_freedom)
wrong_error = 2 * stats.t.sf(1.96, df=degrees_freedom)
explain('\n问题2：t的双侧5%界值，与F的上侧5%界值怎样对应？')
explain(f't界值 = {t_cutoff:.7f}；平方 = {t_cutoff**2:.7f}。')
explain(f'F(1,10)界值 = {f_cutoff:.7f}，与上面的平方相同。')
explain(f'若T实际服从t(10)，错用±1.96：第一类错误率 = {wrong_error:.2%}。')
explain('这个7.84%是零假设下的理论概率，不是一次实验的P值。')

# 多算几种自由度。列表中的每个元素代表一个设定。
# 用循环逐行组织表格，比把许多数组堆在一起更容易逐步检查。
df_values = [1, 5, 10, 30, 100]
critical_rows = []
error_rows = []
for df_value in df_values:
    t_limit = stats.t.ppf(0.975, df=df_value)
    f_limit = stats.f.ppf(0.95, dfn=1, dfd=df_value)
    correct_rate = 2 * stats.t.sf(t_limit, df=df_value)
    wrong_rate = 2 * stats.t.sf(1.96, df=df_value)
    critical_rows.append({
        'df': df_value, 'n_one_sample': df_value + 1,
        't_975': t_limit, 't_975_squared': t_limit**2,
        'f_95': f_limit, 'rejection_at_1_96': wrong_rate,
    })
    error_rows.append({
        'df': df_value, 'n_one_sample': df_value + 1,
        'correct_t': correct_rate, 'wrong_1_96': wrong_rate,
    })
critical = pd.DataFrame(critical_rows)
error_rates = pd.DataFrame(error_rows)
explain(critical.to_string(index=False))
critical.to_csv(out_dir / 'critical_values.csv', index=False, encoding='utf-8-sig')
error_rates.to_csv(out_dir / 'error_rates.csv', index=False, encoding='utf-8-sig')
# allclose允许浮点计算中的极小误差；不是用“约等于”替代统计判断。
assert np.allclose(critical['t_975_squared'], critical['f_95'], rtol=1e-10)
assert np.allclose(error_rates['correct_t'], 0.05, atol=1e-10)

# %% 3. 读入同一份教学数据，先看清每一列
# sample_id：虚构编号；group：A/B组；x：教学自变量；response：教学响应。
# 这是人为编写的数据，没有缺失值，不代表真实实验和效应。
data = pd.read_csv(project_dir / '数据' / 'teaching_data.csv', encoding='utf-8-sig')
explain('\n问题3：两组均数比较，t检验和方差分析会给出同一结果吗？')
explain(data.head().to_string(index=False))
assert not data[['x', 'response']].isna().any().any()
assert data['group'].isin(['A', 'B']).all()

# loc[条件, 列名]：只取符合条件的行与指定的一列。
group_a = data.loc[data['group'] == 'A', 'response']
group_b = data.loc[data['group'] == 'B', 'response']
explain(f'A组n={len(group_a)}，均数={group_a.mean():.4f}。')
explain(f'B组n={len(group_b)}，均数={group_b.mean():.4f}。')

# equal_var=True：明确选择“合并方差”的两独立样本t检验。
# 默认是双侧检验。本例检验两个总体均数相等，前提由教学模型设定。
t_result = stats.ttest_ind(group_a, group_b, equal_var=True)
anova_result = stats.f_oneway(group_a, group_b)

# A编码为0、B编码为1。含截距的组别回归，斜率就是B均数减A均数。
# linregress用于“一个解释变量+截距”的线性回归，不是多元回归函数。
group_code = (data['group'] == 'B').astype(int)
group_regression = stats.linregress(group_code, data['response'])
group_t = group_regression.slope / group_regression.stderr
comparison = pd.DataFrame({
    'method': ['pooled_t_squared', 'one_way_ANOVA_F', 'regression_group_t_squared'],
    'statistic': [t_result.statistic**2, anova_result.statistic, group_t**2],
    'p_value': [t_result.pvalue, anova_result.pvalue, group_regression.pvalue],
})
explain(comparison.to_string(index=False))
explain('解释：t平方与F约为19.7388；双侧t与F的P值相同。t的符号受相减顺序影响。')
comparison.to_csv(out_dir / 't_anova_regression_equivalence.csv', index=False, encoding='utf-8-sig')
assert np.allclose(comparison['statistic'], comparison['statistic'].iloc[0], rtol=1e-10)
assert np.allclose(comparison['p_value'], comparison['p_value'].iloc[0], rtol=1e-10)

# %% 4. 换成连续自变量：一元回归的斜率t与整体F怎样对应？
# 这里单独演示response~x，不是同时控制group和x的多元模型。
regression = stats.linregress(data['x'], data['response'])
slope_t = regression.slope / regression.stderr
fitted_response = regression.intercept + regression.slope * data['x']
residual = data['response'] - fitted_response
# 残差平方和：模型还没有解释的波动。
residual_ss = np.sum(residual**2)
# 回归平方和：拟合值围绕总体样本均数的波动。
regression_ss = np.sum((fitted_response - data['response'].mean())**2)
residual_df = len(data) - 2  # 已估计截距、斜率两个系数。
regression_ms = regression_ss / 1
residual_ms = residual_ss / residual_df
overall_f = regression_ms / residual_ms
regression_table = pd.DataFrame({
    'slope_t': [slope_t], 'slope_t_squared': [slope_t**2],
    'overall_F': [overall_f], 'residual_df': [residual_df],
})
explain('\n问题4：一元回归的t²和整体F是否一致？')
explain(regression_table.to_string(index=False))
explain('答案：两者约为119.3944。多元回归中，单个系数的t²不等于模型整体F。')
regression_table.to_csv(out_dir / 'simple_regression_equivalence.csv', index=False, encoding='utf-8-sig')
assert np.isclose(slope_t**2, overall_f, rtol=1e-10)

# %% 5. 自己“造”出卡方、t和F：随机模拟帮助理解，不代替数学证明
# 一个随机种子方便本语言重复运行；不同语言的同种子不保证同一序列。
rng = np.random.default_rng(20260927)
repetitions = 50000
nu = 10
# 每行10个独立标准正态值。axis=1表示沿着一行求和。
normal_values = rng.normal(size=(repetitions, nu))
chi_values = np.sum(normal_values**2, axis=1)
# 再独立生成分子Z，不能从上面同一行中拿一个值充当“独立”分子。
z_values = rng.normal(size=repetitions)
t_values = z_values / np.sqrt(chi_values / nu)
f_values = z_values**2 / (chi_values / nu)
assert np.allclose(t_values**2, f_values)
probabilities = [0.025, 0.5, 0.975]
simulation_summary = pd.DataFrame({
    'probability': probabilities,
    'empirical_t': np.quantile(t_values, probabilities),
    'theoretical_t': stats.t.ppf(probabilities, df=nu),
})
explain('\n问题5：模拟t的分位点，与理论t(10)接近吗？')
explain(simulation_summary.to_string(index=False))
explain('经验分位点不必逐位等于理论值；重复次数增多通常能减小模拟误差。')
simulation_summary.to_csv(out_dir / 'construction_simulation_quantiles.csv', index=False, encoding='utf-8-sig')
simulation_data = pd.DataFrame({
    'z': z_values, 'u_chisq10': chi_values, 't10': t_values, 'f1_10': f_values,
})
simulation_data.to_csv(out_dir / 'construction_draws.csv.gz', index=False, compression='gzip')

# %% 6. 生成五张图：先跑通上面的统计，再去看配套绘图文件
# 将“代码”目录加入本次Python会话的模块搜索路径，不改系统环境。
code_dir = str(Path(__file__).resolve().parent)
if code_dir not in sys.path:
    sys.path.insert(0, code_dir)
from plot_distributions import draw_figures

draw_figures(project_dir)
explain('\n读图：01构造关系；02 t与正态；03卡方；04 F；05错用界值的错误率。')
explain('练习：改变自由度，看尾部变化；不要把密度高度读成概率。')
explain('全部自动核对通过。理论结果、曲线坐标与模拟数据已分别保存。')
versions = f'Python {sys.version.split()[0]}; NumPy {np.__version__}; pandas {pd.__version__}; SciPy {scipy.__version__}; Matplotlib {matplotlib.__version__}'
(out_dir / '运行记录.txt').write_text('\n'.join(messages) + '\n' + versions, encoding='utf-8')
# Spyder可在Plots中查看。命令行环境若无图形界面，可设置MPLBACKEND=Agg运行。
if matplotlib.get_backend().lower() != "agg":
    plt.show()
else:
    plt.close("all")
