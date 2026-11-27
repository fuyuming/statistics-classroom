# %% 0. 准备：在Spyder以“练习项目”为工作目录，依次运行单元或按F5。
# 共享CSV是教学模拟资料；数据含义与模拟来源见README，文件使用UTF-8。
from pathlib import Path
import os
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt

alpha = 0.05
confidence = 1 - alpha
unit_scale = 1.0  # 第二遍可改为0.133322，将mmHg统一换算成kPa。
save_outputs = os.getenv('SAVE_OUTPUTS', 'false').lower() == 'true'
unit_label = 'mmHg' if unit_scale == 1 else 'kPa' if unit_scale == 0.133322 else '换算单位'
# 整文件可按脚本位置定位；Spyder没有__file__时使用当前项目目录。
project_dir = Path(__file__).resolve().parents[1] if '__file__' in globals() else Path.cwd()
output_dir = project_dir / '运行结果' / 'Python'
figure_dir = project_dir / '图形' / 'Python'

# %% 1. 数据：每行一名受试者，先看字段、人数和各组波动。
h = pd.read_csv(project_dir / 'data' / 'hetero_normal.csv')
# 只在这里换算一次，后续统计量全部用换算后的数据重新计算。
h['value'] = h['value'] * unit_scale
print(h.head())
print(h.groupby('group').size())
names = ['G1', 'G2', 'G3']
# loc按方案选出该组受试者；各组长度可不同。
groups = [h.loc[h['group'] == name, 'value'].to_numpy() for name in names]
ni = np.array([len(x) for x in groups])
mi = np.array([np.mean(x) for x in groups])
# ddof=1让样本方差用n-1作分母；默认ddof=0不适用于这里。
si = np.array([np.std(x, ddof=1) for x in groups])
k = len(groups)
if np.any(ni < 2) or not np.all(np.isfinite(si)) or np.any(si <= 0):
    raise ValueError('每组至少2个完整观测且样本方差必须为正。')

# %% 2. 精度与Welch：从均数方差到加权中心、校正和近似F。
v = si**2 / ni
w = 1 / v
W = np.sum(w)
weighted_mean = np.sum(w * mi) / W
summary_tab = pd.DataFrame({'group': names, 'n': ni, 'mean': mi,
                            'sd': si, 'variance': si**2, 'v': v, 'w': w})
print(summary_tab.to_string(index=False))
# A汇总精度加权的均数偏离；C纳入各组方差估计的不确定性。
A = np.sum(w * (mi - weighted_mean)**2) / (k - 1)
C = np.sum((1 - w / W)**2 / (ni - 1))
f_welch = A / (1 + 2 * (k - 2) * C / (k**2 - 1))
df1 = k - 1
df2 = (k**2 - 1) / (3 * C)
# sf是右尾概率；输入观察到的F及两个自由度。
p_welch = stats.f.sf(f_welch, df1, df2)
# ppf输入累计概率1-alpha，按两个自由度（含小数）返回临界F。
f_critical = stats.f.ppf(confidence, df1, df2)
welch_tab = pd.DataFrame({'weighted_mean': [weighted_mean], 'W': [W],
                          'A': [A], 'C': [C], 'F': [f_welch],
                          'df1': [df1], 'df2': [df2], 'F_critical': [f_critical], 'P': [p_welch]})
print(welch_tab.to_string(index=False))
# SciPy 1.16起，equal_var=False明确采用Welch，返回F及P。
welch_function = stats.f_oneway(*groups, equal_var=False)
print('SciPy Welch对照：', welch_function)
verdict = '拒绝' if p_welch < alpha else '未拒绝'
print(f'在α={alpha}下{verdict}三组总体均数相等；整体结果不指出具体组对。')

# %% 3. GH：i取前组，j取后组，每对使用自己的SE、自由度与q参照。
rows = []
for i in range(k - 1):
    for j in range(i + 1, k):
        difference = mi[j] - mi[i]
        se = np.sqrt(v[i] + v[j])
        # 近似自由度只由当前两个组的均数方差与各自n-1构成。
        nu = (v[i] + v[j])**2 / (v[i]**2 / (ni[i] - 1) + v[j]**2 / (ni[j] - 1))
        # 沿用Welch t的SE和自由度，再将t的绝对值转换到q尺度。
        t_statistic = difference / se
        q = np.sqrt(2) * abs(t_statistic)
        # 两组只有一对，用k=2核对t/q尺度；正式GH仍用全体组数k。
        p_single = 2 * stats.t.sf(abs(t_statistic), nu)
        p_q_two = stats.studentized_range.sf(q, 2, nu)
        # k始终是全部组数；sf给出已调整P，不再追加Bonferroni。
        p_adj = stats.studentized_range.sf(q, k, nu)
        # ppf输入累计概率，返回q临界值。除根号2还原差值SE的尺度。
        q_critical = stats.studentized_range.ppf(confidence, k, nu)
        half_width = q_critical * se / np.sqrt(2)
        rows.append([names[i], names[j], difference, se, nu, q, q_critical,
                     difference - half_width, difference + half_width, p_adj,
                     t_statistic, p_single, p_q_two])
gh_tab = pd.DataFrame(rows, columns=['group1', 'group2', 'difference', 'SE', 'df',
                                    'q', 'q_critical', 'CI_low', 'CI_high', 'P_adj',
                                    't_statistic', 'P_single', 'P_q_k2'])
print('P_single为单对Welch双侧P，P_q_k2用于两组q尺度对照，本题全部三对的调整结果为P_adj。')
print(gh_tab.to_string(index=False))
for row in gh_tab.itertuples():
    verdict = '支持总体均数有差异' if row.P_adj < alpha else '证据不足以确认总体均数有差异'
    print(f'{row.group2}−{row.group1}: 差值={row.difference:.4f} {unit_label}，'
          f'近似{confidence:.0%}同时CI=[{row.CI_low:.4f}, {row.CI_high:.4f}] {unit_label}，'
          f'调整P={row.P_adj:.6g}；{verdict}。')

# %% 4. 完整软件对照：接口名为tukey_hsd，False分支实际执行GH。
gh_function = stats.tukey_hsd(*groups, equal_var=False)
gh_ci = gh_function.confidence_interval(confidence)
# 原生矩阵[j,i]才是后组减前组；打印同方向的三行，便于核对。
function_rows = []
for i in range(k - 1):
    for j in range(i + 1, k):
        function_rows.append([names[j] + '−' + names[i], gh_function.statistic[j, i],
                              gh_ci.low[j, i], gh_ci.high[j, i], gh_function.pvalue[j, i]])
print(pd.DataFrame(function_rows, columns=['comparison', 'difference', 'CI_low', 'CI_high', 'P_adj']))

# %% 5. 图形：圆点为个体，黑色菱形为均数；区间图以0为参照。
from matplotlib import font_manager
available_fonts = {font.name for font in font_manager.fontManager.ttflist}
font_candidates = ['PingFang SC', 'Heiti SC', 'Heiti TC', 'Microsoft YaHei',
                   'Noto Sans CJK SC', 'SimHei', 'Arial Unicode MS']
chinese_fonts = [name for name in font_candidates if name in available_fonts]
if not chinese_fonts:
    raise RuntimeError('请安装Noto Sans CJK SC中文字体，重启Spyder后再运行绘图节。')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = chinese_fonts
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['font.size'] = 12
fig1, ax1 = plt.subplots(figsize=(7.6, 5.2))
for i, x in enumerate(groups):
    ax1.scatter(i + np.linspace(-0.15, 0.15, len(x)), x, alpha=0.75, label=names[i])
ax1.scatter(range(k), mi, color='black', marker='D', label='组均数')
ax1.set_xticks(range(k), names)
ax1.set_xlabel('降压方案')
ax1.set_ylabel(f'收缩压下降值（{unit_label}）')
ax1.set_title('个体下降值与各组均数')
ax1.legend()
fig1.tight_layout()
fig2, ax2 = plt.subplots(figsize=(7.6, 5.2))
positions = np.arange(len(gh_tab))
ax2.errorbar(gh_tab['difference'], positions,
             xerr=[gh_tab['difference'] - gh_tab['CI_low'], gh_tab['CI_high'] - gh_tab['difference']],
             fmt='o', capsize=4)
ax2.set_yticks(positions, [r.group2 + '−' + r.group1 for r in gh_tab.itertuples()])
ax2.invert_yaxis()
ax2.axvline(0, color='grey', linestyle='--')
ax2.set_xlabel(f'平均下降值之差（{unit_label}）')
ax2.set_ylabel('方案比较')
ax2.set_title(f'均数差与近似{confidence:.0%}同时区间')
fig2.tight_layout()

# %% 6. 可选保存：关闭时只在IDE看结果；开启后每种语言各写自己的目录。
if save_outputs:
    output_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    summary_tab.to_csv(output_dir / 'group_summary.csv', index=False)
    welch_tab.to_csv(output_dir / 'welch_intermediate.csv', index=False)
    gh_tab.to_csv(output_dir / 'games_howell.csv', index=False)
    fig1.savefig(figure_dir / 'raw_points_means.png', dpi=180)
    fig2.savefig(figure_dir / 'games_howell_ci.png', dpi=180)
plt.show()
