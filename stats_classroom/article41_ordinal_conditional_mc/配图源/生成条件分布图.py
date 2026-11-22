# 先在本篇项目根运行R脚本并打开保存开关；本图使用同一R模拟流。
from pathlib import Path
import sys
import pandas as pd
article = Path(__file__).resolve().parent.parent
project = article / '配套项目/article41_ordinal_conditional_mc'
if not project.is_dir():
    # GitHub中练习项目与正文位于同一目录；本地完整材料保留配套项目子目录。
    project = article
if not (project / 'results/R_mini_exact_distribution.csv').is_file():
    raise FileNotFoundError('请先在本篇项目运行R脚本并开启保存，生成results中的枚举与模拟结果。')
sys.path.insert(0, str(project / 'scripts'))
from plot_ordinal_conditional_mc import plot_results
exact = pd.read_csv(project / 'results/R_mini_exact_distribution.csv')
observed = pd.read_csv(project / 'results/R_observed_stats.csv', index_col=0).T.to_dict()
sim = pd.read_csv(project / 'results/R_mini_sim_T.csv')['T'].to_numpy()
plot_results(exact, observed, {'mini': {'T_sim': sim}}, article / '文章配图/图1_条件分布.png')
