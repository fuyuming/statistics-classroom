"""统计课堂31：两独立率的大样本近似与 Fisher 精确检验。
题干：患者独立随机分组，各100人；4周同一标准下达标50/30。
目标是总体达标率差；双侧 H0: p_new=p_control，alpha=0.05。
"""
from pathlib import Path
import math
import numpy as np
import pandas as pd
from scipy.stats import norm, chi2_contingency, hypergeom, fisher_exact

SAVE_OUTPUT = False
HERE = Path(__file__).resolve().parent if "__file__" in globals() else Path.cwd()
dat = pd.read_csv(HERE / "data" / "large_sample_table.csv")

# %% 1. 读入数据并按组名取值
new = dat.loc[dat["group"].eq("new")].iloc[0]
control = dat.loc[dat["group"].eq("control")].iloc[0]
x_new = int(new.event)
n_new = int(new.total)
x_control = int(control.event)
n_control = int(control.total)
p_new = x_new/n_new
p_control = x_control/n_control
RD = p_new-p_control
observed_rates = pd.DataFrame({"p_new":[p_new], "p_control":[p_control], "rate_difference":[RD]})
print("观察率与率差：\n", observed_rates.to_string(index=False))

# %% 2. H0 合并率、期望表与 npq 检查
# H0 假定两组共享一个达标率，所以用总达标数/总人数；期望数和 npq 检查近似条件。
p_pool = (x_new+x_control)/(n_new+n_control)
expected_success = np.array([n_new,n_control])*p_pool
expected_failure = np.array([n_new,n_control])*(1-p_pool)
npq = np.array([n_new,n_control])*p_pool*(1-p_pool)
null_checks = pd.DataFrame({"group":["new","control"], "expected_success":expected_success, "expected_failure":expected_failure, "npq":npq})
print(f"H0合并率={p_pool:.15f}；期望表与npq（期望数和npq均满足5/10经验检查）：\n", null_checks.to_string(index=False))

# %% 3. 手算 Z、双侧正态 P
# SE0 使用 H0 合并率；Z 是观察率差相对 H0 标准误，双侧 P 为两侧正态尾部之和。
SE0 = math.sqrt(p_pool*(1-p_pool)*(1/n_new+1/n_control))
Z = RD/SE0
P_normal = 2*norm.sf(abs(Z))
manual_z = pd.DataFrame({"p_pool":[p_pool], "SE0":[SE0], "Z":[Z], "P_normal":[P_normal], "alpha":[0.05]})
print("H0下合并率用于SE0；正态尾函数给双侧P：\n", manual_z.to_string(index=False))

# %% 4. 完整近似函数核对：卡方检验不作连续性校正
# 无校正的 2x2 卡方检验与双侧 Z 检验等价，用函数结果核对手算 P。
table = np.array([[x_new, int(new.no_event)], [x_control, int(control.no_event)]])
chi2, P_chi2, dof, expected = chi2_contingency(table, correction=False)
chi2_check = pd.DataFrame({"P_normal":[P_normal], "P_chi2_no_correction":[P_chi2]})
print("完整近似函数核对（无连续性校正）：\n", chi2_check.to_string(index=False))

# %% 5. 固定边际下超几何分布与概率排序双侧 P
# 固定总达标数80和新方案人数100，A 是新方案达标数；按表概率不大于观察表求双侧P。
N = n_new+n_control
K = x_new+x_control
draw = n_new
A = np.arange(0, K+1)
q_comb = np.array([math.comb(K, int(a))*math.comb(N-K, draw-int(a))/math.comb(N, draw) for a in A])
q = hypergeom.pmf(A, N, K, draw)
q_observed = float(hypergeom.pmf(x_new, N, K, draw))
P_hyper_two_sided = float(q[q <= q_observed*(1+1e-7)].sum())
E_A = float((A*q).sum())
Var_A = float(((A-E_A)**2*q).sum())
hypergeom_check = pd.DataFrame({"sum_probability":[q.sum()], "E_A":[E_A], "Var_A":[Var_A], "P_ordered_two_sided":[P_hyper_two_sided]})
print("固定边际的超几何分布：概率总和应为1；概率排序得到双侧P：\n", hypergeom_check.to_string(index=False))

# %% 6. Fisher 精确检验核对（只报告 P）
oddsratio, fisher_p = fisher_exact(table, alternative="two-sided")
fisher_check = pd.DataFrame({"P_hyper":[P_hyper_two_sided], "P_fisher":[fisher_p]})
print("完整 Fisher 函数核对（只报告P）：\n", fisher_check.to_string(index=False))

# %% 7. 可选保存结果（不绘图）
if SAVE_OUTPUT:
    out = HERE / "results" / "python"
    out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"A":A, "q_combination":q_comb, "q_hypergeom":q}).to_csv(out/"large_sample_probability.csv", index=False)
    pd.DataFrame({"RD":[RD],"SE0":[SE0],"Z":[Z],"P_normal":[P_normal],"P_chi2":[P_chi2],"P_hyper":[P_hyper_two_sided],"P_fisher":[fisher_p],"E_A":[E_A],"Var_A":[Var_A]}).to_csv(out/"large_sample_comparison.csv", index=False)
