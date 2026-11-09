# 三类人数的完整概率表：20人、+1/0/-1三类（UTF-8）。
#%% 读取路径与教学参数
from pathlib import Path
import math
import pandas as pd
from scipy.stats import multinomial

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent
DATA = PROJECT / "数据"
OUT = PROJECT / "运行结果" / "Python"
OUT.mkdir(parents=True, exist_ok=True)
n = 20
theta = -0.45
p_t = 0.44
q_plus = (1 - p_t + theta) / 2
q_zero = p_t
q_minus = (1 - p_t - theta) / 2

#%% 观察组合：先算一种具体安排，再乘以安排种数
r_obs, t_obs, s_obs = 1, 14, 5
observed_arrangements = math.comb(n, r_obs) * math.comb(n-r_obs, t_obs)
single_observation = q_plus**r_obs * q_zero**t_obs * q_minus**s_obs
observed_probability = observed_arrangements * single_observation
observed = pd.DataFrame([dict(n_plus=r_obs, n_zero=t_obs, n_minus=s_obs,
    arrangements=observed_arrangements, single_arrangement_probability=single_observation,
    probability=observed_probability)])
print(observed.to_string(index=False))

#%% 双循环列出全部231个三类人数
tail = pd.read_csv(DATA / "wang_lower_tail_n20.csv")
tail_set = set(zip(tail.n_plus, tail.n_zero, tail.n_minus))
rows = []
for r in range(n + 1):
    for t in range(n - r + 1):
        s = n - r - t
        arrangements = math.comb(n, r) * math.comb(n-r, t)
        single = q_plus**r * q_zero**t * q_minus**s
        probability = float(multinomial.pmf([r, t, s], n, [q_plus, q_zero, q_minus]))
        if not math.isclose(probability, arrangements * single, rel_tol=1e-12, abs_tol=1e-15):
            raise RuntimeError("multinomial与拆解公式不一致")
        rows.append(dict(n_plus=r, n_zero=t, n_minus=s, arrangements=arrangements,
                         single_arrangement_probability=single, probability=probability,
                         in_lower_tail=(r, t, s) in tail_set))
table = pd.DataFrame(rows)
table.to_csv(OUT / "三类人数的231组概率.csv", index=False, encoding="utf-8")
tail_probability = table.loc[table.in_lower_tail, "probability"].sum()
all_probability = table.probability.sum()
summary = pd.DataFrame([dict(n=n, theta=theta, p_t=p_t, q_plus=q_plus, q_zero=q_zero,
    q_minus=q_minus, combinations=len(table), tail_combinations=int(table.in_lower_tail.sum()),
    other_combinations=int((~table.in_lower_tail).sum()), all_probability=all_probability,
    tail_probability=tail_probability, other_probability=all_probability-tail_probability,
    observed_probability=observed_probability, observed_arrangements=observed_arrangements)])
summary.to_csv(OUT / "三类人数概率汇总.csv", index=False, encoding="utf-8")
print(summary.to_string(index=False))

#%% 固定Wang下限尾部的候选theta核对：4001点后在峰值邻域局部优化
from scipy.optimize import minimize_scalar
candidate_theta = [-0.45, -0.44923, -0.44]
cutoff = 0.025
if len(tail) != 141 or ((tail[["n_plus", "n_zero", "n_minus"]] < 0).any().any()
                       or (tail[["n_plus", "n_zero", "n_minus"]] % 1 != 0).any().any()
                       or (tail[["n_plus", "n_zero", "n_minus"]].sum(axis=1) != n).any()):
    raise ValueError("wang_lower_tail_n20.csv必须恰好有141行，且每行是非负整数并且总和为20")
coefficient = math.factorial(n) / (tail.n_plus.map(math.factorial) * tail.n_zero.map(math.factorial) * tail.n_minus.map(math.factorial))

def tail_probability(theta_value, p_value):
    q_plus_value = (1 - p_value + theta_value) / 2
    q_minus_value = (1 - p_value - theta_value) / 2
    if min(q_plus_value, p_value, q_minus_value) < -1e-12:
        return float("nan")
    q_plus_value, p_value, q_minus_value = max(0, q_plus_value), max(0, p_value), max(0, q_minus_value)
    return float((coefficient * q_plus_value**tail.n_plus * p_value**tail.n_zero * q_minus_value**tail.n_minus).sum())

candidate_rows = []
for theta_value in candidate_theta:
    grid = [i * (1 - abs(theta_value)) / 4000 for i in range(4001)]
    values = [tail_probability(theta_value, p) for p in grid]
    peak = max(range(len(values)), key=values.__getitem__)
    left = grid[max(0, peak - 1)]
    right = grid[min(len(grid) - 1, peak + 1)]
    refined = minimize_scalar(lambda p: -tail_probability(theta_value, p), bounds=(left, right), method="bounded",
                              options={"xatol": 1e-15})
    candidate_rows.append(dict(theta=theta_value, p_t_at_max=refined.x,
                               max_tail_probability=-refined.fun, cutoff=cutoff))
candidate_out = pd.DataFrame(candidate_rows)
candidate_out.to_csv(OUT / "候选率差的尾部概率.csv", index=False, encoding="utf-8")
print(candidate_out.to_string(index=False))
