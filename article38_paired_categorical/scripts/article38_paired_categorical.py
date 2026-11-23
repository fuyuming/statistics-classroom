# %% 统计课堂38：完整教学项目（UTF-8）
# 两例均为n=8的假设教学计数：例1为三法二分类，例2为前后多分类。
# 对应问题：例1问三种检测的总体阳性率是否相同及配对差异；例2问前后边际是否同质、反向路径是否对称。
# 首次运行默认只显示；改SAVE_RESULTS=True可保存结果。依赖见requirements.txt。
from pathlib import Path
import os
import itertools
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.contingency_tables import cochrans_q, SquareTable
from statsmodels.stats.multitest import multipletests

SAVE_RESULTS = os.getenv("ARTICLE38_SAVE", "0") == "1"
ROOT = Path.cwd()
# IDE没有__file__时，只要求工作目录为项目根目录；完整运行也可从脚本定位。
if not (ROOT / "data").is_dir() and "__file__" in globals():
    ROOT = Path(__file__).resolve().parents[1]
assert (ROOT / "data").is_dir(), "请先把工作目录设为本项目根目录。"
records = pd.read_csv(ROOT / "data/binary_records.csv")
assert records.columns.tolist() == ["id", "A", "B", "C"]
assert records["id"].is_unique
assert not records.isna().any().any() and records[["A", "B", "C"]].isin([0, 1]).all().all()
# id标识受检者；Y只取三种方案的0/1结果，每人一行。
Y = records[["A", "B", "C"]].to_numpy(dtype=int)
# groupby将相同的三法组合合并，size统计该组合的人数。
patterns = (records.groupby(["A", "B", "C"], sort=False).size().reset_index(name="count"))
patterns["pattern"] = patterns[["A", "B", "C"]].astype(str).agg("".join, axis=1)
patterns = patterns[["pattern", "A", "B", "C", "count"]].sort_values("pattern").reset_index(drop=True)
transition_records = pd.read_csv(ROOT / "data/transition_records.csv")
assert transition_records.columns.tolist() == ["id", "before", "after"]
assert transition_records["id"].is_unique and len(transition_records) == len(records)
assert transition_records[["before", "after"]].isin(labels if "labels" in globals() else ["A", "B", "C"]).all().all()
transitions = pd.read_csv(ROOT / "data/transition_counts.csv")
print("逐人记录前10条（id及A/B/C）：\n", records.head(10).to_string(index=False))
print("由逐人记录重汇总的三检测模式：\n", patterns)
print("前后转移：\n", transitions)
labels = ["A", "B", "C"]
# Y只取三结果列；每人一行、每法8条对应记录。
M = (transition_records.assign(count=1)
     .pivot_table(index="before", columns="after", values="count", aggfunc="sum", fill_value=0)
     .reindex(index=labels, columns=labels, fill_value=0).to_numpy())
assert np.isin(Y, [0, 1]).all()
print("Y前六行（只保留三法结果）：\n", Y[:6])
print("行=前，列=后：\n", M)

# %% 本文自编辅助函数：只整理计算，不改变统计目标
# 普通辅助函数定义在所有调用之前，分节请按顺序运行。
def paired_summary(b, c_count, n):
    discordant = b + c_count
    p_raw = 1.0
    if discordant > 0:
        p_raw = stats.binomtest(c_count, discordant, p=0.5).pvalue
    rd = (c_count - b) / n
    return dict(b=b, c=c_count, rate_diff=rd, p_raw=p_raw)


def sm_formula(tab):
    row_total = tab.sum(axis=1)
    col_total = tab.sum(axis=0)
    difference = row_total - col_total
    covariance = -(tab + tab.T).astype(float)
    np.fill_diagonal(covariance, row_total + col_total - 2*np.diag(tab))
    # 去掉最后一类；所有差的和为0。solve解线性方程，避免手写矩阵求逆。
    reduced = covariance[:-1, :-1]
    d = difference[:-1]
    return float(d @ np.linalg.solve(reduced, d))


def bowker_formula(tab):
    index = np.triu_indices(len(tab), 1)
    totals = (tab + tab.T)[index]
    difference = (tab - tab.T)[index]
    assert (totals > 0).all(), "本文函数只实现非零反向格对。"
    return float(np.sum(difference**2 / totals))


def cochran_q_mc(Y, B, seed):
    """完整配对Y的MC：H0下每行随机重排，保持该人的行和；argsort选择随机位置。
    统计Q右尾超越次数，返回p_mc、MCSE及seed，便于单独复现调用。
    """
    rng = np.random.default_rng(seed)
    n, k = Y.shape
    row_totals = Y.sum(axis=1)
    total = row_totals.sum()
    D = np.sum(row_totals * (k - row_totals))
    observed = k * (k - 1) * np.sum((Y.sum(axis=0) - total/k)**2) / D
    sim_columns = np.zeros((B, k), dtype=int)
    for row_sum in row_totals.astype(int):
        scores = rng.random((B, k))
        sim_columns += np.argsort(scores, axis=1) < row_sum
    sim_statistic = k * (k - 1) * np.sum((sim_columns - total/k)**2, axis=1) / D
    exceedances = int(np.sum(sim_statistic >= observed - 1e-12))
    p_mc = exceedances / B
    return {"test": "Cochran_Q", "B": B, "seed": seed,
            "exceedances": exceedances, "p_mc": p_mc,
            "mcse": np.sqrt(p_mc * (1 - p_mc) / B)}


def bowker_mc(M, B, seed):
    """完整转移表M的MC：H0下固定反向格对总数，方向为独立Binomial(m,.5)。
    统计Bowker右尾超越次数，返回p_mc、MCSE及seed，便于单独复现调用。
    """
    rng = np.random.default_rng(seed)
    assert M.shape == (3, 3)
    index = np.triu_indices(3, 1)
    totals = (M + M.T)[index]
    assert np.all(totals > 0)
    observed = np.sum((M - M.T)[index]**2 / totals)
    sim_counts = rng.binomial(totals.astype(int), 0.5, size=(B, len(totals)))
    sim_statistic = np.sum((2 * sim_counts - totals)**2 / totals, axis=1)
    exceedances = int(np.sum(sim_statistic >= observed - 1e-12))
    p_mc = exceedances / B
    return {"test": "Bowker", "B": B, "seed": seed,
            "exceedances": exceedances, "p_mc": p_mc,
            "mcse": np.sqrt(p_mc * (1 - p_mc) / B)}


# %% 1 Cochran Q
n, k = Y.shape
columns = Y.sum(axis=0)
rows = Y.sum(axis=1)
total = rows.sum()
D = np.sum(rows * (k - rows))
Q = k * (k - 1) * np.sum((columns - total/k)**2) / D
q_native = cochrans_q(Y)                    # statsmodels原生Q，行=区组
print(f"Q节教学输入：n={n}，k={k}；列和C={columns.tolist()}；T={total}，D={D}；"
      f"df={k-1}；渐近P={stats.chi2.sf(Q, k-1):.12g}")
q_function_comparison = pd.DataFrame({
    "test": ["Cochran_Q"], "formula": [Q],
    "function": [q_native.statistic],
    "abs_error": [abs(Q-q_native.statistic)]})
print(q_function_comparison.to_string(index=False))
assert abs(q_native.statistic-Q) < 1e-10
print("Q：", Q)
print("Q回答总体阳性概率是否相同，不评价哪法准确。")

# %% 2 给定每人行和的Cochran Q条件精确枚举
# 本节专用于三法：U数唯一阳性位置，W数唯一阴性位置。
# 两组对象独立；多项配置的概率相乘，再按Q右尾累加。
n1 = int(np.sum(rows == 1))
n2 = int(np.sum(rows == 2))
n3 = int(np.sum(rows == 3))
observed_sum_squares = int(np.sum(columns ** 2))
probability_sum = 0.0
p_exact = 0.0
configurations = 0
for u1 in range(n1 + 1):
    for u2 in range(n1 - u1 + 1):
        U = np.array([u1, u2, n1 - u1 - u2])
        probability_U = stats.multinomial.pmf(U, n1, [1 / 3] * k)
        for v1 in range(n2 + 1):
            for v2 in range(n2 - v1 + 1):
                W = np.array([v1, v2, n2 - v1 - v2])
                probability_W = stats.multinomial.pmf(W, n2, [1 / 3] * k)
                conditional_columns = n3 + n2 + U - W
                probability = probability_U * probability_W
                probability_sum += probability
                configurations += 1
                if int(np.sum(conditional_columns ** 2)) >= observed_sum_squares:
                    p_exact += probability
q_conditional = pd.DataFrame([dict(configurations=configurations,
                                   person_assignments=3**n1,
                                   probability_sum=probability_sum,
                                   p_exact=p_exact)])
assert abs(probability_sum - 1) < 1e-10
print(q_conditional.to_string(index=False))
print("条件于每人行和；U/W分别是行和1/2对象的位置计数。")

# %% 3 三对二分类：条件精确McNemar、Holm、配对率差
pairs = list(itertools.combinations(range(k), 2))
pair_rows = []
for first, second in pairs:
    x = Y[:, first]
    y = Y[:, second]
    b = int(np.sum((x == 1) & (y == 0)))
    c_count = int(np.sum((x == 0) & (y == 1)))
    answer = paired_summary(b, c_count, n)
    answer["comparison"] = labels[second] + "-" + labels[first]
    pair_rows.append(answer)
pairwise = pd.DataFrame(pair_rows)
# 返回对象第2项(pvals_corrected，索引1)是调整P；method指定Holm。
pairwise["holm"] = multipletests(pairwise.p_raw, method="holm")[1]
pairwise = pairwise[["comparison", "b", "c", "rate_diff", "p_raw", "holm"]]
print(pairwise.to_string(index=False))
pairwise_teaching_check = pairwise[["comparison", "b", "c"]].copy()
pairwise_teaching_check["m"] = pairwise_teaching_check["b"] + pairwise_teaching_check["c"]
print(pairwise_teaching_check.to_string(index=False))
assert pairwise_teaching_check["m"].sum() == D
print(f"例1三对不一致数m=b+c之和={pairwise_teaching_check['m'].sum()}，与Q节D相同。")
print("rate_diff乘100为百分点；两例均以8人为分母。")

# %% 4 同一方表：SM问边际，Bowker问对称
SM = sm_formula(M)
Bowker = bowker_formula(M)
# shift_zeros=False禁止软件替零格加0.5；本例无零格。
square = SquareTable(M, shift_zeros=False)
sm_native = square.homogeneity(method="stuart_maxwell")
bow_native = square.symmetry(method="bowker")
assert abs(sm_native.statistic-SM) < 1e-10
assert abs(bow_native.statistic-Bowker) < 1e-10
sm_function_comparison = pd.DataFrame({
    "test": ["Stuart_Maxwell", "Bowker"],
    "formula": [SM, Bowker],
    "function": [sm_native.statistic, bow_native.statistic],
    "abs_error": [abs(SM-sm_native.statistic), abs(Bowker-bow_native.statistic)]})
print(sm_function_comparison.to_string(index=False))
key_results = pd.DataFrame({
    "method": ["Cochran_Q", "Stuart_Maxwell", "Bowker"],
    "statistic": [Q, SM, Bowker],
    "df": [k-1, len(M)-1, len(M)*(len(M)-1)//2]})
# sf是右尾概率，而不是把原始分类资料当正态。
key_results["p_value"] = stats.chi2.sf(key_results.statistic, key_results.df)
print(key_results.to_string(index=False))

# %% 5 逐类别“此类/其他”：分母为8人
category_rows = []
for i, label in enumerate(labels):
    leave = int(M[i, :].sum() - M[i, i])
    enter = int(M[:, i].sum() - M[i, i])
    answer = paired_summary(leave, enter, int(M.sum()))
    answer["category"] = label
    category_rows.append(answer)
categories = pd.DataFrame(category_rows)
categories["holm"] = multipletests(categories.p_raw, method="holm")[1]
print(categories.to_string(index=False))
print("本节比较各类别后减前的边际比例，不比较单条路径。")

# %% 6 反向路径：精确P、Holm及条件方向概率CP区间
path_rows = []
for first, second in pairs:
    forward = int(M[first, second])
    reverse = int(M[second, first])
    fit = stats.binomtest(forward, forward+reverse, p=0.5)
    ci = fit.proportion_ci(confidence_level=0.95, method="exact")
    path_rows.append(dict(path=labels[first]+"->"+labels[second],
                          forward=forward, reverse=reverse, total=forward+reverse,
                          theta=forward/(forward+reverse), p_raw=fit.pvalue,
                          lo95=ci.low, hi95=ci.high))
paths = pd.DataFrame(path_rows)
paths["holm"] = multipletests(paths.p_raw, method="holm")[1]
print(paths.to_string(index=False))
print("theta与CP区间以该路径互转者为分母；区间为逐项95%。")

# %% 7 自编条件枚举：Bowker统计量的右尾排序
sizes = paths.total.to_numpy(dtype=int)
probability_sum = 0.0
p_exact = 0.0
configurations = 0
# 本例三条路径；乘积枚举每条路径从0到total的正向人数。
possible_counts = [range(size+1) for size in sizes]
for counts in itertools.product(*possible_counts):
    counts = np.asarray(counts)
    statistic = np.sum((2*counts-sizes)**2/sizes)
    probability = np.prod(stats.binom.pmf(counts, sizes, p=0.5))
    probability_sum += probability
    configurations += 1
    if statistic >= Bowker-1e-12:
        p_exact += probability
enumeration = pd.DataFrame([dict(configurations=configurations,
                                  probability_sum=probability_sum, p_exact=p_exact)])
assert abs(probability_sum-1) < 1e-10
print(enumeration.to_string(index=False))
print("自编统计量右尾精确枚举；无蒙特卡罗。")

# %% 8 条件蒙特卡罗：固定行和/固定反向路径总数，在H0下模拟
# B=100000、seed=20261123；两个helper各自独立设seed。
# p_mc估计条件精确尾概率，mcse为Monte Carlo标准误，不是另一χ²近似；不模拟SM。
B = 100000
seed = 20261123
q_mc = cochran_q_mc(Y, B, seed)
bowker_mc_fit = bowker_mc(M, B, seed)
monte_carlo = pd.DataFrame({
    "test": [q_mc["test"], bowker_mc_fit["test"]], "B": B, "seed": seed,
    "exceedances": [q_mc["exceedances"], bowker_mc_fit["exceedances"]],
    "p_mc": [q_mc["p_mc"], bowker_mc_fit["p_mc"]],
    "mcse": [q_mc["mcse"], bowker_mc_fit["mcse"]],
    "exact_reference": [q_conditional.p_exact.iloc[0], enumeration.p_exact.iloc[0]]})
print(monte_carlo.to_string(index=False))

# %% 9 循环转移练习：另一组8人的M2
M2 = np.array([[1, 2, 0], [0, 1, 2], [2, 0, 0]], dtype=float)
exercise = pd.DataFrame([dict(SM=sm_formula(M2), Bowker=bowker_formula(M2))])
print("改动后的行和/列和：", M2.sum(axis=1), M2.sum(axis=0))
print(exercise.to_string(index=False))

# %% 10 可选保存：三语输出隔离
if SAVE_RESULTS:
    output_dir = ROOT / "results/Python"
    output_dir.mkdir(parents=True, exist_ok=True)
    result_tables = dict(key_results=key_results, q_conditional=q_conditional,
                         pairwise=pairwise,
                         categories=categories, paths=paths,
                         enumeration=enumeration, exercise=exercise,
                         monte_carlo=monte_carlo)
    for name, table in result_tables.items():
        table.to_csv(output_dir / (name+".csv"), index=False)
