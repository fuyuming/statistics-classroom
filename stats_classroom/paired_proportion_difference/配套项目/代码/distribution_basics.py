"""统计课堂27｜基础分布：两人9种变化、20人多项分布与D的精确分布。

完整例题：20名成年人各测前后两次；Z=后测−前测，D为n个人Z的平均数。
给定教学模型q+=(阴转阳)=.05、q0(不变)=.70、q−(阳转阴)=.25；它们与主例观察比例相同，
但不能因此断言是真实总体参数。可在Spyder按代码分节运行，或从项目根/代码目录运行整文件。
"""
# %% 问题：先把两个人的9种有序结果列全
# 每个人的Z只取−1、0、+1；第一人和第二人的顺序不同就是不同的有序组合。
from pathlib import Path
import math
import pandas as pd
from scipy.stats import multinomial

SAVE_RESULTS = False  # 需要导出时直接改为True


def find_project():
    # __file__在Spyder/脚本运行时可用；控制台粘贴本节时，用当前工作目录和其父目录寻找项目根。
    places = [Path.cwd()]
    if "__file__" in globals():
        places.append(Path(__file__).resolve().parent)
    places += [p.parent for p in list(places)]
    for place in places:
        if (place / "数据" / "paired_20.csv").is_file():
            return place
    raise FileNotFoundError("请从配套项目根目录或代码目录启动")


project = find_project()
output_dir = project / "运行结果" / "Python"

# %% 条件：记录n和三类给定概率
n = 20
q_plus = 0.05  # 阴转阳的给定概率
q_zero = 0.70  # 两次结果相同的给定概率
q_minus = 0.25  # 阳转阴的给定概率
prob_by_z = {-1: q_minus, 0: q_zero, 1: q_plus}

# %% 手动计算：两个人9种组合及其概率
# 一行循环固定第一人的Z，再遍历第二人的Z；p1*p2来自不同人的独立性。
two_rows = []
for z1 in (-1, 0, 1):
    for z2 in (-1, 0, 1):
        p1 = prob_by_z[z1]
        p2 = prob_by_z[z2]
        d = (z1 + z2) / 2
        two_rows.append({"z1": z1, "z2": z2, "D": d, "p1": p1, "p2": p2,
                         "probability": p1 * p2})
two_combinations = pd.DataFrame(two_rows)
two_combinations  # 查看9行：D为两人变化相加除以2，probability为p1×p2。
two_combinations  # Spyder中查看9行：每行是一种有序的(Z1,Z2)结果

# %% 调用函数：按D合并两人的9行
# groupby的D是分组键；sum把互斥的有序安排概率相加，不改变原9行表。
two_rate_distribution = two_combinations.groupby("D", as_index=False)["probability"].sum()
two_rate_distribution  # 预期5个D；概率约为0.0625、0.35、0.515、0.07、0.0025

# %% 结果核对：两人9行与5个D
two_probability_total = two_combinations["probability"].sum()
two_rate_probability_total = two_rate_distribution["probability"].sum()
two_rate_distribution
# D从−1到1对应五项概率：.0625、.35、.515、.07、.0025；总和为1。
assert math.isclose(two_probability_total, 1.0)
assert math.isclose(two_rate_probability_total, 1.0)

# %% 手动计算：20人观察组合的安排数与指定安排概率
# 先选nplus人，再从余下的人选nzero人；每一种具体安排的概率是三类概率的幂乘积。
n_plus_obs = 1
n_zero_obs = 14
n_minus_obs = 5
choose_plus = math.comb(n, n_plus_obs)
remaining_after_plus = n - n_plus_obs
choose_zero = math.comb(remaining_after_plus, n_zero_obs)
observed_arrangements = choose_plus * choose_zero
single_arrangement_probability = q_plus ** n_plus_obs * q_zero ** n_zero_obs * q_minus ** n_minus_obs
observed_probability_by_formula = observed_arrangements * single_arrangement_probability
observed_case = pd.DataFrame([{
    "n_plus": n_plus_obs, "n_zero": n_zero_obs, "n_minus": n_minus_obs,
    "arrangements": observed_arrangements,
    "single_arrangement_probability": single_arrangement_probability,
    "probability": observed_probability_by_formula,
}])
observed_case  # arrangements应为232560，probability约为0.07701541

observed_case
# 看arrangements=232560；probability约.07701541，是这一组人数的概率。

# %% 调用函数：用scipy.stats.multinomial.pmf独立核对指定组合
# pmf直接给一组类别人数的概率；它不是上面的手算，而是对公式结果的独立核对。
observed_probability_by_pmf = float(multinomial.pmf(
    [n_plus_obs, n_zero_obs, n_minus_obs], n, [q_plus, q_zero, q_minus]))

# %% 结果核对：手算与pmf
assert math.isclose(observed_probability_by_formula, observed_probability_by_pmf,
                    rel_tol=1e-12, abs_tol=1e-15)

# %% 手动计算：双循环列出231组三类人数
# 外层固定nplus，内层固定nzero；nminus由n−nplus−nzero确定，因此不会漏项或重复。
rows = []
for n_plus in range(n + 1):
    for n_zero in range(n - n_plus + 1):
        n_minus = n - n_plus - n_zero
        arrangements = math.comb(n, n_plus) * math.comb(n - n_plus, n_zero)
        one_arrangement = q_plus ** n_plus * q_zero ** n_zero * q_minus ** n_minus
        formula_probability = arrangements * one_arrangement
        rows.append({"n_plus": n_plus, "n_zero": n_zero, "n_minus": n_minus,
                     "arrangements": arrangements,
                     "single_arrangement_probability": one_arrangement,
                     "probability": formula_probability,
                     "D": (n_plus - n_minus) / n})
all_counts = pd.DataFrame(rows)
all_counts  # 预期231行，每行是一组(n_plus,n_zero,n_minus)

# %% 调用函数：逐行用multinomial.pmf核对231组
# pmf的x、n、p分别是三类人数、总人数、三类概率，顺序始终为+1、0、−1。
pmf_probabilities = []
for row in all_counts.itertuples(index=False):
    # 一次循环只处理一组人数；人数与概率向量都按+1、0、−1排列。
    counts = [row.n_plus, row.n_zero, row.n_minus]
    probability_from_function = multinomial.pmf(counts, n, [q_plus, q_zero, q_minus])
    pmf_probabilities.append(float(probability_from_function))
all_counts["pmf_probability"] = pmf_probabilities

# %% 结果核对：231组、相同D归类与理论矩
assert len(all_counts) == 231
assert (abs(all_counts["probability"] - all_counts["pmf_probability"]) < 1e-12).all()
rate_difference_distribution = all_counts.groupby("D", as_index=False)["probability"].sum()
same_difference_9_groups = all_counts[all_counts["n_plus"] - all_counts["n_minus"] == -4].copy()
mean_D = float((all_counts["D"] * all_counts["probability"]).sum())
var_D = float(((all_counts["D"] - mean_D) ** 2 * all_counts["probability"]).sum())
summary = pd.DataFrame([{
    "n": n, "q_plus": q_plus, "q_zero": q_zero, "q_minus": q_minus,
    "combinations": len(all_counts), "all_probability": all_counts["probability"].sum(),
    "mean_D": mean_D, "var_D": var_D,
    "observed_arrangements": observed_arrangements,
    "observed_probability": observed_probability_by_formula,
}])
assert math.isclose(all_counts["probability"].sum(), 1.0)
assert len(same_difference_9_groups) == 9

# %% 可选保存：需要保存时把本文件顶部SAVE_RESULTS改为True
if SAVE_RESULTS:
    output_dir.mkdir(parents=True, exist_ok=True)
    two_combinations.to_csv(output_dir / "基础分布_两人9种组合.csv", index=False, encoding="utf-8")
    two_rate_distribution.to_csv(output_dir / "基础分布_两人率差概率.csv", index=False, encoding="utf-8")
    all_counts.drop(columns="pmf_probability").to_csv(output_dir / "基础分布_231组概率.csv", index=False, encoding="utf-8")
    rate_difference_distribution.to_csv(output_dir / "基础分布_率差概率.csv", index=False, encoding="utf-8")
    same_difference_9_groups.drop(columns="pmf_probability").to_csv(output_dir / "基础分布_同率差9组概率.csv", index=False, encoding="utf-8")
    summary.to_csv(output_dir / "基础分布_汇总.csv", index=False, encoding="utf-8")

# %% 思考：改q_plus或q_minus后，先预测均值q_plus−q_minus，再运行本节看D的分布如何移动。
