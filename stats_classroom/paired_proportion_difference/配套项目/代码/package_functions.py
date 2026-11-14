"""统计课堂27｜用现成多项分布函数核对三组已知人数概率。"""
# %% 问题与条件：pmf回答“一组人数出现的概率是多少”
from pathlib import Path
import math
import pandas as pd
from scipy.stats import multinomial

SAVE_RESULTS = False  # 需要导出时直接改为True
q_plus = 0.05  # 阴转阳的给定概率
q_zero = 0.70  # 两次结果相同的给定概率
q_minus = 0.25  # 阳转阴的给定概率
# 每个括号依次是n、阴转阳人数、不变人数、阳转阴人数。
cases = [(2, 0, 2, 0), (2, 1, 0, 1), (20, 1, 14, 5)]

# %% 手动计算：组合数乘以一种指定安排的概率
# 这一步让学生看到pmf背后的公式，不把pmf输出冒充手算结果。
manual_rows = []
for n, n_plus, n_zero, n_minus in cases:
    arrangements = math.factorial(n) / (math.factorial(n_plus) * math.factorial(n_zero) * math.factorial(n_minus))
    one_arrangement = q_plus ** n_plus * q_zero ** n_zero * q_minus ** n_minus
    manual_rows.append({"n": n, "n_plus": n_plus, "n_zero": n_zero, "n_minus": n_minus,
                        "arrangements": arrangements,
                        "single_arrangement_probability": one_arrangement,
                        "manual_probability": arrangements * one_arrangement})
manual_results = pd.DataFrame(manual_rows)
manual_results  # 每行先看安排数，再看一种安排概率，最后看两者相乘。
manual_results  # 先看安排数、单一安排概率和手动结果

# %% 调用函数：scipy.stats.multinomial.pmf
# x是三类人数，n是总人数，p是按+1、0、−1排列的概率；pmf返回一个数。
pmf_rows = []
for row in manual_results.itertuples():
    probability = float(multinomial.pmf([row.n_plus, row.n_zero, row.n_minus],
                                         row.n, [q_plus, q_zero, q_minus]))
    pmf_rows.append({"model": f"n{row.n}_{row.n_plus}{row.n_zero}{row.n_minus}",
                     "n": row.n, "n_plus": row.n_plus, "n_zero": row.n_zero,
                     "n_minus": row.n_minus, "probability": probability})
package_results = pd.DataFrame(pmf_rows)
package_results  # 预期概率约为0.49、0.025、0.07701541

# %% 结果核对：函数与手动公式逐行一致
assert (abs(manual_results["manual_probability"] - package_results["probability"]) < 1e-12).all()

package_results
# 三行probability依次为.49、.025、约.07701541，分别对应上面三组人数。

# %% 可选保存：改SAVE_RESULTS=True后写一张表
if SAVE_RESULTS:
    places = [Path.cwd()]
    if "__file__" in globals():
        places.append(Path(__file__).resolve().parent)
    project = next(p for p in places + [p.parent for p in places]
                   if (p / "数据" / "paired_20.csv").is_file())
    out = project / "运行结果" / "Python"
    out.mkdir(parents=True, exist_ok=True)
    package_results.to_csv(out / "package_functions_probabilities.csv", index=False)

# %% 思考：把n=20的(1,14,5)改成另一组人数，先手算D=(n_plus−n_minus)/n，再核对pmf。
