"""统计课堂37：完整构成、指定类别比例、稀疏表检验（第11课教学例题计数，非真实调查）。
faults是每组80台已失效记录器的故障构成，sparse是每组10台的合格/不合格计数。
依赖NumPy/pandas/SciPy/statsmodels；UTF-8。Spyder/VS Code从根目录按九节运行，默认显示；
设置ARTICLE37_SAVE_RESULTS=1后第9节保存四张表到results/Python。"""
# %% 1. 数据与参数
from pathlib import Path
import itertools
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

ROOT = Path.cwd()
if not (ROOT / "data/faults.csv").is_file():
    raise FileNotFoundError("请先把工作目录设为含data文件夹的本项目根目录。")
sys.path.insert(0, str(ROOT / "helpers"))
sys.dont_write_bytecode = True
from intervals import cramers_v_ci, rate_difference_ci  # 本项目自编函数

faults = pd.read_csv(ROOT / "data/faults.csv")
sparse = pd.read_csv(ROOT / "data/sparse.csv")
# 先看字段和全部三行：faults每行一组、分母为80台已失效记录器；sparse每行一组、分母为10台。
print("faults字段：", faults.columns.tolist())
print(faults.to_string(index=False))
print("sparse字段：", sparse.columns.tolist())
print(sparse.to_string(index=False))
# 后续矩阵保留计数列，分别回答完整构成、单类比例和稀疏表检验问题。
rows = faults["row"].tolist()
cols = ["open", "electrode", "board"]
tab = faults[cols].to_numpy()
sparse_tab = sparse[["pass", "fail"]].to_numpy()
conf_level = 0.95
alpha = 0.05
# Python索引从0开始：0=F，1=G，2=H；每个元组为一个组对。
pair_index = list(itertools.combinations(range(len(rows)), 2))

# %% 2. 整体检验与关联效应量
# 问题：三组的完整故障构成是否相同？同时报告整体Q、P和V。
# correction=False明确关闭连续性校正；返回统计量、P、自由度、预期矩阵。
q, p, df, expected = stats.chi2_contingency(tab, correction=False)
# 自编cramers_v_ci(table, conf_level)返回(V, V_low, V_high)。
v, v_low, v_high = cramers_v_ci(tab, conf_level)
main_overall = pd.DataFrame([{
    "q": q, "df": df, "p": p, "n": tab.sum(),
    "V": v, "V_low": v_low, "V_high": v_high,
}])
print(main_overall.to_string(index=False))
print("Q/P回答整体差异，V及区间概括关联强度；这里尚未进行组对Holm调整。")
# 实跑：Q=44.79167，P=4.392904e-09，V=0.3054766。

# %% 3. 完整分类构成的三项组对
# 思路：每次保留三类故障，重算2×3子表，再对三对原始P作Holm调整。
full_records = []
for i, j in pair_index:
    sub_tab = tab[[i, j], :]  # 两行、全部三列
    q, p, df, sub_expected = stats.chi2_contingency(sub_tab, correction=False)
    # 自编函数返回V及区间三元组；conf_level只控制区间水平。
    v, v_low, v_high = cramers_v_ci(sub_tab, conf_level)
    full_records.append({
        "pair": f"{rows[i]}-{rows[j]}", "q": q, "df": df, "p_raw": p,
        "min_expected": sub_expected.min(), "V": v, "V_low": v_low, "V_high": v_high,
    })
full_pairs = pd.DataFrame(full_records)
# multipletests返回多个对象，第[1]项才是调整后的P值向量。
full_pairs["p_holm"] = multipletests(full_pairs["p_raw"], method="holm")[1]
# loc先按条件取行、再取pair列；tolist把组名转成列表，供下方if判断。
detected_full = full_pairs.loc[full_pairs["p_holm"] < alpha, "pair"].tolist()
print(f"完整构成比较（alpha={alpha:g}）：")
if detected_full:
    print("检出差异的组对：", ", ".join(detected_full))
else:
    print("无检出差异的组对。")
full_pairs = full_pairs[[
    "pair", "q", "df", "p_raw", "p_holm", "min_expected", "V", "V_low", "V_high",
]]
print(full_pairs.to_string(index=False))
# p_raw是单个2×3检验，p_holm是F-G、F-H、G-H这一比较族的调整P。
# 实跑：三个Holm P约0.002365、4.298e-09、0.004714。

# %% 4. 单类比例：电极脱落与其他
# 问题：只比较电极脱落比例时，把其余两类合并为other。
electrode = tab[:, 1]  # 第2列，索引1
# axis=1按行求和：每组总数减电极脱落数得到other。
electrode_tab = np.column_stack([electrode, tab.sum(axis=1) - electrode])
category_records = []
for i, j in pair_index:
    sub_tab = electrode_tab[[i, j], :]
    # correction=True在这个2×2表中使用Yates连续性校正。
    q, p, df, _ = stats.chi2_contingency(sub_tab, correction=True)
    # 自编函数按(success,total,conf_level)返回p1,p2,diff,low,high。
    p1, p2, difference, low, high = rate_difference_ci(
        sub_tab[:, 0], sub_tab.sum(axis=1), conf_level
    )
    category_records.append({
        "pair": f"{rows[i]}-{rows[j]}", "q_yates": q, "df": df, "p_raw": p,
        "p1": p1, "p2": p2, "diff": difference, "diff_low": low, "diff_high": high,
    })
category_pairs = pd.DataFrame(category_records)
category_pairs["p_holm"] = multipletests(category_pairs["p_raw"], method="holm")[1]
detected_category = category_pairs.loc[category_pairs["p_holm"] < alpha, "pair"].tolist()
print(f"电极脱落比例比较（alpha={alpha:g}）：")
if detected_category:
    print("检出差异的组对：", ", ".join(detected_category))
else:
    print("无检出差异的组对。")
category_pairs = category_pairs[[
    "pair", "q_yates", "df", "p_raw", "p_holm", "p1", "p2", "diff", "diff_low", "diff_high",
]]
print(category_pairs.to_string(index=False))
print("p_holm仍按三对调整；p1/p2为前后组比例，diff及区间均为前组减后组。")
# 实跑：Holm P=1、1、0.4249341；G-H差0.125，单项CI约[-0.03627,0.28627]。

# %% 5. 稀疏表的Fisher组对
# 思路：每对保留合格/不合格两列；小预期频数时用双侧Fisher，再调整三对P。
sparse_records = []
for i, j in pair_index:
    sub_tab = sparse_tab[[i, j], :]
    # fisher_exact返回statistic,pvalue；conditional OR固定行列边际。
    exact = stats.fisher_exact(sub_tab, alternative="two-sided")
    or_result = stats.contingency.odds_ratio(sub_tab, kind="conditional")
    # kind="conditional"固定行列边际；confidence_interval返回low/high。
    or_ci = or_result.confidence_interval(confidence_level=conf_level)
    # 精确P值不依赖卡方近似；另算预期频数供检查子表。
    sub_expected = np.outer(sub_tab.sum(axis=1), sub_tab.sum(axis=0)) / sub_tab.sum()
    sparse_records.append({
        "pair": f"{sparse['row'][i]}-{sparse['row'][j]}",
        "p_raw": exact.pvalue, "min_expected": sub_expected.min(),
        "OR": or_result.statistic, "OR_low": or_ci.low, "OR_high": or_ci.high,
    })
sparse_pairs = pd.DataFrame(sparse_records)
sparse_pairs["p_holm"] = multipletests(sparse_pairs["p_raw"], method="holm")[1]
detected_sparse = sparse_pairs.loc[sparse_pairs["p_holm"] < alpha, "pair"].tolist()
print(f"小样本合格率比较（alpha={alpha:g}）：")
if detected_sparse:
    print("检出差异的组对：", ", ".join(detected_sparse))
else:
    print("无检出差异的组对。")
sparse_pairs = sparse_pairs[["pair", "p_raw", "p_holm", "min_expected",
                             "OR", "OR_low", "OR_high"]]
print(sparse_pairs.to_string(index=False))
# p_raw来自精确检验，p_holm是三对合格率比较族的调整P；min_expected只作稀疏性检查。
# 实跑：Holm P=1、0.05929983、0.1395570。

# %% 6. 两类问题的调整P并列
# 问题：同一组对在完整构成与指定类别比例问题下，结论是否相同？
comparison = pd.DataFrame({
    "pair": full_pairs["pair"],
    "full_p_holm": full_pairs["p_holm"],
    "category_p_holm": category_pairs["p_holm"],
})
print(comparison.to_string(index=False))
print("两列Holm P分别属于完整构成三对族与指定类别比例三对族。")

# %% 7. 独立检验练习：Hochberg
# 先用两项独立检验的模型演算Hochberg的FWER：各项H0为真时P~Uniform(0,1)。
alpha_proof = 0.05
half_alpha = alpha_proof / 2
prob_any_small = 1 - (1 - half_alpha) ** 2
prob_both_middle = (alpha_proof - half_alpha) ** 2
fwer_two_null = prob_any_small + prob_both_middle
event_table = pd.DataFrame({
    "event": ["至少一个P≤alpha/2", "两P都在(alpha/2, alpha]", "两项H0真时FWER"],
    "probability": [prob_any_small, prob_both_middle, fwer_two_null],
})
print(event_table.to_string(index=False))
print("这是模型的理论误报概率，不是模拟结果或某一次观测的P值。")
print("若只有1项H0为真，误拒绝该项必须满足P≤alpha，因此上界为alpha。")

# 这是单独给定的3项相互独立检验，不来自前面F-G、F-H、G-H共享数据的原例。
# 原始P为test1/test2/test3=.020/.023/.040，按FWER=.05作判断。
p_demo = np.array([0.020, 0.023, 0.040])
practice_labels = ["test1", "test2", "test3"]
practice_holm = multipletests(p_demo, method="holm")[1]
# statsmodels的simes-hochberg对应R的p.adjust(method="hochberg").
practice_hochberg = multipletests(p_demo, method="simes-hochberg")[1]
practice_table = pd.DataFrame({
    "test": practice_labels,
    "raw": p_demo,
    "holm": practice_holm,
    "hochberg": practice_hochberg,
})
print(practice_table.to_string(index=False))
practice_alpha = 0.05
practice_reject_holm = practice_table.loc[practice_table["holm"] <= practice_alpha, "test"].tolist()
practice_reject_hochberg = practice_table.loc[practice_table["hochberg"] <= practice_alpha, "test"].tolist()
print("test1/test2/test3分别对应检验甲/乙/丙。")
print("Hochberg目标是控制FWER；适用于相互独立或符合特定正依赖条件的检验。")
print("本练习采用相互独立前提，原三组对分析继续使用Holm。")
print(f"Holm拒绝项数：{len(practice_reject_holm)}；Hochberg拒绝项数：{len(practice_reject_hochberg)}。")

# 对照：最大P首关通过时全部拒绝，但各项调整P仍可不同。
p_compare = np.array([0.001, 0.020, 0.040])
compare_hochberg = multipletests(p_compare, method="simes-hochberg")[1]
compare_table = pd.DataFrame({
    "test": ["test1", "test2", "test3"],
    "raw": p_compare,
    "hochberg": compare_hochberg,
})
compare_reject_hochberg = compare_table.loc[
    compare_table["hochberg"] <= practice_alpha, "test"
].tolist()
print(compare_table.to_string(index=False))
if p_compare.max() <= practice_alpha:
    print(f"最大P={p_compare.max():g}≤{practice_alpha:g}，Hochberg首关通过，拒绝全部三项。")
else:
    print(f"最大P={p_compare.max():g}>{practice_alpha:g}，首关未通过，继续按Hochberg规则判断。")
print("各项调整P由整组P从后往前取累计最小计算，不一定相同。")
print(f"对照例Hochberg拒绝项数：{len(compare_reject_hochberg)}。")

# %% 8. G-H效应量与区间报告
# 只报告一个示例组对，区间为单项区间，P值仍按三对家族调整。
gh = category_pairs[category_pairs["pair"] == "G-H"].iloc[0]
print(
    f"G-H电极脱落比例差={gh['diff']:.6f}（前组减后组），"
    f"单项{conf_level * 100:g}% CI [{gh['diff_low']:.6f}, {gh['diff_high']:.6f}]，"
    f"Holm P={gh['p_holm']:.6f}。"
)
print("区间没有作三项同时覆盖调整；Holm只用于同一问题的三对P值。")

# %% 9. 可选保存
# 默认不写文件。需保存时先运行 os.environ["ARTICLE37_SAVE_RESULTS"] = "1"，再运行本节。
save_results = os.environ.get("ARTICLE37_SAVE_RESULTS") == "1"
if save_results:
    output_dir = ROOT / "results/Python"
    output_dir.mkdir(parents=True, exist_ok=True)
    results = {
        "main_overall": main_overall, "full_pairs": full_pairs,
        "category_pairs": category_pairs, "sparse_pairs": sparse_pairs,
    }
    for name, result in results.items():
        result.to_csv(output_dir / f"{name}.csv", index=False)
    print("已保存四张表到", output_dir)
