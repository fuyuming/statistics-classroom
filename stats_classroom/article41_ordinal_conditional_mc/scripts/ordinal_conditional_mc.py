# 统计课堂41：双向有序表的条件模拟（按节运行，UTF-8）
# 目标：固定行列边际的完整枚举与Monte Carlo，理解双向有序表。
# 数据：医学统计学第5版p116–117例7-9及另设8例假设小例；依赖见requirements.txt。
# 默认显示摘要、首轮模拟表和图；ARTICLE41_SAVE_RESULTS=1才写results/figures。
# ARTICLE41_SHOW_DETAILS=1才显示完整频数表、逐格演算和尾部归并。

# %% 01.1 读取数据与集中参数
from pathlib import Path
import os, sys, math
import numpy as np
import pandas as pd
from fractions import Fraction

candidates = []
if "__file__" in globals():
    script_dir = Path(__file__).resolve().parent
    candidates.extend([script_dir.parent, script_dir])
cwd = Path.cwd()
candidates.extend([cwd, cwd.parent])
ROOT = next((candidate for candidate in candidates
             if (candidate / "data/ordinal_tables.csv").is_file()), None)
if ROOT is None:
    raise FileNotFoundError("未找到data/ordinal_tables.csv；请将工作目录设为项目根目录或scripts目录。")
sys.path.insert(0, str(ROOT / "scripts"))
from helpers.ordinal_helpers import ordinal_stats, formula_steps, enumerate_tables, table_probability, conditional_mc
from plot_ordinal_conditional_mc import plot_results
raw = pd.read_csv(ROOT / "data/ordinal_tables.csv")
ROW_SCORES = np.arange(1, 5)
COL_SCORES = np.arange(1, 5)
B = 19999
SEED = 20261122
SAVE_RESULTS = os.environ.get("ARTICLE41_SAVE_RESULTS") == "1"
SHOW_DETAILS = os.environ.get("ARTICLE41_SHOW_DETAILS") == "1"
print(f"参数：B={B}, SEED={SEED}, SAVE_RESULTS={int(SAVE_RESULTS)}, SHOW_DETAILS={int(SHOW_DETAILS)}")
print("输入字段示例：")
print(raw.head(3).to_string(index=False))

def get_table(example):
    part = raw.loc[raw.example.eq(example)].sort_values(["row_score", "col_score"])
    return part["count"].to_numpy(int).reshape(4, 4)

main = get_table("main")
mini = get_table("mini")
if SHOW_DETAILS:
    print("主例频数表：")
    print(main)
    print("小例频数表：")
    print(mini)

# %% 01.2 观察摘要：chi2_p为渐近P，不是条件尾部P
observed = {}
for name in ("main", "mini"):
    observed[name] = ordinal_stats(get_table(name), ROW_SCORES, COL_SCORES)
summary_rows = []
for name in ("main", "mini"):
    row = {"example": name}
    for field in ("N", "T", "r", "Q", "chi2_p"):
        row[field] = observed[name][field]
    summary_rows.append(row)
summary = pd.DataFrame(summary_rows)
print("观察摘要（chi2_p为χ²1渐近P；不是完整条件尾部P）")
print(summary.to_string(index=False))
for ans in observed.values():
    assert abs(ans["T"] ** 2 / ans["VarT"] - ans["Q"]) < 1e-10

# %% 01.3 演算中间量（默认只计算并保存；SHOW_DETAILS=1才打印）
formula_rows = []
for name in ("main", "mini"):
    steps = formula_steps(
        get_table(name), observed[name], name,
        ROW_SCORES, COL_SCORES, SHOW_DETAILS,
    )
    formula_rows.extend(steps)

# %% 01.4 枚举小例的全部固定边际表（282张）
row_sums = mini.sum(axis=1).tolist()
col_sums = mini.sum(axis=0).tolist()
enum_rows = []
exact_probability_sum = Fraction(0)
exact_tail_sum = Fraction(0)
for tab_tuple in enumerate_tables(row_sums, col_sums):
    tab = np.array(tab_tuple, dtype=int)
    st = ordinal_stats(tab, ROW_SCORES, COL_SCORES)
    pr = table_probability(tab_tuple, row_sums, col_sums)
    extreme = abs(st["T"]) >= abs(observed["mini"]["T"]) - 1e-10
    exact_probability_sum += pr
    if extreme:
        exact_tail_sum += pr
    enum_rows.append({
        "T": st["T"],
        "probability": float(pr),
        "probability_fraction": str(pr),
        "extreme": extreme,
    })
enum_df = pd.DataFrame(enum_rows)
assert len(enum_df) == 282 and abs(enum_df.probability.sum() - 1) < 1e-12
exact_p = enum_df.loc[enum_df.extreme, "probability"].sum()
assert abs(exact_p - 4 / 105) < 1e-12
print("完整枚举摘要：")
print(f"可能表数={len(enum_df)}，概率和={exact_probability_sum}，极端表数={int(enum_df.extreme.sum())}，加权精确P={exact_p:.7f}。")
tail_summary = (
    enum_df[enum_df.extreme]
    .groupby("T", sort=True)
    .agg(table_count=("T", "size"), probability_sum=("probability", "sum"))
    .reset_index()
)
total_label_sequences = math.factorial(int(mini.sum()))
for col_total in mini.sum(axis=0):
    total_label_sequences //= math.factorial(int(col_total))
tail_summary["label_sequence_count"] = tail_summary["probability_sum"] * total_label_sequences
if SHOW_DETAILS:
    print("enumerated tables", len(enum_df), "probability sum", enum_df.probability.sum())
    print("极端表数", int(enum_df.extreme.sum()), "加权精确P", exact_tail_sum)
    print("\n[mini] 极端T分组（由实际enum_df归并）\n", tail_summary.to_string(index=False), sep="")
    print(f"极端标签序列={int(round(tail_summary['label_sequence_count'].sum()))}/总标签序列={total_label_sequences}（核对常数2520）；理论尾部期望超越次数(B×精确p)={B * exact_p:.12g}，不等于某次实际k。")

# %% 01.5 固定边际Monte Carlo
mc = {}
for name in ("mini", "main"):
    mc[name] = conditional_mc(get_table(name), B, SEED, ROW_SCORES, COL_SCORES)
mc_summary = pd.DataFrame({name: ans["summary"] for name, ans in mc.items()}).T
print("Monte Carlo摘要（k_over_B=k/B；plus1_over_plus1=(k+1)/(B+1)；CP95为尾概率区间，不是r区间）")
print(mc_summary[["B", "k", "k_over_B", "plus1_over_plus1", "MCSE", "CP95_low", "CP95_high"]].to_string())
mini_first = mc["mini"]["first_table"]
print("小例第一轮模拟表（同一批B轮中的第一张）：\n", mini_first, sep="")
print("行合计：", mini_first.sum(axis=1))
print("列合计：", mini_first.sum(axis=0))
first_T = ordinal_stats(mini_first, ROW_SCORES, COL_SCORES)["T"]
assert np.isclose(first_T, mc["mini"]["T_sim"][0])
print(f"第一张表T={first_T:.12g}，与mini_sim_T首行相同。")
main_ans = observed["main"]
main_mc = mc["main"]["summary"]
report = (f"{int(main_ans['N'])}例资料，按当前预设评分，r={main_ans['r']:.6g}，Q={main_ans['Q']:.6g}；"
          f"固定边际双侧条件模拟B={int(main_mc['B'])}，极端{int(main_mc['k'])}次，加一P={main_mc['plus1_over_plus1']:.8g}；")
if int(main_mc["k"]) == 0:
    report += f"条件尾概率的单侧95%二项上限={main_mc['zero_exceedance_one_sided_95_upper']:.8g}。"
print("主例报告：" + report)

# %% 01.6 画图与可选保存
plot_mc = {}
for name, ans in mc.items():
    plot_mc[name] = dict(ans["summary"])
    plot_mc[name]["T_sim"] = ans["T_sim"]
plot_results(enum_df, observed, plot_mc, ROOT / "figures/ordinal_conditional_mc.png" if SAVE_RESULTS else None)
if SAVE_RESULTS:
    out = ROOT / "results"
    out.mkdir(exist_ok=True)
    pd.DataFrame(observed).T.to_csv(out / "Python_observed_stats.csv")
    enum_df.to_csv(out / "Python_mini_exact_distribution.csv", index=False)
    pd.DataFrame(formula_rows).to_csv(out / "Python_formula_steps.csv", index=False)
    tail_summary.to_csv(out / "Python_mini_tail_summary.csv", index=False)
    mc_summary.to_csv(out / "Python_mc_results.csv")
    pd.DataFrame({"T": mc["mini"]["T_sim"]}).to_csv(out / "Python_mini_sim_T.csv", index=False)
    pd.DataFrame(mini_first).to_csv(out / "Python_mini_first_table.csv", index=False, header=False)

# 练习：B改为4999/79999并从首节重跑，比较MCSE；保持观察资料不变。
