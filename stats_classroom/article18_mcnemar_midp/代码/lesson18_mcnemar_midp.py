# 目标：沿用第11课58人例题，比较四P并计算条件第一类错误率。
# 编码UTF-8；从项目根目录运行；共享输入在数据/，结果在运行结果/对应语言。
# 依赖、安装、对象解释和预期结果见README；按四部分顺序执行。
# 统计课堂18：McNemar精确检验与 mid-P；UTF-8；四部分，适合 Spyder 的 # %% 单元。
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import binom, binomtest, chi2

project_dir = Path.cwd()
if not (project_dir / "数据").is_dir():
    project_dir = Path(__file__).resolve().parents[1]
data = pd.read_csv(project_dir / "数据" / "sle_mcnemar.csv")
data
# a、b、c、d按“荧光行、乳胶列”读取，避免把方向读反。
a = int(data.loc[0, "latex_positive"])
b = int(data.loc[0, "latex_negative"])
c = int(data.loc[1, "latex_positive"])
d = int(data.loc[1, "latex_negative"])
N = a + b + c + d
data_summary = pd.DataFrame({"N": [N], "a": [a], "b": [b], "c": [c], "d": [d], "m": [b + c]})
data_summary

# %% 第一部分：读题数据与条件二项模型
# 只有不一致对进入方向检验；在H0下，b~Binomial(m,.5)，m=b+c固定。
m = b + c
if m == 0:
    raise ValueError("m=0时停止检验，不生成新的检验结果。")
k_values = np.arange(m + 1)
# pmf计算P(B=k)：第二个参数是试验次数m，第三个是方向概率0.5。
probabilities = binom.pmf(k_values, m, 0.5)
# T是距中心m/2的距离乘2；正比例缩放不改变极端程度排序。
t_values = np.abs(2 * k_values - m)
t_observed = abs(b - c)
# 方括号筛选所有满足条件的k，sum将对应概率相加。
strict_probability = float(probabilities[t_values > t_observed].sum())
boundary_probability = float(probabilities[t_values == t_observed].sum())
exact_P = strict_probability + boundary_probability
midP = strict_probability + 0.5 * boundary_probability
probability_table = pd.DataFrame({"k": k_values, "prob": probabilities, "T": t_values})

# %% 第二部分：逐项边界，精确法与 mid-P 核对
# exact把边界整块计入；mid-P只计边界的一半。b=c时边界是中心概率，不能用2*leftmid代替。
ordinary_Q = (b - c) ** 2 / m
continuity_Q = max(0, abs(b - c) - 1) ** 2 / m
# chi2.sf(Q, 1)计算自由度1的卡方右尾概率。
observed_four_methods = pd.DataFrame({
    "method": ["ordinary", "continuity_correction", "exact", "midP"],
    "Q": [ordinary_Q, continuity_Q, np.nan, np.nan], "df": [1, 1, np.nan, np.nan],
    "P": [chi2.sf(ordinary_Q, 1), chi2.sf(continuity_Q, 1), exact_P, midP],
    "strict_probability": [np.nan, np.nan, strict_probability, strict_probability],
    "boundary_probability": [np.nan, np.nan, boundary_probability, boundary_probability]})
binomtest_P = float(binomtest(b, m, p=0.5, alternative="two-sided").pvalue)
boundary_components = pd.DataFrame({"b": [b], "c": [c], "m": [m], "t_observed": [t_observed],
    "strict_probability": [strict_probability], "boundary_probability": [boundary_probability],
    "exact_P": [exact_P], "midP": [midP], "binomtest_P": [binomtest_P]})
observed_four_methods

# %% 第三部分：原例与 b=11,c=3 的四方法比较
def four_methods(a_value, b_value, c_value, d_value, label):
    """返回一组四格数的四种P值；m=0时明确返回NA，不伪造检验。"""
    m_value = b_value + c_value
    if m_value == 0:
        return {"exercise": label, "a": a_value, "b": b_value, "c": c_value,
             "d": d_value, "N": a_value + b_value + c_value + d_value,
            
                "ordinary_P": np.nan, "continuity_P": np.nan, "exact_P": np.nan, "midP": np.nan}
    k = np.arange(m_value + 1)
    p = binom.pmf(k, m_value, 0.5)
    T = np.abs(2 * k - m_value)
    t = abs(b_value - c_value)
    strict = p[T > t].sum()
    boundary = p[T == t].sum()
    Q = (b_value - c_value) ** 2 / m_value
    Qc = max(0, abs(b_value - c_value) - 1) ** 2 / m_value
    return {"exercise": label, "a": a_value, "b": b_value, "c": c_value,
         "d": d_value, "N": a_value + b_value + c_value + d_value,
        
            "ordinary_P": chi2.sf(Q, 1), "continuity_P": chi2.sf(Qc,
                 1), "exact_P": strict + boundary, "midP": strict + 0.5 * boundary}

exercise_11_3 = pd.DataFrame([four_methods(a, b, c, d, "original_12_2"),
     four_methods(11, 11, 3, 33, "teaching_11_3")])
exercise_11_3

# %% 第四部分：枚举跳跃、条件第一类错误率与改参数练习
# 固定m=14时，k=0:m是15个可能结果；错误率按各自binomial概率加权，不除以15。
outcomes = []
for k in k_values:
    Q = (2 * k - m) ** 2 / m
    Qc = max(0, abs(2 * k - m) - 1) ** 2 / m
    strict = probabilities[t_values > t_values[k]].sum()
    boundary = probabilities[t_values == t_values[k]].sum()
    outcomes.append({"m": m, "k": int(k), "prob": probabilities[k], "ordinary_P": chi2.sf(Q, 1),
        "continuity_P": chi2.sf(Qc, 1), "exact_P": strict + boundary, "midP": strict + 0.5 * boundary,
        "ordinary_reject_alpha05": chi2.sf(Q, 1) <= .05, "continuity_reject_alpha05": chi2.sf(Qc,
             1) <= .05,
        "exact_reject_alpha05": strict + boundary <= .05, "midP_reject_alpha05": strict + .5 * boundary <= .05})
all_outcomes = pd.DataFrame(outcomes)
conditional_size = pd.DataFrame({"m": [m], "alpha": [.05], "ordinary_size": [(all_outcomes.prob * all_outcomes.ordinary_reject_alpha05).sum()],
    
    "continuity_size": [(all_outcomes.prob * all_outcomes.continuity_reject_alpha05).sum()],
    "exact_size": [(all_outcomes.prob * all_outcomes.exact_reject_alpha05).sum()], "midP_size": [(all_outcomes.prob * all_outcomes.midP_reject_alpha05).sum()],
})
# 人工模型练习：把alpha改为0.01或把m改为10；以下代码会重新枚举，不把答案写死。
alpha_exercise = []
for m_exercise, alpha_exercise_value in [(14, .01), (10, .05)]:
    k_exercise = np.arange(m_exercise + 1)
    p_exercise = binom.pmf(k_exercise, m_exercise, .5)
    T_exercise = abs(2 * k_exercise - m_exercise)
    rows = []
    for index, k_exercise_value in enumerate(k_exercise):
        strict_exercise = p_exercise[T_exercise > T_exercise[index]].sum()
        boundary_exercise = p_exercise[T_exercise == T_exercise[index]].sum()
        exact_value = strict_exercise + boundary_exercise
        mid_value = strict_exercise + .5 * boundary_exercise
        rows.append({"m": m_exercise, "alpha": alpha_exercise_value,
             "k": int(k_exercise_value), "prob": p_exercise[index],
            
            "ordinary_P": chi2.sf((2 * k_exercise_value - m_exercise) ** 2 / m_exercise, 1),
            "continuity_P": chi2.sf(max(0, abs(2 * k_exercise_value - m_exercise) - 1) ** 2 / m_exercise,
                 1),
            "exact_P": exact_value, "midP": mid_value,
            "ordinary_reject": chi2.sf((2 * k_exercise_value - m_exercise) ** 2 / m_exercise,
                 1) <= alpha_exercise_value,
            "continuity_reject": chi2.sf(max(0, abs(2 * k_exercise_value - m_exercise) - 1) ** 2 / m_exercise,
                 1) <= alpha_exercise_value,
            "exact_reject": exact_value <= alpha_exercise_value, "midP_reject": mid_value <= alpha_exercise_value})
    alpha_exercise.extend(rows)
alpha_exercise = pd.DataFrame(alpha_exercise)
summary_rows = []
for m_summary in sorted(alpha_exercise.m.unique()):
    summary_part = alpha_exercise[alpha_exercise.m == m_summary]
    for alpha_summary in sorted(summary_part.alpha.unique()):
        group = summary_part[summary_part.alpha == alpha_summary]
        summary_rows.append({"m": m_summary, "alpha": alpha_summary,
            "ordinary_size": (group.prob * group.ordinary_reject).sum(),
            "continuity_size": (group.prob * group.continuity_reject).sum(),
            "exact_size": (group.prob * group.exact_reject).sum(),
            "midP_size": (group.prob * group.midP_reject).sum()})
alpha_exercise_summary = pd.DataFrame(summary_rows)
# 同时列出拒绝域中的k，便于与概率汇总逐项核对。
for row_index in alpha_exercise_summary.index:
    summary_m = alpha_exercise_summary.loc[row_index, "m"]
    summary_alpha = alpha_exercise_summary.loc[row_index, "alpha"]
    group = alpha_exercise[
        (alpha_exercise.m == summary_m) & (alpha_exercise.alpha == summary_alpha)
    ]
    for method in ["ordinary", "continuity", "exact", "midP"]:
        rejected_k = group.loc[group[method + "_reject"], "k"]
        rejection_text = " ".join(str(value) for value in rejected_k)
        alpha_exercise_summary.loc[row_index, method + "_reject_k"] = rejection_text

# 明确的命名对象和CSV方便整文件运行后回看；不依赖print讲课。
output_dir = project_dir / "运行结果" / "Python"
output_dir.mkdir(parents=True, exist_ok=True)
for name, frame in {"observed_four_methods": observed_four_methods,
     "boundary_components": boundary_components,
                    "all_outcomes": all_outcomes, "conditional_size": conditional_size,
                    "exercise_11_3": exercise_11_3, "alpha_exercise": alpha_exercise,
                    "alpha_exercise_summary": alpha_exercise_summary}.items():
    frame.to_csv(output_dir / f"{name}.csv", index=False, encoding="utf-8")
