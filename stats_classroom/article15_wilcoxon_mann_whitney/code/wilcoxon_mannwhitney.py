# %% 第一部分：看数据→明问题→立假设；显式导入，支持项目根或 code/ 工作目录。
from itertools import combinations
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, rankdata, wilcoxon
def find_data(name):
    candidates = [Path("data") / name, Path("..") / "data" / name]
    for path in candidates:
        if path.exists():
            return path
    raise FileNotFoundError(f"找不到data/{name}；请将工作目录设为配套练习根目录或code目录。")
textbook = pd.read_csv(find_data("textbook_example.csv"))
rd_data = pd.read_csv(find_data("independent_RD.csv"))
# 教材例5-8：甲、乙各5个独立小池，单位万尾；预定方向是甲高于乙。
# H0：两组总体分布相同；前提：两组独立、数据可排序；无并列值。
A = textbook.loc[textbook.group == "A", "value"].to_numpy(float)
B = textbook.loc[textbook.group == "B", "value"].to_numpy(float)
m = len(A)
n = len(B)
# %% 第二部分：选方法→查条件；合并编秩，解释原始秩和与U。
# 本例原始秩和R甲=37、R乙=18，扣除组内秩和15后U甲=22、U乙=3。
pooled_rank = rankdata(np.r_[A, B], method="average")
RA = pooled_rank[:m].sum()
RB = pooled_rank[m:].sum()
UA = RA - m * (m + 1) / 2
UB = RB - n * (n + 1) / 2
win_A = np.greater.outer(A, B).sum()
win_B = np.less.outer(A, B).sum()
win_tie = np.equal.outer(A, B).sum()
effect_A = UA / (m * n)
textbook_summary = pd.Series({"n_A":m,"n_B":n,"R_A":RA,"R_B":RB,"U_A":UA,"U_B":UB,"wins_A":win_A,"wins_B":win_B,"ties":win_tie,"effect_A":effect_A}, name="textbook_example")
textbook_summary
# %% 第三部分：做检验→判显著→看效应；枚举10选5的252个标签分配，U甲≥22有7种。
all_values = np.r_[A, B]
N = len(all_values)
assignment = list(combinations(range(N), m))
u_all = np.array(
    [sum(index + 1 for index in item) - m * (m + 1) / 2 for item in assignment]
)
count_assignments = len(u_all)
p_upper = np.mean(u_all >= UA)
p_lower = np.mean(u_all <= UA)
p_two_sided = min(1.0, 2 * min(p_upper, p_lower))
u_critical_lower = 4
p_critical_lower = np.mean(u_all <= u_critical_lower)
lower_count = np.sum(u_all <= u_critical_lower)
# HL=5；共同形状位置平移模型下，U临界值4对应第(u临界值+1)=5个有序差值2，实际覆盖20/21=.952381。
cross_diff = np.sort((A[:, None] - B[None, :]).ravel())
HL = float(np.median(cross_diff))
lower_index = u_critical_lower
lower_95_nominal = float(cross_diff[lower_index])
actual_confidence = 1 - p_critical_lower
exact_summary = pd.Series({"assignments":count_assignments,"upper_count":np.sum(u_all >= UA),"p_upper":p_upper,"p_lower":p_lower,"p_two_sided":p_two_sided,"u_critical_lower":u_critical_lower,"lower_count":lower_count,"p_critical_lower":p_critical_lower,"HL":HL,"ci_lower":lower_95_nominal,"ci_upper":np.inf,"actual_confidence":actual_confidence}, name="exact_enumeration")
exact_summary
# SciPy返回U=22；与R的W=22同为甲组U型统计量，method=exact对应无并列精确法。
python_result = mannwhitneyu(A, B, alternative="greater", method="exact")
python_result
# 交换组回归检查：交换输入后改用B<A的等价方向，单侧P应为0.0277778。
swap_result = mannwhitneyu(B, A, alternative="less", method="exact")
swap_u = float(swap_result.statistic)
swap_p_upper = np.mean(u_all >= swap_u)
swap_p_lower = np.mean(u_all <= swap_u)
swap_p_two_sided = min(1.0, 2 * min(swap_p_upper, swap_p_lower))
swap_summary = pd.Series({"swapped_U":swap_u,"swapped_p_less":swap_result.pvalue,"swapped_p_two_sided":swap_p_two_sided}, name="swapped_groups")
swap_summary
# %% 第四部分：讲结果→画图表→写结论；RD课件对照；两组n=10、12，有并列值，秩取平均，近似+连续性校正。
RD1 = rd_data.loc[rd_data.group == 1, "RD"].to_numpy(float)
RD2 = rd_data.loc[rd_data.group == 2, "RD"].to_numpy(float)
rd_rank = rankdata(np.r_[RD1, RD2], method="average")
rd_n1 = len(RD1)
rd_n2 = len(RD2)
RA_RD = rd_rank[:rd_n1].sum()
RB_RD = rd_rank[rd_n1:].sum()
UA_RD = RA_RD - rd_n1 * (rd_n1 + 1) / 2
UB_RD = RB_RD - rd_n2 * (rd_n2 + 1) / 2
rd_win_tie = np.equal.outer(RD1, RD2).sum()
rd_effect = UA_RD / (rd_n1 * rd_n2)
rd_result = mannwhitneyu(RD1,RD2,alternative="greater",method="asymptotic",use_continuity=True)
rd_summary = pd.Series({"n_1":rd_n1,"n_2":rd_n2,"R_1":RA_RD,"U_1":UA_RD,"U_2":UB_RD,"ties":rd_win_tie,"effect_1":rd_effect,"p_approx_continuity":rd_result.pvalue}, name="rd_example")
rd_summary
# %% 第五部分：配对资料；同一问题可写成差值的单样本符号秩或paired调用。
# 来源为颜艳、王彤教材课件的8对大鼠维生素A数据，X与Y按pair对应。
# H0：差值分布关于0对称；d=round(X-Y,2)，零差剔除；本例n=8、无绝对差并列。
paired_data = pd.read_csv(find_data("paired_vitaminA.csv"))
X = paired_data["X"].to_numpy(float)
Y = paired_data["Y"].to_numpy(float)
d = np.round(X - Y, 2)
d_nonzero = d[d != 0]
d_rank = rankdata(np.abs(d_nonzero), method="average")
positive_rank_sum = d_rank[d_nonzero > 0].sum()
negative_rank_sum = d_rank[d_nonzero < 0].sum()
signed_rank_summary = pd.Series({"n":len(d_nonzero),"positive_rank_sum":positive_rank_sum,"negative_rank_sum":negative_rank_sum}, name="paired_signed_rank")
signed_rank_summary
one_sample_result = wilcoxon(d, zero_method="wilcox", alternative="two-sided", method="exact")
paired_result = wilcoxon(X, Y, zero_method="wilcox", alternative="two-sided", method="exact")
one_sample_result
paired_result
# %% 可选第六部分：保存可分享的关键结果CSV；默认不写文件。
result_export = pd.Series({"language":"Python","R_A":RA,"R_B":RB,"U_A":UA,"U_B":UB,"p_upper":p_upper,"p_two_sided":p_two_sided,"HL":HL,"ci_lower":lower_95_nominal,"ci_upper":np.inf,"actual_confidence":actual_confidence,"RD_R1":RA_RD,"RD_U1":UA_RD,"RD_U2":UB_RD,"RD_effect":rd_effect,"RD_p":rd_result.pvalue}, name="results")
result_export
save_csv = False
if save_csv:
    project_root = find_data("textbook_example.csv").resolve().parent.parent
    output_dir = project_root / "results"
    output_dir.mkdir(exist_ok=True)
    result_export.to_frame().T.to_csv(output_dir / "Python_results.csv", index=False)
