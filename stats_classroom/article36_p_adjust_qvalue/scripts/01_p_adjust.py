# 统计课堂36：从入选门槛理解四种P值调整。
# 共享输入为1000个教学P值：800项无差异、200项有差异；标签只用于核对。
# is_null只是模拟标签，算法只读取p。在Spyder中以项目根目录为工作目录。

# %% 01.1 准备完整比较族
import os
from pathlib import Path
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
from scipy import stats

if not Path("data/screening_pvalues.csv").is_file():
    raise FileNotFoundError("请将工作目录设为项目根目录，再从01.1开始。")
raw = pd.read_csv("data/screening_pvalues.csv")
print(raw.head())
p = raw["p"].to_numpy(float)
m = len(p)
alpha = 0.05

# %% 01.2 手算Bonferroni与Holm
# order保存排序前的位置；按它把调整值放回原指标顺序。
order = np.argsort(p)
sorted_p = p[order]
rank_p = np.arange(1, m + 1)
bonf_manual = np.minimum(1.0, m * p)
holm_candidate = (m - rank_p + 1) * sorted_p
# 每项之前的检验都要通过，所以从前往后取累计最大值。
holm_sorted = np.minimum(1.0, np.maximum.accumulate(holm_candidate))
holm_manual = np.empty(m)
holm_manual[order] = holm_sorted
holm_steps = pd.DataFrame({"rank": rank_p, "p": sorted_p,
                           "candidate": holm_candidate, "adjusted": holm_sorted})
print(holm_steps.head())

# %% 01.3 手算BH与BY
bh_candidate = m * sorted_p / rank_p
# [::-1]先倒序，minimum.accumulate取累计最小，再倒回升序位置。
# 后面的某个位置通过时，前面的项一起入选，因此查看整个后段。
bh_sorted = np.minimum(1.0, np.minimum.accumulate(bh_candidate[::-1])[::-1])
bh_manual = np.empty(m)
bh_manual[order] = bh_sorted
harmonic_m = np.sum(1.0 / rank_p)
by_manual = np.minimum(1.0, harmonic_m * bh_manual)
bh_steps = pd.DataFrame({"rank": rank_p, "p": sorted_p,
                         "candidate": bh_candidate, "adjusted": bh_sorted})
print(bh_steps.head())

# %% 01.4 完整函数核对
# statsmodels返回四部分；Python 0起始索引[1]取第二个返回对象，即按输入顺序的调整P值。
function_results = pd.DataFrame({
    "bonf": multipletests(p, method="bonferroni")[1],
    "holm": multipletests(p, method="holm")[1],
    "BH": multipletests(p, method="fdr_bh")[1],
    "BY": multipletests(p, method="fdr_by")[1],
})
manual_results = pd.DataFrame({"bonf": bonf_manual, "holm": holm_manual,
                               "BH": bh_manual, "BY": by_manual})
function_check = (manual_results - function_results).abs().max()
print(function_check)
assert function_check.max() < 1e-10

# %% 01.5 两个各自单独构成比较族的三项Holm小例
# 每组三项单独构成m=3的比较族，不改变前面1000项的m。
mini_p_stop = np.array([.020, .023, .040])
mini_p_increasing = np.array([.020, .040, .100])
mini_rank = np.arange(1, 4)
mini_stop_candidate = (3 - mini_rank + 1) * mini_p_stop  # 乘数依次为3、2、1。
mini_stop_holm = np.maximum.accumulate(mini_stop_candidate)  # 前向累计最大值。
mini_increasing_candidate = (3 - mini_rank + 1) * mini_p_increasing
mini_increasing_holm = np.maximum.accumulate(mini_increasing_candidate)
mini_example = pd.DataFrame({
    "example": ["stop"] * 3 + ["increasing"] * 3,
    "p": np.r_[mini_p_stop, mini_p_increasing],
    "holm_candidate": np.r_[mini_stop_candidate, mini_increasing_candidate],
    "holm": np.r_[mini_stop_holm, mini_increasing_holm],
})
print(mini_example)
print("stop例第一关未通过，后续步骤停止；调整值仍返回完整三项。")
expected_holm_values = np.array([.060, .060, .060, .060, .080, .100])  # 仅作核对。
expected_holm = np.r_[multipletests(mini_p_stop, method="holm")[1],
                      multipletests(mini_p_increasing, method="holm")[1]]
np.testing.assert_allclose(mini_example["holm"], expected_holm_values, atol=1e-10)
np.testing.assert_allclose(mini_example["holm"], expected_holm, atol=1e-10)

# %% 01.6 原题三疗法复算：每个2x2比较使用Yates连续性校正
# therapy前缀避免污染前面1000项的p、m和调整结果。
therapy_raw = pd.read_csv("data/therapies.csv")
therapy_table = therapy_raw[["effective", "ineffective"]].to_numpy()
therapy_pairs = [(0, 1), (0, 2), (1, 2)]
therapy_names = ["physical-drug", "physical-topical", "drug-topical"]
therapy_raw_p = np.empty(3)
therapy_expected_min = np.empty(3)
for therapy_i, therapy_rows in enumerate(therapy_pairs):
    # correction=True：明确开启2×2卡方检验的Yates连续性校正。
    therapy_fit = stats.chi2_contingency(
        therapy_table[list(therapy_rows), :], correction=True)
    therapy_raw_p[therapy_i] = therapy_fit.pvalue
    therapy_expected_min[therapy_i] = therapy_fit.expected_freq.min()
therapy_bonf = np.minimum(1.0, 3 * therapy_raw_p)
therapy_order = np.argsort(therapy_raw_p)
therapy_sorted = therapy_raw_p[therapy_order]
therapy_holm_sorted = np.minimum(1.0, np.maximum.accumulate(
    (3 - np.arange(3)) * therapy_sorted))
therapy_holm = np.empty(3)
therapy_holm[therapy_order] = therapy_holm_sorted
therapy_results = pd.DataFrame({"pair": therapy_names, "raw_p_yates": therapy_raw_p,
                                "bonferroni_p": therapy_bonf, "holm_p": therapy_holm,
                                "expected_min": therapy_expected_min})
print(therapy_results)
np.testing.assert_allclose(therapy_raw_p,
                           [.016754306570385496, .000009344615515985061,
                            .047763774036015214], atol=1e-10)
assert np.sum(therapy_bonf <= .05) == 1
assert np.all(therapy_holm <= .05)

# %% 01.7 查看完整结果与发现数
results = pd.concat([raw, manual_results], axis=1)
counts = [np.sum(p <= alpha)]
counts.extend((manual_results <= alpha).sum().tolist())
summary_result = pd.DataFrame({"method": ["raw", "bonf", "holm", "BH", "BY"],
                               "discoveries": counts})
print(summary_result)
# 这些是当前输入的发现数；FWER和FDR则描述重复研究的错误性质。

# %% 01.8 可选保存
# 需要保存时先设置os.environ["ARTICLE36_SAVE_RESULTS"] = "1"。
if os.environ.get("ARTICLE36_SAVE_RESULTS") == "1":
    Path("results").mkdir(exist_ok=True)
    results.to_csv("results/Python_01_p_adjust.csv", index=False)
    summary_result.to_csv("results/Python_01_summary.csv", index=False)
    mini_example.to_csv("results/Python_01_mini_example.csv", index=False)
    therapy_results.to_csv("results/Python_01_therapies.csv", index=False)
