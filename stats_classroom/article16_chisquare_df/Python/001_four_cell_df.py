# %% 统计课堂16：四格表为什么自由度是1
# 直接运行对象查看结果；不依赖裸 __file__。
from pathlib import Path
import json
import numpy as np
from scipy import __version__ as scipy_version
from scipy.stats import chi2, chi2_contingency, norm

# %% ① 看数据与立假设：数据与边际合计
O = np.array([[225, 275], [128, 272]], dtype=float)
row_total = O.sum(axis=1)
col_total = O.sum(axis=0)
N = O.sum()
O
row_total
col_total
N

# %% ② 选方法与查条件：自由度与期望数
df = (O.shape[0] - 1) * (O.shape[1] - 1)
E = np.outer(row_total, col_total) / N
E
df

# %% ③ 做检验与看效应：未校正 Pearson χ²
cell_contribution = (O - E) ** 2 / E
chi_uncorrected = cell_contribution.sum()
p_uncorrected = chi2.sf(chi_uncorrected, df)
cell_contribution
chi_uncorrected
p_uncorrected

# %% ③ 做检验与看效应：scipy 官方函数与 Yates 校正
chi_yates, p_yates, df_scipy, E_scipy = chi2_contingency(O, correction=True)
chi_yates
p_yates
df_scipy
E_scipy

# %% ③ 做检验与看效应：pooled Z、Z²、率差与区间
x = O[:, 0]
p_hat = x / row_total
p_pool = x.sum() / N
rate_difference = p_hat[0] - p_hat[1]
se_null = np.sqrt(p_pool * (1 - p_pool) * np.sum(1 / row_total))
z_pooled = rate_difference / se_null
z_squared = z_pooled ** 2
se_estimate = np.sqrt(np.sum(p_hat * (1 - p_hat) / row_total))
ci_uncorrected = rate_difference + np.array([-1, 1]) * norm.ppf(.975) * se_estimate
# Approximate interval matching this case and R prop.test default correction; not exact.
half_width_corrected = norm.ppf(.975) * se_estimate + 0.5 * np.sum(1 / row_total)
ci_corrected = rate_difference + np.array([-1, 1]) * half_width_corrected
result = {
    "rate_difference": rate_difference,
    "rate_difference_pp": 100 * rate_difference,
    "ci_uncorrected": ci_uncorrected,
    "ci_corrected": ci_corrected,
    "z_pooled": z_pooled,
    "z_squared": z_squared,
    "chi_uncorrected": chi_uncorrected,
    "p_uncorrected": p_uncorrected,
    "chi_yates": chi_yates,
    "p_yates": p_yates,
    "df": df,
}
result

# %% ④ 讲结果：版本与可复核结果文件（可选保存）
project_dir = Path.cwd()
output_dir = project_dir / "output" / "Python"
output_dir.mkdir(parents=True, exist_ok=True)
(output_dir / "001_result.txt").write_text(
    "scipy_version=" + scipy_version + "\n" + repr(result) + "\n",
    encoding="utf-8",
)
(output_dir / "001_result.json").write_text(
    json.dumps({k: (v.tolist() if hasattr(v, "tolist") else float(v))
                for k, v in result.items()}, ensure_ascii=False, indent=2) + "\n",
    encoding="utf-8",
)
