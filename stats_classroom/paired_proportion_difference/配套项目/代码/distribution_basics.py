# 基础分布：固定已知的三类概率，展示率差D的精确分布。
# %% 定位项目目录，支持Spyder F5和从项目根运行
from pathlib import Path
import math
import pandas as pd
from scipy.stats import multinomial
script_dir = Path(__file__).resolve().parent
candidates = [Path.cwd(), script_dir, script_dir.parent]
project = next(p for p in candidates if (p / "数据").is_dir() and (p / "运行结果").is_dir())
out_dir = project / "运行结果" / "Python"
out_dir.mkdir(parents=True, exist_ok=True)

# %% 固定模型参数与n=2的D=0小演示
n = 20; q_plus, q_zero, q_minus = 0.05, 0.70, 0.25
n2_demo = pd.DataFrame({"n_plus": [0, 1], "n_zero": [2, 0], "n_minus": [0, 1], "probability": [q_zero**2, 2*q_plus*q_minus]})
print(n2_demo)
print(pd.DataFrame({"D": [-1, -0.5, 0, 0.5, 1], "probability": [q_minus**2, 2*q_minus*q_zero, q_zero**2+2*q_plus*q_minus, 2*q_plus*q_zero, q_plus**2]}))

# %% 观察人数：安排数、指定人员安排概率和总概率
n_plus_obs, n_zero_obs, n_minus_obs = 1, 14, 5
observed_arrangements = math.comb(n, n_plus_obs) * math.comb(n-n_plus_obs, n_zero_obs)
single_arrangement_probability = q_plus**n_plus_obs * q_zero**n_zero_obs * q_minus**n_minus_obs
observed_probability = observed_arrangements * single_arrangement_probability
print(pd.DataFrame([{ "n_plus": n_plus_obs, "n_zero": n_zero_obs, "n_minus": n_minus_obs, "arrangements": observed_arrangements, "single_arrangement_probability": single_arrangement_probability, "probability": observed_probability }]))

# %% 双循环枚举231组三类人数，并以multinomial核对拆解公式
rows = []
for n_plus in range(n+1):
    for n_zero in range(n-n_plus+1):
        n_minus = n-n_plus-n_zero
        arrangements = math.comb(n, n_plus) * math.comb(n-n_plus, n_zero)
        single = q_plus**n_plus * q_zero**n_zero * q_minus**n_minus
        probability = float(multinomial.pmf([n_plus, n_zero, n_minus], n, [q_plus, q_zero, q_minus]))
        if not math.isclose(probability, arrangements*single, rel_tol=1e-12, abs_tol=1e-15): raise RuntimeError("multinomial与拆解公式不一致")
        rows.append({"n_plus": n_plus, "n_zero": n_zero, "n_minus": n_minus, "arrangements": arrangements, "single_arrangement_probability": single, "probability": probability, "D": (n_plus-n_minus)/n})
table = pd.DataFrame(rows)
table.to_csv(out_dir / "基础分布_231组概率.csv", index=False, encoding="utf-8")

# %% 按相同D相加，并汇总全概率、均值和方差
rate_diff = table.groupby("D", as_index=False)["probability"].sum()
rate_diff.to_csv(out_dir / "基础分布_率差概率.csv", index=False, encoding="utf-8")
mean_d = float((table["D"]*table["probability"]).sum())
summary = pd.DataFrame([{ "n": n, "q_plus": q_plus, "q_zero": q_zero, "q_minus": q_minus, "combinations": len(table), "all_probability": table["probability"].sum(), "mean_D": mean_d, "var_D": ((table["D"]-mean_d)**2*table["probability"]).sum(), "observed_arrangements": observed_arrangements, "observed_probability": observed_probability }])
summary.to_csv(out_dir / "基础分布_汇总.csv", index=False, encoding="utf-8")
print(summary)
