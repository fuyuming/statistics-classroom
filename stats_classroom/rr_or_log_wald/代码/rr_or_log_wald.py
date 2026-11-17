# 统计课堂：RR、OR与log-Wald区间（Python / Spyder，UTF-8）
# 数据为人为设定的教学数据，不是真实患者资料。按 # %% 从上到下运行。
from pathlib import Path
import os
import numpy as np
import pandas as pd
from scipy import stats

# 改为True可保存；环境变量用于教师批量核对。
SAVE_RESULTS = False
if os.environ.get("RR_OR_SAVE_RESULTS") == "1":
    SAVE_RESULTS = True

# 分节运行时 __file__ 可能不存在；此时请把项目目录设为当前工作目录。
if "__file__" in globals():
    project_dir = Path(__file__).resolve().parents[1]
else:
    project_dir = Path.cwd()
if not (project_dir / "数据" / "rehab_rates.csv").exists():
    raise FileNotFoundError("请把当前工作目录设为项目根目录。")

# %% 0 读取与检查：先看字段，再确认每组唯一且计数有效
data = pd.read_csv(project_dir / "数据" / "rehab_rates.csv", encoding="utf-8")
print("数据预览：\n", data)
if list(data.columns) != ["group", "success", "failure"]:
    raise ValueError("输入字段必须为 group, success, failure。")
if len(data) != 2 or data["group"].nunique() != 2 or set(data["group"]) != {"new", "control"}:
    raise ValueError("group 必须恰好有唯一的 new 和 control 两行。")
counts = data[["success", "failure"]]
if not np.isfinite(counts.to_numpy(dtype=float)).all() or counts.isna().any().any():
    raise ValueError("success/failure 必须是有限且非缺失的数值。")
if (counts <= 0).any().any() or not np.equal(counts, np.floor(counts)).all().all():
    raise ValueError("success/failure 必须是非负整数；本练习同时计算RR和OR，请使用四格正整数；零计数需另选方法。")
new = data.loc[data.group == "new"].iloc[0]
control = data.loc[data.group == "control"].iloc[0]
a = int(new.success)
b = int(new.failure)
c = int(control.success)
d = int(control.failure)
n1 = a + b
n0 = c + d
p1 = a / n1
p0 = c / n0
print(f"两组人数：new={n1}, control={n0}；成功率：{p1:.6f}, {p0:.6f}")

# %% 1 RD、RR、odds与OR：先算点估计，再说明它们比较的对象
rd = p1 - p0
rr = p1 / p0
odds1 = p1 / (1 - p1)
odds0 = p0 / (1 - p0)
or_value = odds1 / odds0
print(f"RD={rd:.6f}（百分点差为{rd*100:.1f}）；RR={rr:.6f}（概率比）；OR={or_value:.6f}（两组odds之比）。")
print(f"odds(new)={odds1:.6f}, odds(control)={odds0:.6f}；RD不是RR或OR。")

# %% 2 从二项方差与变换斜率得到两个log尺度SE
# p(1-p)/n 是样本率方差；log(p)斜率为1/p，logit(p)斜率为1/[p(1-p)]。
var_p1 = p1 * (1 - p1) / n1
var_p0 = p0 * (1 - p0) / n0
slope_log_p1 = 1 / p1
slope_log_p0 = 1 / p0
slope_logit_p1 = 1 / (p1 * (1 - p1))
slope_logit_p0 = 1 / (p0 * (1 - p0))
var_log_p1 = slope_log_p1**2 * var_p1
var_log_p0 = slope_log_p0**2 * var_p0
var_logit_p1 = slope_logit_p1**2 * var_p1
var_logit_p0 = slope_logit_p0**2 * var_p0
se_log_rr = np.sqrt(var_log_p1 + var_log_p0)
se_log_or = np.sqrt(var_logit_p1 + var_logit_p0)
se_log_rr_cells = np.sqrt(1 / a - 1 / n1 + 1 / c - 1 / n0)
se_log_or_cells = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
print(f"var(p1)={var_p1:.8f}, var(p0)={var_p0:.8f}")
print(f"log斜率={slope_log_p1:.6f}, {slope_log_p0:.6f}；logit斜率={slope_logit_p1:.6f}, {slope_logit_p0:.6f}")
print(f"变换后方差：log p1={var_log_p1:.8f}, log p0={var_log_p0:.8f}；logit p1={var_logit_p1:.8f}, logit p0={var_logit_p0:.8f}")
print(f"SE_log_RR={se_log_rr:.15f}（四格公式={se_log_rr_cells:.15f}）；SE_log_OR={se_log_or:.15f}（四格公式={se_log_or_cells:.15f}）。")

# %% 3 log尺度95%区间，再用exp回到RR/OR尺度
# ppf输入左侧累计概率；0.975对应双侧5%区间右侧界值。
z = stats.norm.ppf(0.975)
log_rr = np.log(rr)
log_or = np.log(or_value)
log_rr_lower = log_rr - z * se_log_rr
log_rr_upper = log_rr + z * se_log_rr
log_or_lower = log_or - z * se_log_or
log_or_upper = log_or + z * se_log_or
# exp把log尺度端点逐一还原为比值尺度。
rr_lower = np.exp(log_rr_lower)
rr_upper = np.exp(log_rr_upper)
or_lower = np.exp(log_or_lower)
or_upper = np.exp(log_or_upper)
print(f"z_.975={z:.12f}；RR log区间=[{log_rr_lower:.8f},{log_rr_upper:.8f}]，原尺度=[{rr_lower:.8f},{rr_upper:.8f}]")
print(f"OR log区间=[{log_or_lower:.8f},{log_or_upper:.8f}]，原尺度=[{or_lower:.8f},{or_upper:.8f}]")

# %% 4 交换两组：比值倒数，SE不变，端点变为[1/U,1/L]
swapped_rr = p0 / p1
swapped_or = odds0 / odds1
swapped_se_log_rr = np.sqrt(1 / c - 1 / n0 + 1 / a - 1 / n1)
swapped_se_log_or = np.sqrt(1 / c + 1 / d + 1 / a + 1 / b)
swapped_log_rr = np.log(swapped_rr)
swapped_log_or = np.log(swapped_or)
swapped_log_rr_lower = swapped_log_rr - z * swapped_se_log_rr
swapped_log_rr_upper = swapped_log_rr + z * swapped_se_log_rr
swapped_log_or_lower = swapped_log_or - z * swapped_se_log_or
swapped_log_or_upper = swapped_log_or + z * swapped_se_log_or
swapped_rr_lower = np.exp(swapped_log_rr_lower)
swapped_rr_upper = np.exp(swapped_log_rr_upper)
swapped_or_lower = np.exp(swapped_log_or_lower)
swapped_or_upper = np.exp(swapped_log_or_upper)
print(f"交换后RR={swapped_rr:.8f}≈1/RR；OR={swapped_or:.8f}≈1/OR；实际重算SE={swapped_se_log_rr:.15f}, {swapped_se_log_or:.15f}。")
print(f"交换后RR log区间=[{swapped_log_rr_lower:.8f},{swapped_log_rr_upper:.8f}]，与倒数端点误差={max(abs(swapped_rr_lower-1/rr_upper), abs(swapped_rr_upper-1/rr_lower)):.3e}。")
print(f"交换后OR log区间=[{swapped_log_or_lower:.8f},{swapped_log_or_upper:.8f}]，与倒数端点误差={max(abs(swapped_or_lower-1/or_upper), abs(swapped_or_upper-1/or_lower)):.3e}。")

# %% 5 人数同时乘4：率与点估计不变，SE减半；自主练习可改倍数
scale = 4
a4 = scale * a
b4 = scale * b
c4 = scale * c
d4 = scale * d
se_rr4 = np.sqrt(1/a4 - 1/(scale*n1) + 1/c4 - 1/(scale*n0))
se_or4 = np.sqrt(1/a4 + 1/b4 + 1/c4 + 1/d4)
print(f"人数乘{scale}后点估计不变；SE_log_RR={se_rr4:.15f}、SE_log_OR={se_or4:.15f}，均为原来的{1 / np.sqrt(scale):.4f}倍。")

# %% 可选保存：同一schema，默认关闭
metrics = ["p1", "p0", "RD", "RR", "OR", "z", "se_log_RR", "se_log_OR", "log_RR_lower", "log_RR_upper", "log_OR_lower", "log_OR_upper", "RR_lower", "RR_upper", "OR_lower", "OR_upper"]
values = [p1, p0, rd, rr, or_value, z, se_log_rr, se_log_or, log_rr_lower, log_rr_upper, log_or_lower, log_or_upper, rr_lower, rr_upper, or_lower, or_upper]
results = pd.DataFrame({"metric": metrics, "value": values})
if SAVE_RESULTS:
    out = project_dir / "运行结果" / "Python"
    out.mkdir(parents=True, exist_ok=True)
    results.to_csv(out / "results.csv", index=False, encoding="utf-8")
    print("已保存运行结果/Python/results.csv")
