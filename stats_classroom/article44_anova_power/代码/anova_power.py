"""统计课堂44：从大/小效应理解功效与样本量。UTF-8。
先按# %%从01节顺序运行；共享数据是规划参数，不是已测样本。
依赖numpy/scipy；绘图另需pandas/matplotlib。默认不覆盖随包输出。
"""
# %% 01 先看两种均数组合：个体波动相同，均数分散程度不同
from pathlib import Path
import csv
import os
import runpy
import numpy as np
from scipy import stats, optimize

try:
    root = Path(__file__).resolve().parents[1]
except NameError:
    root = Path.cwd()
if root.name == "代码":
    root = root.parent
with open(root / "数据/planning_scenarios.csv", encoding="utf-8") as handle:
    rows = list(csv.DictReader(handle))
with open(root / "数据/settings.csv", encoding="utf-8") as handle:
    settings = {row["setting"]: float(row["value"]) for row in csv.DictReader(handle)}
print("输入字段及示例行：", rows[:2])
alpha = settings["alpha"]
targets = [settings["target_power"], settings["high_target_power"]]
n_check = int(settings["n_check"])
save_results = os.getenv("ARTICLE44_SAVE_RESULTS", "0") == "1"
cases = {}
for name in ["large", "small"]:
    part = sorted([row for row in rows if row["scenario"] == name], key=lambda row: row["group"])
    means = np.array([float(row["mean_MPa"]) for row in part])
    sds = np.array([float(row["sd_MPa"]) for row in part])
    if len(means) != 4 or not np.all(sds == sds[0]) or sds[0] <= 0:
        raise ValueError("本练习每种情形应有4组，且共用一个正的SD。")
    cases[name] = (means, sds[0])
    print(name, "规划均数/MPa", means, "共同SD", sds[0])

# %% 02 效应量先决定要找的信号；再在目标功效下求每组n
# f把组均数散布除以个体组内SD；这里按k个规划均数取平方偏离平均。
def cohens_f(means, sd):
    return np.sqrt(np.mean((means - means.mean())**2)) / sd


def exact_power(n, means, sd, level):
    k = len(means)
    effect = cohens_f(means, sd)
    noncentrality = k * n * effect**2
    cutoff = stats.f.ppf(1-level, k-1, k*(n-1))
    if noncentrality == 0:
        power = stats.f.sf(cutoff, k-1, k*(n-1))
    else:
        power = stats.ncf.sf(cutoff, k-1, k*(n-1), noncentrality)
    return dict(n_per_group=n, f=effect, power=power)


def plan_n(target, means, sd, level):
    effect = cohens_f(means, sd)
    if effect == 0 and target > level:
        raise ValueError("总体均数全相同，拒绝概率为alpha，无法达到更高功效目标。")
    if exact_power(2, means, sd, level)["power"] >= target:
        continuous_n = np.nan
        n = 2
    else:
        upper = 4
        while exact_power(upper, means, sd, level)["power"] < target:
            upper *= 2
        # 连续根用于定位；整数n还要逐一回代。
        continuous_n = optimize.brentq(
            lambda value: exact_power(value, means, sd, level)["power"]-target, 2, upper)
        n = int(np.ceil(continuous_n))
        while exact_power(n, means, sd, level)["power"] < target:
            n += 1
        while n > 2 and exact_power(n-1, means, sd, level)["power"] >= target:
            n -= 1
    previous = exact_power(n-1, means, sd, level)["power"] if n > 2 else np.nan
    return dict(target=target, f=effect, continuous_n=continuous_n,
                n_per_group=n, total_n=n*len(means),
                power=exact_power(n, means, sd, level)["power"], previous_power=previous)

planning = []
for name, (means, sd) in cases.items():
    for target in targets:
        answer = plan_n(target, means, sd, alpha)
        answer["scenario"] = name
        planning.append(answer)
        print("02", name, "目标", target, "每组需", answer["n_per_group"],
              "实际功效", answer["power"], "前一整数功效", answer["previous_power"])
small_mu, small_sd = cases["small"]
small_plan_n = plan_n(targets[0], small_mu, small_sd, alpha)["n_per_group"]

# %% 03 再回头问：固定每组5个时，两种效应的检出机会相同吗？
fixed_n = []
for name, (means, sd) in cases.items():
    answer = exact_power(n_check, means, sd, alpha)
    fixed_n.append({"scenario": name, "n_per_group": answer["n_per_group"],
                    "f": answer["f"], "power": answer["power"]})
print("03 固定样本量的理论功效：", fixed_n)

# %% 04 补充：SD敏感性；不替换主例的规划参数
sensitivity = []
for sd in [2.5, 5., 10.]:
    answer = plan_n(targets[0], small_mu, sd, alpha)
    answer["sigma"] = sd
    sensitivity.append(answer)
print("04 保留小差异均数组合，改变共同SD：", sensitivity)
# %% 05 可选保存；三张核心结果表
result_tables = dict(planning=planning, fixed_n=fixed_n,
                     sensitivity=sensitivity)
if save_results:
    result_dir = root / "运行结果/Python"
    result_dir.mkdir(parents=True, exist_ok=True)
    for name, records in result_tables.items():
        with open(result_dir / f"{name}.csv", "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(records[0]))
            writer.writeheader()
            writer.writerows(records)
if os.getenv("ARTICLE44_PLOT", "1") == "1":
    runpy.run_path(str(root / "代码/plot_power.py"))
