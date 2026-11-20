# 统计课堂35｜关联强度与区间。输入为第11课的假设工程试验汇总计数。
# 300台独立记录器随机分3组各100台：A硅胶、B环氧、C热缩封装。
# 弯折100次后接模拟器记录10min，无持续≥1s数据缺口判合格（课堂阈值）。
# 主表按A/B/C行、合格/不合格列排列为[[90,10],[80,20],[65,35]]。
# 在项目根目录按节顺序运行；默认不写文件，ARTICLE35_SAVE_RESULTS=1才保存。

# %% 02.1 准备输入
import os
import sys
from pathlib import Path
import numpy as np
import pandas as pd

# 从项目根目录加载既有自编函数；不使用__file__，Spyder分节可独立运行。
sys.path.insert(0, str(Path.cwd() / "helpers"))
from association_ci import cramers_v_ci, interpret_v_funder

raw = pd.read_csv("data/device_pass.csv")
print(raw)
tab = raw[["pass", "fail"]].to_numpy(dtype=float)
assert np.all(tab.sum(0) > 0) and np.all(tab.sum(1) > 0)

# %% 02.2 调用函数计算V及双侧95%区间
# 这是本项目配套自编函数，内部使用SciPy数值接口且没有R后端。
fit = cramers_v_ci(tab, level=.95)

# %% 02.3 读取结果
# V为整体关联强度，V_L/V_U为区间下限/上限，confidence记录置信水平。
V = float(fit["V"])
V_L = float(fit["V_L"])
V_U = float(fit["V_U"])
ci = pd.DataFrame([{ "V": V, "CI_low": V_L, "CI_high": V_U, "confidence": .95 }])
print(ci.to_string(index=False))

# %% 02.4 效应分级
interpretation = pd.DataFrame({
    "component": ["lower", "estimate", "upper"],
    "V": [V_L, V, V_U],
    "label": interpret_v_funder([V_L, V, V_U]),
    "rule": ["funder2019"] * 3,
})
assert np.all((interpretation["V"] >= 0) & (interpretation["V"] <= 1))
print(interpretation.to_string(index=False))

# %% 02.5 本例解释与报告
# 本段用Python格式化字符串组织报告句。
print(f"本例Cramér V={V:.6f}，双侧95%区间=[{V_L:.6f}, {V_U:.6f}]；区间表示估计不确定性。")
# 第0、1、2行分别是下限、点估计、上限的等级。
labels = interpretation["label"].to_numpy()
print(f"按Funder–Ozer规则，点估计为{labels[1]}；区间两端等级为{labels[0]}/{labels[2]}。")

# %% 02.6 可选保存
if os.environ.get("ARTICLE35_SAVE_RESULTS") == "1":
    Path("results").mkdir(exist_ok=True)
    ci.to_csv("results/Python_ci.csv", index=False)
    interpretation.to_csv("results/Python_interpretation.csv", index=False)
