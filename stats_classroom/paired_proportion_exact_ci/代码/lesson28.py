"""统计课堂28：20人配对率差的 Wald 复算（UTF-8）。
数据：数据/paired_20.csv；依赖：numpy、matplotlib；运行：Spyder分节或整文件；
输出：默认显示一幅图，SAVE_RESULTS=True时写入运行结果/Python。
完整四方法比较请运行 lesson28.R；本入口只用原生 Python 复算 Wald。
"""
# %% 准备：定位共享数据
from pathlib import Path
import csv
import math

SAVE_RESULTS = False
places = [Path.cwd()]
if "__file__" in globals():
    places.append(Path(__file__).resolve().parent)
places.extend([p.parent for p in list(places)])
root = next((p for p in places if (p / "数据/paired_20.csv").is_file()), None)
if root is None:
    raise FileNotFoundError("请将工作目录设为项目根或代码目录")

# %% 问题1：这20人的后测率减前测率是多少？
# 先显示字段和前3行，再把同一个人的后测减前测并计数。
with (root / "数据/paired_20.csv").open(encoding="utf-8", newline="") as stream:
    records = list(csv.DictReader(stream))
print("字段：", list(records[0]))
print("前3行：", records[:3])
before = [int(row["before"]) for row in records]
after = [int(row["after"]) for row in records]
changes = [y - x for x, y in zip(before, after)]
n_plus = changes.count(1)
n_zero = changes.count(0)
n_minus = changes.count(-1)
n = len(changes)
D = (n_plus - n_minus) / n
print(f"(N+,T,N-)=({n_plus},{n_zero},{n_minus}), n={n}, D={D:.2f}")

# %% 问题2：当前confidence下的Wald区间是多少？
# SE=sqrt(((b+c)/n-D^2)/n)；z为双侧正态分位点，给出当前区间。
from statistics import NormalDist
import numpy as np
import matplotlib.pyplot as plt

confidence = 0.95
se = math.sqrt(((n_plus + n_minus) / n - D**2) / n)
z = NormalDist().inv_cdf((1 + confidence) / 2)
half_width = z * se
current_row = (confidence, D, se, D - half_width, D + half_width)
print(f"Wald标准误={se:.10f}；当前confidence={confidence:.2f}；z={z:.8f}")
print("置信水平      estimate       SE       lower       upper")
print(f"{confidence:.2f}       {D:.6f}  {se:.10f}  {current_row[3]:.7f}  {current_row[4]:.7f}")

# %% 问题3：改变置信水平时区间怎样变化？
# 固定这三档便于比较；i是行索引，levels[i]是该行的置信水平。
levels = (0.90, 0.95, 0.99)
rows = []
for level in levels:
    level_z = NormalDist().inv_cdf((1 + level) / 2)
    level_half_width = level_z * se
    rows.append((level, D, se, D - level_half_width, D + level_half_width))
print("固定比较的置信水平：", levels)
for level, estimate, std_err, lower, upper in rows:
    print(f"{level:.2f}       {estimate:.6f}  {std_err:.10f}  {lower:.7f}  {upper:.7f}")

# 横向图展示固定三档；当前confidence单独打印，不改变图中标签。
fig, ax = plt.subplots(figsize=(7, 3.2))
y = np.arange(len(rows))
for i, (level, estimate, _, lower, upper) in enumerate(rows):
    ax.plot([lower, upper], [i, i], color="#2f6f8f", linewidth=3)
    ax.plot(estimate, i, "o", color="#d95f02")
ax.axvline(D, color="0.45", linestyle="--", linewidth=1)
ax.set_yticks(y, [f"{level:.0%}" for level, *_ in rows])
ax.set_xlabel("Paired rate difference (after − before)")
ax.set_title(f"Wald intervals for paired difference (D={D:.2f})")
fig.tight_layout()
plt.show()

if SAVE_RESULTS:
    output = root / "运行结果/Python"
    output.mkdir(parents=True, exist_ok=True)
    with (output / "Wald区间.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["confidence", "estimate", "standard_error", "lower", "upper"])
        writer.writerows(rows)
# 修改练习：把问题2的confidence=.95改为.90，只重跑问题2，观察当前区间半宽变窄。
