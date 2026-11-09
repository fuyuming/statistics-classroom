# Wang排序前三步：仅演示第三位候选的选择（Python标准库版）
from pathlib import Path
import csv
import math

#%% 相对路径与固定教学参数
HERE = Path(__file__).resolve().parent
OUT = HERE.parent / "运行结果" / "Python"
OUT.mkdir(parents=True, exist_ok=True)
n = 20
alpha_tail = 0.025
candidates = {"A": (19, 0, 1), "B": (18, 2, 0)}

#%% 给定theta、pT时，临时加入一个候选后的尾概率
def tail_probability(theta, p_t, candidate):
    a = (1 - p_t + theta) / 2
    b = (1 - p_t - theta) / 2
    base = a**20 + 20 * a**19 * p_t
    return base + (20 * a**19 * b if candidate == "A" else 190 * a**18 * p_t**2)

#%% 全域网格后，在网格峰附近作黄金分割，并比较两个端点
def peak_over_pt(theta, candidate):
    right = 1 - abs(theta)
    grid = [right * i / 400 for i in range(401)]
    values = [tail_probability(theta, p, candidate) for p in grid]
    k = max(range(401), key=lambda i: values[i])
    best_p, best_q = grid[k], values[k]
    left = grid[max(0, k - 1)]
    right_local = grid[min(400, k + 1)]
    golden = (math.sqrt(5) - 1) / 2
    x1 = right_local - golden * (right_local - left)
    x2 = left + golden * (right_local - left)
    for _ in range(100):
        if tail_probability(theta, x1, candidate) < tail_probability(theta, x2, candidate):
            left, x1 = x1, x2
            x2 = left + golden * (right_local - left)
        else:
            right_local, x2 = x2, x1
            x1 = right_local - golden * (right_local - left)
    for p in (0.0, 1 - abs(theta), (left + right_local) / 2):
        q = tail_probability(theta, p, candidate)
        if q > best_q:
            best_p, best_q = p, q
    return best_p, best_q

#%% 对每个候选单独求maxQ(theta)=.025的最小根
rows = []
for candidate, counts in candidates.items():
    lower = 0.0
    upper = 1.0
    for _ in range(70):
        middle = (lower + upper) / 2
        if peak_over_pt(middle, candidate)[1] >= alpha_tail:
            upper = middle
        else:
            lower = middle
    peak_p_t, tail_q = peak_over_pt(upper, candidate)
    rows.append({"candidate": candidate, "n_plus": counts[0], "n_zero": counts[1], "n_minus": counts[2],
                 "D": (counts[0] - counts[2]) / n, "trial_lower": upper, "peak_p_t": peak_p_t,
                 "tail_probability": tail_q, "selected_for_third": False})

#%% 比较计算出的临时下限，相同最大值一起选入
best_lower = max(row["trial_lower"] for row in rows)
for row in rows:
    row["selected_for_third"] = abs(row["trial_lower"] - best_lower) <= 1e-10

#%% A候选在pT=0时的解析式复核
a_at_a = (1 + rows[0]["trial_lower"]) / 2
print("A_analytic_pT0", a_at_a**20 + 20 * a_at_a**19 * (1 - a_at_a))

#%% 保存两类教学结果
fields = ["candidate", "n_plus", "n_zero", "n_minus", "D", "trial_lower", "peak_p_t", "tail_probability", "selected_for_third"]
with (OUT / "Wang排序第三位试算.csv").open("w", newline="", encoding="utf-8-sig") as f:
    writer = csv.DictWriter(f, fieldnames=fields)
    writer.writeheader()
    writer.writerows(rows)
tail_rows = []
for candidate, counts in candidates.items():
    for triple in ((20, 0, 0), (19, 1, 0), counts):
        tail_rows.append({"candidate": candidate, "n_plus": triple[0], "n_zero": triple[1], "n_minus": triple[2]})
with (OUT / "Wang排序候选尾部.csv").open("w", newline="", encoding="utf-8-sig") as f:
    writer = csv.DictWriter(f, fieldnames=["candidate", "n_plus", "n_zero", "n_minus"])
    writer.writeheader()
    writer.writerows(tail_rows)
assert abs(rows[0]["tail_probability"] - alpha_tail) < 1e-10
assert abs(rows[1]["tail_probability"] - alpha_tail) < 1e-10
assert rows[1]["selected_for_third"] is True
assert all(sum(row[k] for k in ("n_plus", "n_zero", "n_minus")) == n for row in tail_rows)
print(rows)
