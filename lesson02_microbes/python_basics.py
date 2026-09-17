"""Python基础练习：类型、结构、选择、条件、循环和函数。"""
import numpy as np
import pandas as pd
# 1. 基础类型：float、int、str、bool；None表示没有值。
weight, count, label, is_valid = 2.3, 7, "control", True
for value in (weight, count, label, is_valid):
    print(type(value).__name__)
# 2. list是可变序列；tuple是不可变序列；dict按键查值。
values = [1.8, 2.0, 2.4]
shape = (2, 3)
metadata = {"group": "control", "unit": "g"}
print(values[0], shape, metadata["unit"])
# ndarray适合数值运算；DataFrame是带列名的二维表。
a = np.array(values)
print(a * 1000) # 每个数乘1000；list * 2则是重复序列
print(values * 2)
tab = pd.DataFrame({"id": ["C01", "C02"], "weight": [1.8, 2.0]})
print(tab["weight"])
# 3. if后加冒号，代码块用四个空格缩进。
if len(values) >= 2:
    print(np.std(values, ddof=1))
else:
    print("有效观测不足两个")
# 4. 循环逐个取值；range(3)给出0、1、2，不包含3。
for value in values:
    print(value * 1000)
# 5. 函数：输入数值序列，按pandas缺失规则排除缺失。
def mean_valid(values):
    x = pd.Series(values, dtype="float64").dropna()
    if x.empty:
        return float("nan")
    return x.mean()
print(mean_valid([1.8, 2.0, None])) # 1.9
assert abs(mean_valid([1.8, 2.0, None]) - 1.9) < 1e-12
