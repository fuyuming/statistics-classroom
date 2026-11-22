"""输入4x4计数与同序一致权重；输出[Po, Pe, κ, 一般delta SE, Wald下限, 上限]。
梯度是一格联合概率改变对κ的一阶系数；不是sklearn的区间接口。"""
import numpy as np

def kappa_parts(counts, weight):
    counts = np.asarray(counts, dtype=float)
    weight = np.asarray(weight, dtype=float)
    if counts.shape != (4, 4):
        raise ValueError("counts必须是4×4计数表。")
    if weight.shape != (4, 4):
        raise ValueError("weight必须是4×4一致权重矩阵。")
    if not np.all(np.isfinite(counts)) or np.any(counts < 0) or np.any(counts != np.floor(counts)):
        raise ValueError("counts必须为有限的非负整数。")
    if not np.all(np.isfinite(weight)) or np.any(weight < 0) or np.any(weight > 1):
        raise ValueError("weight必须为0到1之间的有限数。")
    if not np.allclose(np.diag(weight), 1.0, atol=1e-12, rtol=0):
        raise ValueError("weight对角线必须全为1。")
    n = counts.sum()
    if not np.isfinite(n) or n <= 0:
        raise ValueError("配对表总人数必须为正数。")
    # 每格人数/n是经验联合比例；axis=1按行求和，axis=0按列求和。
    p = counts / n
    row = p.sum(axis=1)
    col = p.sum(axis=0)
    po = np.sum(weight * p)
    # outer把行比例与列比例逐一相乘，构造同一顺序的机会参照表。
    pe = np.sum(weight * np.outer(row, col))
    if 1 - pe <= 1e-14:
        return np.array([po, pe, np.nan, np.nan, np.nan, np.nan])
    kappa = (po - pe) / (1 - pe)
    # 梯度同时包含真实联合概率与两边际的改变，保留同人评级的关联。
    row_score = weight @ col  # 按行等级汇总机会一致分。
    col_score = row @ weight  # 按列等级汇总机会一致分。
    gradient = (weight - (1 - kappa) * (row_score[:, None] + col_score[None, :])) / (1 - pe)
    variance = (np.sum(p * gradient ** 2) - np.sum(p * gradient) ** 2) / n
    se = np.sqrt(max(0.0, variance))
    # 标准正态0.975分位数，用于双侧95%一般Wald区间。
    z = 1.959963984540054
    return np.array([po, pe, kappa, se, kappa - z * se, kappa + z * se])
