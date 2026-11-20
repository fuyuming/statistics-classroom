"""本项目自编的非中心χ²区间函数，不调用R后端。"""
import numpy as np
from scipy import stats, optimize


def lambda_limit(q, df, target):
    """固定观察值q和自由度，寻找CDF等于target的非负lambda。"""
    at_zero = stats.ncx2.cdf(q, df, 0.0)
    # 非负参数空间内没有正根时，端点取0并保留边界标记。
    if q == 0 or at_zero <= target:
        return 0.0, True
    upper = max(1.0, q)
    while stats.ncx2.cdf(q, df, upper) > target:
        upper = upper * 2.0
    objective = lambda value: stats.ncx2.cdf(q, df, value) - target
    answer = optimize.brentq(objective, 0.0, upper, xtol=1e-12)
    return answer, False


def cramers_v_ci(counts, level=0.95):
    """完整步骤：联合计数→Pearson Q→lambda区间→总体V区间。"""
    counts = np.asarray(counts, dtype=float)
    if counts.ndim != 2 or min(counts.shape) < 2:
        raise ValueError("至少需要两行、两列")
    if not np.all(np.isfinite(counts)) or np.any(counts < 0):
        raise ValueError("计数必须有限且非负")
    if np.any(counts.sum(axis=0) <= 0) or np.any(counts.sum(axis=1) <= 0):
        raise ValueError("各行、列合计必须为正")
    if not 0 < level < 1:
        raise ValueError("置信水平须在0与1之间")
    n = counts.sum()
    expected = np.outer(counts.sum(axis=1), counts.sum(axis=0)) / n
    q = np.sum((counts - expected) ** 2 / expected)
    df = (counts.shape[0] - 1) * (counts.shape[1] - 1)
    m = min(counts.shape[0] - 1, counts.shape[1] - 1)
    lower, boundary_l = lambda_limit(q, df, (1 + level) / 2)
    upper, boundary_u = lambda_limit(q, df, (1 - level) / 2)
    return dict(n=n, Q=q, df=df, m=m, V=np.sqrt(q/(n*m)),
                lambda_L=lower, lambda_U=upper,
                V_L=min(1.0, np.sqrt(lower/(n*m))),
                V_U=min(1.0, np.sqrt(upper/(n*m))),
                boundary_L=boundary_l, boundary_U=boundary_u)


def interpret_v_funder(values):
    """按Funder 2019阈值给V标注；返回与输入同形状的字符串数组。"""
    values = np.asarray(values, dtype=float)
    if np.any(~np.isfinite(values)) or np.any((values < 0) | (values > 1)):
        raise ValueError("Cramer's V必须位于[0,1]")
    # 左闭右开，只有最后一档包含0.40及以上，避免边界歧义。
    bins = np.array([0.05, 0.10, 0.20, 0.30, 0.40])
    labels = np.array(["tiny", "very small", "small", "medium", "large", "very large"])
    return labels[np.digitize(values, bins, right=False)]
