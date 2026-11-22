"""本篇自编区间函数；主脚本调用即可，想研究算法时再读这里。"""
import numpy as np
from scipy import stats
from scipy.optimize import brentq


def cramers_v_ci(table, conf_level=0.95):
    """未作偏差修正的V及非中心χ²反演双侧区间，与R effectsize对齐。"""
    q, _, df, _ = stats.chi2_contingency(table, correction=False)
    scale = np.sum(table) * min(table.shape[0] - 1, table.shape[1] - 1)
    alpha = 1 - conf_level

    def lambda_limit(target):
        # λ≥0；若目标在边界外，取0。向右扩展区间，直到根被括住。
        at_zero = stats.ncx2.cdf(q, df, 0)
        if q == 0 or at_zero <= target:
            return 0.0
        upper = max(1.0, q)
        while stats.ncx2.cdf(q, df, upper) > target:
            upper *= 2
        return brentq(
            lambda value: stats.ncx2.cdf(q, df, value) - target,
            0, upper, xtol=1e-12,
        )

    low_lambda = lambda_limit(1 - alpha / 2)
    high_lambda = lambda_limit(alpha / 2)
    return (
        np.sqrt(q / scale),
        min(1.0, np.sqrt(low_lambda / scale)),
        min(1.0, np.sqrt(high_lambda / scale)),
    )


def rate_difference_ci(success, total, conf_level=0.95):
    """前组减后组：R prop.test(correct=TRUE)的单项近似区间，返回p1,p2,差,下限,上限。"""
    p1, p2 = np.asarray(success) / np.asarray(total)
    n1, n2 = total
    difference = p1 - p2
    reciprocal_sum = 1 / n1 + 1 / n2
    correction = min(0.5, abs(difference) / reciprocal_sum)
    se = np.sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2)
    width = stats.norm.ppf((1 + conf_level) / 2) * se
    width += correction * reciprocal_sum
    return p1, p2, difference, max(-1.0, difference - width), min(1.0, difference + width)
