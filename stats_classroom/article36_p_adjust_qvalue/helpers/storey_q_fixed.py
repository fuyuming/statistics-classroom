"""本篇固定lambda、常规估计模式的Storey q自编函数。"""
import numpy as np


def storey_q_fixed(p, lambda_value=0.5):
    p = np.asarray(p, dtype=float)
    if p.ndim != 1 or p.size == 0 or np.any(~np.isfinite(p)):
        raise ValueError("p必须是一维非空、有限数值向量。")
    if np.any((p < 0) | (p > 1)) or not 0 <= lambda_value < 1:
        raise ValueError("P值须在[0,1]内，lambda须在[0,1)内。")
    pi0 = min(1.0, np.mean(p >= lambda_value) / (1 - lambda_value))
    if pi0 <= 0:
        raise ValueError("所选lambda右侧无P值，无法使用本篇pi0估计；请核查设置。")
    # 后段最小值让每项取得能够将自己纳入的最小估计水平。
    m = p.size
    order = np.argsort(p)
    sorted_p = p[order]
    rank_p = np.arange(1, m + 1)
    candidate = m * sorted_p / rank_p
    bh_sorted = np.minimum(1.0, np.minimum.accumulate(candidate[::-1])[::-1])
    q = np.empty(m)
    q[order] = pi0 * bh_sorted
    return pi0, q
