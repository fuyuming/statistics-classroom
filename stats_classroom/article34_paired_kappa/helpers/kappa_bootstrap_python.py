"""第34篇 Bootstrap helper；端点严格仿照 boot:::norm.inter。"""
import numpy as np
from scipy.stats import norm

def kappa_from_counts(v):
    v = np.asarray(v, dtype=np.int64)
    if v.size != 4 or np.any(v < 0) or v.sum() == 0:
        return np.nan
    a, b, c, d = v
    den = (a+b)*(b+d) + (c+d)*(a+c)
    return np.nan if den == 0 else 2.0*(a*d-b*c)/den

def norm_inter_r(t, alpha):
    """R boot::norm.inter: r=(B+1)alpha，正态分位尺度插值。"""
    t = np.sort(np.asarray(t, dtype=float)[np.isfinite(t)])
    R = len(t); rk = (R+1)*np.asarray(alpha, dtype=float)
    out = np.empty_like(rk)
    for i, (rki, ai) in enumerate(zip(rk, np.asarray(alpha))):
        k = int(np.trunc(rki))
        if k <= 0: out[i] = t[0]
        elif k >= R: out[i] = t[-1]
        elif k == rki: out[i] = t[k-1]
        else:
            q0, q1 = norm.ppf(k/(R+1)), norm.ppf((k+1)/(R+1))
            out[i] = t[k-1] + (norm.ppf(ai)-q0)/(q1-q0)*(t[k]-t[k-1])
    return out

def bca_levels(z0, acc, alpha=(.025, .975)):
    za = norm.ppf(alpha)
    return norm.cdf(z0 + (z0+za)/(1-acc*(z0+za)))
