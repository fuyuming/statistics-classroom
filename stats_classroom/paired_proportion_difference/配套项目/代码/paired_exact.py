# 配对率差的三项无条件精确法（UTF-8，Python原生，不调用R）。
# 本文件的区间按D排序，与PairedCI的Wang归纳排序区间分别报告。
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.special import gammaln, xlogy


def paired_exact(nplus, nzero, nminus, conf_level=0.95, grid_points=201):
    """方向为后测减前测。返回按D排序的central区间及零差双侧P值。"""
    counts = np.asarray([nplus, nzero, nminus], dtype=float)
    if not np.all(np.isfinite(counts)) or np.any(counts != np.floor(counts)):
        raise ValueError("三类人数须为有限整数。")
    if np.any(counts < 0) or counts.sum() == 0:
        raise ValueError("人数须非负，总人数须大于0。")
    nplus, nzero, nminus = map(int, counts)
    n = nplus + nzero + nminus
    if not 0 < conf_level < 1:
        raise ValueError("置信水平须在0与1之间。")
    if not np.isfinite(grid_points) or int(grid_points) != grid_points or grid_points < 3:
        raise ValueError("grid_points须为至少3的整数。")
    grid_points = int(grid_points)
    observed = (nplus - nminus) / n
    tail_alpha = (1 - conf_level) / 2
    # 列出所有合法三类计数。每一行人数相加等于n。
    cells = np.array([(r, n-r-m, m) for r in range(n+1) for m in range(n-r+1)])
    differences = cells[:, 0] - cells[:, 2]
    observed_integer = nplus - nminus
    upper_mask = (differences >= observed_integer).astype(float)
    lower_mask = (differences <= observed_integer).astype(float)
    log_coefficient = gammaln(n+1) - gammaln(cells+1).sum(axis=1)

    def maximum_tail(delta, upper):
        # 给定theta=delta，未知相同结果概率p0在[0,1-|theta|]中。
        left = 0.0
        right = max(0.0, 1.0 - abs(delta))
        mask = upper_mask if upper else lower_mask

        def probability(p0):
            p0 = np.atleast_1d(p0)
            qp = np.clip((1-p0+delta)/2, 0, 1)
            qm = np.clip((1-p0-delta)/2, 0, 1)
            probabilities = np.column_stack([qp, p0, qm])
            # xlogy(0,0)=0，对应概率幂0^0=1，不把边界计数丢掉。
            log_mass = log_coefficient[None, :] + xlogy(
                cells[None, :, :], probabilities[:, None, :]
            ).sum(axis=2)
            mass = np.exp(log_mass)
            return np.clip(mass @ mask, 0, 1)

        if right <= left:
            return float(probability(left)[0])
        grid = np.linspace(left, right, grid_points)
        values = probability(grid)
        best = float(values.max())
        if best >= 1 - 1e-14:
            return 1.0
        # 保留端点、网格最大值，再细化所有检出的局部峰。
        peaks = np.flatnonzero(
            (values[1:-1] >= values[:-2])
            & (values[1:-1] >= values[2:])
            & ((values[1:-1] > values[:-2]) | (values[1:-1] > values[2:]))
        ) + 1
        peaks = np.unique(np.append(peaks, np.argmax(values)))
        for index in peaks:
            if 0 < index < grid_points - 1:
                fit = minimize_scalar(
                    lambda q: -float(probability(q)[0]),
                    bounds=(grid[index - 1], grid[index + 1]),
                    method="bounded", options={"xatol": 1e-12},
                )
                if not fit.success:
                    raise RuntimeError("概率最大化未收敛，请检查输入或增加搜索网格。")
                best = max(best, float(-fit.fun))
        return best

    upper_p = maximum_tail(0.0, True)
    lower_p = maximum_tail(0.0, False)
    p_value = min(1.0, 2 * min(upper_p, lower_p))

    def endpoint(left, right, upper_tail):
        # 下限使用随delta递增的右尾P，上限使用随delta递减的左尾P。
        for _ in range(60):
            middle = (left + right) / 2
            value = maximum_tail(middle, upper_tail)
            if (value < tail_alpha) == upper_tail:
                left = middle
            else:
                right = middle
            if right - left <= 1e-7:
                break
        # 取外侧括号，避免末位舍入把区间缩窄。
        return left if upper_tail else right

    lower = -1.0 if observed == -1 else endpoint(-1.0, observed, True)
    upper = 1.0 if observed == 1 else endpoint(observed, 1.0, False)
    return dict(method="paired_fixed_D_central", estimate=observed,
                lower=lower, upper=upper, p_value=p_value, grid_points=grid_points)


# %% 默认教学例题。交换nplus与nminus，可检查相减方向。
if __name__ == "__main__":
    result = paired_exact(nplus=1, nzero=14, nminus=5)
    print("20人：1个+1、14个0、5个-1，方向为后测减前测。")
    print(f"率差={result['estimate']:.6f}，双侧P={result['p_value']:.8f}")
    print(f"95%区间=[{result['lower']:.8f}, {result['upper']:.8f}]")
    print("这是按D排序的三项无条件精确区间，README另列PairedCI区间。")
    output = Path(__file__).resolve().parent.parent / "运行结果" / "Python"
    output.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([result]).to_csv(output/"共同精确区间.csv", index=False)
