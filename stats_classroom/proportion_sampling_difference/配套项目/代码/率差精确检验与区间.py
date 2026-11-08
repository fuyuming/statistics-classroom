# %% 两个独立率之差的无条件精确检验和区间（UTF-8）
# 教学构造例题：甲社区2/10，乙社区0/10，方向为甲减乙。
# 在Spyder打开本文件并F5。只用numpy、scipy、pandas，不调用R。
# 依据：exact2x2文档中的固定排序和central双侧构造，详见项目README。
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.stats import binom


def exact_risk_difference(a, n1, b, n2, conf_level=0.95, grid_points=501):
    """按率差排序的无条件精确法。数值搜索精度由grid_points控制。

    a/n1为第一组，b/n2为第二组。区间通过反演单侧检验计算。
    统计定义依赖二项模型，未知概率的最大化仍采用数值搜索。
    """
    counts = np.asarray([a, n1, b, n2], dtype=float)
    if not np.all(np.isfinite(counts)) or np.any(counts != np.floor(counts)):
        raise ValueError("阳性人数和样本量须为有限整数。")
    a, n1, b, n2 = map(int, counts)
    if n1 <= 0 or n2 <= 0 or not (0 <= a <= n1 and 0 <= b <= n2):
        raise ValueError("样本量须大于0，阳性人数须在0至样本量之间。")
    if not 0 < conf_level < 1:
        raise ValueError("置信水平须在0与1之间。")
    if not np.isfinite(grid_points) or int(grid_points) != grid_points or grid_points < 3:
        raise ValueError("grid_points须为至少3的整数。")
    grid_points = int(grid_points)
    observed = a / n1 - b / n2
    tail_alpha = (1 - conf_level) / 2

    # 两组计数交叉相乘，相同率差有相同整数，避免浮点减法打破并列。
    first_counts = np.arange(n1 + 1)
    second_counts = np.arange(n2 + 1)
    differences = first_counts[:, None] * n2 - second_counts[None, :] * n1
    observed_integer = a * n2 - b * n1
    upper_mask = (differences >= observed_integer).astype(float)
    lower_mask = (differences <= observed_integer).astype(float)

    def maximum_tail(delta, upper):
        # 给定delta=p1-p2，令q=p2，q的合法区间由两个率均在[0,1]确定。
        left = max(0.0, -delta)
        right = min(1.0, 1.0 - delta)
        mask = upper_mask if upper else lower_mask

        def probability(q):
            q = np.atleast_1d(q)
            p1 = np.clip(q + delta, 0, 1)
            p2 = np.clip(q, 0, 1)
            first = binom.pmf(first_counts[None, :], n1, p1[:, None])
            second = binom.pmf(second_counts[None, :], n2, p2[:, None])
            # 每个q都计算独立联合概率之和，mask选中这一侧的计数组合。
            result = np.sum((first @ mask) * second, axis=1)
            return np.clip(result, 0, 1)

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
    return dict(method="unconditional_difference_central", estimate=observed,
                lower=lower, upper=upper, p_value=p_value, grid_points=grid_points)


# %% 先运行例题，再尝试交换两组
if __name__ == "__main__":
    first_positive = 2
    first_total = 10
    second_positive = 0
    second_total = 10
    result = exact_risk_difference(first_positive, first_total, second_positive, second_total)
    print("第一组2/10，第二组0/10，计算第一组减第二组。")
    print(f"率差={result['estimate']:.6f}，双侧P={result['p_value']:.8f}")
    print(f"95%区间=[{result['lower']:.8f}, {result['upper']:.8f}]")
    print("区间包含0，在双侧α=0.05下，差异证据不足。")
    print("本方法按率差排序，对应R扩展的整数排序核对项，详见README。")
    print("可把grid_points从501改为1001，比较数值搜索的末位变化。")
    output = Path(__file__).resolve().parent.parent / "运行结果" / "Python"
    output.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([result]).to_csv(output / "率差精确_三语言共同方法.csv", index=False)
