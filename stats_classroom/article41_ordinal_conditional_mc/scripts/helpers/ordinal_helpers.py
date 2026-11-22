"""Shared calculations for the ordinal conditional simulation lesson."""
from fractions import Fraction
import math
import numpy as np
from scipy import stats


def ordinal_stats(tab, row_scores, col_scores):
    tab = np.asarray(tab, dtype=int)
    row_scores = np.asarray(row_scores, dtype=float)
    col_scores = np.asarray(col_scores, dtype=float)
    n = int(tab.sum())
    ii, jj = np.meshgrid(row_scores, col_scores, indexing="ij")
    mr = float((tab * ii).sum() / n)
    mc = float((tab * jj).sum() / n)
    num = float((tab * (ii - mr) * (jj - mc)).sum())
    lxx = float((tab * (ii - mr) ** 2).sum())
    lyy = float((tab * (jj - mc) ** 2).sum())
    r = num / math.sqrt(lxx * lyy)
    q = (n - 1) * r * r
    return {"N": n, "T": num, "Lxx": lxx, "Lyy": lyy, "r": r,
            "VarT": lxx * lyy / (n - 1), "Q": q,
            "chi2_p": stats.chi2.sf(q, 1)}


def formula_steps(tab, ans, name, row_scores, col_scores, show_details=False):
    n = int(np.sum(tab))
    row_scores = np.asarray(row_scores)
    col_scores = np.asarray(col_scores)
    row_total = float((tab * row_scores[:, None]).sum())
    col_total = float((tab * col_scores[None, :]).sum())
    row_mean, col_mean = row_total / n, col_total / n
    row1_lxx = float((tab[0, :] * (row_scores[0] - row_mean) ** 2).sum())
    rows = [{"example": name, "quantity": q, "value": v}
            for q, v in [("row_score_total", row_total), ("col_score_total", col_total),
                         ("row_mean", row_mean), ("col_mean", col_mean),
                         ("sqrt(Lxx*Lyy)", math.sqrt(ans["Lxx"] * ans["Lyy"])),
                         ("VarT", ans["VarT"]), ("SDT", math.sqrt(ans["VarT"])),
                         ("Z", ans["T"] / math.sqrt(ans["VarT"])),
                         ("row_1_Lxx_contribution", row1_lxx)]]
    for col in (0, 3):
        product = (row_scores[0] - row_mean) * (col_scores[col] - col_mean)
        # 单元格行的value留空；具体中心化积与T贡献各有专列。
        rows.append({"example": name, "quantity": f"(1,{col + 1})", "value": np.nan,
                     "count": int(tab[0, col]), "centered_product": product,
                     "T_contribution": int(tab[0, col]) * product})
    if show_details:
        print(f"\n[{name}] 演算中间量")
        print(f"两评分总和=({row_total}, {col_total})；均值=({row_mean:.12g}, {col_mean:.12g})")
        print(f"第一行（评分{row_scores[0]:g}）对Lxx贡献={row1_lxx:.12g}")
        for row in rows[-2:]:
            print(f"格{row['quantity']}：一人的中心化积={row['centered_product']:.12g}；乘人数后的T贡献={row['T_contribution']:.12g}")
        print(f"sqrt(Lxx*Lyy)={math.sqrt(ans['Lxx']*ans['Lyy']):.12g}；VarT={ans['VarT']:.12g}；SDT={math.sqrt(ans['VarT']):.12g}；Z={ans['T']/math.sqrt(ans['VarT']):.12g}")
    return rows


def compositions(total, caps):
    # caps是各列剩余名额；当前行每格不能超过对应名额。
    if len(caps) == 1:
        if total <= caps[0]:
            yield (total,)
        return
    for x in range(min(total, caps[0]) + 1):
        for tail in compositions(total - x, caps[1:]):
            yield (x,) + tail


def enumerate_tables(row_sums, col_sums):
    # 每填完一行就扣除remain，只有所有列名额用完才保留完整表。
    def rec(i, remain, rows):
        if i == len(row_sums):
            if not any(remain):
                yield tuple(tuple(x) for x in rows)
            return
        for row in compositions(row_sums[i], remain):
            next_remain = tuple(remain[j] - row[j] for j in range(len(remain)))
            yield from rec(i + 1, next_remain, rows + [row])
    yield from rec(0, tuple(col_sums), [])


def table_probability(tab, rows, cols):
    # 固定边际下按标签序列计数加权；不同表不能等概率处理。
    numerator = math.prod(math.factorial(int(x)) for x in rows) * math.prod(math.factorial(int(x)) for x in cols)
    denominator = math.factorial(int(sum(rows))) * math.prod(math.factorial(int(x)) for row in tab for x in row)
    return Fraction(numerator, denominator)


def shuffle_table(tab, rng):
    # 保留每个人的行标签，只随机打乱固定的列标签，得到同一条件模型。
    rows = np.repeat(np.arange(tab.shape[0]), tab.sum(axis=1))
    cols = np.repeat(np.arange(tab.shape[1]), tab.sum(axis=0))
    shuffled = rng.permutation(cols)
    result = np.zeros_like(tab, dtype=int)
    np.add.at(result, (rows, shuffled), 1)
    return result


def conditional_mc(tab, B, seed, row_scores, col_scores):
    rng = np.random.default_rng(seed)
    obs = ordinal_stats(tab, row_scores, col_scores)
    values = np.empty(B)
    first_table = None
    for b in range(B):
        sampled = shuffle_table(tab, rng)
        if b == 0:  # 首表来自同一批B轮，不额外抽样。
            first_table = sampled.copy()
        values[b] = ordinal_stats(sampled, row_scores, col_scores)["T"]
    k = int(np.sum(np.abs(values) >= abs(obs["T"]) - 1e-10))
    p_hat = k / B
    ci = stats.binomtest(k, B).proportion_ci(confidence_level=.95, method="exact")
    summary = {"B": B, "seed": seed, "k": k, "k_over_B": p_hat,
               "plus1_over_plus1": (k + 1) / (B + 1),
               "MCSE": math.sqrt(p_hat * (1 - p_hat) / B),
               "CP95_low": ci.low, "CP95_high": ci.high,
               "zero_exceedance_one_sided_95_upper": 1 - .05 ** (1 / B) if k == 0 else np.nan,
               "Tobs": obs["T"]}
    return {"summary": summary, "T_sim": values, "first_table": first_table}
