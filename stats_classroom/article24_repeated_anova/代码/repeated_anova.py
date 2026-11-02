"""第24篇：8人×4条件重复测量方差分析（Spyder 可直接运行）。"""
from pathlib import Path
import itertools
import numpy as np
import pandas as pd
from scipy import stats, optimize
from statsmodels.stats.anova import AnovaRM

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "数据" / "四种声音条件反应率.csv"
OUT = ROOT / "运行结果" / "Python"
OUT.mkdir(parents=True, exist_ok=True)
CONDS = ["A", "B", "C", "D"]

def analysis(response, label):
    # ① 每一行仍是同一名受试者，将宽表展开成长表交给标准函数。
    response = np.asarray(response, dtype=float)
    n, k = response.shape
    long = pd.DataFrame({"subject": np.repeat(np.arange(1, n + 1), k),
                         "condition": np.tile(CONDS, n),
                         "response": response.ravel()})
    rm = AnovaRM(long, "response", "subject", within=["condition"]).fit()
    long.to_csv(OUT / f"{label}_long_data.csv", index=False)
    rm.anova_table.to_csv(OUT / f"{label}_AnovaRM.csv")
    # ② 先算均值，再分别计算总/条件/受试者/误差平方和。
    grand, cm, sm = response.mean(), response.mean(0), response.mean(1)
    ss_total = np.sum((response - grand) ** 2)
    ss_condition = n * np.sum((cm - grand) ** 2)
    ss_subject = k * np.sum((sm - grand) ** 2)
    ss_error = ss_total - ss_condition - ss_subject
    df1, df_subject, df2 = k - 1, n - 1, (n - 1) * (k - 1)
    ms_c, ms_e = ss_condition / df1, ss_error / df2
    fval, pval = ms_c / ms_e, stats.f.sf(ms_c / ms_e, df1, df2)
    # 手算F/P必须与标准函数实际输出一致。
    assert np.isclose(fval, rm.anova_table.loc["condition", "F Value"], rtol=1e-10)
    assert np.isclose(pval, rm.anova_table.loc["condition", "Pr > F"], rtol=1e-10)
    # 每个单元格相对行均值、列均值和grand mean的中心化剩余量
    residual = response - sm[:, None] - cm[None, :] + grand
    assert np.isclose(np.sum(residual**2), ss_error, atol=1e-10)
    pd.DataFrame(response-sm[:, None], columns=CONDS).to_csv(OUT / f"{label}_person_centered.csv", index=False)
    # ③ AnovaRM不提供GG校正；从协方差矩阵单独手算epsilon。
    D = np.eye(k) - np.ones((k, k)) / k
    S = np.cov(response, rowvar=False, ddof=1)
    eps = np.trace(D @ S) ** 2 / ((k - 1) * np.trace(D @ S @ D @ S))
    gg1, gg2 = eps * df1, eps * df2
    pgg = stats.f.sf(fval, gg1, gg2)
    eta = ss_condition / (ss_condition + ss_error)
    # ④ 经典非中心F区间：使用未GG调整的df，不是GG稳健区间。
    ncp = []
    for target in (.975, .025):
        hi = 1.0
        while stats.ncf.cdf(fval, df1, df2, hi) > target:
            hi *= 2
        ncp.append(optimize.brentq(lambda z: stats.ncf.cdf(fval, df1, df2, z) - target, 0, hi))
    eta_low, eta_high = (z / (z + df2) for z in ncp)
    rows = [{"scenario": label, "measure": "repeated_rm", "term": "condition",
             "SS": ss_condition, "df": df1, "MS": ms_c, "F": fval, "P": pval,
             "partial_eta2": eta, "CI_low": eta_low, "CI_high": eta_high},
            {"scenario": label, "measure": "repeated_rm", "term": "subject",
             "SS": ss_subject, "df": df_subject, "MS": ss_subject / df_subject},
            {"scenario": label, "measure": "repeated_rm", "term": "error",
             "SS": ss_error, "df": df2, "MS": ms_e},
            {"scenario": label, "measure": "repeated_rm", "term": "total", "SS": ss_total, "df": n*k-1},
            {"scenario": label, "measure": "GG", "term": "condition", "epsilon_GG": eps,
             "df1_GG": gg1, "df2_GG": gg2, "F": fval, "P_GG": pgg},
            {"scenario": label, "measure": "independent_wrong_demo", "term": "condition",
             "F": stats.f_oneway(*[response[:, j] for j in range(k)]).statistic,
             "P": stats.f_oneway(*[response[:, j] for j in range(k)]).pvalue,
             "df1": k-1, "df2": n*k-k}]
    # ⑤ 六对独立受试者差值：逐对t检验，Bonferroni P与同时区间。
    for i, j in itertools.combinations(range(k), 2):
        d = response[:, j] - response[:, i]  # 后条件−前条件
        t = stats.ttest_rel(response[:, j], response[:, i])
        rows.append({"scenario": label, "measure": "paired_t", "term": f"{CONDS[j]}-{CONDS[i]}",
                     "mean_diff_percentage_points": d.mean(), "SD": d.std(ddof=1),
                     "SE": d.std(ddof=1) / np.sqrt(n), "df": n-1, "t": t.statistic,
                     "rawP": t.pvalue, "bonferroniP": min(6*t.pvalue, 1),
                     "CI_level": 1-.05/6,
                     "CI_low": d.mean()-stats.t.ppf(1-.05/12, n-1)*d.std(ddof=1)/np.sqrt(n),
                     "CI_high": d.mean()+stats.t.ppf(1-.05/12, n-1)*d.std(ddof=1)/np.sqrt(n),
                     "Shapiro_P": stats.shapiro(d).pvalue})
    pd.DataFrame(rows).to_csv(OUT / f"{label}_summary.csv", index=False)
    pd.DataFrame(residual, columns=CONDS).to_csv(OUT / f"{label}_centered_residual.csv", index=False)
    pd.DataFrame({"condition": CONDS, "condition_mean": cm}).to_csv(OUT / f"{label}_condition_means.csv", index=False)
    pd.DataFrame({"subject": np.arange(1, n+1), "subject_mean": sm}).to_csv(OUT / f"{label}_subject_means.csv", index=False)
    result = pd.DataFrame(rows)
    result[result["measure"] == "paired_t"].to_csv(OUT / f"{label}_paired_t.csv", index=False)
    return np.array([fval, pval, eps, pgg, eta])

if __name__ == "__main__":
    x = pd.read_csv(DATA).set_index("subject")[CONDS].to_numpy(float)
    original = analysis(x, "original")
    offsets = np.arange(x.shape[0], dtype=float) * 2
    changed = analysis(x + offsets[:, None], "offsets_0_14")
    assert np.allclose(original, changed, rtol=1e-10, atol=1e-10)
    # ⑥ 生成正文图：绘图模块也可单独运行。
    import runpy
    runpy.run_path(str(Path(__file__).with_name("plot_repeated.py")), run_name="__main__")
    print("Python completed; outputs:", OUT)
