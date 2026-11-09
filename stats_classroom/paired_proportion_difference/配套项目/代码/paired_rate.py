# %% 配对率差三项模型：Python 主入口
from pathlib import Path
import math
import numpy as np
import pandas as pd
from scipy.stats import norm

ROOT = Path(__file__).resolve().parents[1] if "__file__" in globals() else Path.cwd()
if ROOT.name == "代码":
    ROOT = ROOT.parent
DATA = ROOT / "数据"
OUT = ROOT / "运行结果" / "Python"
OUT.mkdir(parents=True, exist_ok=True)

def distribution(n, qp, qm):
    """不同人独立；每人变化取 -1,0,+1，逐人卷积而不是模拟。"""
    q0 = 1 - qp - qm
    probability = np.array([1.0])
    kernel = np.array([qm, q0, qp])
    for _ in range(n):
        probability = np.convolve(probability, kernel)
    k = np.arange(-n, n + 1)
    d = k / n
    theta = qp - qm
    sd = math.sqrt((qp + qm - theta ** 2) / n)
    edge = (d[:-1] + d[1:]) / 2
    normal = np.empty(2 * n + 1)
    normal[0] = norm.cdf(edge[0], theta, sd)
    normal[-1] = 1 - norm.cdf(edge[-1], theta, sd)
    normal[1:-1] = norm.cdf(edge[1:], theta, sd) - norm.cdf(edge[:-1], theta, sd)
    return pd.DataFrame({'k': k, 'D': d, 'exact_probability': probability,
                         'normal_cell_probability': normal})

def main():
    scenarios = [('normal_n2',2,.05,.25), ('normal_n20',20,.05,.25),
                 ('normal_n1000',1000,.05,.25), ('rare_n1000',1000,.001,.002)]
    moments = []
    for name,n,qp,qm in scenarios:
        table = distribution(n,qp,qm)
        table.to_csv(OUT/f'distribution_{name}.csv',index=False)
        theta = qp-qm
        variance = qp+qm-theta**2
        moments.append({'scenario':name,'n':n,'qp':qp,'q0':1-qp-qm,'qm':qm,
                        'theta':theta,'mean_D':theta,'var_D':variance/n,
                        'max_abs_cdf_error':np.max(np.abs(table.exact_probability.cumsum()-table.normal_cell_probability.cumsum()))})
    pd.DataFrame(moments).to_csv(OUT/'moments.csv',index=False)
    people = pd.read_csv(DATA/'paired_20.csv')
    assert np.isin(people[["before","after"]].to_numpy(), [0,1]).all()
    changes = people.after - people.before
    assert np.array_equal(changes.to_numpy(), people.z.to_numpy())
    examples = pd.read_csv(DATA/'scenarios.csv')
    rows=[]
    for row in examples.itertuples(index=False):
        if row.scenario == 'main_20':
            n = len(people)
            plus = int((changes == 1).sum())
            zero = int((changes == 0).sum())
            minus = int((changes == -1).sum())
            both_positive = int(((people.before == 1) & (people.after == 1)).sum())
        else:
            n = int(row.n)
            plus = int(row.nplus)
            zero = int(row.nzero)
            minus = int(row.nminus)
            both_positive = int(row.both_positive)
        D = (plus-minus)/n
        se = math.sqrt(((plus+minus)/n-D**2)/n)
        before = (both_positive+minus)/n
        after = (both_positive+plus)/n
        h = 2*math.asin(math.sqrt(after))-2*math.asin(math.sqrt(before))
        rows.append({'scenario':row.scenario,'n':n,'nplus':plus,'nzero':zero,'nminus':minus,'D':D,
                     'SE_Wald':se,'CI_low_Wald':D-norm.ppf(.975)*se,'CI_high_Wald':D+norm.ppf(.975)*se,
                     'before_rate':before,'after_rate':after,'Cohen_h_marginal':h})
    pd.DataFrame(rows).to_csv(OUT/'paired_examples.csv',index=False)
    print(pd.DataFrame(rows).to_string(index=False))

    # 本语言的绘图入口依次生成三张图，保存在Python自己的图像目录。
    import runpy
    runpy.run_path(str(ROOT / "代码" / "plot_figures.py"), run_name="__main__")

    # 逐组展示231种人数结果的安排数、概率与尾部归属。
    runpy.run_path(str(ROOT / "代码" / "three_category_probabilities.py"), run_name="__main__")

if __name__=='__main__': main()
