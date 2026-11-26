# %% 1. 读入三种支架设计的教学数据，单位MPa
# 每组5个独立试件；组内编号无配对关系。完整题干和数据见README。
import os
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf

if "__file__" in globals():
    project_root = Path(__file__).resolve().parents[1]
else:
    project_root = Path.cwd()
    if project_root.name == "代码":
        project_root = project_root.parent
sys.path.insert(0, str(project_root/"代码"))
from residplot import residplot

data = pd.read_csv(project_root/"数据"/"anova_example.csv")
print(data.head())
SAVE_OUTPUTS = os.getenv("RESIDPLOT_SAVE_OUTPUTS", "0") == "1"
output_dir = project_root/"运行结果"/"Python"
if SAVE_OUTPUTS:
    output_dir.mkdir(parents=True, exist_ok=True)


def run_case(case_data, case_name):
    # C(group)表示把组别作为分类变量；模型拟合的是各组均数。
    fit = smf.ols("value ~ C(group)", data=case_data).fit()
    result = residplot(fit)
    z = result["studentized_residuals"]
    summary = pd.DataFrame([{
        "n": len(z), "df_resid": fit.df_resid, "bandwidth": result["bandwidth"],
        "mean_studentized": z.mean(), "sd_studentized": z.std(ddof=1),
        "max_abs_studentized": np.abs(z).max()}])
    print(case_name)
    print(summary.to_string(index=False))
    if SAVE_OUTPUTS:
        residuals = case_data.copy()
        residuals["fitted"] = np.asarray(fit.fittedvalues)
        residuals["residual"] = np.asarray(fit.resid)
        residuals["leverage"] = fit.get_influence().hat_matrix_diag
        residuals["studentized"] = z
        tables = {"residuals": residuals, "histogram": pd.DataFrame(result["histogram"]),
                  "curves": pd.DataFrame(result["curves"]), "summary": summary}
        for name, table in tables.items():
            table.to_csv(output_dir/f"{case_name}_{name}.csv", index=False, float_format="%.17g")
        result["figure"].savefig(output_dir/f"{case_name}_residplot.png", dpi=160)
        plt.close(result["figure"])
    return fit, result

# %% 2. 原例：调用函数，观察3条分布图层和底部残差位置
base_fit, base_plot = run_case(data, "base")

# %% 3. 改动练习：只改A3第5个值，再拟合一个新模型
outlier_data = data.copy()
selected = (outlier_data["group"] == "A3") & (outlier_data["replicate"] == 5)
outlier_data.loc[selected, "value"] = 44.0
outlier_fit, outlier_plot = run_case(outlier_data, "outlier")
print("比较两张图：改动后哪一处残差位置最突出？最大绝对值由数据重新计算。")
if not SAVE_OUTPUTS:
    plt.show()
