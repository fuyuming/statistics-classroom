# %% 1. 两份输入：20个独立教学观测用于正态性；三组模量用于方差齐性。
# 不把相互关联的拟合残差当作独立原始样本来套正态性检验。
import os
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from statsmodels.stats.diagnostic import normal_ad
import statsmodels.formula.api as smf

if "__file__" in globals():
    root = Path(__file__).resolve().parents[1]
else:
    root = Path.cwd()
    if root.name == "代码":
        root = root.parent
sys.path.insert(0,str(root/"代码"))
from studentized_qq import studentized_qq
save_outputs = os.getenv("RESIDPLOT_SAVE_OUTPUTS","0") == "1"
outdir = root/"运行结果"/"Python"
if save_outputs:
    outdir.mkdir(parents=True,exist_ok=True)
x = pd.read_csv(root/"数据"/"normality_example.csv")["value"].to_numpy()
d = pd.read_csv(root/"数据"/"anova_example.csv")
print("Normality sample:",x)
print(d.head())
alpha = 0.05

# %% 2. 正态性：两个H0都为来自某个正态总体，均数和方差未知
sw = stats.shapiro(x)
a2, ad_p = normal_ad(x)
print("Shapiro-Wilk W, P:",sw.statistic,sw.pvalue)
print("Anderson-Darling A2, P:",a2,ad_p)
normal_figure, normal_axes = plt.subplots(figsize=(7,5))
stats.probplot(x,dist="norm",plot=normal_axes)
normal_axes.set_title("Independent sample: normal Q-Q plot")
normal_figure.tight_layout()
if save_outputs:
    normal_figure.savefig(outdir/"normal_qq.png",dpi=160)
    plt.close(normal_figure)

# %% 3. 方差齐性：两者H0均为各组总体方差相等
samples = [d.loc[d["group"]==g,"value"].to_numpy() for g in ["A1","A2","A3"]]
bartlett = stats.bartlett(*samples)
# 明确median，与R car::leveneTest(center=median)的Brown-Forsythe修正对应。
levene = stats.levene(*samples,center="median")
methods = ["Shapiro-Wilk","Anderson-Darling","Bartlett","Levene-median"]
statistics = [sw.statistic,a2,bartlett.statistic,levene.statistic]
pvalues = [sw.pvalue,ad_p,bartlett.pvalue,levene.pvalue]
tests = pd.DataFrame({"method":methods,"sample":["independent20","independent20","groups15","groups15"],
    "statistic":statistics,"p_value":pvalues,"alpha":alpha,
    "reject":[p<alpha for p in pvalues]})
print(tests.to_string(index=False))

# %% 4. 模型Q-Q：对应car模型方法，使用外部学生化残差和t(df_resid-1)
fit = smf.ols("value ~ C(group)",data=d).fit()
qq = studentized_qq(fit)
fig, ax = plt.subplots(figsize=(7,5))
ax.boxplot(samples,tick_labels=["A1","A2","A3"],showfliers=False)
for i,values in enumerate(samples,start=1):
    offsets = np.linspace(-0.07,0.07,len(values))
    ax.scatter(i+offsets,values,s=25,color="0.25",zorder=3)
ax.set(xlabel="Group",ylabel="Compressive modulus (MPa)",title="Group spread and observations")
fig.tight_layout()
if save_outputs:
    tests.to_csv(outdir/"diagnostics_tests.csv",index=False,float_format="%.17g")
    qq["points"].to_csv(outdir/"studentized_qq.csv",index=False,float_format="%.17g")
    qq["figure"].savefig(outdir/"studentized_qq.png",dpi=160)
    fig.savefig(outdir/"group_boxplot.png",dpi=160)
    plt.close("all")
else:
    plt.show()
