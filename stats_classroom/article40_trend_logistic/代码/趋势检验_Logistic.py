# %% 01.1 准备与输入
# 统计课堂40：假设心电干扰案例，不是实测数据。UTF-8。
# 目标：得分/趋势检验、Logistic拟合、OR及区间、Wald/LRT和预测。
# 原题：假设200名不同患者各提供一段心电记录，每段有100个已标注QRS。
# 将记录随机分为4组，每段仅接受一种干扰强度，各患者记录相互独立。
# 固定算法处理全部记录；漏检超过5个QRS判为失败，不超过5个判为成功。
# 4组各50人，50 Hz干扰峰峰值为0/.1/.2/.3 mV，预评分为0:3，失败数为5/10/18/27。
# 问题：是否有有序趋势？趋势统计量怎样对应Logistic模型？
# 每0.1 mV的失败胜算比及双侧95% Wald区间是多少，各组预测失败概率是多少？
# 为什么趋势检验与模型默认P值不同？采用双侧alpha=.05。
# 在Spyder打开本项目文件夹，按节运行；依赖和保存开关见README。
import os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm
import matplotlib.pyplot as plt

try:
    CODE_DIR = Path(__file__).resolve().parent
except NameError:
    CODE_DIR = Path.cwd()
ROOT_DIR = CODE_DIR.parent
if not (ROOT_DIR / "数据" / "共享教学数据.csv").exists():
    ROOT_DIR = Path.cwd()
DATA_PATH = ROOT_DIR / "数据" / "共享教学数据.csv"
RESULT_DIR = ROOT_DIR / "运行结果" / "Python"
try:
    save_outputs
except NameError:
    save_outputs = os.environ.get("ARTICLE40_SAVE_OUTPUTS") == "1"

d = pd.read_csv(DATA_PATH)
print(d)
print("题目：检查共享数据，并准备失败事件、总数和观察失败率。")
required = ["score", "amplitude_mV", "failed", "n"]
assert all(col in d.columns for col in required)
assert np.isfinite(d[required].to_numpy(dtype=float)).all()
assert np.all(d["n"] > 0) and np.all(d["n"] == np.floor(d["n"]))
assert np.all(d["failed"] == np.floor(d["failed"]))
assert np.all((d["failed"] >= 0) & (d["failed"] <= d["n"]))
mV_per_score = 0.1
assert np.isfinite(mV_per_score) and mV_per_score > 0
# 教学单位参数：改为1即可把score按mV使用；括号必留，默认换算因子为1。
d["score"] = d["score"] * (0.1 / mV_per_score)
d["passed"] = d["n"] - d["failed"]
d["observed_rate"] = d["failed"] / d["n"]
N = int(d["n"].sum())

# %% 01.2 二项Logistic拟合、coef、OR及95% Wald区间
print("问题：先看模型结果：失败概率如何随评分变化，斜率的OR及95% Wald区间是多少？")
# statsmodels的GLM用每组两列计数作为endog；这与R的cbind(failed,passed)同义。
# add_constant为评分列增加截距列，四行仍对应四组。
# 第一列failed是建模事件，第二列passed是非事件。
X_fit = sm.add_constant(d[["score"]], has_constant="add")
fit = sm.GLM(d[["failed", "passed"]], X_fit, family=sm.families.Binomial()).fit(tol=1e-12, maxiter=100)
# 返回对象fit保存系数(params)、标准误(bse)、检验量和预测功能。
beta = fit.params["score"]
beta_se = fit.bse["score"]
z_wald = fit.tvalues["score"]
p_wald = fit.pvalues["score"]
OR = np.exp(beta)
# ppf输入累计概率，0.975使双侧区间两端各留2.5%。
# [-1,1]分别构造下限、上限；索引[0]和[1]读取这两个端点。
OR_CI = np.exp(beta + np.array([-1, 1]) * stats.norm.ppf(0.975) * beta_se)
# norm.ppf(0.975)是双侧95% Wald区间的正侧临界值；OR是每增加1分的优势比。
OR_per_0_1mV = np.exp(beta * 0.1 / mV_per_score)
print(fit.summary2().tables[1])
print(f"beta={beta:.10f}, SE={beta_se:.10f}, OR(每{mV_per_score:.6g} mV)={OR:.10f}, 95% Wald CI=[{OR_CI[0]:.10f}, {OR_CI[1]:.10f}]")
print(f"换算为每0.1 mV OR={OR_per_0_1mV:.10f}。")

# %% 01.3 四组预测与观察率+拟合曲线图
print("问题：用刚才的模型预测四组失败概率，和观察失败率如何对应？")
d["fitted_rate"] = fit.predict(X_fit)
print(d)
x_curve = np.linspace(d["amplitude_mV"].min(), d["amplitude_mV"].max(), 200)
# 模型用score拟合；图横轴换成mV时，预测输入要除以教学单位还原为score。
X_curve = sm.add_constant(pd.DataFrame({"score": x_curve / mV_per_score}), has_constant="add")
plt.rcParams.update({"font.size": 14})
plt.figure(figsize=(7.5, 5.2))
plt.scatter(d["amplitude_mV"], d["observed_rate"], color="black", label="Observed rate")
plt.plot(x_curve, fit.predict(X_curve), color="#1f77b4", label="Logistic fit")
plt.ylim(0, 0.65)
plt.xlim(-0.015, 0.315)
plt.xticks([0, 0.1, 0.2, 0.3])
plt.xlabel("50 Hz interference amplitude (mV, peak-to-peak)", fontsize=12)
plt.ylabel("Failure probability")
plt.legend(frameon=False)
plt.tight_layout()

# %% 01.4 趋势U/I_eff/Q/P与Wald、Rao、LRT比较
print("问题：趋势Q如何核对模型检验？这里I_eff作为得分U的方差估计使用。")
p0 = d["failed"].sum() / N
score_mean = (d["n"] * d["score"]).sum() / N
U = ((d["score"] - score_mean) * (d["failed"] - d["n"] * p0)).sum()
I_eff = p0 * (1-p0) * (d["n"] * (d["score"] - score_mean)**2).sum()
Q_manual = U**2 / I_eff
P_trend = stats.chi2.sf(Q_manual, 1)
# 全1列表示只有截距的共同率模型。
X_null = np.ones((len(d), 1))
null = sm.GLM(d[["failed", "passed"]], X_null, family=sm.families.Binomial()).fit(tol=1e-12, maxiter=100)
Q_Rao = Q_manual
P_Rao = stats.chi2.sf(Q_Rao, 1)
Q_LRT = null.deviance - fit.deviance
P_LRT = stats.chi2.sf(Q_LRT, 1)
logLik_null = null.llf
logLik_fit = fit.llf
Q_from_logLik = 2 * (logLik_fit - logLik_null)
assert np.isclose(Q_from_logLik, Q_LRT, rtol=1e-10)
print(pd.DataFrame([dict(N=N, p0=p0, score_mean=score_mean, U=U, I_eff=I_eff, Q=Q_manual, P=P_trend)]))
print(f"Wald P={p_wald:.6g}；手算score/Rao P={P_Rao:.6g}；LRT P={P_LRT:.6g}。")
print(f"LogLik_null={logLik_null:.10f}, LogLik_fit={logLik_fit:.10f}, 2倍差={Q_from_logLik:.10f}")
print("Python按相同公式手算score/Rao；R脚本另示范lmtest接口。三者均检验beta=0，有限样本P不必相同。")

# %% 01.5 整理结果与可选保存
print("报告：以下对象给出样本量、趋势检验、单位OR及拟合结果。")
result = pd.DataFrame([dict(N=N, p0=p0, score_mean=score_mean, U=U, I_eff=I_eff,
    Q_manual=Q_manual, Q_trend=Q_manual, Q_Rao=Q_Rao, P_trend=P_trend, P_Rao=P_Rao,
    Q_LRT=Q_LRT, P_LRT=P_LRT, logLik_null=logLik_null, logLik_fit=logLik_fit,
    Q_from_logLik=Q_from_logLik, beta=beta, beta_se=beta_se, OR=OR,
    OR_low=OR_CI[0], OR_high=OR_CI[1], Z_Wald=z_wald, P_Wald=p_wald,
    mV_per_score=mV_per_score)])
print(result.to_string(index=False))
print(f"可读报告：本次{N}段记录；未校正趋势Q={Q_manual:.10f}（df=1），双侧P={P_trend:.6g}；二项Logistic每{mV_per_score:.6g} mV的OR={OR:.10f}，双侧95% Wald区间=[{OR_CI[0]:.10f}, {OR_CI[1]:.10f}]。")
if save_outputs:
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    result.to_csv(RESULT_DIR / "Python_summary.csv", index=False)
    d.to_csv(RESULT_DIR / "Python_group_predictions.csv", index=False)
    plt.savefig(RESULT_DIR / "Python_observed_fitted.png", dpi=120)
plt.show()
