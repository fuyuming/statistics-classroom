# 统计课堂43｜方差分析的正态性、方差齐性：指什么，一定要检验吗？

本项目是三语言（R/Python/MATLAB）的单一主例场景A教学资料。A的总体均值 μ=(10,20,30)、共同扰动标准差 σ=2，单位为 MPa；数据由 `set.seed(20261124)` 生成，15行、三组各5个观测。CSV保留 `scenario, group, replicate, value, mu, z` 字段，`mu` 与 `z` 用于展示模拟真值，不是实测数据。

进入“配套项目”目录后，可运行 `代码/anova_residuals.R`、`anova_residuals.py` 或 `anova_residuals.m`。RStudio打开 `ANOVA条件与残差.Rproj`；Python在Spyder将工作目录设为项目根目录，按 `# %%` 分节或整文件运行；MATLAB将Current Folder设为项目根目录，按 `%%` 分节或整文件运行。主入口均只计算场景A：组摘要、手算ANOVA及原生核对、残差、Q-Q、按组Shapiro/AD状态、Levene（中位数中心）与Bartlett，共6类结果。

设置 `SAVE_RESULTS=true` 可将6个统一CSV写入对应 `结果/R`、`结果/Python` 或 `结果/MATLAB`；默认不保存。`代码/生成模拟资料.R` 仅按固定种子重建同一份A数据，不会自动覆盖主入口数据。

依赖：R需 `car`、`ggplot2`、`nortest`（可运行 `install.packages(c("car", "ggplot2", "nortest"))`）；Python依赖 `requirements.txt` 中的 pandas、NumPy、SciPy、statsmodels、Matplotlib、seaborn；MATLAB需 Statistics and Machine Learning Toolbox。场景A每组n=5，Anderson–Darling按本练习约定记录为 `not_run_n_lt_8`。三张诊断图为残差-拟合值、原始残差Q-Q和residplot分布图，在RStudio/Spyder/MATLAB中可直接查看；无桌面环境设置 `SHOW_PLOTS=false`。辅助 `residplot.*` 与 `shapiro_wilk.m` 为配套函数。

R/Python/MATLAB结果字段与行序一致；Q-Q表为A场景15行，理论位置使用 `(rank-0.5)/15`。是否需要正式前提检验，要结合已有依据、样本量和诊断目的判断；P值不等于证明前提成立。
