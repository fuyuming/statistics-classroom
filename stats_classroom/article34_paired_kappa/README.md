# 统计课堂34：配对率差与一致性

58名可疑系统性红斑狼疮患者，各接受免疫荧光法A与乳胶凝集法B检测抗核抗体。本项目使用颜艳、王彤《医学统计学》第5版例7-3、表7-3的四格联合汇总计数：`[[11,12],[2,33]]`。行是A，列是B，两轴均按阳性、阴性排序；每名患者恰好落入一个格子，不把CSV的四行当成四个人。

同一例题依次计算率差、观察一致率、机会一致率、κ、一般多项delta标准误与95%Wald区间、κ=0的右尾近似检验。配对McNemar的条件二项精确P用于衔接前文。

## 跟着正文逐步推导：03脚本

第三份脚本按正文的顺序计算：每人的X/Y/C三个0、1记录→三个样本均数→κ的三个系数→四格W值→m₁、m₂→方差/58→SE→1.96 Wald区间。完整推导已经写在公众号正文，脚本用于逐步代入及函数核对。

- RStudio：打开`scripts/03_kappa_stepwise.R`，从03.1到03.6依次运行。
- Spyder：打开`scripts/03_kappa_stepwise.py`，从第一个`# %%`单元顺序运行；不依赖`__file__`。
- MATLAB：打开`scripts/m03_kappa_stepwise.m`，从第一个`%%`节顺序运行。

三个入口都可在新会话独立执行，工作目录设为本项目根目录。依赖沿用下文：R为vcd，Python为NumPy/pandas/statsmodels，MATLAB调用现有helpers/kappa_summary.m及其所需工具箱。

正文的p、q、r是未知总体概率；代码用`p_hat/q_hat/r_hat/e_hat`表示由四格人数得到的代入值。手算的Wald临界值为1.96；完整函数核对κ与一般SE。输出`stepwise_result`依次含18列：`p_hat,q_hat,r_hat,e_hat,kappa,cr,cp,cq,Wa,Wb,Wc,Wd,m1,m2,variance,se,ci_low,ci_high`。

先核对三个系数约2.257718、0.678831、0.254562，m₁≈2.039000045、m₂≈4.929051347；方差≈0.0133022442、SE≈0.1153353554，95%Wald区间四位小数为[0.2290,0.6811]。概率矩的分母为58；它是经验概率代入计算，不乘n/(n−1)。

默认只显示结果。按下文设置`ARTICLE34_SAVE_RESULTS=1`时，另保存`results/R_stepwise.csv`、`Python_stepwise.csv`或`MATLAB_stepwise.csv`。01保留另一条等价的四格公式及边际枚举，02计算成对Bootstrap区间。

## 先打开项目，再选择一种语言

保留`data/`、`scripts/`、`helpers/`的相对位置。首次建议从第一个准备节开始，逐节运行并查看命名对象；中间小节使用前面已经创建的变量。

### R：RStudio

1. 打开`article34_paired_kappa.Rproj`，工作目录就是项目根目录。
2. 如未安装vcd，在Console运行`install.packages("vcd")`。
3. 打开`scripts/01_paired_kappa.R`，从“输入”节起选择完整语句或代码块，按Run顺序执行。循环要连同大括号一起运行。
4. 手算后查看`point_estimates`、`interval_result`、`function_comparison`和`null_test`。vcd提供κ和一般SE/区间；κ=0的Rosner右尾检验另用零假设SE计算。

整文件命令行入口：`Rscript scripts/01_paired_kappa.R`。

### Python：Spyder

1. 将项目文件夹设为Spyder工作目录；解释器需有requirements.txt列出的NumPy、pandas、SciPy及statsmodels。
2. 如需安装，在对应解释器环境运行`python -m pip install -r requirements.txt`。
3. 打开`scripts/01_paired_kappa.py`，从第一个`# %%`单元起，使用Run current cell顺序运行；也可用F5运行整文件。控制台分节不依赖`__file__`。
4. `fit.std_kappa`是一般SE，`fit.std_kappa0`是零假设SE；按事先的κ>0方向读取`pvalue_one_sided`。`function_comparison`把手算与函数字段并排。

整文件命令行入口：`python scripts/01_paired_kappa.py`。

### MATLAB

1. Current Folder设为本项目根目录，需要Statistics and Machine Learning Toolbox提供`normcdf`、`norminv`、`binocdf`。
2. 在编辑器打开`scripts/m01_paired_kappa.m`，从“输入”节开始使用Run Section顺序运行；循环/条件块要完整执行。
3. “调用函数”节加载`helpers/`，调用本项目自编的`kappa_summary.m`。它是同一统计方法的完整公式实现，采用等价的展开方差，不称作MATLAB原生κ接口。
4. 查看`function_comparison`中manual与function_value的对应关系；当前数据两列数值应相符。

整文件可在项目根目录运行`run('scripts/m01_paired_kappa.m')`。该脚本在必要时定位回项目根目录，按节运行时保留已正确设置的目录。

## 02：Bootstrap运行顺序与方法

在完成01后，打开对应的02脚本，从`02.1`开始按节顺序运行；02也可在新会话独立运行。三语入口分别是`scripts/02_kappa_bootstrap.R`、`02_kappa_bootstrap.py`和`m02_kappa_bootstrap.m`。R需安装`boot`（脚本实际调用`boot::boot`和`boot::boot.ci`）；Python沿用`requirements.txt`中的NumPy、pandas、SciPy；MATLAB需Statistics and Machine Learning Toolbox。

`data/bootstrap_counts.csv`是按四种联合结果类别重抽样得到的50000行四格计数，每行合计58；它来自R `boot::boot`、`RNGkind("Mersenne-Twister","Inversion","Rejection")`及seed=20261120的真实生成，不是虚构的逐人病历。R脚本会再次真实调用`boot::boot`并核对共享计数。Python/MATLAB读取这份共享表，核对相同输入下的算法；其中seed记录R的生成条件，不控制本版新一轮随机抽样。百分位与BCa的覆盖属于近似推断；有限重抽次数另外带来Monte Carlo计算误差。BCa使用严格小于的`z0`、58个delete-one jackknife影响值及加速量；Python/MATLAB按`boot::norm.inter`在正态分位尺度插值，不使用NumPy默认线性分位插值。

`z0`由“小于原κ的重抽比例”换成正态分位点；`acceleration`由58个删一κ的变化计算；二者合起来把2.5%和97.5%调整成`alpha_low/high`，再读取对应重抽端点。R的`boot.ci()`负责整个区间计算；Python的`norm_inter_r()`、MATLAB的`kappa_bootstrap_norm_inter()`只负责给定累计概率的端点插值，BCa分位水平由脚本前一步计算。

02最终结果对象统一包含17列：`B`重复次数；`seed`随机种子；`invalid`不可定义κ样本数；`kappa`原样本κ；`bootstrap_mean`κ*均值；`bootstrap_bias`均值−原样本κ；`bootstrap_variance`κ*样本方差（B−1）；`bootstrap_se`κ*样本标准差；`bias_probability`严格小于原样本κ的比例；`z0`偏差校正分位点；`acceleration`jackknife加速量；`alpha_low/alpha_high`BCa调整后的累计概率；`percentile_low/percentile_high`百分位区间端点；`bca_low/bca_high`BCa区间端点。目标值约为百分位[0.221476510067, 0.675776397516]、BCa[0.219258992519, 0.672521957340]，z0=0.013135109852，acceleration=−0.012863109946。

02默认只显示结果、不写文件；需保存时，按后文设置`ARTICLE34_SAVE_RESULTS=1`，输出`results/R_bootstrap.csv`、`Python_bootstrap.csv`或`MATLAB_bootstrap.csv`。

## 先核对这些结果

- n=58；A阳性23/58、B阳性13/58；B−A率差−0.172414，即−17.24个百分点。
- Pₒ≈0.758621、Pₑ≈0.557075、κ≈0.455034。
- 一般SE≈0.115335；双侧95%Wald区间约[0.228980,0.681087]。
- κ=0的SE₀≈0.120950，Z≈3.762165，κ>0的单侧近似P≈0.00008422。
- 不一致人数m=14；McNemar条件二项双侧精确P≈0.012939。

一般标准误推导见[method_notes.md](method_notes.md)。总体κ区间允许新样本的边际变化；固定观察边际的枚举范围不是已知的总体参数界限，不能据此裁剪一般区间。

## 经验分级

Landis–Koch原表六档为：<0 poor、0.00—0.20 slight、0.21—0.40 fair、0.41—0.60 moderate、0.61—0.80 substantial、0.81—1.00 almost perfect；原表按两位小数展示。Rosner三级为0≤κ<0.4 marginal、0.4≤κ≤0.75 good、κ>0.75 excellent。本例分别为moderate与good；这是经验描述，应连同区间、用途解读。

01中另外保留了固定边际枚举作为可选补充，不属于本篇正文练习。正文的推导与手算按03运行，Bootstrap按02运行。

## 可选保存结果

首次练习默认只显示结果，不写文件。需要保存时，在相应IDE中设置一次，再运行最后的保存小节：

- R：`Sys.setenv(ARTICLE34_SAVE_RESULTS="1")`
- Python：`os.environ["ARTICLE34_SAVE_RESULTS"]="1"`
- MATLAB：`setenv('ARTICLE34_SAVE_RESULTS','1')`

`results/`中的三语summary.csv保存相同的18个统计字段，bounds.csv保存`a,b,c,d,rd,po,kappa`枚举表。参考CSV可以与自己的计算逐项比较；三种语言分别写入各自前缀，互不覆盖。
