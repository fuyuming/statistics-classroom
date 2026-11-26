# 统计课堂42：诊断检验与残差图的三语练习

本项目提供R、Python、MATLAB的残差分布图、Q–Q图和四种常见检验。用自己熟悉的一种语言跑通即可。

## 先从哪里开始

下载ZIP并解压，进入`统计课堂42_诊断检验与图形`目录，先读本文件。

本篇GitHub目录：https://github.com/fuyuming/statistics-classroom/tree/main/stats_classroom/article42_diagnostic_tests 。如果下载整个仓库ZIP，解压后进入`statistics-classroom-main/stats_classroom/article42_diagnostic_tests`，再按本文件运行。

- **R/RStudio**：打开`统计课堂42.Rproj`，先按节运行`代码/demo_residplot.R`，再运行`代码/demo_diagnostics.R`。示例依赖`car`、`nortest`；初次使用可运行`install.packages(c("car", "nortest"))`。
- **Python/Spyder**：将解压目录作为项目文件夹，打开`代码/demo_residplot.py`，按`# %%`分节运行或按F5；再打开`代码/demo_diagnostics.py`。依赖见`requirements.txt`，在Spyder实际使用的解释器环境安装：`python -m pip install -r requirements.txt`。
- **MATLAB**：把解压目录设为Current Folder，打开`代码/demo_residplot.m`按节运行，再运行`代码/demo_diagnostics.m`。脚本会加入同目录函数路径，需要Statistics and Machine Learning Toolbox。

默认在IDE显示图，不写出文件。需要导出数表和图片时，先把环境变量`RESIDPLOT_SAVE_OUTPUTS`设为`1`；结果分别写到`运行结果/R`、`运行结果/Python`、`运行结果/MATLAB`。示例不改变原CSV。

命令行备选（从项目根目录运行）：

```text
Rscript 代码/demo_residplot.R
Rscript 代码/demo_diagnostics.R
python 代码/demo_residplot.py
python 代码/demo_diagnostics.py
matlab -batch "run('代码/demo_residplot.m'); run('代码/demo_diagnostics.m')"
```

## 例题一：三种支架的残差分布和组间离散程度

假设比较A1、A2、A3三种组织工程支架的压缩模量。每种设计独立制备5个试件，每件只测量一次，单位MPa；各组是不同试件，同一行只是组内编号，没有配对关系。

这是已有课堂正态模型生成、保留两位小数的教学数据，不是真实实验记录。完整输入如下，也保存在`数据/anova_example.csv`。

| 组内编号 | A1 | A2 | A3 |
|---|---|---|---|
| 1 | 8.97 | 20.02 | 27.56 |
| 2 | 9.95 | 19.59 | 28.03 |
| 3 | 13.78 | 21.47 | 28.66 |
| 4 | 8.62 | 19.94 | 30.83 |
| 5 | 5.39 | 19.43 | 31.48 |

本例要做三件事：拟合`value ~ group`的组均数模型；观察外部学生化残差的分布及Q–Q图；在α＝0.05下，用Bartlett和中位数中心Levene检查各组总体方差是否相等。这里不开展组间均数显著性比较。

`residplot`返回学生化残差、直方图数表、两条密度曲线及带宽。原例最大绝对外部学生化残差约为3.19778209。练习只将A3第5个值由31.48改为44.00，另拟合一个模型；该值约增至6.71794128。观察图中哪一处位置最突出，以及蓝色、红色曲线怎样变化。示例始终按输入重新计算这些数。

## 例题二：20个独立观测的正态性检验

假设同一类支架独立制备20个试件，每件只测一次模量（MPa）。下面的固定教学样本来自N(10, 2²)模型，经NumPy default_rng(20261126)生成并保留两位小数；见`数据/normality_example.csv`。

9.25、9.49、11.88、12.78、7.17、10.14、13.76、11.77、6.58、12.22、8.77、8.11、10.80、10.17、9.32、12.80、11.32、9.58、10.38、8.13。

在α＝0.05下分别进行Shapiro–Wilk与Anderson–Darling正态性检验，并画正态Q–Q图。两个H₀都为总体属于某个均数、方差未知的正态分布。

正态检验使用这组独立原始观测；模型残差另用于诊断绘图。R的`nortest::ad.test`要求至少8个观测，因此没有对例题一每组仅5个数强行调用AD。

## 函数对应

先看各自的H₀和拒绝方向，再核对统计量；软件采用不同近似校准时，P值可能略有差异。

| 用途 | R | Python | MATLAB |
|---|---|---|---|
| Shapiro–Wilk | `shapiro.test(x)` | `scipy.stats.shapiro(x)` | 本课`shapiro_wilk(x)` |
| AD正态性 | `nortest::ad.test(x)` | `statsmodels.stats.diagnostic.normal_ad(x)` | `adtest(x,'Distribution','norm')` |
| 中位数中心Levene | `car::leveneTest(value ~ group, center=median)` | `scipy.stats.levene(..., center="median")` | `vartestn(...,'TestType','BrownForsythe')` |
| Bartlett | `bartlett.test(value ~ group)` | `scipy.stats.bartlett(...)` | `vartestn(...,'TestType','Bartlett')` |
| 原始样本正态Q–Q图 | `car::qqPlot(x, distribution="norm")` | `scipy.stats.probplot(x, dist="norm")` | `qqplot(x)` |
| 模型学生化残差Q–Q图 | 本课`studentized_qq(fit)`调用`car::qqPlot` | 本课`studentized_qq(fit)` | 本课`studentized_qq(fit)` |
| 残差分布图 | 本课`residplot(fit)` | 本课`residplot(fit)` | 本课`residplot(fit)` |

MATLAB的`shapiro_wilk.m`为本课按Royston/SciPy算法移植的函数，接受3—5000个完整有限观测；相应BSD许可随代码提供。它不是用Shapiro–Francia替换Shapiro–Wilk。

**不要把AD的P值当作三语必须逐位相同的值。** R的nortest和Python的normal_ad使用同一类Stephens修正近似；MATLAB原生adtest在参数未知时按临界值表插值。本次A²约为0.1449997844，三语一致；P分别约为R/Python 0.96148283、MATLAB 0.97342092。三者均未拒绝H₀。

各统计量怎样构造、参照什么分布，见[统计量与参照分布](统计量与参照分布.md)。

## 残差图到底画的是什么

`residplot`面向无权重、满秩的普通OLS模型。R使用`rstudent(fit)`，Python使用`fit.get_influence().resid_studentized_external`，MATLAB使用`fit.Residuals.Studentized`；都是删除当前观测后估计误差方差的外部学生化残差。

- R：`result <- residplot(fit, nbreaks=10, draw=TRUE)`。
- Python：`result = residplot(fit, nbreaks=10, draw=True)`。
- MATLAB：`result = residplot(fit, 10, true)`。

`draw=FALSE/False/false`只计算，不创建图。返回字段为`studentized_residuals`、`histogram`、`curves`、`bandwidth`；Python在绘图时另返回`figure`和`axes`。

为便于核对，三语统一为10个等宽区间、n−1样本SD、type=7四分位数、nrd0规则带宽、512个相同网格点和显式Gaussian核密度；`nbreaks`可改。灰色柱是密度直方图，蓝线是按学生化残差均数、SD拟合的正态形状参照，红虚线是核密度，棕色短线是残差位置。

这是按《R in Action》思路整理的课堂改写：使用确定性位置短线及统一分箱，没有原书的随机jitter；KDE也统一了计算细节。蓝线作形状对照，不表示学生化残差具有精确正态分布。

`car::qqPlot(fit)`的模型方法默认采用学生化残差的t参照；本例误差自由度12，参照自由度为11。三语模型Q–Q图统一采用四分位参考线、关闭模拟包络，从而可逐点比较坐标。原始向量的正态Q–Q图仍以正态分布为参照。

原始样本Q–Q图和箱线图调用各软件的原生函数。它们的默认绘图位置、参考线或四分位算法可能不同，小样本箱线图的箱体和须尤其可能不同；叠加的观测点始终来自同一份数据。逐点统一的是模型学生化Q–Q图，以及`residplot`的残差、分箱和密度曲线数表。

## 已保存的结果

- `base/outlier_residplot.png`及相应的残差、直方图、曲线、摘要CSV。
- `diagnostics_tests.csv`：四种检验的统计量、P值、α和拒绝与否。
- `normal_qq.png`：20个独立观测的正态Q–Q图。
- `studentized_qq.png/.csv`：组均数模型的学生化残差Q–Q图及坐标。
- `group_boxplot.png`：按组箱线图与全部观测点。

改变单个观测、分箱数或数据单位后，先运行模型拟合，再调用绘图函数；输出将随输入重算。运行核对范围见[VALIDATION.md](VALIDATION.md)。

## 原始资料

- [Kabacoff第二版第8章源码，Listing8.6](https://github.com/kabacoff/RiA2/blob/master/Ch08%20Regression.R)
- [R rstudent说明](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/influence.measures.html)
- [statsmodels外部学生化残差](https://www.statsmodels.org/stable/generated/statsmodels.stats.outliers_influence.OLSInfluence.resid_studentized_external.html)
- [MATLAB LinearModel.Residuals](https://www.mathworks.com/help/stats/linearmodel.html)
- [R nortest AD说明](https://search.r-project.org/CRAN/refmans/nortest/html/ad.test.html)
- [MATLAB AD校准说明](https://www.mathworks.com/help/stats/adtest.html)
- [car Q–Q图说明](https://cran.r-project.org/web/packages/car/refman/car.html)
