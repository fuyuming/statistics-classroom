# 统计课堂24：同一批人测了4次，方差分析怎样保留配对关系？

本项目对应8名受试者、A/B/C/D四种声音条件的反应率（%）例题。课堂脚本注明数据来自颜艳等《医学统计学》第5版例8-9。本项目沿用课堂表格，不包含真实受试者身份信息。四列表示声音条件，不能自行解释为四个时间点。

## 先运行哪一个文件？

三种语言任选一种即可，不需要依次运行三个。脚本都会读取原始CSV、计算原表和练习表，并把结果写入各自的目录。

- **RStudio：**打开根目录的`统计课堂24.Rproj`，再打开`代码/repeated_anova.R`，点击Source。需要R以及`rstatix`、`effectsize`扩展包；如缺包，先在Console运行`install.packages(c("rstatix", "effectsize"))`。项目脚本不会自行安装。R会生成QQ图PDF及数值表。
- **Spyder：**打开`代码/repeated_anova.py`，按F5。需要NumPy、pandas、SciPy、statsmodels、matplotlib；可在自己的Python环境中用`python -m pip install -r requirements.txt`安装。主脚本使用自身路径找数据，运行后调用`plot_repeated.py`生成四幅中文图。中文图需要系统中有PingFang SC、Noto Sans CJK SC或Arial Unicode MS字体。
- **MATLAB：**需要Statistics and Machine Learning Toolbox。打开`代码/repeated_anova.m`并运行；脚本根据自己的位置找数据，使用`fitrm`、`ranova`、`epsilon`、`mauchly`、`ttest`并核对手算结果。R2025a已实跑。MATLAB不提供本项目所用的原生Shapiro–Wilk函数，本版输出差值QQ图PDF。

命令行用户在项目根目录运行：

```text
Rscript 代码/repeated_anova.R
python 代码/repeated_anova.py
```

数据保留在`数据/四种声音条件反应率.csv`，不要移动或只复制入口脚本。根目录Rproj用于RStudio；`文章配图/`包含预先生成的图，尚未安装绘图库时也能查看。

## 按什么顺序看计算？

1. 看原始表和个体轨迹，确认同一人四次测量对应。
2. 把8×4宽表展开成长表，用标准函数拟合重复测量模型。
3. 手算全体、个人和条件均值，再算SS、MS、F与P。
4. 从协方差矩阵计算GG系数、调整自由度，核对标准函数。
5. 报告偏η²及经典非中心F双侧95%区间。
6. 对六对条件分别算同一个人的差值、配对t、Bonferroni P和同时区间。
7. 对给每人加不同常数的练习表重复分析，理解保留配对关系的意义。

基本条件：不同受试者相互独立；原表完整且平衡；经典模型关注受试者内条件对比的正态性及球形性。六组差值QQ图和Shapiro检验是辅助诊断，不能证明联合正态性。Mauchly在小样本下判断能力有限。实际研究须考虑条件顺序、缺失和携带效应；有缺失时不能直接照搬本例平方和算法。

## 核心结果

| 项目 | 原表结果 |
|---|---:|
| 总平方和 | 57.555 |
| 条件平方和 | 14.3875 |
| 受试者平方和 | 32.38 |
| 剩余平方和 | 10.7875 |
| 条件/误差自由度 | 3 / 21 |
| F | 9.33603708 |
| P | 0.0004033041 |
| GG ε | 0.6774634161 |
| GG自由度 | 2.032390248 / 14.226731738 |
| GG P | 0.0024937909 |
| 偏η² | 0.5714995035 |
| 经典模型95%CI | [0.2175051, 0.7350484] |

偏η²是条件SS除以“条件SS＋误差SS”，分母排除了受试者SS。区间用未GG校正的F模型反演非中心参数，不能称为GG校正区间。

配对比较方向统一为B−A、C−A、D−A、C−B、D−B、D−C；差值单位是百分点。每对用自己的差值SD，标准误为SD/√8，自由度7。六个原始P值乘6并截断到1；R再用`p.adjust()`核对。每对区间置信水平为1−0.05/6≈99.1667%，相应模型下整组六个区间同时覆盖概率至少95%。例如D−A平均差1.7625，Bonferroni P≈0.0287546，区间[0.18525,3.33975]。

## 输出文件怎样读？

`运行结果/R`、`Python`、`MATLAB`分别存放结果。`original_`表示原表；`offsets_0_14_`表示练习表。

- `*_summary.csv`：平方和、自由度、MS/F/P、GG与偏η²。部分语言同时包含两两比较；完整六对统一见`*_paired_t.csv`。
- `*_paired_t.csv`：差值均值、SD/SE、t、原始P、Bonferroni P、区间和实际置信水平。
- `*_long_data.csv`：32行长表；同一个受试者连续四行。
- `*_condition_means.csv`、`*_subject_means.csv`：条件和个人均值。
- `*_person_centered.csv`：每条记录减去本人的四次均值，图1下半部分使用这些数。
- `*_centered_residual.csv`：再减条件均值、加回总均值的逐格剩余量；其平方和等于10.7875。
- R的`*aov.txt`、`*_rstatix_GG.csv`、`*_mauchly.csv`和`*_effectsize.csv`是标准函数实际输出。rstatix表部分数值取整显示，完整精度见summary。
- Python的`*_AnovaRM.csv`是标准函数实际输出；AnovaRM不直接提供GG，脚本另行计算。
- MATLAB的`*_ranova.csv`、`*_epsilon.csv`、`*_mauchly.csv`是标准函数实际输出。
- R与MATLAB的`*_difference_QQ.pdf`：六对差值诊断图。正文四幅中文版数据图由Python绘图模块生成。

## 练习答案

把第i个人的四个数都加上2(i−1)，即依次加0、2、4、6、8、10、12、14。这是人为构造的教学变换，不是新增观察数据。

每一对条件的差值相减时，所加常数会抵消。因此条件间平均差、配对t、Bonferroni P/区间，以及重复测量F/P、GG和偏η²都不变。受试者SS从32.38变成549.18，总SS从57.555变成574.355。

把四列错误当作独立组时，分母仍混有人际差异。原表F≈3.11075、P≈0.042195；变换后F≈0.239806、P≈0.867828。两张表的条件差值完全相同，这种错误算法却因人际基线变化而改变结果。

## 官方方法说明

- R aov：https://stat.ethz.ch/R-manual/R-devel/library/stats/html/aov.html
- rstatix anova_test：https://rpkgs.datanovia.com/rstatix/reference/anova_test.html
- effectsize eta_squared：https://easystats.github.io/effectsize/reference/eta_squared.html
- Python AnovaRM：https://www.statsmodels.org/stable/generated/statsmodels.stats.anova.AnovaRM.html
- MATLAB ranova：https://www.mathworks.com/help/stats/repeatedmeasuresmodel.ranova.html

本项目提供重复测量方差分析的三语言脚本、示例数据与预先生成的运行结果，供读者下载后复现。
