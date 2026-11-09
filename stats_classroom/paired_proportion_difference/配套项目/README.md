# 配对样本率差：配套练习

本项目与文章、配图放在同一个主文件夹。三套主脚本与三套绘图入口均保留，按相同教学问题对应。数据均为教学构造。选择自己熟悉的一种语言，先核对20人例题，再比较不同样本量的分布。

## 先从基础分布开始

配合基础正文v2，先运行`distribution_basics`，只看“给定三种变化概率后，人数和率差怎样分布”。它不引入未知参数、尾部或排序算法。

- R：打开`配对率差.Rproj`，Source `代码/distribution_basics.R`，仅需基础R。
- Python：用Spyder打开项目文件夹，F5运行`代码/distribution_basics.py`。
- MATLAB：将Current Folder设为“代码”，运行`distribution_basics`，无需统计工具箱。

三套统一使用qplus=.05、qzero=.70、qminus=.25。先显示2人的五个D概率，再保存`基础分布_231组概率.csv`与`基础分布_率差概率.csv`，最后用`基础分布_汇总.csv`核对全概率、均值−.20和方差.013。20人观察组合(1,14,5)的概率为.0770154091。

正文的`PairedCI()`调用用于下一步：从观察人数估计未知总体率差。它与前面固定已知q的分布计算是不同问题。完整分析和内部排序演示按下面的入口继续运行。

## 继续运行完整分析


- **R**：双击`配对率差.Rproj`，打开`代码/paired_rate.R`并Source，得到精确分布、标准误及Wald区间，并自动运行`代码/plot_figures.R`生成正文三张统计图。随后安装`ExactCIdiff`、`PropCIs`，Source `代码/paired_intervals.R`，复现正文PairedCI与Tango区间，并运行三语言共同的精确区间示例。
- **Python**：用Spyder打开本项目文件夹，打开`代码/paired_rate.py`并按F5。主入口会自动调用Python绘图生成三张对照图；然后运行`代码/paired_exact.py`。如提示缺包，在Spyder使用的同一解释器环境中执行`python -m pip install -r requirements.txt`。
- **MATLAB**：将Current Folder设为`代码`，运行`paired_rate`，会自动调用MATLAB绘图生成三张对照图；随后运行`paired_exact`。使用基础MATLAB函数，无需统计工具箱。

输出分别位于`运行结果/R`、`运行结果/Python`和`运行结果/MATLAB`。`paired_exact.R`是R区间入口自动调用的辅助函数，保留在同一代码目录。**正文三张统计图统一由`代码/plot_figures.R`生成**，也可单独Source该文件；PNG与PDF放在项目上级的“文章配图”中。Python绘图脚本保留作对照，图像单独写入`运行结果/Python/图像`，不覆盖正文R图。

R包首次安装命令：

```r
install.packages(c("ExactCIdiff", "PropCIs"))
```

终端也可以在项目根目录执行：

```text
Rscript 代码/paired_rate.R
Rscript 代码/paired_intervals.R
Rscript 代码/plot_figures.R
python 代码/paired_rate.py
python 代码/paired_exact.py
python 代码/plot_figures.py
```

## 三套脚本怎样对应？

| 计算内容 | R | Python | MATLAB |
| --- | --- | --- | --- |
| 基础正文：同一已知概率模型的精确分布 | `代码/distribution_basics.R` | `代码/distribution_basics.py` | `代码/distribution_basics.m` |
| 个人记录、精确分布、均值方差、标准误与效应量 | `代码/paired_rate.R` | `代码/paired_rate.py` | `代码/paired_rate.m` |
| 三语言共同的固定D排序精确区间 | `代码/paired_exact.R`（由区间入口调用） | `代码/paired_exact.py` | `代码/paired_exact.m` |
| 图1、图2、图3及网格试算 | `代码/plot_figures.R` | `代码/plot_figures.py` | `代码/plot_figures.m` |
| 20人三类人数231组的逐项概率 | `代码/three_category_probabilities.R` | `代码/three_category_probabilities.py` | `代码/three_category_probabilities.m` |
| Wang排序第三位的试算机制 | `代码/wang_order_first_steps.R` | `代码/wang_order_first_steps.py` | `代码/wang_order_first_steps.m` |

三套使用同一批输入和同一统计定义，各自在本语言完成计算。Python的网格部分拆在`plot_grid_search.py`，由Python绘图入口自动调用。正文采用R图；Python和MATLAB图放在各语言运行结果的“图像”子目录。R的`paired_intervals.R`还调用原包`PairedCI()`，其Wang构造与共同练习的固定D排序方法分别说明。

## 输入与计算顺序

`数据/paired_20.csv`逐人记录id、before、after和z，z为后测减前测。默认20人有8人两次阳性、6人两次阴性、1人阴转阳、5人阳转阴。前测率13/20，后测率9/20，D=-0.20。三语言主入口读取这份个人数据计算主例，不将结果另行写死。

`数据/scenarios.csv`另列1000人的教学例题。两次阳性400人、两次阴性300人、阴转阳50人、阳转阴250人，边际率与主例相同。

给定概率的分布练习有四个场景：qplus=.05、qminus=.25、qzero=.70时n=2、20、1000，以及n=1000、qplus=.001、qminus=.002、qzero=.997的变化稀少模型。三类人数按`(Nplus,Nzero,Nminus)`排列时，联合概率为`(qplus,qzero,qminus)`。逐人累加时，数组从负总差值排到正总差值，因此卷积核按`[qminus,qzero,qplus]`排列。

输出`distribution_*.csv`给出k、D=k/n、精确概率及对应正态单格概率；`moments.csv`记录给定模型的均值、方差与累积概率最大差；`paired_examples.csv`给出配对Wald标准误、区间和边际Cohen's h。已知参数模型下的近似误差比较，不等于未知参数实际研究中任何区间覆盖率的保证。

## 两种精确区间如何对应？

**正文R结果采用Wang的归纳区间构造。** `ExactCIdiff::PairedCI(n12=1,t=14,n21=5,grid.one=60,grid.two=40)`得到点估计−0.20、95%区间[−0.44923,0.06303]。完整输出见`运行结果/R/PairedCI_实际打印.txt`。

**三语言共同练习采用固定D排序的无条件精确构造。** 给定候选总体率差theta，令相同结果概率s在[0,1-|theta|]内变化，另外两个概率为(1-s+theta)/2与(1-s-theta)/2。枚举合法三项计数，对D不小于或不大于观察值的组合累加概率，再对s取最大尾概率。双侧P取两倍较小单侧P，上限为1；95%区间反演两侧各2.5%的检验。

共同练习只按D排序；PairedCI的归纳构造还能区分相同D下的其他计数信息。因此两种区间分别报告，不能当作同一函数的三语言移植。共同方法默认输出P≈0.50344467，95%区间约[−0.61761990,0.27891486]，保存在各语言的`共同精确区间.csv`。该练习用于学习三项概率求和及区间反演，正文的小样本区间解释对应PairedCI。

Python与MATLAB均在本语言中计算，不调用R。统计方法采用有限三项分布，软件仍需数值最大化及端点求解。默认网格201点，局部峰细化并包含端点；可改为401点比较末位结果。R共同算法采用相同设置。

## 大样本与效应量

主例插件SE约0.11402，1000人例约0.01612。后者Wald近似95%区间约[−0.2316035,−0.1683965]。`PropCIs::scoreci.mp`当前源代码返回(c-b)/n方向，故本文后测减前测使用b=250、c=50，Tango近似区间约[−0.2319192,−0.1685890]。

效应量以净率差及百分点解释，并同时报告区间。Cohen's h仅使用前后边际率，不能替代配对标准误。两个人群大小例的边际率相同，因此h同为约−0.40486。

先尝试交换前后记录或精确函数的nplus与nminus，检查点估计反号、双侧P相同、区间反号后交换端点。计算规模增大时，精确区间耗时也随算法和观察结果变化。

## 方法来源

- ExactCIdiff作者说明（PairedCI的三项模型与归纳区间）：https://journal.r-project.org/articles/RJ-2013-026/index.html
- Tango score区间的R函数文档：https://search.r-project.org/CRAN/refmans/PropCIs/html/scoreci.mp.html

实际运行版本、核对范围与数值差异见`VALIDATION.md`。

## 图3：把网格搜索逐点算出来

运行`代码/plot_figures.R`。输入`数据/wang_lower_tail_n20.csv`为正文20人例在Wang归纳排序下、包含观察结果（1，14，5）的实际下限尾部，共141组三类计数，由当前ExactCIdiff的PairedCIone构造提取；按n_plus、n_zero、n_minus列存放。此图沿用该固定尾部，未替换为按D排序的共同精确算法。

脚本固定候选theta=−0.45，pT在[0,0.55]内变化，每试一个pT就将141个三项概率相加。先在全区间试11点，最高点为0.495；再在相邻点[0.44,0.55]之间试21点，最高点为0.4895。R的试算结果输出到`运行结果/R/wang_grid_*.csv`。Python对照脚本`代码/plot_grid_search.py`的CSV和图像放在`运行结果/Python`中。图保存在上级文章配图目录，同时提供PNG与PDF。

11与21是便于观察的教学网格，不是正文PairedCI调用的60与40，也不是函数每一步内部运算的完整跟踪。这里找的是固定theta下的最大尾部概率，不是把pT估计为网格最高点；theta的区间端点仍由外层搜索获得。

## 231组人数怎样计算概率？

三个主入口会自动调用本语言的`three_category_probabilities`脚本，也可以单独运行它。固定theta=−.45、pT=.44，三种概率为.055、.44、.505；先从观察人数(1,14,5)算起，分别显示安排种数232560、一种指定安排的概率1.8414377874e−8和人数结果的总概率.0042824477，再枚举全部231组。

打开各语言结果目录的`三类人数的231组概率.csv`：`n_plus/n_zero/n_minus`是三类人数；`arrangements`是能产生这组人数的人员安排种数；`single_arrangement_probability`是一种安排的概率；`probability`为相乘后的总概率；`in_lower_tail`说明该组是否属于主例的Wang下限尾部。`三类人数概率汇总.csv`给出总组数231、尾部141组、其余90组及其概率和。全部概率之和为1，尾部为.02375331196746，其余为.97624668803254。

## 下限尾部的计入标准

本例用97.5%单侧下限构造配合双侧95%区间。`数据/wang_lower_tail_n20.csv`按Wang次序存放观察组及其之前的141组，末行为观察组(1,14,5)。前141组无并列，因此本例第1至141位计入，第142位以后不计入；一般情形下与观察组同一层级的并列结果需一并计入。`数据/下限尾部判定示例.csv`另给三组实际位置135、141、152及判定。概率在给定theta与pT后另算，不能用单组概率是否大于观察组概率代替Wang次序。

## Wang排序怎样选第三位？

`wang_order_first_steps.R/.py/.m`只演示排序开头的真实比较：已排好的两组为(20,0,0)、(19,1,0)，待排A=(19,0,1)、B=(18,2,0)。各候选单独加入前两组，对这三个组的概率求和，再在合法pT范围寻找最大值；改变theta求最大尾概率=.025的最小根。A与B的临时下限约.50253447、.60859627，因此B排第三。尚未选中的A在下一轮需重新试算。

三个主入口会自动执行这一步，也可以单独运行相应脚本。输出`Wang排序第三位试算.csv`（两行候选及临时下限等）和`Wang排序候选尾部.csv`（每个候选各三组）。这三份脚本演示同一步选择，不是完整PairedCI移植；完整Wang区间仍由R包函数计算。
