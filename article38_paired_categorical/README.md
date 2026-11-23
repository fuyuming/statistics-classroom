# 统计课堂38：怎样比较重复测量的分类计数资料？

本项目对应McNemar、Cochran Q、Stuart–Maxwell、Bowker四方法及本文实际使用的事后比较、区间和条件枚举。两例均为8人的假设教学计数，无缺失；三套语言任选一种。

## 下载、打开、逐节运行

完整下载解压，进入本目录，再选一种语言。首次只显示结果，默认不写文件；如需保存，把脚本顶部`SAVE_RESULTS`改为真（R为`TRUE`、Python为`True`、MATLAB为`true`）。自动整文件验证也支持环境变量`ARTICLE38_SAVE=1`。各语言各自写入`results/R`、`results/Python`、`results/MATLAB`。

- **RStudio**：打开`article38_paired_categorical.Rproj`，工作目录即本目录。首次在Console执行`install.packages(c("DescTools", "ANSM5"))`。打开`scripts/article38_paired_categorical.R`，从输入与辅助函数开始按节选中运行，理解后可Source整文件。R不自动恢复/保存旧工作空间。Q与Stuart–Maxwell调用DescTools，Q的MC调用ANSM5，Bowker的χ²近似调用基础R，条件精确枚举与MC调用脚本自编函数。
- **Python / Spyder**：在Spyder建立/打开本文件夹项目，把工作目录设为本目录，打开`scripts/article38_paired_categorical.py`。从第一节依次用“运行当前单元”执行，最后可F5整文件运行。在Spyder实际使用的解释器环境安装`requirements.txt`列出的numpy、pandas、scipy、statsmodels。可以在对应解释器的终端执行`python -m pip install -r requirements.txt`。命令行备选：在本目录运行`python scripts/article38_paired_categorical.py`。分节模式不要求`__file__`存在。
- **MATLAB**：把本目录设为Current Folder，打开`scripts/article38_paired_categorical.m`，从第一节依次运行，或整文件Run。Current Folder决定相对数据位置，搜索路径用于定位函数。脚本末尾的本地辅助函数与前面的各节一同保存在完整脚本中，运行时请打开该脚本。本篇使用基础MATLAB的特殊函数，无需Statistics and Machine Learning Toolbox。

## 数据结构与目标

`data/binary_records.csv`是例1完整8人的逐人教学记录，包含`id,A,B,C`四列。`id`为受检者编号，一人一行；A/B/C为同一人在三种培养方案下的结果。三语主脚本直接读取这8行，并只取A/B/C构造Y；每种方案均有8份对应结果。

`data/binary_patterns.csv`为这8行按A/B/C组合汇总的频数对照表，字段为`pattern,A,B,C,count`。pattern保留8种模式文本，count依次为1,6,0,0,0,0,0,1。主脚本每次从逐人记录重新汇总，在Console和工作区显示频数表；统计输入来自逐人CSV。

正文表1展示全部8人的记录，表2汇总实际出现的3种模式（000、001、111）。资料为本篇假设教学数据。
例1分为两问：先比较三法总体检出率；有总体证据后，再对预设三对方案作两两比较，报告精确P与配对阳性率差。总体结果见`key_results`，条件精确核对见`q_conditional`，两两比较见`pairwise`；每对率差均以全部8人为分母。例2同样基于8人，按各节依次完成总体、逐类和路径问题。

`data/transition_records.csv`是例2完整8人的逐人转移记录，字段为`id,before,after`；三语主脚本直接读取并按A/B/C顺序汇总M。`data/transition_counts.csv`为完整9格频数对照表，不作为统计输入。这里A/B/C是无序菌群类别，行=前、列=后；每个人只进入一个联合格。

例2的四问依次比较总体构成、反向转移是否平衡，并估计逐类比例差和互转者中的方向概率。前两问见`key_results`，后两问见`categories`与`paths`；率差统一按第28天减第0天。

题干规定比较族与家族错误率控制目标。Holm用于调整三项检验P值；路径表保留Clopper–Pearson区间。

## 按问题核对输出

| 输出表 | 回答什么 |
|---|---|
| key_results | 例1 Q=12，df=2，χ²近似P=e⁻⁶；例2结果不变 |
| q_conditional | 固定每人行和，用U/W归并28个汇总配置；另标明729种逐人安排，概率和约1，`p_exact`=1/243 |
| pairwise | 三对精确McNemar+Holm；报告后一法减前一法的精确P与率差（分母8） |
| categories | 各类别“该类/其余类”的前后配对比较；分母8，报告率差与精确P |
| paths | 三条反向路径的条件方向概率、精确二项P+Holm、CP逐项95%区间；路径总数2/1/2 |
| enumeration | 固定反向路径总数2/1/2，18配置；概率和约1，`p_exact`=0.25 |
| exercise | 另一组8人转移表，SM=0、Bowker=6 |
| monte_carlo | 固定行和/反向路径总数的H0条件蒙特卡罗尾概率估计 |

### 字段与分母速查

| 字段/表 | 含义与分母/单位 |
|---|---|
| `pairwise`、`categories` 的 `rate_diff` | 后减前/后一法减前一法；两例分母均为8；保存为概率，乘100为百分点 |
| `categories` 的 `rate_diff` | 后减前的率差；分母为8，保存为概率，乘100为百分点 |
| `p_raw`、`holm` | 条件精确二项/McNemar原始P及Holm调整结果 |
| `paths` 的 `theta`、`lo95/hi95` | 反向路径中正向人数/互转总人数；路径总数2/1/2，CP区间估计互转者中的条件方向概率 |
| 脚本中的 `U/W` | `U`为行和1对象的唯一阳性位置计数，`W`为行和2对象的唯一阴性位置计数；各有A/B/C三个分量 |
| `q_conditional` 的 `probability_sum`、`p_exact` | 全部配置的概率和应约1；`p_exact`为Q至少达到观测值的条件右尾概率 |
| `enumeration` 的 `probability_sum`、`p_exact` | 固定反向路径总数2/1/2的配置概率和与Bowker统计量右尾P |
| `monte_carlo` 的字段 | `test,B,seed,exceedances,p_mc,mcse,exact_reference`；本轮B=100000、seed=20261123，`p_mc`为条件精确尾概率的模拟估计，`mcse=sqrt(p_mc*(1-p_mc)/B)` |

R和Python会并排显示Q公式值、函数值与绝对差。三语还显示各对`m=b+c`与`D`的关系；首次在Console查看，需要时再保存结果。

三语主入口默认执行一节条件蒙特卡罗，并显示/可保存`monte_carlo.csv`。Python/MATLAB提供同名自编函数`cochran_q_mc(Y,B,seed)`和`bowker_mc(M,B,seed)`；R的Q调用`ANSM5::cochran.q`，Bowker使用自编`bowker_mc(M,B,seed)`。以`B=100000`、`seed=20261123`调用后，将两行结果合并并补入`exact_reference`。Cochran Q模拟在每人内部随机重排0/1，保持每人行和；Bowker模拟固定三条反向路径总数，从独立Binomial(m, 1/2)生成方向。不同语言的随机序列不要求一致。`p_mc`是条件精确尾概率的Monte Carlo估计，`mcse`是其标准误，不是另一种χ²近似；不提供SM的MC。分类逐人记录与原有七张确定性结果表不变。

## 原生接口与自编计算

R调用`DescTools::CochranQTest`、`DescTools::StuartMaxwellTest`、`stats::mcnemar.test`、`binom.test`、`p.adjust`；控制台核对表覆盖Q和SM/Bowker的公式与原生函数。3×3的Bowker对称性检验由基础R `mcnemar.test`完成，连续性校正只作用于2×2。

Python调用`statsmodels.stats.contingency_tables.cochrans_q`、`SquareTable.homogeneity(method="stuart_maxwell")`、`symmetry(method="bowker")`、SciPy二项检验及statsmodels Holm；控制台核对表覆盖Q和SM/Bowker的公式与原生函数。`shift_zeros=False`保留输入的原始计数；本例保留零计数格，各反向格对合计均大于0。

MATLAB按正文公式自编全部检验及Holm，分布概率/分位点由`gammainc`、`gammaln`、`betaincinv`、`erfinv`计算。三语Q条件计算均为自编的三法多项计数枚举，针对给定行和的可交换模型。R的`bowker_exact3(M)`仅接受3×3且三对反向合计均>0的表，按固定反向合计的二项条件分布对Bowker统计量作右尾精确枚举，不使用蒙特卡罗；本例结果为18个配置、概率和1、`p_exact`=0.25。条件枚举对应Bowker的转移对称性零假设；Stuart–Maxwell的边际同质问题见正文。

基础R的`mcnemar.test`在反向两格合计为0时，计算中会出现0/0；`DescTools::StuartMaxwellTest`提供渐近检验接口。含零反向格对或奇异协方差的表，需要另选具有相应处理的实现。

本篇自编SM/Bowker函数用于给定非退化教学数据。反向格对合计为0、协方差奇异时会停止。此时需根据有效格对和矩阵秩确定自由度，再选择适合该表的计算方法。

R脚本末尾附有可选的ANSM5精确核对代码。`max.exact.perms`设为300000可运行8×3例1；默认上限100000，而本例穷举门槛为6^7=279936。若超过上限而输出`pval.mc`，不能视作穷举精确值。已实跑的`fit_exact$pval.exact`为0.00411522633744865。

## 第二次运行的练习

变式表描述另一组8人的前后分类资料。其三问是：各类人数及反向路径人数多少、SM和Bowker统计量多少、两个数为何不同。

脚本以独立M2构造这张表：M2=[[1,2,0],[0,1,2],[2,0,0]]。先预测，再查看行列边际3/3/2、SM=0与Bowker=6；该练习用来体会边际汇总与反向路径的区别。

主入口读取例1逐人记录与例2逐人转移记录，两个频数CSV仅供对照。文章配图另在本地交付的“文章配图”目录；运行本项目三语主入口无需绘图库或TeX。
