# 统计课堂38：怎样比较重复测量的分类计数资料？

本项目对应McNemar、Cochran Q、Stuart–Maxwell、Bowker四方法及本文实际使用的事后比较、区间和条件枚举。两例都是假设教学计数，无缺失；三套语言任选一种。

## 下载、打开、逐节运行

完整下载解压，进入本目录，再选一种语言。首次只显示结果，默认不写文件；如需保存，把脚本顶部`SAVE_RESULTS`改为真（R为`TRUE`、Python为`True`、MATLAB为`true`）。自动整文件验证也支持环境变量`ARTICLE38_SAVE=1`。各语言各自写入`results/R`、`results/Python`、`results/MATLAB`。

- **RStudio**：打开`article38_paired_categorical.Rproj`，工作目录即本目录。首次在Console执行`install.packages("DescTools")`。打开`scripts/article38_paired_categorical.R`，从输入与辅助函数开始按节选中运行，理解后可Source整文件。R不自动恢复/保存旧工作空间。Q与Stuart–Maxwell调用DescTools，Bowker调用基础R。
- **Python / Spyder**：在Spyder建立/打开本文件夹项目，把工作目录设为本目录，打开`scripts/article38_paired_categorical.py`。从第一节依次用“运行当前单元”执行，最后可F5整文件运行。在Spyder实际使用的解释器环境安装`requirements.txt`列出的numpy、pandas、scipy、statsmodels。可以在对应解释器的终端执行`python -m pip install -r requirements.txt`。命令行备选：在本目录运行`python scripts/article38_paired_categorical.py`。分节模式不要求`__file__`存在。
- **MATLAB**：把本目录设为Current Folder，打开`scripts/article38_paired_categorical.m`，从第一节依次运行，或整文件Run。Current Folder决定相对数据位置，搜索路径用于定位函数。脚本末尾的本地辅助函数与前面的各节一同保存在完整脚本中，运行时请打开该脚本。本篇使用基础MATLAB的特殊函数，无需Statistics and Machine Learning Toolbox。

## 数据结构与目标

`data/binary_records.csv`是例1完整8人的逐人教学记录，包含`id,A,B,C`四列。`id`为受检者编号，一人一行；A/B/C为同一人在三种培养方案下的结果。三语主脚本直接读取这8行，并只取A/B/C构造Y；每种方案均有8份对应结果。

`data/binary_patterns.csv`为这8行按A/B/C组合汇总的频数对照表，字段为`pattern,A,B,C,count`。pattern保留8种模式文本，count依次为1,6,0,0,0,0,0,1。主脚本每次从逐人记录重新汇总，在Console和工作区显示频数表；统计输入来自逐人CSV。

正文表1展示全部8人的记录，表2汇总实际出现的3种模式（000、001、111）。资料为本篇假设教学数据。
例1分为两问：先比较三法总体检出率；有总体证据后，再对预设三对方案作两两比较，报告精确P与配对阳性率差。总体结果见`key_results`，条件精确核对见`q_conditional`，两两比较见`pairwise`；每对率差均以全部8人为分母。区间练习集中在例2的`categories`表。

`data/transition_counts.csv`包含`before,after,count`三列。例2假设随访144名成年居民，第0天和第28天各检测一次。为本例约定：甲、乙、丙三组指示菌群中，相对丰度最高者分别对应A/B/C型；并列按甲、乙、丙顺序归类，前后规则相同。这里A/B/C是无序菌群类别，行=第0天、列=第28天。每个人只进入一个联合格，前后两次记录共同组成配对输入。

例2的四问依次比较总体构成、反向转移是否平衡，并估计逐类比例差和互转者中的方向概率。前两问见`key_results`，后两问见`categories`与`paths`；率差统一按第28天减第0天。

题干规定比较族、家族错误率控制目标与区间覆盖目标。Holm用于调整三项检验P值；区间练习集中在例2 `categories`表。

## 按问题核对输出

| 输出表 | 回答什么 |
|---|---|
| key_results | 例1 Q=12，df=2，χ²近似P=e⁻⁶；例2结果不变 |
| q_conditional | 固定每人行和，用U/W归并28个汇总配置；另标明729种逐人安排，概率和约1，`p_exact`=1/243 |
| pairwise | 三对精确McNemar+Holm；报告后一法减前一法的精确P与率差（分母8） |
| categories | 各类别“该类/其余类”的前后配对比较；分母144，后减前，含逐项及同时区间列 |
| paths | 三条反向路径的条件方向概率、精确二项P+Holm、CP逐项95%区间；分母为路径互转者 |
| enumeration | 固定反向路径总数21/18/18，7942配置；概率和约1，`p_exact`为Bowker统计量条件右尾P≈1.00269005809e−6 |
| exercise | 更改循环转移例后SM=0、Bowker=24 |

### 字段与分母速查

| 字段/表 | 含义与分母/单位 |
|---|---|
| `pairwise`、`categories` 的 `rate_diff` | 后减前/后一法减前一法；分母分别为全体配对对象（例1=8、例2=144）；保存为概率，乘100为百分点 |
| `categories` 的 `se` | 由144份个人差值计算的率差标准误；单位与率差相同 |
| `categories` 的 `lo95/hi95`、`lo_bonf/hi_bonf` | 例2率差区间；前者为逐项95%，后者为三项Bonferroni同时区间，仍用全体分母 |
| `p_raw`、`holm` | 条件精确二项/McNemar原始P及本项目选用的Holm调整结果；率差估计和区间按相应公式另算 |
| `paths` 的 `theta`、`lo95/hi95` | 反向路径中正向人数/互转总人数；分母为三条路径各自互转者21/18/18；CP区间估计互转者中的条件方向概率 |
| 脚本中的 `U/W` | `U`为行和1对象的唯一阳性位置计数，`W`为行和2对象的唯一阴性位置计数；各有A/B/C三个分量 |
| `q_conditional` 的 `probability_sum`、`p_exact` | 全部配置的概率和应约1；`p_exact`为Q至少达到观测值的条件右尾概率 |
| `enumeration` 的 `probability_sum`、`p_exact` | 固定反向路径总数21/18/18的配置概率和与Bowker统计量右尾P |

R和Python会并排显示Q公式值、函数值与绝对差。三语还显示各对`m=b+c`与`D`的关系；首次在Console查看，需要时再保存结果。

## 原生接口与自编计算

R调用`DescTools::CochranQTest`、`DescTools::StuartMaxwellTest`、`stats::mcnemar.test`、`binom.test`、`p.adjust`；控制台核对表覆盖Q和SM/Bowker的公式与原生函数。3×3的Bowker对称性检验由基础R `mcnemar.test`完成，连续性校正只作用于2×2。

Python调用`statsmodels.stats.contingency_tables.cochrans_q`、`SquareTable.homogeneity(method="stuart_maxwell")`、`symmetry(method="bowker")`、SciPy二项检验及statsmodels Holm；控制台核对表覆盖Q和SM/Bowker的公式与原生函数。`shift_zeros=False`保留输入的原始计数；本例所有格子都有观测。

MATLAB按正文公式自编全部检验及Holm，分布概率/分位点由`gammainc`、`gammaln`、`betaincinv`、`erfinv`计算。三语Q条件计算均为自编的三法多项计数枚举，针对给定行和的可交换模型。三语条件Bowker枚举均为自编统计量右尾排序，与配置概率排序的精确检验不必同P。条件枚举对应Bowker的转移对称性零假设；Stuart–Maxwell的边际同质问题见正文。

基础R的`mcnemar.test`在反向两格合计为0时，计算中会出现0/0；`DescTools::StuartMaxwellTest`提供渐近检验接口。含零反向格对或奇异协方差的表，需要另选具有相应处理的实现。

本篇自编SM/Bowker函数用于给定非退化教学数据。反向格对合计为0、协方差奇异时会停止。此时需根据有效格对和矩阵秩确定自由度，再选择适合该表的计算方法。

## 第二次运行的练习

表5描述同一4周设计下另一组144人的假设资料，分别在第0/28天按相同规则分类。其三问是：各类人数及反向路径人数多少、SM和Bowker统计量多少、两个数为何不同。

脚本以M2副本构造这张表：把A→B设为15、C→C设为30，其余计数取M中的值。先预测，再查看行列边际48/48/48、SM=0与Bowker=24。表5三对反向转移各为15比3，每对对Bowker贡献8，合计24；该练习用来体会边际汇总与反向路径的区别。

主入口读取逐人记录与前后转移计数，频数CSV供对照。文章配图另在本地交付的“文章配图”目录，由独立`make_assets.py`和相同的转移计数CSV生成；运行本项目三语主入口无需绘图库或TeX。
