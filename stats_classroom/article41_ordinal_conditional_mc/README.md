# 统计课堂41：双向有序表的条件模拟

完整下载并解压本项目后，先选一种熟悉的语言。三套程序使用相同输入、评分、双侧尾部，任选一种都能完成正文任务，不需跨语言调用。

## 第一次按节运行

- **R / RStudio**：双击 `article41_ordinal_conditional_mc.Rproj`，打开 `scripts/ordinal_conditional_mc.R`。从“01 数据”起逐节选中运行（Windows Ctrl+Enter / macOS Cmd+Enter），Console看结果、Plots看图。整文件可Source。首节自动载入辅助函数。
- **Python / Spyder**：建立或打开文件夹项目，选择本目录作为项目根，工作目录也设为本目录。打开 `scripts/ordinal_conditional_mc.py`，从第一个 `# %%` 开始用分节运行（Ctrl+Enter），F5可执行整文件。脚本支持有或无`__file__`的入口；按节运行时先执行输入和模拟节，再运行绘图节。
- **MATLAB**：Current Folder设为项目根，在命令窗口执行 `addpath('scripts')`，打开 `scripts/ordinal_conditional_mc.m`，从第一个`%%`开始按Run Section依次运行；命令窗口输入 `ordinal_conditional_mc` 可整文件运行。Current Folder决定相对数据路径，搜索路径只决定去哪里找函数，两者作用不同。首节自动加入 `scripts/helpers`。

首节显示参数状态和输入的前三行字段示例。观察摘要给出N/T/r/Q/chi2_p，chi2_p是χ²1渐近P；枚举摘要给出可能表数、概率和、极端表数及加权精确P。MC完成后显示首张小表及实际行列合计、k、k/B、加一P、MCSE与二项95%区间，并打印由当前结果组成的主例报告。应核对：主例T≈149.0324、r≈0.478375、Q≈63.3895；小例T=8、r=.8、Q=4.48；枚举282张表、概率和1、20张极端表的加权精确P=4/105≈.0380952。三语随机序列允许不同。

## 安装实际依赖

R需要 `DescTools`、`coin`，可在RStudio执行 `install.packages(c("DescTools", "coin"))`。Python需要requirements.txt中的numpy、pandas、scipy、matplotlib；请在Spyder实际使用的解释器环境安装（备选命令：`python -m pip install -r requirements.txt`）。MATLAB需要 Statistics and Machine Learning Toolbox（chi2cdf与binofit）。

命令行备选入口（均在项目根）：`Rscript scripts/ordinal_conditional_mc.R`、`python scripts/ordinal_conditional_mc.py`；MATLAB用 `addpath('scripts'); ordinal_conditional_mc`。

## 输入、来源与统计定义

共享文件 `data/ordinal_tables.csv` 的每行对应一个格子，字段为 `example`（main主例/mini假设小例）、`row_score`、`col_score`（各1到4）、`count`（例数）。每表16行，没有缺失值或缺失处理；0表示该格没有记录。脚本先展示字段及示例行，再依评分排序形成频数表。

主例原题来自颜艳、王彤《医学统计学》第5版p116–117例7-9、表7-11，278例尸解资料：年龄组20～/30～/40～/≥50岁，对应硬化等级−/+/++/+++。不同个体按题设独立，两边事先评分1:4。mini另设8例假设尸解资料，记录相同的两项指标、采用相同组序和评分，各年龄组与各硬化等级均为2例。它是独立构造的教学计数，不是从278例中抽出的样本。

目标：独立H0下给定行列边际，以 `abs(T)>=abs(Tobs)`（含等号）定义双侧尾部；T是中心化评分积和，r=T/sqrt(Lxx Lyy)，Var(T)=Lxx Lyy/(N−1)，Q=(N−1)r²，渐近P为χ²1右尾。少量浮点容差只防止相等端点因舍入漏计。边际概率按 `prod(row!)*prod(col!)/(N!*prod(cell!))`，不能把282张不同表等概率处理。

小表的完整枚举是三语自编算法。R通过原生 `r2dtable` 做正确条件抽样，`DescTools::MHChisqTest`和`coin::lbl_test(approximate(...))`核对完整检验；coin的exact接口对本4×4报两样本问题限制。Python/MATLAB通过自编程序置换个体列标签，属于相同条件模型的另一实现。χ²近似、完整枚举和Monte Carlo的来源分别打印。

## 输出与保存

首次默认只显示结果和图。需要写CSV/PNG时，将环境变量 `ARTICLE41_SAVE_RESULTS` 设为1，再从首节重跑：R用 `Sys.setenv(ARTICLE41_SAVE_RESULTS="1")`，Python在导入os后用 `os.environ["ARTICLE41_SAVE_RESULTS"]="1"`，MATLAB用 `setenv('ARTICLE41_SAVE_RESULTS','1')`。恢复默认可设置为0。结果写入 `results/`，图写入 `figures/`，按语言命名互不覆盖；随包结果是本轮实际计算。

`*_mini_exact_distribution.csv`给出每张可能表的T、概率、是否极端；按T归并即完整统计量分布。`*_mini_sim_T.csv`（Python为Python_mini_sim_T.csv）保存每次小表模拟的T。主文件用本次内存中的概率与模拟T自动调用绘图函数，无需旧CSV。上图蓝色为中心、橙色为含±Tobs尾部；下图比较各T值的精确概率与模拟频率。英文图标对应中文正文图注。

三套主脚本仍计算“演算中间量”和尾部归并，但默认不打印；设置`ARTICLE41_SHOW_DETAILS=1`后显示两评分总和/均值、第一行对Lxx贡献、两格中心化积/T贡献、sqrt(Lxx·Lyy)、VarT、SDT、Z及极端表分组。保存时写出`*_formula_steps.csv`、`*_mini_tail_summary.csv`和`*_mini_first_table.csv`；最后一个是无表头、无行名的4×4纯数值矩阵，按当前行列顺序存储同批模拟首张小表。其边际与输入相同，T等于`*_mini_sim_T.csv`首行；Python读取时用`pd.read_csv(..., header=None)`。默认只显示结果，不写新文件。底层Python算法在`scripts/helpers/ordinal_helpers.py`；先理解主入口，再按需阅读辅助函数。

## 显示详细演算与从外部目录运行

默认`ARTICLE41_SHOW_DETAILS`未开启，适合第一次按节运行；需要完整频数矩阵、逐格演算和六组极端T尾部归并时，先设置环境变量为1，再从首节运行。设置方法与保存开关相同，例如Python用`os.environ["ARTICLE41_SHOW_DETAILS"]="1"`。Python脚本优先按自身路径定位项目根，Spyder无`__file__`时可将工作目录设为项目根或`scripts`；MATLAB首节会在当前目录找不到数据时按脚本路径回到项目根。三语均在首节显示B、SEED、保存开关和详细开关状态。
R的4×4 `coin` exact接口探针（该接口会提示“不属于two-sample problem”）也仅在详细模式显示；默认首轮只显示正常核对结果。

二项区间/MCSE描述条件尾概率的计算误差，不是r的置信区间。k=0时单侧95%上限为1−.05^(1/B)；单侧上限与双侧95%区间的上端不同。上限列仅在零超越时填写，其余留空。正文主要采用R的B=19999、seed=20261122小例k=765；其他语言随机结果以其实际输出为准。

## 第二次改参数

先将B改为4999，再改为79999，均从首节重跑：MCSE怎样变化？观察单位数是否改变？再试两边评分同时倒序（4:1），比较r、Q及双侧P；仅倒序一边看r符号。改评分是理解统计量，不是在多次尝试后挑最小P；实际分析应预先确定评分。R如改评分，coin段显式scores参数会同步采用。

核对提示：同一小例中，B约增至4倍时，MCSE通常约减半，观察资料仍为8例；某次模拟P值与精确值的距离不一定随B增大而单调缩小。仅将一边的评分倒序，r变号，Q与精确双侧P保持不变；两边都倒序，r也保持不变。若同时固定B、种子及抽表实现，本项目同一批表的双侧极端次数也保持不变。
