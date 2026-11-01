# 统计课堂18：McNemar精确检验与mid-P

本篇合并讲McNemar配对四格表的样本量要求、−1连续性校正、条件精确检验与mid-P。比较四种P值，拆开严格尾部与边界概率，再枚举条件第一类错误率。本目录为公众号文章配套的公开练习代码。

公开项目目录：[GitHub 项目目录](https://github.com/fuyuming/statistics-classroom/tree/main/stats_classroom/article18_mcnemar_midp)；[仓库 ZIP](https://github.com/fuyuming/statistics-classroom/archive/refs/heads/main.zip)。

## 下载解压后怎样开始

完整解压，保留“数据、代码、运行结果”文件夹。用自己熟悉的一种语言跑通即可，三语言读取同一份CSV。

- **RStudio**：打开`统计课堂18-McNemar精确检验与midP.Rproj`，再打开`代码/lesson18_mcnemar_midp.R`并Source。只用基础R，不需安装扩展包。
- **Spyder**：把工作目录设为本项目文件夹，打开`代码/lesson18_mcnemar_midp.py`，F5运行全文，或依次运行四个`# %%`单元。依赖numpy、pandas、scipy；复现配图另需matplotlib。缺包时在Spyder当前解释器对应的终端执行`python -m pip install -r requirements.txt`；先核对Spyder使用的Python解释器，避免装到另一个环境。
- **MATLAB**：将Current Folder设为本项目文件夹，打开`代码/lesson18_mcnemar_midp.m`，Run运行全文或依次Run Section。仅加入搜索路径不能代替设置Current Folder。脚本用`gammaln`计算二项概率、`erfc`计算自由度1的卡方右尾，无需Statistics and Machine Learning Toolbox。末尾局部函数需支持脚本内局部函数的MATLAB版本（本次实际使用R2025a）。

先运行全部四部分。R的Environment、Spyder的Variable Explorer、MATLAB的Workspace均可查看命名结果；选中对象名执行也可显示。整文件运行不依赖打印长段文字讲课，各部分含义和预期结果见下表。结果自动另存于`运行结果/R`、`运行结果/Python`、`运行结果/MATLAB`。

## 题目、数据与四个问题

共享文件`数据/sle_mcnemar.csv`来自第11课和去年分类资料课件中的58人课堂例题，不是新增临床数据。行是免疫荧光阳性/阴性；`latex_positive`、`latex_negative`是乳胶阳性/阴性人数。四格为`[[11,12],[2,33]]`，所以a=11、b=12、c=2、d=33，N=58，m=b+c=14。数据为完整的非负整数计数，无缺失；不能把缺失值当0，也不要把百分比填进人数格。

| 步骤 | 回答什么 | 预期结果或查看对象 |
| --- | --- | --- |
| 1 | H₀下哪些方向计数可能出现？ | k=0到14；`probabilities`总和为1；结果不等可能 |
| 2 | 边界如何改变P值？ | 严格尾部0.0018310546875；边界0.0111083984375；`observed_four_methods`列出四种P |
| 3 | 11比3的教学变更有什么不同？ | `exercise_11_3`两行比较；仅第一行为原题，第二行为人工教学变更 |
| 4 | 固定m时整套拒绝规则多常误报？ | `conditional_size`；改变α或m后看`alpha_exercise_summary`及拒绝取值 |

原例普通近似、连续性校正、精确法、mid-P的P值依次约为0.0075263152、0.0161569313、0.0129394531、0.0073852539。11比3的教学变更依次约为0.0325094446、0.0613688291、0.0573730469、0.0351562500。它们是同一份数据的不同检验处理，不是四份独立证据。

## 统计定义与边界

前提是独立配对单位、明确的二分类和配对关系。H₀为总体边际阳性率相等；给定m后，方向计数B服从Binomial(m,0.5)。本课双侧统计量`T=abs(2*k-m)`是k距中心m/2的距离乘2。严格尾部为`P(T>tobs)`，边界为`P(T==tobs)`；精确P等于二者之和，mid-P等于严格尾部加半个边界。

`b=c`时中心边界只出现一次，mid-P仍按该定义计算，不能机械地把左尾mid-P乘2。主入口遇到m=0会停止，不生成新结果；R/Python的第三部分辅助函数遇m=0返回NA，便于批量比较时明确标记不可计算。脚本停止后，旧CSV不会自动删除，不能把先前运行的输出当成这次的结果。

拒绝规则统一为`P<=alpha`。给定m=14、α=0.05，普通近似和mid-P的拒绝范围为k≤3或k≥11，条件错误率0.057373046875；校正近似和精确法为k≤2或k≥12，条件错误率0.012939453125。按每个k的二项概率加权，不能除以15。这些是H₀下整套规则的条件第一类错误率，既不是一次样本P值，也不是对随机m再平均的无条件错误率。

## 改参数练习与结果文件

先改α，再改m，观察拒绝域与误报率如何跳跃；不要为追求显著性改数据。脚本已包含两个可核对的练习：m=14、α=0.01时普通/mid-P错误率0.012939453125，校正/精确为0.0018310546875；m=10、α=0.05时四种方法均为0.021484375。详细拒绝取值由代码重新枚举。

| CSV | 内容 |
| --- | --- |
| observed_four_methods.csv | 原题四方法统计量与P值 |
| boundary_components.csv | 严格尾部、完整边界、精确P和mid-P |
| all_outcomes.csv | 固定m=14的全部15个结果、概率、P值及拒绝标记 |
| conditional_size.csv | α=0.05时四方法条件错误率 |
| exercise_11_3.csv | 原题与人工11比3的对照 |
| alpha_exercise.csv | 两种改参数练习的全部取值与拒绝标记 |
| alpha_exercise_summary.csv | 四种错误率和各自拒绝取值 |

## 复现三幅科学配图

项目含`文章配图/`的三张图与对应坐标CSV。先运行R主脚本，再在项目根目录执行`python 制作配图.py`；它读取R枚举输出，导出二项概率图、边界减半图和P值跳跃图。这里不是随机模拟，无需设置种子。图脚本使用matplotlib并需要中文字体（如PingFang SC或Noto Sans CJK SC）；缺字体会给出明确提示。它只负责导出图片，课堂主脚本的数值计算不依赖绘图环境。三个语言的确定性结果已交叉核验。

半格边界刻度图也可单独复现：在项目根目录执行`python 绘制半格校正.py`，导出`文章配图/04-半格到减一.png`及其`04-刻度参数.json`。

来源与定义：R官方https://stat.ethz.ch/R-manual/R-devel/library/stats/html/binom.test.html；Fagerland等2013，https://link.springer.com/article/10.1186/1471-2288-13-91。完整验证见`VALIDATION.md`。
