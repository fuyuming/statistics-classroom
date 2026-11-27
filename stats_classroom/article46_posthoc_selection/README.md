# 统计课堂46｜多组均数怎样比较？

本项目配合《多组均数怎样比较？LSD、Tukey、Dunnett和SNK怎么选》。先确定研究目标与比较族，再读差值、区间与P值。任选一种语言跑通即可。

## 先下载解压，再按节运行

拿到本篇ZIP后完整解压，进入`article46_posthoc_selection`，先读本README。保留`data`、`代码`与项目文件的相对位置。若下载整个系列仓库ZIP，解压后进入`stats_classroom/article46_posthoc_selection`。

- **R / RStudio**：打开`统计课堂46_多组均数方法选择.Rproj`，打开`代码/posthoc_selection.R`，从第1节开始，按节选中运行；确认数据和MSE后，再继续各比较。需要整文件时点Source。
- **Python / Spyder**：打开本文件夹项目，把工作目录设为含`data`的项目根，打开`代码/posthoc_selection.py`，按`# %%`依次运行。第1节支持没有`__file__`的分节方式；F5可运行整文件。
- **MATLAB**：把Current Folder设为本项目根，打开`代码/posthoc_selection.m`，按`%%`依次运行。把代码文件夹添加到搜索路径只帮助MATLAB找到函数，不会把Current Folder改为数据所在目录。

三语均分7节：设置与数据 → 整体模型/MSE → LSD → Tukey → Dunnett → SNK → 绘图及可选保存。每节依赖前面已经运行的节；单独跳到后面会缺少变量。

## 依赖

R使用基础`stats`和`multcomp`包。在当前RStudio中执行`install.packages("multcomp")`。

Python需要NumPy、pandas、SciPy≥1.15（使用`dunnett`的`rng`参数）、Matplotlib。在Spyder实际使用的解释器环境安装`requirements.txt`；命令行可用`python -m pip install -r requirements.txt`。

MATLAB需要Statistics and Machine Learning Toolbox。已在R2025a运行。SNK部分的学生化极差概率与分位数由项目代码目录中的studentized_range_cdf、studentized_range_inv自编数值函数计算；函数只调用公开MATLAB函数，不依赖R/Python后端。较长算法供核对使用，新手可先按主脚本分节运行。

命令行备选入口，均从项目根运行：

```text
Rscript 代码/posthoc_selection.R
python 代码/posthoc_selection.py
matlab -batch "run('代码/posthoc_selection.m')"
```

## 数据和分析计划

`data/scaffold.csv`含`group`和`modulus`两列。A1–A4各5个独立试件，每件一个有效压缩模量，单位MPa；无缺失。数字沿用李春喜等《生物统计学》第6版例7-1，支架情境与单位为教学改编，并非真实支架实验。课堂原件保持不变。

模型假定观察独立、各组总体正态且具有共同方差。双侧α=.05；区间为95%。LSD使用全部四组共同MSE和16个误差自由度；不是只抽两组重新估计方差。Tukey保护全部六对。Dunnett只保护A2−A1、A3−A1、A4−A1，A1须事先指定为参照。这里同时计算四种方法便于学习，实际研究须事先确定计划。

SNK为等组量的**教学自编实现**：按均数降序排列，跨度r包含两个端点；对每对的所有包含区间取局部极差P的最大值，执行“大区间未拒绝，其内部比较不能拒绝”的规则。`raw_range_P`仅为局部极差P；`step_P`是逐步判定值，不是强FWER调整P，且不提供Tukey式的同时95%区间。

## 第一次核对这些结果

- MSE=.4285 MPa²，误差df=16；整体F=12.01828，P≈.00022669。
- A1−A3=1.10 MPa：LSD普通P≈.01722，普通95%区间[.22235,1.97765] MPa。
- 同一差值Tukey调整P≈.07352，同时95%区间[−.08448,2.28448] MPa。
- Dunnett的A3−A1=−1.10 MPa，调整P约.044，同时95%区间约[−2.17,−.03] MPa。正文演示可直接用R的`summary(dun_model)`/`confint(dun_model)`默认接口；配套主入口为便于三语核对提高数值积分精度，因此正文统一报告差值2位、区间端点2位、调整P约`.022`、`.044`、`<.001`（默认积分末位可能略异；例如A2−A1约.0224、A3−A1约.0440、A4−A1约.00010）。
- SNK排序A1>A3>A2>A4，A1−A3局部P≈.01722、step_P≈.02207；只有A3−A2不拒绝，其余五对拒绝。

Dunnett段还会打印共同参照的尺度关系：本例`MSE/n_A1=.0857 MPa²`，每个处理−A1差值方差为`.1714 MPa²`，A2−A1与A3−A1协方差为`.0857 MPa²`，相关系数为`.5`。这些数值由脚本根据MSE、样本量和`UNIT_SCALE`动态计算，用于说明联合t参照的相关结构；MSE代入只是演示共同尺度。SNK段会按变量筛选A1−A3的三个包含区间（r=2、3、4），并打印其`raw_range_P`及最大值`step_P`，与正文逐步判定相互核对。

脚本打印的结论由当前结果生成。LSD表的`reject`按普通P判定；前面单独打印整体F门是否通过，不能把该列当作六比较强FWER结果。未拒绝不等于等效。

`difference`为有方向差值，`SE`为均数差标准误，`CI_low/CI_high`为区间端点。R的Tukey原输出为后组减前组，Python/MATLAB保留前组减后组；比较跨语结果时应先统一方向。MATLAB的组号1–4分别代表A1–A4；SNK的i/j是排序位置，high/low为原组标识。

## 图与保存

第7节绘制全部实际观测点和组均数。首次只显示图；代码顶部SAVE_OUTPUTS改为TRUE/True/true后，再从第1节运行，更新运行结果/语言名/中的5张CSV及图片。教学核验也支持环境变量SAVE_OUTPUTS=1。MATLAB的LSD节同时保留手算t/SE表和公开multcompare(...,'CType','lsd')核对表。随包已保留三语实际输出，代码会从共享输入重新计算，不读取这些旧结果。

RStudio的图出现在Plots窗；非交互式Rscript可能创建默认`Rplots.pdf`，这不是本篇必需材料。Python默认不强制Agg，方便Spyder查看；无图形界面自动运行可设置`MPLBACKEND=Agg`。

## 第二次：只改单位

把第1节`UNIT_SCALE`从1改成1000，再从第1节运行：MPa转为kPa，均数、差值、SE和区间都乘1000，MSE乘1000000；F、t、q、P和显著判断保持不变。图轴自动改成kPa。不要再同时改CSV，否则会重复缩放。也可修改一个观察值观察共同MSE如何改变，但不以追求显著为目标。

本篇项目目录：
https://github.com/fuyuming/statistics-classroom/tree/main/stats_classroom/article46_posthoc_selection
