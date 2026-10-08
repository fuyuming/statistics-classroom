# 精读05｜多控制一个变量，为什么反而会错？

先完整下载并解压项目，不要只下载单个脚本。数据已附，每种语言可以独立运行。

## 打开哪里、结果在哪里

- RStudio：打开 `精读05_v2.Rproj`，打开 `代码/lesson05_v2.R`，Source整个文件。需作者的rethinking及其依赖、systemfonts、ragg。运行后Plots面板有两张中文图；输出在 `运行结果/v2/R/`。
- R tidyverse版：打开 `代码/lesson05_tidyverse_v2.R` 并Source。除原R依赖外，需要dplyr、tidyr、readr、ggplot2，可用 `install.packages(c("dplyr","tidyr","readr","ggplot2"))` 安装。仍使用rethinking::quap；两张ggplot图在Plots显示，输出到 `运行结果/v2/R_tidyverse/`。不需要先运行R原版。
- Spyder：在所用Python环境安装 `requirements_v2.txt` 所列numpy、scipy、matplotlib。打开 `代码/lesson05_v2.py`，F5运行。输出在 `运行结果/v2/Python/`，共四张图。使用环境默认图形后端；可在Plots面板或自动/Qt后端窗口查看。无显示服务时只能保存图。
- MATLAB：打开 `代码/lesson05_v2.m`，按Run。通过脚本位置定位数据，不必手工设置路径；基础MATLAB即可。两张Figure保持打开，结果在 `运行结果/v2/MATLAB/`。

脚本采用UTF-8。需要PingFang SC、Microsoft YaHei、Noto Sans CJK SC或Heiti SC之一，缺字体会提示。中文字体与文件编码是两回事。

因果图和公式：单独运行 `代码/lesson05_figures_v2.py`，导出到 `文章配图_v2/`。封面由内置imagegen制作，不能由该绘图脚本再生。

## 计算范围

1. WaffleDivorce婚姻模型：严格采用PPT31页的标准化、Normal(0,0.2)截距、Normal(0,0.5)斜率、Exponential(1)残差标准差先验。Normal第二项为标准差。三语言均保留sigma不确定性和联合协方差。R调用quap，Python/MATLAB优化后验并计算解析Hessian作二次近似，不是网格、MCMC或调用R。
2. 四种结构：按PPT13、43、55、69页给定概率，独立精确枚举所有0/1组合。前三种里的A为独立无关占位变量；第四种A是中间变量Z的后代。joint.csv概率总和逐模型为1，非真实观测数据。
3. plant_expectation.csv：原书171页生成规则给出的平均增长3.5、4.7及总效应1.2；没有声称重跑原书全部植物后验模型。
4. Python额外绘制先验条件均值线与干预对比。干预预测差图用联合后验抽样混合两个独立残差的差；分布比beta_M后验宽，不是已识别的个体反事实效应分布。
5. 基金、年龄幸福、祖辈例子为文章中的原课机制讲解，不声称本项目复现其全部随机模拟或后验拟合。

## 输入与输出字段

数据/v2/WaffleDivorce.csv从作者rethinking包导出三列，50个地区，每行一个地区。Divorce离婚率（每千名成年人）、Marriage结婚率（教材每千名成年人率）、MedianAgeMarriage结婚年龄中位数（岁）。按原课忽略测量误差建模；没有缺失值。不涉及华夫饼店数，也没有混用总店数与人均店数。

posterior.csv行序：alpha、beta_M、beta_A、sigma；列是近似后验均值、标准差、89%区间端点。associations.csv的model：1叉、2管、3对撞、4中间变量的后代；group：0全体、1固定Z=0、2固定Z=1、3固定A=0、4固定A=1。correlation是理论相关。

先验图的M固定在0，两边用同一批标准正态数缩放；不是完整先验预测数据。干预图使用固定种子20261008。后验数值与验证范围见VALIDATION_v2.md。

## 自己改一次

只改后代A的准确程度，把最后一个伯努利概率从0.1+0.8Z改成0.5。先猜相关，再运行：A不再带有Z的信息，固定A应与全体一样。不要把X到Z或Z到Y的规则一并修改。

原课：https://speakerdeck.com/rmcelreath/statistical-rethinking-2023-lecture-05
作者包与安装说明：https://github.com/rmcelreath/rethinking
项目：https://github.com/fuyuming/statistics-classroom/tree/main/rethinking/article05_elemental_confounds

数据口径补充：Marriage与Divorce为2009年ACS数据，MedianAgeMarriage为2005—2010年指标；共49州及哥伦比亚特区，缺内华达州。`数据/v2/WaffleDivorce_地区说明.csv`保留地名供读者查阅，不改变计算脚本使用的三列输入。

## 批注与阅读顺序

先读数据字段与标准化，再看似然、先验和参数顺序；随后读取posterior.csv中的均值、标准差及89%区间。R调用quap，Python/MATLAB的目标函数逐项标明负对数似然与先验，逆Hessian保留联合后验协方差。四种结构随后独立枚举，注释说明每行联合概率、分组归一化与相关计算。

R tidyverse版使用readr/dplyr/tidyr整理数据、计算条件相关，以ggplot2绘图；posterior.csv与原版字段一致，另有带参数名的posterior_named.csv便于阅读。两个R版本均输出相同模型和枚举数值。
