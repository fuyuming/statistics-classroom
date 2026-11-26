# 统计课堂44：ANOVA功效与样本量

本项目用四组独立、正态、同方差总体演示研究前规划。三语入口都按01—05五节顺序运行：01读CSV，02求目标功效下的每组n，03固定n功效，04共同SD敏感性，05可选保存/绘图。均数与SD是规划参数，不是一次小样本的观察结果。

## 题干与结果

比较4种支架设计的压缩模量（MPa），预定经典单因素ANOVA总体F检验，α=.05，目标功效80%并另看90%；默认固定每组n=5。large均数30、35、40、45、SD=5；small均数36、37、38、39、SD=5。预期large每组n=4、small每组n=56（90%目标为72）；固定n=5时功效约97.1%与10.1%。把small的SD改为10，80%目标应为每组220。样本数指独立实验单位。

完整解压后保留“代码、数据、运行结果、文章配图”的目录层级，再选一种语言开始。

## 运行

R需要`pwr`，使用`pwr::pwr.anova.test()`；Python需要requirements.txt中的numpy、scipy及绘图依赖；MATLAB需要Statistics and Machine Learning Toolbox。入口文件分别为`代码/anova_power.R`、`代码/anova_power.py`、`代码/anova_power.m`。Python分节环境无`__file__`时按当前项目目录定位输入。默认不保存结果；设置`ARTICLE44_SAVE_RESULTS=1`可写入每语种planning、fixed_n、sensitivity三表。设置`ARTICLE44_PLOT=0`关闭绘图。绘图脚本仅生成文章图1、图3、图4。

### R / RStudio

打开`统计课堂44_ANOVA功效与样本量.Rproj`，再打开R入口；首次安装`install.packages("pwr")`。按01—05节顺序运行（Ctrl/Cmd+Enter），也可Source整文件。

### Python / Spyder

在Spyder实际使用的解释器环境安装`requirements.txt`所列依赖，命令行可用`python -m pip install -r requirements.txt`。在Spyder打开项目文件夹和`代码/anova_power.py`，从第一个`# %%`单元顺序运行或F5整文件。绘图需要Noto Sans CJK SC、PingFang SC、SimHei或Arial Unicode MS等中文字体；也可命令行运行`python 代码/anova_power.py`。

### MATLAB

需要Statistics and Machine Learning Toolbox。将MATLAB Current Folder设为项目根目录，打开`代码/anova_power.m`，按01—05节Run Section或整文件运行；入口会加入`代码/helpers`路径。

在运行第01节之前按需设置保存开关；默认关闭，不覆盖随包结果。

```r
Sys.setenv(ARTICLE44_SAVE_RESULTS = "1")
```

```python
import os
os.environ["ARTICLE44_SAVE_RESULTS"] = "1"
```

```matlab
setenv('ARTICLE44_SAVE_RESULTS','1');
```

结果分别写入`运行结果/R`、`Python`或`MATLAB`；Python同时保存数值图和坐标。恢复关闭可把`ARTICLE44_SAVE_RESULTS`设为`0`；关闭绘图可将`ARTICLE44_PLOT`设为`0`。

## 看哪张结果表

| 表 | 要核对什么 |
| --- | --- |
| planning.csv | 两种效应、80%/90%目标的连续计算值、整数每组n、实际功效与前一整数功效 |
| fixed_n.csv | 当前固定n_check下，两个情形的f及功效 |
| sensitivity.csv | 小差异均数不变，SD取2.5/5/10时的样本量 |

`continuous_n`是求根所得的计算值；`n_per_group`为最终整数每组数量，`total_n`为总数量。`power`是将整数n回代后的功效；`previous_power`帮助核对前一个整数是否尚未达标。

## 共享输入与解释

`数据/planning_scenarios.csv`含每种情形四组总体均数和共同SD，`数据/settings.csv`含α、目标功效和固定n。f是组均数散布相对共同SD的标准化效应量。正式方案应说明有实际意义的均数组合、SD依据、α、目标功效、主要检验及损耗预留；整体F检验功效不等于两两比较功效。三语入口中Python/MATLAB保留非中心F作为内部计算实现，公开输出只呈现n、f、power三项。

## SD敏感性练习

保留原文件副本后，将small四行的`sd_MPa`从5改为10，从01节运行到03节，核对small 80%规划为每组220；不要把一次观察η²自动当成未来真实效应。

主要来源：Navidi《Statistics for Engineers and Scientists》第6版§9.1的功效定义及R pwr官方文档。
