# 统计课堂47｜η²、ω²与偏η²

下载并完整解压本篇ZIP后，先进入本目录。任选一种语言，先按1–6节运行，再试整文件运行。三套计算都在各自语言内完成，不依赖另一语言后端。

## 题目与数据

假设团队比较A1–A4四种支架孔结构，材料与外形尺寸相同，各独立制备5个试件，共20个；制备、测试顺序随机，每个试件测一次有效压缩模量MPa。模型假定独立、正态、共同方差。求整体检验、η²/ω²/偏η²、带截距R²与调整R²、A1−A3差值及全部6对的Tukey同时95%区间。

`data/scaffold.csv`中`group`为分类设计、`modulus`为MPa。20数字原样复制课堂现用CSV，源于李春喜等《生物统计学》第6版例7-1，只改情境和单位，非实际支架实验。无缺失，脚本遇到模量缺失会报错，避免默默改变样本。

`data/factorial_demo.csv`为完整人工2×2辅助例：孔结构design=D1/D2，温度temperature=low/high，各格2个独立试件，每件一次MPa。D1low[9,11]、D2low[11,13]、D1high[13,15]、D2high[15,17]；独立、正态、共同误差方差，拟合含交互作用完整模型，只理解孔结构两种分母。

## R与RStudio

1. 打开本目录`统计课堂47.Rproj`，不恢复旧工作空间，UTF-8编码。
2. 打开`R/anova_effect_sizes.R`。首次在控制台执行`install.packages("effectsize")`，本文实际核对版本1.0.3。
3. 从第1节选中代码运行，按1–6节顺序继续；RStudio用Ctrl/Cmd+Enter运行所选代码。整文件可点Source，或在项目根运行`source("R/anova_effect_sizes.R")`。
4. `aov`返回拟合对象；`anova`提取平方和表；`eta_squared/omega_squared`从模型计算效应量表。明确`partial=FALSE, ci=NULL`只求普通效应量点估计。普通单因素的`partial=TRUE`会提示与η²相等并仍返回`Eta2`列。

## Python与Spyder

1. 用自己的Python环境安装`requirements.txt`列出的numpy、pandas、scipy、statsmodels、matplotlib。例如终端运行`python -m pip install -r requirements.txt`，该python须与Spyder内核一致。
2. Spyder中新建文件夹项目或将工作目录设为本目录，打开`Python/anova_effect_sizes.py`，从第1个`# %%`单元开始逐节执行；也可F5整文件运行。
3. 脚本支持Spyder没有`__file__`的分节方式；此时工作目录必须含`data`。整文件运行能从脚本路径定位。
4. `statsmodels.ols`拟合带截距分类模型；`C(group)`表明分类。`scipy.stats.tukey_hsd(*groups)`接收四组数组，`*`把列表展开为四个参数，返回组对差值和调整P矩阵；置信区间对象的`low/high`按同一[i,j]方向读取（第i组减第j组）。

## MATLAB

1. 将本目录设为Current Folder，打开`MATLAB/anova_effect_sizes.m`。需要Statistics and Machine Learning Toolbox。
2. 从第1个`%%`开始按节运行；整文件点Run，或在项目根执行`run('MATLAB/anova_effect_sizes.m')`。
3. Current Folder决定相对数据路径；添加搜索路径主要帮助查找函数，不能代替正确工作目录。整文件还能按脚本位置定位项目根。
4. `anova1`返回P、ANOVA表和stats结构；`multcompare`接收该stats，明确Tukey与Alpha=.05。输出列为组i、组j、下界、i−j差值、上界、调整P。`anovan`在第5节提取完整二因素模型的来源与误差平方和。

## 先核对什么，再改什么

六节顺序：读数据与设计 → 手算平方和/效应量 → 模型/R²/包核对 → A1−A3差值与区间 → 2×2分母 → 图形与单位练习。

- 主例SSB=15.4495，SSW=6.856，SST=22.3055，MSE=.4285，F(3,16)≈12.0182808246，P≈.000226687542。
- η²=偏η²=R²≈.6926318621；ω²原值≈.6230315827；调整R²≈.6350003362。
- A1−A3=1.10 MPa；Tukey同时95%CI约[−.084476，2.284476] MPa，调整P≈.0735237。
- 辅助例SSA=8、温度SSB=32、交互SSAB=0、SSE=8、SST=48；孔结构η²=1/6，偏η²=.5。
- 第6节全数×1000，从MPa换为kPa后重新拟合与重算：差值/区间端点×1000；η²/ω²/R²不变。脚本检查的是这些恒等关系，没有用原例常数阻止读者改动。

η²是样本变异比例。ω²用于偏差修正，仍非总体真实值；`omega2_raw`保留公式值，`omega2_report=max(0,raw)`采用与effectsize相同的非负截断约定。这里没有计算η²或ω²区间；所报区间针对A1−A3的均数差。

Field规则在本篇固定写`rules="field2013"`，界值.01/.06/.14，只作粗略参照。官方较新文档另列`cohen1992r7`；本机1.0.3不据此假设新别名可用。`cohen1992`指第8行多重相关语境，转换界值.02/.13/.26。专业意义应回到MPa和研究需求。

## 保存、结果和整文件备用入口

默认只打印、看图，不写结果。需另存时，R先执行`Sys.setenv(SAVE_RESULTS="1")`，Python在首节设置环境变量前执行`import os; os.environ["SAVE_RESULTS"]="1"`，MATLAB先执行`setenv('SAVE_RESULTS','1')`，再从第1节重跑。恢复默认改为`"0"`。图形在各自IDE图窗显示；Python绘图窗口需关闭后整文件才结束。

项目根的命令行备用入口：

```sh
Rscript R/anova_effect_sizes.R
python Python/anova_effect_sizes.py
matlab -batch "run('MATLAB/anova_effect_sizes.m')"
```

`运行结果/R`、`Python`、`MATLAB`各有5个实际数值CSV和1张实际PNG，不互相覆盖：主例结果、A1−A3区间、二因素分母、kPa效应量、kPa区间及四组散点图。修改输入并开启保存后会覆盖对应语言结果。

## 方法文档

- effectsize：https://easystats.github.io/effectsize/reference/eta_squared.html
- 解释规则：https://easystats.github.io/effectsize/reference/interpret_omega_squared.html
- SciPy Tukey：https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.tukey_hsd.html
- MATLAB比较：https://www.mathworks.com/help/stats/multcompare.html
- R Tukey：https://stat.ethz.ch/R-manual/R-devel/library/stats/html/TukeyHSD.html

本篇公开项目：https://github.com/fuyuming/statistics-classroom/tree/main/stats_classroom/article47_anova_effect_sizes

从仓库目录学习时进入`练习项目`；下载本篇独立项目ZIP时，完整解压进入`article47_anova_effect_sizes`。简短运行依据见`VALIDATION.md`。
