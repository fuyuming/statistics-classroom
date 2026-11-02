# 统计课堂23：Kruskal–Wallis之后的Dunn检验

完整下载并解压本项目，先打开本篇文件夹，再选熟悉的一种语言运行。核心数值分析三套语言独立可用。

## 数据与本篇问题

`数据/小鼠存活天数.csv`来自2025年多组非参数检验课堂，并由2026年第12课沿用。按原有记录转存33个完整存活时间：A组11只、B组10只、C组12只。字段`id`是记录编号，`group`是组别，`days`是存活天数。每只小鼠为独立观察单位；相同天数保留为并列值，数据无缺失。这里按全部存活时间已观察到分析；有删失时需使用适合的生存分析方法。

问题依次为：各组原始分布与正态性如何；三组是否存在分布差异；Dunn怎样比较均秩；平均秩差的标准误怎样得到；三对如何做Holm；为什么两两Wilcoxon的P值不完全相同。

## 选一种语言运行

### R与RStudio

双击`统计课堂23.Rproj`，在Files打开`代码/dunn_rank.R`，点击Source。项目不恢复或保存旧工作空间，编码为UTF-8。首次先在Console安装：

```r
install.packages("rstatix")
```

主脚本使用R自带stats与rstatix；工作目录应为本篇项目根目录，能看到`数据/`与`代码/`。终端运行也从此目录执行`Rscript 代码/dunn_rank.R`。结果写入`运行结果/R/`。主脚本先调用`代码/check_distribution.R`，逐组做Shapiro–Wilk检验，并保存原始分布与正态Q–Q图；参考线与`qqline()`默认的四分位数算法一致。该R图使用Noto Sans CJK SC字体；若本机未安装，可把模块中的`family`改为本机中文字体，例如Windows的Microsoft YaHei。只重画这张图时，在读入数据后运行`source("代码/check_distribution.R")`和`check_distribution(dat)`即可。

### Python与Spyder

在Spyder把本篇文件夹打开为现有文件夹项目，打开`代码/dunn_rank.py`，按F5。主脚本从自身位置定位数据，不依赖临时变量。先在Spyder所用的Python环境安装NumPy、pandas、SciPy、Matplotlib：

```bash
python -m pip install -r requirements.txt
```

结果写入`运行结果/Python/`。脚本末尾调用同目录`plot_dunn.py`，生成三张文章科学图。图1、3来自共享数据，图2**明确读取随包附带的R置换结果**，使图中数字与正文一致；Python自己的置换结果保留在Python目录。想重新生成图2所用数据，先运行R脚本；也可只用本包已保存的R CSV绘图。终端可单独运行`python 代码/plot_dunn.py`。

中文字体优先选PingFang SC、Noto Sans CJK SC、Microsoft YaHei或SimHei；缺少时安装其中一种。数字分析不受字体影响。

### MATLAB

需要Statistics and Machine Learning Toolbox。将Current Folder设为本篇文件夹，打开并运行`代码/dunn_rank.m`。脚本用自身位置定位数据，结果写入`运行结果/MATLAB/`。`kruskalwallis`与`tiedrank`用于标准整体检验和并列秩，Dunn–Holm按正文公式分步计算。

MATLAB的`ranksum(...,'method','approximate')`默认包含连续性校正。输出另列这个默认P值，并用相同校正公式核对；`P_no_continuity`列则按两组重编秩的标准误计算无校正双侧正态P，用于与本文R/Python相同设定比较。

## 脚本和输出怎样看

每个主文件按“数据与完整函数→手算Dunn→两两Wilcoxon→置换→Bootstrap→改单位练习”排列。R入口在完整检验前增加各组分布检查。打印的表格都保存为CSV：

| 文件 | 内容 |
| --- | --- |
| `normality.csv`（R） | 每组样本量、Shapiro–Wilk W与P |
| `qq_coordinates.csv`、`qq_lines.csv`（R） | QQ图坐标及参考线截距、斜率 |
| `ranks.csv` | 每只小鼠的天数、组别和全体统一秩 |
| `groups.csv` | 各组只数、秩和、平均秩 |
| `ties.csv` | 每个天数出现次数，用于检查并列 |
| `kw.csv` | H、近似P、秩方差、并列修正合计、整体秩ε² |
| `rank_effect_definition.csv`（R） | 秩的组间/总平方和、秩ε²与秩ANOVA的η²核对、所选经验分级 |
| `dunn.csv` | B−A、C−A、C−B的均秩差、SE、Z、原始P和Holm P |
| `holm_steps.csv` | 排序后的原始P、3/2/1乘数、乘积、累计最大值和调整P |
| `pairwise_wilcoxon.csv` | 三对重新编秩的无连续性校正P及对应Holm P |
| `ab_rerank.csv` | 同一批A/B观察在全体排序和两组排序中的秩 |
| `permutation_differences.csv` | 10000轮固定秩随机分组后的B−A均秩差 |
| `permutation_summary.csv` | 置换均值、模拟标准差与理论SE |
| `bootstrap_epsilon.csv` | 20000轮组内有放回重抽的整体秩ε² |
| `bootstrap_summary.csv` | 原始效应量、2.5%与97.5%百分位点 |
| `unit_check.csv` | 将天数乘10后重新分析的结果 |
| `versions.txt` | 实际运行的软件版本 |

`P_raw`为双侧正态近似P。`P_holm`已对三对比较调整，无需再把它交给p.adjust。R主脚本把手算的Z、原始P、调整P与rstatix逐项核对；Python/MATLAB的确定性结果与R交叉核对，记录见VALIDATION。

## 先核对这些数，再做练习

- 三组秩和105、216、240，平均秩约9.54545、21.6、20。
- 全体秩方差92.65625，Σ(t³−t)为324。
- 秩的组间平方和930.8727273、总平方和2965，比值0.3139537等于H/(N−1)，也等于对全体统一秩做单因素ANOVA的η²；约31.4%是样本秩尺度上的平方和比例。
- 本篇借用effectsize文档Field（2013）的η²分界：0.01、0.06、0.14分别为小、中等、大效应起点；0.314属大效应，Bootstrap区间约[0.078,0.664]仍较宽。分级仅作粗略参照。
- KW H约10.04652，近似P约0.00658304；秩ε²约0.3139537。
- B−A均秩差12.05455，SE约4.20582，Z约2.86616。
- 三对Dunn原始P约0.00415486、0.00927085、0.69786382；Holm P约0.01246459、0.01854170、0.69786382。
- 4天在全部33个观察中取秩8，在A/B的21个观察中重新取秩7。
- A/B无连续性校正的Wilcoxon近似P约0.00415805。

**练习1：天数全部乘10。** 程序确实重新计算秩、Dunn与Holm，核对P和秩效应保持不变；原始天数的中位数乘10。

**练习2：看均秩差的波动。** 置换模拟每轮把相同的33个秩随机分成A11、B10、C12，记录B−A均秩差。许多轮的标准差接近4.20582；有限轮模拟会有波动。这是无差异参照下的演示。

**练习3：区分置换与Bootstrap。** Bootstrap在A、B、C各自组内有放回抽原只数，重新合并编秩；它保留样本中观察到的组别差异，描述整体秩效应的估计不确定性。每种语言seed固定，但随机数算法不同，区间不要求逐位相同。本篇正文采用R的20000轮百分位区间约[0.078,0.664]。秩ε²有0边界，不用Bootstrap区间是否跨0替代KW检验；这里也没有将整体效应区间当作三对差值的同时区间。

## 配图与来源

`文章配图/`含四张PNG及绘图坐标、参数。原有统一秩、置换及重编秩三张图由`代码/plot_dunn.py`生成；新增`04-各组分布与QQ图.png`由R模块`代码/check_distribution.R`生成。正文图号随阅读顺序排列，文件名前缀保持不变。课堂封面与文章稿由上层文章目录管理，不属于数据计算入口。

R Dunn说明：https://rpkgs.datanovia.com/rstatix/reference/dunn_test.html

原始方法：Dunn OJ (1964), Multiple Comparisons Using Rank Sums, Technometrics 6(3):241–252。
https://www.stat.cmu.edu/technometrics/59-69/VOL-06-03/v0603241.pdf

整体秩效应与Bootstrap限制：https://easystats.github.io/effectsize/reference/rank_epsilon_squared.html

经验分级来源（本篇借用于秩尺度）：https://easystats.github.io/effectsize/reference/interpret_omega_squared.html

Dunn＋Holm与两两Wilcoxon＋Holm均为相应条件下的独立组秩比较方案，统一编秩范围不同。
