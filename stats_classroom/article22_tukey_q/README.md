# 统计课堂22配套练习：Tukey 检验与学生化极差 q

下载完整项目并解压后，先进入本篇文件夹，再选择以下一种语言运行。

本练习用三个小麦品种的地块产量说明：方差分析后用 Tukey 同时比较定位差异，并理解组数增加时 q 分布临界值为何提高。

## 数据与研究设计

`数据/小麦产量.csv` 有 3 列、每列 10 个观测：每列代表一个品种，每个数来自一个独立地块。数据源于第12课例9-1的课堂配套文件，按原数值转存。原题未给产量单位，因此练习中保留原题数值，不补写单位；数据无缺失。同行的三个数不是同一地块的配对观测。

## 运行方法

### R

用 RStudio 打开 `统计课堂22.Rproj`，在右下角 Files 中打开 `代码/tukey_q.R`，点击 **Source**（或运行 Source）。使用R自带`stats`与额外的`car`包。首次在RStudio Console运行`install.packages("car")`；之后无需重复安装。输出写入 `运行结果/R/`，并生成诊断图 `qq_diagnostics.png`、`q_tail.png`。

### Python

用 Spyder 打开本“配套练习”文件夹（现有文件夹项目），打开 `代码/tukey_q.py`，按 **F5** 运行主脚本；主脚本会自动调用 `plot_tukey_q.py` 绘制文章配图。先按 `requirements.txt` 安装 NumPy、pandas、SciPy、Matplotlib：

```bash
python -m pip install -r requirements.txt
```

请在Spyder当前使用的Python环境安装上述依赖。若中文字体不可用，可安装Noto Sans CJK SC；数值计算不受字体影响。

输出写入 `运行结果/Python/`，配图写入 `文章配图/`。无图形界面时可使用 Matplotlib 的 Agg 后端。

### MATLAB

需要 Statistics and Machine Learning Toolbox。在 MATLAB 中把 Current Folder 设为本“配套练习”目录，运行 `代码/tukey_q.m`。正式的方差分析和多重比较使用公共函数 `anova1`、`multcompare`（Tukey-Kramer）；`internal.stats.stdrinv` 与 `internal.stats.stdrcdf` 只用于 q 分布演示，属于内部接口，可能随 MATLAB 版本变化。输出写入 `运行结果/MATLAB/`。

## 分析顺序与前提

三套脚本均先给正态性诊断与方差齐性检验，再执行ANOVA和Tukey。R/Python使用逐组Shapiro–Wilk，三种语言均绘制Q–Q图。Levene中位数版（Brown–Forsythe）与均值版分别比较到组中心的绝对偏差；Bartlett直接使用各组方差，对非正态更敏感。诊断结果需结合设计和图形判断，P大于0.05不证明前提成立。

本例Shapiro–Wilk P为0.7517、0.5992、0.0493；品种3需结合Q–Q图进一步查看。Levene中位数版、均值版、Bartlett P依次为0.4071、0.3747、0.8601。MATLAB没有原生Shapiro–Wilk，脚本保留Q–Q图并明确这一差别。

Tukey HSD使用学生化极差分布控制全部两两比较；本例HSD临界差值约47.59。输出P已经调整，对于同一组三对比较不再调用`p.adjust()`。

## 输出 CSV 怎么读

三种语言的输出目录结构相同（MATLAB 少数诊断表略有不同）：

- `variance_tests.csv`：Levene中位数版、均值版和Bartlett的统计量及P值；三语言定义一致。
- `groups.csv`：各组样本量、均值、组内方差、组内平方和。
- `anova.csv`：SSA、SSE、分子/分母自由度、MSE、F 值和总体检验 P 值。
- `pairs.csv`：品种 `2-1`、`3-1`、`3-2` 的均值差、q、Tukey 调整 P 值及同时 95% 区间。
- `k2_check.csv`：当组数 `k=2` 时，q 右尾概率与双侧 t 概率的核对。
- `q_samples.csv`：模拟中 `k=2/3/5` 的每次 q 值（共 `B=100000` 次）。
- `simulation.csv`：每个组数的理论 `q95`、模拟 95% 分位数和超过理论临界值的比例。
- `unit_check.csv`：所有数据乘 10 后的均值差、区间和 P；用于观察量纲变化。
- `diagnostics.csv`（R/Python）：各组 Shapiro-Wilk 检验 P，只有诊断提示，不替代设计与图形判断。

Python的`qq_diagnostics_before_anova.png`是在主分析前产生的诊断图；正文竖排Q–Q图另存于`文章配图/04-正态QQ诊断.png`，其坐标与R的qqnorm/qqline定义一致。

R/MATLAB 的 `q_tail.png`、Python/R/MATLAB 的 Q–Q 图是可视化诊断或分布展示，不是新的检验结论。

## 练习答案与应看到的数值

本数据的组内均方误差为 `MSE=1842`，残差自由度为 `df=27`。品种 `2-1` 的均值差为 `-90`，同时 95% 区间约为 `[-137.59, -42.41]`，调整后 `P≈0.00020175`。三组的 q 分布练习固定 `df=27`：

| 组数 k | q 的 95% 临界值 |
|---:|---:|
| 2 | 2.9017 |
| 3 | 3.5064 |
| 5 | 4.1305 |

把全部数据乘以 10 后，均值差和置信区间乘以 10；q 值与 P 值不变。

模拟采用

$$
Q=\frac{\max_i Z_i-\min_i Z_i}{\sqrt{U/\nu}}.
$$

其中ν＝27；每个Zᵢ表示标准化后的组均值误差，彼此独立且服从标准正态分布。U与所有Zᵢ独立，服从χ²(27)，U/ν表示共同方差估计相对于真实方差的波动。分子和分母分别展示均值与方差估计的抽样波动。`B=100000` 时超越理论 95% 临界值的比例应在 0.05 附近。为直观看组数影响，三种 k 共用同一批 Z 的前 2、3、5 列；这是演示设计，不是三组独立模拟。即使相同 seed，R、Python、MATLAB 的随机数算法和序列也可能不同，因此三语言模拟逐行不必相同。

## 参考文档

- [R `TukeyHSD`](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/TukeyHSD.html)
- [R `qtukey`/`ptukey`](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/Tukey.html)
- [SciPy `studentized_range`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.studentized_range.html)
- [MATLAB `multcompare`](https://www.mathworks.com/help/stats/multcompare.html)

三种语言均已实际运行，确定性结果与模拟特征的核对见`VALIDATION.md`。本篇项目拟上传目录为`stats_classroom/article22_tukey_q`，当前尚未上传。


## 为什么q分布实现多重比较调整？

三种语言主脚本均输出`family_control.csv`。本例MSE＝1842，残差自由度27，每组10个地块。预先只比较一对时，双侧t检验的q尺度门槛为2.9017265，对应本次均值差39.3823。把这个单对门槛同时用于三对，在三个总体均值相等且模型成立时，整组误报概率约为11.9253%。Tukey使用k＝3的q分布95%分位点3.5064261，对应本次均值差47.5893，使整组误报概率为5%。

三对绝对均值差的最大值就是极差，因此至少一对越界等价于最大差越界。q分布控制最大差越界的概率，就控制了整组误报概率。这里使用理论分布，重复抽样时MSE也重新估计；39.38与47.59是本次数据对应的门槛。三对比较共享组数据，不能按独立检验公式计算其联合误报概率。

品种3与2相差72，用同一个MSE计算，单对P≈0.0008514676，Tukey调整P≈0.0023769529。


## 整轮误报模拟与不等样本量SSA

三种语言主脚本新增`null_family_simulation.csv`：每轮生成3组各10个独立正态观测，总体均值均为720、总体方差均为1842，seed设20261102，共10000轮。每轮重算MSE与三对差异；任一对超过门槛，就把该轮记为有误报。分别比较逐对未调整t门槛与Tukey门槛。两个比例应在理论值约11.9253%与5%附近，有限次模拟会随机波动；三种语言的随机数算法不同，不要求逐轮相同。正文图05由R主脚本直接绘制，采用R的模拟结果；Python运行结果与独立配图另存，避免覆盖正文R图。

`ssa_demo.csv`逐步按组均值和人数计算：均值固定756、666、738；人数10/10/10时合并均值720、SSA45360；演示人数5/10/15时合并均值717、SSA40230。第二种为教学构造的汇总例子，演示先求加权总均值，再对均值差平方、乘本组人数并求和。


## R怎样数t规则的误报？

R脚本对每轮先算三对t与双侧P，`rowSums(sim_p < .05) > 0`将每轮标成有没有至少一对误报，再求和。结果1163轮（11.63%）；用最大q等价门槛得到同一计数，脚本已核对。Tukey规则有518轮（5.18%）。`null_trial_flags.csv`保留逐轮标记，`null_example.csv`列出第4轮的三对结果，`数据/模拟示例.csv`保存该轮30个值，供正文R代码复现。三组真实均值在模拟中同为720，故任一对被拒绝都属于误报。t使用三组共同MSE、df=27。

正文图05直接由R生成；Python绘制自己的`05-整组误报模拟-Python.png`，不覆盖R图。针对同一比较族，Tukey输出已调整P，无需再用p.adjust重复校正。
