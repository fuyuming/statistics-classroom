# 统计课堂35：关联强度与Cramér V区间

本项目沿用第11课的假设工程试验：300台独立心电记录器随机分成三种封装各100台，A为硅胶、B为环氧、C为热缩封装。弯折100次后接模拟器记录10分钟，没有持续达到1秒的数据缺口判合格。这些次数、时长和阈值均为课堂设定，每台设备只采用一种封装、接受一次判定。

共享数据在`data/device_pass.csv`，三行是A/B/C三个组，`pass`、`fail`两列是合格、不合格台数：`[[90,10],[80,20],[65,35]]`，不是三台设备。保留`data/`、`scripts/`、`helpers/`的相对位置，选择熟悉的一种语言即可。

## 从哪份脚本开始？

完整下载解压后，把本项目文件夹设为工作目录。

| 主题 | R | Python | MATLAB |
| --- | --- | --- | --- |
| Q、P、C、V及四格φ | `scripts/01_association_strength.R` | `scripts/01_association_strength.py` | `scripts/m01_association_strength.m` |
| 调用函数计算V区间、分级和报告 | `scripts/02_cramer_v_ci.R` | `scripts/02_cramer_v_ci.py` | `scripts/m02_cramer_v_ci.m` |

01按01.1—01.5运行；02可独立开始，不用先跑01。三语言02均按同样顺序分为六节：

1. 准备输入，查看原表。
2. 调用区间函数。
3. 读取点估计及上下限。
4. 用Funder–Ozer（2019）规则分级。
5. 结合原题解释结果、形成报告。
6. 按需保存。

首次从准备节顺序执行；后面的小节使用前面创建的对象。主脚本先调用完整函数，底层数值计算可在`helpers/`中选读。

### RStudio

1. 打开`article35_association_strength.Rproj`。
2. 首次需要时安装：`install.packages(c("effectsize", "rstatix", "report"))`。
3. 打开01，从01.1顺序运行，查看原表、`expected`、`point`、`phi_result`、`scaling`与`function_comparison`。
4. 打开02，从02.1顺序运行，查看`fit_v`、`ci`、`interpretation`及`report_result`。

02直接使用`effectsize::cramers_v(tab, adjust=FALSE, ci=.95, alternative="two.sided")`计算，`report::report()`生成文字报告。R主脚本不加载自编区间函数。`adjust=FALSE`选择与正文一致、未经小样本偏差修正的V，`alternative="two.sided"`明确采用双侧区间。

整文件备选：`Rscript scripts/01_association_strength.R`或`Rscript scripts/02_cramer_v_ci.R`。

### Python：Spyder

1. 把解压后的项目根目录设为Spyder工作目录，解释器需要NumPy、pandas、SciPy；需要时用该解释器执行`python -m pip install -r requirements.txt`。
2. 01按`# %% 01.1`至01.5运行，02按02.1至02.6运行，也可F5执行整文件。
3. 02调用配套自编`cramers_v_ci(tab, level=.95)`。完整函数在`helpers/association_ci.py`，返回字典；主脚本将其中的`V`、`V_L`、`V_U`整理为课堂结果表。

该函数使用SciPy的非中心χ²累计概率与求根接口，不调用R。02从项目根目录加载`helpers/`，分节运行不依赖`__file__`。

整文件备选：`python scripts/01_association_strength.py`或`python scripts/02_cramer_v_ci.py`。

### MATLAB

1. 把Current Folder设为项目根目录，需要Statistics and Machine Learning Toolbox。
2. 打开01或02，从第一个`%%`节开始Run Section；循环、条件块连同`end`完整执行。
3. 02调用配套自编`cramers_v_ci(tab,.95)`，函数在`helpers/cramers_v_ci.m`。主脚本从结果的`V`、`V_L`、`V_U`读取点估计、下限、上限。

该函数使用`ncx2cdf`与`fzero`完成数值计算，是本项目自编封装。准备节加载`helpers/`；若完整`run`暂时进入`scripts/`，脚本会在找不到输入时回到项目根目录。

整文件备选：`run('scripts/m01_association_strength.m')`或`run('scripts/m02_cramer_v_ci.m')`。

## 先核对什么？

- 原表n＝300，Q≈18.657938，df＝2，P≈0.0000888。
- Pearson C≈0.241974，Cramér V≈0.249385；A/C两组200台的有符号φ≈0.299342。
- V的双侧95%近似区间约为[0.1269,0.3570]。
- 按Funder–Ozer经验规则，下限、点估计、上限分别为small（小）、medium（中等）、large（大）。

02统一的`ci`表有四列：`V`、`CI_low`、`CI_high`、`confidence`；最后一列为0.95。R和自编函数的求根精度略有差异，本例按2e-7容差核对区间，四位小数相同。三组按设计各占1/3，本例V描述这一等权比较下的整体关联。

原题跑通后，可只在02的`tab`中重新排列A/B/C行次，再从02.2顺序运行，观察V和区间是否改变。

## 函数的分工

| 任务 | R | Python | MATLAB |
| --- | --- | --- | --- |
| 01点估计核对 | stats、effectsize与rstatix | SciPy的列联表检验及association | 手算与配套association_summary |
| 02点估计和区间 | effectsize现成接口 | 配套自编cramers_v_ci | 配套自编cramers_v_ci |
| 02效应分级 | effectsize::interpret_cramers_v | 配套interpret_v_funder | 配套interpret_v_funder |
| 02报告 | report包自动报告及结果解释 | 配套脚本组织中文报告句 | 配套脚本组织中文报告句 |

Python/MATLAB的函数使用相同的非中心χ²区间定义；`level=.95`表示双侧95%区间，点估计不作小样本偏差修正。`helpers/`中保留求根和边界处理，想了解内部算法时再阅读；`association_ci.R`是同一算法的选读版本，R主入口无需加载。这里使用的是近似区间。

## 可选保存

默认只显示结果。需要保存时，在相应IDE设置一次，再顺序运行脚本：

- R：`Sys.setenv(ARTICLE35_SAVE_RESULTS="1")`
- Python：`os.environ["ARTICLE35_SAVE_RESULTS"] = "1"`
- MATLAB：`setenv('ARTICLE35_SAVE_RESULTS','1')`

各语言只写自己的`R_`、`Python_`或`MATLAB_`前缀。01保存前三类结果，02保存后两类，`results/`附本例的15份参考CSV。

| 后缀 | 字段 |
| --- | --- |
| `point.csv` | n、r、s、df、m、Q、p、C、Cmax、V |
| `phi.csv` | a/b/c/d、n_sub、Q_sub、phi、phi_pearson、phi_row_flipped、V_sub |
| `scaling.csv` | multiplier、n、Q、p、Q_over_n、C、V |
| `ci.csv` | V、CI_low、CI_high、confidence |
| `interpretation.csv` | component（lower/estimate/upper）、V、label、rule |

官方说明：[effectsize](https://easystats.github.io/effectsize/reference/phi.html)、[report](https://easystats.github.io/report/reference/report.htest.html)、[rstatix](https://rpkgs.datanovia.com/rstatix/reference/cramer_v.html)、[SciPy association](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.contingency.association.html)、[SciPy ncx2](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ncx2.html)、[MATLAB ncx2cdf](https://www.mathworks.com/help/stats/ncx2cdf.html)。实际核验见[VALIDATION.md](VALIDATION.md)。
