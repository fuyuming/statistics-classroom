# 统计课堂37：分类计数数据怎样两两比较？——从原始检验到Hochberg调整

本项目对应“明哥的微生物世界”统计课堂37。数据是第11课P31/P33的教学例题计数，不是现实调查；比较对象先确定，再计算原始P和整批调整结果。只需选择一种语言，从第1节顺序运行到第9节。

## 先打开项目

- **RStudio**：双击根目录的 Rproj，再打开 scripts/article37_categorical_posthoc.R。用 Run 逐行运行，或在章节右上角运行本节。章节标题以四个短横线结束，可在文档大纲中导航。
- **Python（Spyder）**：在Spyder中通过“Projects→New Project”选择已有的本项目根目录，确认工作目录也指向这里；打开唯一主入口 scripts/article37_categorical_posthoc.py；Spyder中按# %%逐节用菜单“Run→Run cell”（macOS可用Cmd+Return，Windows/Linux可用Ctrl+Enter；需要移到下一节可用Shift+Enter），需要整文件时按F5。请使用安装依赖时的同一个Python环境。
- **MATLAB**：将 Current Folder 设为本项目根目录，打开 scripts/article37_categorical_posthoc.m；每个 %% 为一节，点击 Run Section。添加搜索路径并不等于改变 Current Folder。

三种入口均先检查 data/faults.csv 是否存在，避免在错误目录里读错数据。本篇每种语言只有一个主入口和4张结果表；项目合计6个helper，由相应主脚本调用，helpers是选读的自编函数。

## 一次性安装依赖

R控制台：

~~~r
install.packages(c("effectsize", "rcompanion", "rstatix", "report"))
~~~

Python在当前解释器环境的终端：

~~~text
python -m pip install -r requirements.txt
~~~

本项目使用`scipy.stats.contingency.odds_ratio`，Python环境请使用SciPy 1.10或更高版本。

MATLAB需要 **Statistics and Machine Learning Toolbox**，用于 fishertest、chi2cdf、ncx2cdf、normcdf、norminv；fzero 是基础MATLAB函数。本项目不使用其他工具箱。

## 九节练习与文章对应

| 节 | 操作 | 主要对象 / 应读结果 |
| --- | --- | --- |
| 1 | 读入两张独立试验表，设定95%区间与0.05判断水平 | tab为F/G/H故障构成；sparse_tab为A/B/C合格/不合格 |
| 2 | 整体Pearson检验与V，并用report组织文字结果 | main_overall：Q约44.792，V约0.3055；report用于组织整体检验的文字结果 |
| 3 | 三个2×3子表，各自重算预期频数 | full_pairs：三对Holm P均小于0.05 |
| 4 | 电极脱落与其他，两组比例比较 | category_pairs：Holm P为1、1、0.424934；G-H差12.5个百分点 |
| 5 | 各10台的独立试验，用双侧Fisher | sparse_pairs：三对最小E为1.5、4、4.5 |
| 6 | 并列比较两类问题的调整P | comparison：完整构成与指定类别比例的三对Holm P |
| 7 | 独立检验练习：先演算再比较Holm与Hochberg | 先用两项独立检验理论演算FWER（.049375、.000625、合计.05），再以`p_demo`给定test1/test2/test3的原始P=.020/.023/.040输出raw/holm/hochberg；另以`p_compare`的`.001/.020/.040`演示最大P首关通过时三项均拒绝，调整P仍从后往前取累计最小；两组都不来自F-G/F-H/G-H计数表 |
| 8 | 报告G-H比例差效应量 | 比例差、单项区间与同一三对族的Holm P |
| 9 | 可选保存四张结果表 | 默认不保存，启用后写到results/对应语言 |

faults.csv中 open、electrode、board 分别为接头开路、电极脱落、主板连接故障。每行80台都是已失效记录器，每台只有一个主要故障；分析故障构成，不是全部投入使用设备的失效率。sparse.csv是另取30台独立记录器的试验，A硅胶灌封、B环氧灌封、C热缩套管各10台；每台弯折100次后连接心电模拟器记录10 min，期间无持续≥1 s的数据缺口判为合格。合格/不合格计数分别为1/9、2/8、7/3。三对合格率作为同一比较族，采用双侧检验、Holm调整及单项95%区间。

## 函数、参数和字段怎么读

- R使用原生 chisq.test、prop.test、fisher.test、p.adjust；第2节另用 report::report 组织整体检验文字；第3—5节还实际调用正文的 rcompanion / rstatix 入口，与逐对子表结果核对。第7节先用两项独立检验的P~Uniform(0,1)模型演算至少一个P≤.025、两P都在(.025,.05]及两项H0真时FWER，再用 `stats::p.adjust(method = "hochberg")` 与 Holm 并列；Hochberg目标是控制FWER，适用于相互独立或符合特定正依赖条件的检验，本练习采用相互独立前提，原三组对分析继续使用Holm。effectsize::cramers_v 的 adjust=FALSE 对应未作小样本偏差修正的课堂V，alternative="two.sided"请求双侧区间。
- 课堂单次两率比较与批量入口相接：pairwise_prop_test 经 pairwise.prop.test 逐对使用 prop.test；pairwise_fisher_test 经 rstatix::fisher_test，底层使用 stats::fisher.test。若分析方案选用无条件精确率差，可逐对用 exact2x2::uncondExact2x2（parmtype="difference", nullparm=0），再汇总原始P作p.adjust。该函数的率差方向为第二组减第一组，conf.int=TRUE请求率差区间；fisher.test的conf.int则针对优势比OR。本项目执行第4节的近似率比较和第5节的Fisher；无条件精确率差在正文作为课堂函数入口提示。
- Python使用SciPy的chi2_contingency、fisher_exact、contingency.odds_ratio(kind="conditional")，statsmodels的multipletests；第7节用 `multipletests(method="simes-hochberg")` 对应R的Hochberg；intervals.py提供比例差区间函数。
- MATLAB使用原生fishertest得到精确P；第7节在本节内分步排序、乘m,m−1,…,1、反向累计最小并还原原序，不新增helper；fisher_or_ci、rate_difference_ci等是自编函数，其名称与文件一一对应；true/false决定2×2是否采用Yates校正。

所有三对结果的顺序为F-G、F-H、G-H；稀疏表为A-B、A-C、B-C。

| 字段 | 含义 |
| --- | --- |
| q / q_yates、df | Pearson Q / Yates校正Q，以及自由度 |
| p_raw、p_holm | 单项原始P、该分支完整比较族的Holm调整P |
| min_expected | 各对子表自己的最小预期频数 |
| V、V_low、V_high | 未修正Cramér V及单项双侧95%非中心χ²近似区间 |
| p1、p2、diff | 前组、后组的电极脱落比例，diff=p1-p2 |
| diff_low、diff_high | 带连续性校正的正态近似比例差单项95%区间，与R prop.test(correct=TRUE)一致 |
| OR、OR_low、OR_high | 固定边际条件Fisher检验的条件优势比及单项双侧区间 |

比例差乘100，单位是百分点。V、比例差和OR区间均为**单项区间**，没有构造三项同时区间；P值另按同一问题的三对结果作Holm判断。这里的条件OR采用前一组相对于后一组的合格优势，OR>1表示前一组的合格优势更高；区间跨1时不能据此明确方向。固定行、列合计后的条件OR估计，与直接用两组样本优势之比估计OR的数值可以不同：本例A–C的样本OR为0.047619，Fisher条件OR约为0.0579，结果表使用后者。

## 显示与保存

默认整文件或逐节运行仅显示结果，不修改项目文件。需要导出时，在各自控制台启用开关，再运行第9节：

~~~r
Sys.setenv(ARTICLE37_SAVE_RESULTS = "1")
~~~

~~~python
os.environ["ARTICLE37_SAVE_RESULTS"] = "1"
~~~

~~~matlab
setenv('ARTICLE37_SAVE_RESULTS', '1')
~~~

第9节只保存 main_overall、full_pairs、category_pairs、sparse_pairs 四张CSV，每种语言的列名、标签、顺序一致；第7节练习只在控制台显示，不新增CSV。恢复默认分别用 Sys.unsetenv("ARTICLE37_SAVE_RESULTS")、os.environ.pop("ARTICLE37_SAVE_RESULTS", None)、setenv('ARTICLE37_SAVE_RESULTS','')。

终端整文件运行也要从项目根目录开始：

~~~text
Rscript scripts/article37_categorical_posthoc.R
python scripts/article37_categorical_posthoc.py
matlab -batch "addpath('scripts'); article37_categorical_posthoc"
~~~

再跑一次时，可以做以下练习，先预测再核对：

- **Hochberg练习**：在第7节把第三个原始P从`0.040`改为`0.060`，先预测再运行。预期Holm为`.060/.060/.060`，Hochberg为`.046/.046/.060`；这是独立教学向量的演算练习，不改变原三组对分析。
- **Hochberg首关对照**：第7节另有`p_compare = .001/.020/.040`的独立试验向量。最大P `.040≤.05`，首关通过意味着拒绝全部三项；从后往前取累计最小后，Hochberg调整P预期为`.003/.040/.040`，三项均拒绝。不要把它与前面的`p_demo`或原三组对分析合并。
- **区间水平**：将第1节的 conf_level 从0.95改为0.99，从第1节顺序重跑到第8节。预期区间变宽，点估计与检验P值保持不变；控制台的置信水平标签也应显示99%。
- **作差方向**：将G–H单次比例比较的组序换成H–G。预期差为−12.5个百分点，95%区间变为［−28.63，3.63］个百分点，双侧P值不变。

后两项练习分别改变区间水平与报告方向，保留原来的观察计数。

### 单次H–G比例比较（完成前面主脚本后执行）

先完成主脚本，再在下方单次比较中把`conf_level`设为0.95。输入组序均为先H、后G，因此比例差是H−G；这里的P是单次比较的原始P。预期比例差`-0.125`，95%区间`[-0.286273, 0.036273]`，Yates双侧P=`0.1416447`。

R：

~~~r
conf_level <- 0.95
hg <- stats::prop.test(c(25, 35), c(80, 80), correct = TRUE,
                        conf.level = conf_level)
c(diff = unname(hg$estimate[1] - hg$estimate[2]),
  low = hg$conf.int[1], high = hg$conf.int[2], p = hg$p.value)
~~~

Python：

~~~python
conf_level = 0.95
hg_tab = np.array([[25, 55], [35, 45]])  # 行顺序H、G；列为电极脱落、其他
q_hg, p_hg, df_hg, _ = stats.chi2_contingency(hg_tab, correction=True)
hg_p1, hg_p2, hg_diff, hg_low, hg_high = rate_difference_ci(
    np.array([25, 35]), np.array([80, 80]), conf_level)
print(hg_diff, hg_low, hg_high, p_hg)
~~~

MATLAB：

~~~matlab
conf_level = 0.95;
hg_tab = [25 55; 35 45];  % 行顺序H、G
[q_hg, df_hg, p_hg] = pearson_test(hg_tab, true);
[hg_p1, hg_p2, hg_diff, hg_low, hg_high] = ...
    rate_difference_ci(hg_tab(:, 1), sum(hg_tab, 2), conf_level);
disp([hg_diff, hg_low, hg_high, p_hg]);
~~~

## 来源与核对

主例来自第11次课P31，辅助小样本例来自P33；两者均为课堂教学计数，用于练习本篇比较方法。方法和函数说明：
[R卡方](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/chisq.test.html)、
[R prop.test区间](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/prop.test.html)、
[effectsize的V](https://easystats.github.io/effectsize/reference/phi.html)、
[report检验报告](https://easystats.github.io/report/reference/report.htest.html)、
[rcompanion手册](https://cran.r-project.org/web/packages/rcompanion/rcompanion.pdf)、
[rstatix比例检验](https://rpkgs.datanovia.com/rstatix/reference/prop_test.html)、
[rstatix Fisher检验](https://rpkgs.datanovia.com/rstatix/reference/fisher_test.html)、
[SciPy条件OR](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.contingency.odds_ratio.html)。

随附结果来自本篇脚本的实际运行。运行环境与三语对齐结果见文章交付核对记录；不以事先填入的数值代替计算。

本次实跑：R 4.6.1，effectsize 1.0.3、rcompanion 2.5.4、rstatix 1.1.0、report 0.6.4；Python 3.12.12，NumPy 2.4.3、pandas 3.0.1、SciPy 1.17.1、statsmodels 0.14.6；MATLAB R2025a。P值等统计量及比例差区间跨语言误差小于1e-10。R包的V区间采用优化分位点的实现，另两语直接求累计概率方程的根；本次端点最大差约1.35e-6，四位小数一致。条件OR及精确区间由各语言实际求解；R原生函数的默认求根精度使OR最大相差约0.000023、区间端点约0.000234，表中三位小数一致。
