# 统计课堂39：有序多分类的加权κ

对应文章：《轻度判成中度，与轻度判成重度，能算同样的不一致吗？》。本地配套项目；三种语言任选一种。

## 先下载、打开，再逐节运行

本次交付含独立ZIP。完整解压到一个文件夹，不要在压缩包内直接打开脚本。

- **R/RStudio**：双击`统计课堂39.Rproj`，工作目录自动为项目根。先打开`代码/weighted_kappa.R`，从第1节开始用Ctrl/Cmd+Enter按顺序运行，再用Source执行全文件。首次可在Console安装`install.packages(c("DescTools", "vcd"))`。
- **Python/Spyder**：打开解压后的项目文件夹（可用Projects → New Project → Existing directory），把工作目录设为项目根。打开`代码/weighted_kappa.py`，按`# %%`顺序用Ctrl+Enter逐节运行，再按F5全文件运行。先在Spyder实际使用的解释器环境安装`python -m pip install -r requirements.txt`；requirements列出分析和绘图依赖，不假定Spyder已装齐。脚本也兼容分节时没有`__file__`的运行方式。
- **MATLAB**：把**Current Folder**设为项目根，打开`代码/weighted_kappa.m`，按`%%`从第1节顺序Run Section，再Run运行全文件。第1节会把`代码/`加入函数搜索路径：搜索路径用于找到`kappa_parts.m`，Current Folder用于定位数据，两者作用不同。不需要额外工具箱。

新会话按这个顺序上手：第1节看原表→第2节看权重→第3.1节跟正文手算主例线性κ→第3.2节看六行完整比较和CI→第4节看Bootstrap→第5节按需保存。后续各节使用前面生成的变量；不要跳过数据准备。第一次先看配套自编函数的输入和输出，一般SE的内部推导可在第二遍学习。

命令行备选（均在项目根运行）：

```text
Rscript 代码/weighted_kappa.R
python 代码/weighted_kappa.py
matlab -batch "run('代码/weighted_kappa.m')"
```

## 这道题研究什么？

假设某医院从同期需要疼痛评估的就诊患者中随机抽取100人，每人仅纳入一次就诊，患者之间相互独立。两位评定者分别阅读每人的同次就诊标准化资料，互相盲评过去24小时疼痛影响日常活动的程度：无为没有疼痛；轻为有疼痛但不影响活动；中为影响活动但仍能完成；重为导致部分活动无法完成。每位评定者对每人给出一个等级，全部100人都有完整A/B评级。独立单位为患者，同一患者的两个结果必须配对保留。主例为教学假设。

本题计算完全同级比例、未加权κ、线性与二次加权κ，并给出两种加权κ的双侧95%Wald及成对Bootstrap百分位区间。

正文第七节介绍Nugent评分方法比较、ISUP病理分级和KL影像阅片三项真实研究的目标与资料结构；本项目的计算练习对应前六节的100人假设主例和对照，原始研究见正文参考文献7—9。

主例中有68人完全同级；另设一个对照，仅移动分歧，保持对角人数及全部行列合计相同。对照用于解释分歧位置，不是第二批实际患者或临床验证结果。

`数据/paired_tables.csv`是三语及数据图共享的唯一输入，每张表包含16格，含零格。脚本检查case为main/control、格子不重复、计数为非负整数、每张表总人数为100；读错文件或误改总人数会明确报错：

| 字段 | 意义 |
| --- | --- |
| case | main为主例，control为对照 |
| A | 行，评定者A的等级位置：0无、1轻、2中、3重 |
| B | 列，评定者B，使用相同定义和顺序 |
| count | 该联合格的人数；每张表合计100 |

本数据无缺失。实际研究需先定义完整配对分析集，交代缺失情况。不要把没有出现的等级自动删除；权重依赖完整的预设类别顺序。

## 核对什么结果？

先看第1节输入表，再看第2节的线性与二次权重。第3.1节逐步显示n、完全同级比例、每格比例、边际、独立参照、Po、Pe、κ；第3.2节显示六行完整比较、一般SE、95%Wald区间及Lo、Le、Lo/Le。第4节成对Bootstrap重抽10000轮；每轮重新计算边际，线性和二次使用同一批重抽表。每一节应整块运行；第一遍按节运行，不要把多行语句逐行拆开。

本节看哪个对象、先核对什么：

| 节 | IDE中可查看的对象 | 先核对 |
| --- | --- | --- |
| 1 | `tables` | 两张16格原表、总人数100 |
| 2 | `weights` | 三套同序一致权重 |
| 3.1 | `linear_steps`、`distance_steps` | n、68%、Po、Pe、κ；对应式4/5的距离0—3逐项人数、机会概率与贡献 |
| 3.2 | `distance_comparison`、`control_checks`、`results`、`function_comparison`、`report_text` | 对照四种距离计数及对角/行列合计不变检查；六行结果、原生/自编身份、完整报告句 |
| 4 | `bootstrap_results` | B、有效/无定义轮次、百分位区间 |

| 表 | 未加权κ | 线性κ | 二次κ |
| --- | ---: | ---: | ---: |
| 主例 | 0.573333 | 0.680000 | 0.776000 |
| 对照 | 0.573333 | 0.552000 | 0.520000 |

主例线性Po=0.866667、Pe=0.583333、SE≈0.052535、95%Wald区间[0.577034,0.782966]。二次Po≈0.937778、Pe≈0.722222、SE≈0.047104、区间[0.683678,0.868322]。完全同级比例仍为0.68。

正文还用`(1-Po)/(1-Pe)`解释对照表中κ的高低：线性规则的比值为0.448，二次为0.480；用1减去该比值，分别得到κ=0.552和0.520。第3节会直接计算、显示这一关系；三语的`summary.csv`也都保留Po、Pe、Lo、Le、loss_ratio和exact_agreement字段。

## 权重、函数与区间的身份

一致权重：未加权为对角1其余0，线性为`1-|i-j|/3`，二次为`1-(|i-j|/3)^2`；分歧损失为`1-w`。i/j取0—3。线性作为本例预设主分析；二次考察不同损失判断，不能按κ高低倒选规则。

- R原生`DescTools::CohenKappa()`显式指定`Unweighted`、`Equal-Spacing`、`Fleiss-Cohen`及`conf.level=.95`；另外用`vcd::Kappa()`的Weighted/ASE和`confint()`核对。DescTools默认未加权，vcd默认的加权行用线性，不能比较不同默认值。
- Python原生`sklearn.metrics.cohen_kappa_score()`计算点估计，显式`labels=[0,1,2,3]`和`weights=None/linear/quadratic`。它不提供本文的SE/区间；这些量由`kappa_parts.py`自编定义函数计算。
- MATLAB全部用`kappa_parts.m`自编实现，不声称存在内置κ接口。
- 三语`kappa_parts.*`使用相同一般多项delta方差，由联合格概率估计抽样波动。Wald区间采用标准正态0.975分位点，不是κ=0检验的SE₀。本篇不另外做假设检验。
- Bootstrap百分位取2.5%与97.5%，R用type=7，Python用method=linear，MATLAB自编位置`1+(有效轮次-1)q`处的线性插值。种子均为20261122，各语言随机流不同，区间端点不强求完全相同。
- 无定义轮次会计数报告；若全部无定义则停止。出现无定义不能悄悄把剩余轮次当作可靠推断。本次三语两种权重均10000轮有效、0轮无定义。

## 默认不写文件，怎样重建结果与图？

第一次只在Console/Command Window查看即可。第1节准备数据后已设Python `False`、R `FALSE`、MATLAB `false`，第5节不会覆盖它；在IDE日常保存时只改这一处为对应的`True/TRUE/true`。批量运行可设置`ARTICLE39_SAVE_RESULTS=1`覆盖为打开，结果写入`运行结果/对应语言/`，三语互不覆盖。

**三语分析主脚本以数表为结果，任选一种即可完成正文计算。** `draw_article.py`是独立的文章配图重建工具：公式图对应正文固定例题，数据图按共享CSV绘制，标题中的最远分歧人数也由输入算出。它不作为R或MATLAB分析的运行依赖。

文章图和7张公式可独立重建，不依赖旧运行结果：从项目根执行`python 代码/draw_article.py`，写入`文章配图/`。科学图从同一共享CSV读取计数。图中深蓝框为对角，橙框为相差三级；字体优先使用可用中文字体，缺少中文字体时需安装适合的字体再生成。

## 第二遍做一个有意义的改动

先保持三种预设规则，对照主例与control：哪些格子变化，距离各增加多少，为什么未加权κ不变？预测后再运行。

若要另拟损失规则，先说明依据，并列出按无、轻、中、重排列的4×4一致权重矩阵；矩阵须有限、各元素在0—1之间且对角线全为1（不强制对称）。第3.2节运行后，可把主例和该矩阵传给配套自编函数：R `kappa_parts(tables$main, new_weight)`、Python `kappa_parts(tables["main"], new_weight)`、MATLAB `kappa_parts(tables{1}, new_weight)`。输入是4×4计数和同序权重，返回Po、Pe、κ、一般SE、Wald下限、上限；主脚本的原生核对仍只针对预设三种规则。三语任选一种练习即可。

## 来源与范围

正文衔接第34篇，完整题干、解释及参考见同批交付的文章。Rosner 8版§10.8 pp.434—435仅给加权κ阅读入口；线性/二次公式与原生函数接口核对官方文档：

- https://andrisignorell.github.io/DescTools/reference/CohenKappa.html
- https://cran.r-project.org/web/packages/vcd/refman/vcd.html#Kappa
- https://scikit-learn.org/stable/modules/generated/sklearn.metrics.cohen_kappa_score.html
- https://support.sas.com/documentation/onlinedoc/stat/142/freq.pdf

拟上传目录为`stats_classroom/article39_weighted_kappa`，当前本地制作，未新增线上提交。经验分档不能代替权重依据、配对原表、区间和研究用途。
