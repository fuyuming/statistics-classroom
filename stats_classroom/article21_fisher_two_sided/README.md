# 统计课堂21：Fisher精确检验的双侧P值

沿用第11课P3—11、2025分类资料例7-4：33名HBsAg阳性孕妇分为预防注射组与非预防组，比较新生儿感染结局。本项目只有课堂汇总计数，没有个人记录。本篇提供完整R、Python、MATLAB代码与复现结果。

## 下载后怎样开始

完整解压并保留“数据、代码、运行结果”目录，用熟悉的一种语言跑通即可。

- **RStudio**：打开`统计课堂21-Fisher双侧P值.Rproj`，再打开`代码/lesson21_fisher_two_sided.R`，Source全文。基础R即可，不需要扩展包。
- **Spyder**：打开本文件夹，打开`代码/lesson21_fisher_two_sided.py`，F5运行，或按六个`# %%`单元顺序执行。需要numpy、pandas、scipy、matplotlib；在Spyder实际使用的解释器环境安装`requirements.txt`中的依赖。命令行为`python -m pip install -r requirements.txt`。
- **MATLAB**：将Current Folder设为本项目目录，打开`代码/lesson21_fisher_two_sided.m`，Run全文，或依次执行六节。只加入搜索路径不等于设置Current Folder。脚本使用基础MATLAB，不依赖统计工具箱；本次运行版本为R2025a。

每份主脚本分为六部分：读表与检验、逐表概率、单双侧P、OR及区间、交换行列练习、保存结果。先从头运行，避免沿用工作区中的旧变量。R在Environment、Python在Variable Explorer、MATLAB在Workspace中可查看中间对象。

## 数据字段和检验定义

`数据/hbv_prevention.csv`有两行，`group`为Prevention（预防注射）和Control（非预防），`positive`和`negative`为感染阳性/阴性人数。原表为`[[4,18],[5,6]]`。每个观察单位属于一行一列；使用独立母婴单位、二分类结局的模型。计数必须为非负整数，零格合法，缺失值不能填作0；本项目要求每行列合计均大于0。

原假设为两组总体感染率相等。给定行合计22、11和列合计9、24，第一格X可取0至9。每张表为`[[x,22-x],[9-x,2+x]]`，其超几何概率为`choose(9,x)*choose(24,22-x)/choose(33,22)`。

概率排序双侧：累加概率不大于观察表概率的表，包含相等的概率。程序采用1e-7相对容差以包含浮点运算中的并列值。另一种central双侧定义为`min(1,2*min(left_P,right_P))`，与概率排序定义分别保存。

## 先核对哪些结果

| 对象 | 原题结果 |
| --- | --- |
| 全部10张表的概率和 | 1 |
| 观察表概率q₄ | 0.0876272781723 |
| 概率排序计入的X | 0、1、2、3、4、9 |
| 左单侧P，X≤4 | 0.108147229115 |
| 右单侧P，X≥4 | 0.979480049057 |
| 双倍较小单侧P（central） | 0.216294458230 |
| 概率排序双侧P | 0.121044750578 |
| 样本OR | 0.266666666667 |
| 条件OR、等尾95%区间 | 约0.279，[0.040,1.767] |

R与SciPy都用标准Fisher函数核对概率枚举。SciPy的`fisher_exact().statistic`为样本OR；另外用`contingency.odds_ratio(kind="conditional")`取得条件OR。R的`fisher.test()`直接给条件OR。MATLAB用非中心超几何权重计算条件MLE，并由上下尾概率分别等于0.025求等尾区间。其条件MLE使模型期望X等于观察值4；这部分辅助函数供进一步学习。

**P值与区间方法分别说明。** 本项目主P值是概率排序双侧；默认OR区间是等尾区间，与central双侧规则对应，不能把它称为概率排序P值的严格反演区间。R默认求根精度使条件OR及端点末位略异于SciPy/MATLAB，差别不影响文中三位小数的报告；详见VALIDATION。

## 输出文件与练习答案

各语言分别保存`运行结果/R`、`运行结果/Python`、`运行结果/MATLAB`，每种语言三份CSV。

| 文件 | 内容 |
| --- | --- |
| enumeration.csv | 10张表四格计数、概率和是否计入双侧P；X为第一格 |
| summary.csv | 边际、单双侧P、两种OR、条件OR区间、两组观察率 |
| swap_exercise.csv | 原表、交换两行、交换两列后的P和OR结果 |

练习：只交换两行，数据内容保持不变。双侧P保持0.12104475；左、右单侧P互换。样本OR从0.2666667变为3.75，条件OR约变为3.583；区间按`[1/原上限,1/原下限]`反转。交换列也有相同的倒数关系；解释对象随阳性/阴性方向一起改变。不能改变检验方向来追求较小的P值。

## 两幅图怎样复现

Python主脚本末尾自动调用`代码/plot_fisher_two_sided.py`。它从共享计数重新计算确定性超几何概率，输出到`文章配图/`：

- `01-超几何概率与双侧范围.png`：完整概率分布与4、8、9三根柱的单独对照；CSV保留全部概率、是否计入、观察概率门槛。
- `02-观察表与右侧两表.png`：三张具体计数表及相同边际；CSV保留四格计数、概率和纳入状态。

绘图需要matplotlib和中文字体（如Arial Unicode MS或Noto Sans CJK SC）。仅运行R/MATLAB数值计算不需要Python。图没有模拟误差，不需要随机种子。图中的预防/非预防标签对应原CSV行顺序；交换练习在结果表内完成，不改共享输入标签。

本篇完整项目：https://github.com/fuyuming/statistics-classroom/tree/main/stats_classroom/article21_fisher_two_sided

完整仓库：https://github.com/fuyuming/statistics-classroom

仓库ZIP：https://github.com/fuyuming/statistics-classroom/archive/refs/heads/main.zip

本篇上传后完整解压仓库ZIP，进入`statistics-classroom-main/stats_classroom/article21_fisher_two_sided/`，先读本README。

方法来源：R fisher.test文档 https://stat.ethz.ch/R-manual/R-devel/library/stats/html/fisher.test.html ；R超几何文档 https://stat.ethz.ch/R-manual/R-devel/library/stats/html/Hypergeometric.html ；SciPy fisher_exact https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.fisher_exact.html ；双侧及匹配区间的讨论 https://cran.r-project.org/web/packages/exact2x2/vignettes/exact2x2.pdf 。
