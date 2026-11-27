# 统计课堂45｜Welch与Games–Howell新手练习

先完整解压本项目，再选一种熟悉的语言运行。三套脚本使用相同CSV、α=0.05、全部三对比较与近似95%同时区间，差值统一为后组减前组。

## 第一次从哪里开始

1. **R/RStudio**：双击`统计课堂45.Rproj`，打开`R/welch_gh.R`。从“0.准备”开始，在同一新会话中选中一节并运行，依次到最后一节；先看Environment中的中间量、Console结果，再到Plots看图。整文件可点Source，或在Console输入`source("R/welch_gh.R")`。仅用R自带函数，无需安装扩展包。
2. **Python/Spyder**：将工作目录设为解压后的`练习项目`，打开`Python/welch_gh.py`，按`# %%`单元从0节向下运行。在Variable Explorer查看`summary_tab`、`welch_tab`和`gh_tab`，在Plots查看两图；也可按F5执行整文件。新控制台逐节运行不依赖`__file__`。命令行备选为`python Python/welch_gh.py`。
3. **MATLAB**：将Current Folder设为`练习项目`，打开`MATLAB/welch_gh.m`，从0节依次点击Run Section；整文件点Run。把目录加入搜索路径与设置Current Folder是不同操作。需要Statistics and Machine Learning Toolbox；当前验证版本为R2025a。

后面的节使用前面算出的变量。若刚启动新会话，要从第0节开始；若改了数据或`unit_scale`，也从第0节重新向下运行。

Python实际依赖是NumPy、pandas、SciPy≥1.16、matplotlib。请在Spyder实际使用的解释器环境中安装，备选命令为`python -m pip install -r requirements.txt`。SciPy两个`equal_var=False`接口需要1.16及以后版本。

## 数据和题目

`data/hetero_normal.csv`每行一名受试者，`group`为G1/G2/G3，`value`为基线收缩压减去8周收缩压，单位mmHg；各组不同受试者12人，共36人。不同受试者独立，每人只有一个用于分析的下降值，跨组没有配对关系。CSV无缺失值；脚本不默默填补缺失，要求每组至少2个完整观测且方差为正。

这份数据来自课堂现用教学模拟，原模拟种子为20261130，三组生成均数10/12/16、标准差1/3/6，保留两位小数。它不是实际临床资料。推断只使用CSV中的给定样本，不把生成参数当成已知总体参数，也不以三语言相同种子重生成数据。

题目设各组近似正态，允许总体方差不同。先用Welch检验三个总体均数是否相等，再用Games–Howell比较全部三对。独立性由设计判断；放宽共同方差不取消组内分布、极端值等条件。

## 各节算什么，先核对哪些结果

- 读入后显示字段、样本量和均数/SD。样本SD采用n−1分母。
- `v=s²/n`是均数方差估计，`w=1/v`是精度权重，按权重求共同均数的估计；随后算A、C、Welch F、近似自由度、临界F与右尾P。
- GH逐对沿用两独立样本Welch t的SE和Welch–Satterthwaite自由度，先算t，再由q＝√2×|t|转换尺度，计算调整P和同时区间。R以`ptukey/qtukey`计算，Python以`studentized_range`计算，并对照完整GH接口的相同方向矩阵元素。
- 本题预期Welch F≈11.868577、分母自由度≈15.562533、P≈0.0007405；右尾α=0.05、自由度(2,15.562533)的临界F≈3.654113。G2−G1、G3−G1、G3−G2差值分别约2.4942、8.1567、5.6625 mmHg，区间分别约[−0.0997,5.0880]、[3.0864,13.2269]、[0.2561,11.0689] mmHg。

脚本中的结论由当前结果生成。P和区间均已纳入GH调整，不再追加Bonferroni。R原生`oneway.test`返回整体F/df/P；Python的`tukey_hsd(..., equal_var=False)`执行GH。

**MATLAB实现限制**：Welch与GH按正文定义自编，学生化极差使用R2025a的`internal.stats.stdrinv/stdrcdf`。这是非公开工具箱接口，不是MATLAB公开原生GH函数；其他版本可能不可用，需复核接口和数值。项目不调用R/Python作为隐藏后端。

## 第二次：只改变单位

保持全部原始数据、分组、α和比较族不变，将第0节`unit_scale`从1改成0.133322，统一把mmHg换成kPa，重新运行全部节。

均数、SD、差值、SE和区间端点乘0.133322，方差乘其平方，权重除以其平方；A、C、F、自由度、q和P保持不变。换算发生在读入之后，后续每项计算都重算，输出与图的单位也随之改变。这个练习用于理解尺度，不用于寻找显著结果。

## 输出与保存

默认`save_outputs`关闭，只在IDE查看。R命令行整文件运行时只打印结果；开启保存后输出PNG，不调用默认PDF设备。需要保存时，可把第0节该变量直接设为`TRUE`/`True`/`true`（对应R/Python/MATLAB），或在启动前设环境变量`SAVE_OUTPUTS=true`。

- `运行结果/R`、`运行结果/Python`、`运行结果/MATLAB`：每种语言的组摘要、Welch中间量和全部GH结果CSV，均带表头。
- `图形/`各语言子目录：个体点与均数图、三对同时区间图。图题、坐标轴和图例统一使用中文；方案代码G1/G2/G3及单位mmHg/kPa保留。
- 图默认在IDE显示；自动化验证用无界面后端，不强制学生使用它。

交付结果是原始mmHg输入的实跑结果，单位练习的结果不覆盖本例交付CSV。运行记录范围见`VALIDATION.md`。

## 中文绘图字体

三语脚本的图题、坐标轴和图例均使用中文，并处理中文字体。macOS优先使用苹方/黑体，Windows使用微软雅黑，Linux可安装Noto Sans CJK SC；缺少中文字体时应先安装，重启相应IDE后重跑绘图节。Python/MATLAB会从候选字体中选择本机可用字体；R按系统设定`plot_font`，可改为本机已安装的中文字体。UTF-8保证源文本编码，中文字体保证图中字形显示，两者都要满足。

## 用R重画公众号正文数据图

公众号正文的两幅数据图统一使用`R/plot_article.R`绘制，采用共享CSV并在R中重新计算均数和GH区间。打开R项目后，先运行`source("R/plot_article.R")`，再运行`draw_article()`；两图与图2区间坐标CSV保存到项目内`文章配图`。函数的三个输入分别是数据文件、图像目录、坐标文件，均有可直接使用的默认值。命令行也可运行`Rscript R/plot_article.R`。

`R/welch_gh.R`及Python/MATLAB主入口用于逐步计算和练习；公众号使用的正式数据图以这份R绘图实现为准。

## q分布参照的同题对照

程序中的`P_single`是单对Welch t双侧P；`P_q_k2`是在相同SE和自由度下、以k=2的q参照计算的尾概率，用于看清t与q的尺度对应；`P_adj`才是本题全部三对采用k=3的GH调整P。G3−G1的前两项均约0.001128，第三项约0.002974。三语均从输入计算这些列；单位换算后，三项P均不变。

R正文PNG使用R自带Cairo设备；运行`capabilities("cairo")`应返回`TRUE`，并保证已安装中文字体。本轮已在上述R环境验证。
