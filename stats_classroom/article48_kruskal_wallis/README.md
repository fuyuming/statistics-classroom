# 统计课堂48：Kruskal–Wallis H与零分布

本项目配合《H统计量从哪里来，精确P值与χ²近似怎样算？》。用自己熟悉的一种语言完成即可。

## 下载解压后，从哪里开始

完整解压项目，保留“代码、数据、文章配图、运行结果”相对位置。项目根目录是包含本README和Rproj的文件夹。

- **RStudio**：双击根目录`统计课堂48_KW.Rproj`，打开`代码/kw_practice.R`。按第1—5节顺序选中一节、Ctrl/Cmd+Enter执行。整文件可点Source，或在Console运行`source("代码/kw_practice.R")`。统计计算仅需base R；绘图需要ggplot2，可先运行`install.packages("ggplot2")`。
- **Python / Spyder**：以解压后的根目录建立文件夹项目，确认工作目录为该根目录，打开`代码/kw_practice.py`。先运行准备节，再按`# %%`顺序Ctrl+Enter运行当前单元；F5运行整文件。安装依赖须使用Spyder当前解释器，可在相同环境运行`python -m pip install -r requirements.txt`。命令行备选：在根目录运行`python 代码/kw_practice.py`。顺序分节不依赖`__file__`。
- **MATLAB**：将Current Folder设为项目根，打开`代码/kw_practice.m`，从第1节依次Run Section，再按需Run整文件。Current Folder决定相对数据路径；添加代码到搜索路径不能代替设置Current Folder。需要Statistics and Machine Learning Toolbox。命令行备选：根目录执行`matlab -batch "run('代码/kw_practice.m')"`。

第1节读取数据和设置参数；第2节推算兔例H；第3节枚举六瓶条件零分布；第4节进行兔例Monte Carlo；第5节绘图；第6节给出修改练习。后节使用前节中间量，不能跳过准备节。三种语言分别保存结果，不互相覆盖。

## 数据与独立单位

`数据/rabbit.csv`每行一只兔，24行，四组各6只。`drug`为药物名称，`score`为改善指标，单位为分。来源Rosner, *Fundamentals of Biostatistics*, 8e, 表12.17，p597。每兔一眼药物、一眼盐水，10分钟后两眼接受花生四烯酸钠挑战。本项目取挑战后15分钟资料，指标为“(基线−15分钟评分)药物眼−(基线−15分钟评分)盐水眼”，越大表示相对改善越多。独立单位是兔。输入保留当前课程CSV的原始24个值；原数据无缺失，项目要求完整资料。

`数据/six_bottles.csv`为假设教学数据，每行一独立培养瓶。`condition`为培养条件A/B/C，各2瓶；`concentration_mg_L`为同一培养时点的产物浓度（mg/L），依次1、2、3、4、5、6。用于展示完整枚举，不是真实培养实验。

## 先核对哪些结果

兔例四组平均秩依次16.25、14.1666667、15.25、4.3333333；秩和97.5、85、91.5、26。B=546.5833333，总秩平方和1065，实际秩方差46.30434783，C=0.926086956522，H*=10.93166667，H=11.80414710。三语原生KW函数采用χ²近似，df=3，P=0.008085087714。

六瓶原分组A={1,2}、B={3,4}、C={5,6}得到H=32/7，完整身份划分90个，尾部含相等有6个；精确P=0.066666666667，χ²(2)近似P=0.101701392304。完整H分布存入每语的`aux_exact_distribution.csv`。

兔例199999次独立均匀置换，固定每组大小，各次可抽到重复划分。每次重排有身份单位的秩；并列频数固定。用`(b+1)/(模拟次数+1)`报告Monte Carlo P，用`sqrt(P*(1-P)/模拟次数)`近似MC标准误。seed=20261130仅便于在同一语言环境复算，不保证跨语言相同序列。

### 并列秩的方差链条与三个对比坐标

第2节还打印并列校正的完整方差链条：兔例 `Vbase=N(N+1)/12=50`，无并列总秩平方和 `Tbase=(N-1)Vbase=1150`；评分频数 `1,5,5,4,9` 的并列损失 `(t³−t)/12` 分别为 `0,10,10,5,60`，所以 `Tactual=1065`、`Vactual=Tactual/(N-1)=46.3043478261`，且 `C=Vactual/Vbase`。这说明并列值减少了可用的秩变异，C不是凭空加入的除数。

当四组样本量相等（本兔例均为6）时，代码另把四个平均秩写成三个逐步对比坐标：`m1−m2`、`(m1+m2)/2−m3`、`(m1+m2+m3)/3−m4`。各自按 `sqrt(2V/n)`、`sqrt(3V/(2n))`、`sqrt(4V/(3n))` 标准化，输出差值、SE、Z及`Z²`贡献到每语`运行结果/contrast_components.csv`；三项`Z²`之和应与H相同。这是理解H分解的坐标演示，不是三次独立假设检验，也不应分别报告P值。若输入不是四组等样本量，代码会跳过该表而保留其余计算。

本次Python尾部679次，P=0.003400000，MCSE=0.000130163；R与MATLAB各644次，P=0.003225000，MCSE=0.000126780。比较模拟输出时，差异应按合并MC标准误衡量；确定性的H、C及完整枚举分布应在浮点容差内相同。有限次模拟不等于完整精确值。

## 图、保存与输出

正文主图采用R＋ggplot2生成。三种语言默认显示平均秩图及六瓶精确CDF/χ² CDF对照，不自动保存图片或CSV。把第1节save_figures设为真可保存PNG，把save_results设为真可保存CSV；也可运行前设置环境变量SAVE_FIGURES=1、SAVE_RESULTS=1。自动检查使用KW_HEADLESS=1关闭显示，Python另用MPLBACKEND=Agg；学生在IDE中无需设置这些变量。

开启save_results后，数表写入运行结果/语言/，便于逐项核对；默认只打印，不创建结果目录。附带的console.txt是本次实际整文件运行输出。summary.csv汇总H、C、近似P及MC结果，rabbit_ranks.csv为组秩汇总。图片的英文标签便于三环境移植，正文提供中文解释。CDF纵轴均为累计概率，精确分布用阶梯，观察H用竖线；不是把概率质量与密度直接比较。

## 再改一次分组

先备份六瓶CSV。保留浓度和每组2瓶，将浓度1至6对应的`condition`依次改为`A,B,C,C,B,A`，然后在新会话运行第1、3节（或整文件）。三组平均秩均为3.5，应得H=0、精确P=1，脚本判断随当前输入改变。重跑第5节须先运行第2节。完成后恢复原CSV。改变分组是学习统计量的实验，真实研究应按实际设计分析。

## 实际环境与验证

本次使用Python 3.12.12、NumPy 2.4.3、SciPy 1.17.1、pandas 3.0.1、Matplotlib 3.10.8；R 4.6.1、ggplot2 4.0.3；MATLAB R2025a及Statistics and Machine Learning Toolbox。具体数字证据见附带实际结果，核验范围见随稿核对记录。较旧环境的输出排版可能不同，先核对数值与统计定义。


## 可选：直接调用R包核对

原`代码/kw_practice.R`的统计计算仍只需基础R。若要对照现成函数，先在自己的R环境运行`install.packages(c("coin", "kSamples", "rstatix"))`，再从项目根目录打开`代码/kw_R包核对.R`，按第1—4节顺序运行，或`source("代码/kw_R包核对.R")`。默认打印结果，不保存CSV；运行前设`SAVE_RESULTS=1`可保存比较表。

`stats::kruskal.test()`和`rstatix::kruskal_test()`采用χ²近似，不会像`wilcox.test(exact=NULL)`那样自动切换精确/渐近方法。当前核对环境为R4.6.1，wilcox.test的新版本已支持并列秩的条件精确推断；不要将旧版本的并列限制直接套用到新版本。

`kSamples::qn.test(test="KW", method="exact", Nsim=90)`完整枚举六瓶90种划分，返回H=32/7、χ²近似P=0.101701392304及精确P=1/15；必须检查返回的method实际为exact。Nsim不足合法划分数时，它会改用模拟。

兔例用coin的approximate(nresample=199999)和kSamples的simulated/Nsim=199999，分别从seed=20261130开始。本次coin尾部615次、P=615/199999=0.003075015375、MCSE约0.000123806；kSamples尾部649次、P=649/199999=0.003245016225、MCSE约0.000127171。当前包返回b/B；同一批抽样按正文加1规则换算，分别为0.003080和0.003250。包的Monte Carlo区间衡量计算波动；不是药物效应的置信区间。

本次包版本coin1.4.5、kSamples1.2.12、rstatix1.1.0；kSamples依赖SuppDists1.1.9.9。其他版本先核对算法、实际method与MC次数，再按计算误差比较结果。包函数使用同一份输入，既有三语教学入口及其结果保留。
