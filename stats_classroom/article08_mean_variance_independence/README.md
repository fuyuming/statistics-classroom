# 均值与方差为什么独立？

明哥的微生物世界 · 统计课堂。本课从计算公式出发，解释抽样独立性、联合正态、协方差和t统计量。目录接续article07；正文暂不标系列编号，待实际发表顺序确认。

## 下载、打开、运行

本篇：https://github.com/fuyuming/statistics-classroom/tree/main/stats_classroom/article08_mean_variance_independence

完整ZIP：https://github.com/fuyuming/statistics-classroom/archive/refs/heads/main.zip

下载并完整解压，进入 `stats_classroom/article08_mean_variance_independence`。普通练习不用重新生成数据，三种语言直接读取附带CSV。

- **R/RStudio**：双击 `均值与方差的独立性.Rproj`，打开 `代码/mean_variance.R`，点击Source。仅基础R，无需额外包。工作目录应是本篇目录；命令行为 `Rscript 代码/mean_variance.R`。
- **Python/Spyder**：Projects → New Project → Existing directory，选择本篇目录；打开 `代码/mean_variance.py`，按F5。需NumPy、pandas、SciPy、Matplotlib；在实际解释器环境运行 `python -m pip install -r requirements.txt`。脚本自动调用同目录 `plot_mean_variance.py`。命令行为 `python 代码/mean_variance.py`。
- **MATLAB**：把Current Folder切到本篇目录，打开 `代码/mean_variance.m` 点击Run，或 `run('代码/mean_variance.m')`。`tinv`需要Statistics and Machine Learning Toolbox。Current Folder不是全局搜索路径，无需全局添加本课目录。

输出在 `运行结果/Python`、`运行结果/R`、`运行结果/MATLAB`，重复运行覆盖对应结果，不修改输入。Python生成正文的两张图；R、MATLAB在 `文章配图` 下各自子目录生成对应英文标签图，统计内容相同。Python检测常用中文字体，缺少时自动使用英文标题。文件使用UTF-8；字体缺字与源文件编码错误是两类问题。

## 按什么顺序学习

1. 人工数据：均值应为10、20、10，样本方差为4、4、16。
2. n=2的M、D：核对S²=D²/2，以及经验Cov(M,D)=[Var(X₁)−Var(X₂)]/2。有限样本的两列经验方差不同，所以经验Cov(M,D)不必为0；恒等式本身仍成立。
3. n=10重复抽样：对比正态总体与平移指数总体的均值—方差关系；按固定门槛筛选均值，查看样本方差的条件摘要。
4. t构造：逐行核对两种代数表达相等。正态模型才有这里使用的精确t(9)分布保证；代数相等本身不能保证分布。
5. 图2：比较独立正态、Y=X²、Y=随机正负X。理论协方差都为0，只有第一种独立。

每节打印答案与解释。Python用 `# %%`、MATLAB用 `%%` 分节；首次请完整运行，之后再逐节学习。

## 数据字典与生成方式

全部数据为无量纲教学数据，无缺失；脚本遇到缺失/非有限值会停止，不静默删除。没有真实观测、个人信息或实验结论。

| 输入文件 | 行数 | 字段与来源 |
| --- | --- | --- |
| `数据/position_examples.csv` | 3 | group为A/B/C；x1、x2、x3为人为编写的三个观测 |
| `数据/normal_samples.csv` | 20,000 | sample_id为模拟编号；x1至x10为每行10个独立N(10,4)观测 |
| `数据/shifted_exponential_samples.csv` | 20,000 | 同上；每个观测为8+Exponential(scale=2)，总体均值10、方差4 |
| `数据/joint_pairs.csv` | 1,200 | pair_id为编号；x与independent_y独立标准正态；sign独立等概率取−1、1 |
| `数据/generation_parameters.json` | — | 种子、NumPy版本、模型、单位、抽样次数等 |

`代码/generate_data.py`按NumPy default_rng/PCG64、种子20260930重新生成上述文件。仅主动改变模拟设计时运行；运行后要重跑三套分析和图，并同步正文数值。不同语言相同种子不保证同序列，所以这里统一生成CSV再跨语言读取。原始模拟输入保留17位有效数字。

## 结果字典

三种语言保存同名CSV：

- `position_results.csv`：原始人工数据、mean、variance。
- `covariance_identity.csv`：cov_MD、half_var_difference及两者的误差；max_pair_variance_error验证n=2的方差恒等式。这里的Cov是跨重复抽样结果的经验协方差，分母为重复次数减1。
- `normal_statistics.csv`、`shifted_exponential_statistics.csv`：每份样本的mean、variance；M_n2、D_n2、S2_n2只用该行前两个观测；Z、Q、T则用全行10个观测及已知总体μ=10、σ=2。
- `sampling_summary.csv`：均值和方差的模拟平均、经验协方差/相关、理论协方差；双侧t(9)界值2.26215716285410；t_rejection_fraction是20,000次模拟的拒绝比例，不是单次P值。两模型都在真实μ=10下计算，指数模型也套同一t界值作比较。rejection_mcse为二项比例的蒙特卡洛标准误√[p̂(1−p̂)/B]。
- `conditional_summary.csv`：all、mean_lt_9_5、mean_gt_10_5三种筛选下的入选数及S²均值、25%、50%、75%分位数。分位数统一为R type=7 / NumPy linear；MATLAB显式线性插值。
- `figure2_coordinates.csv`：图2完整坐标；independent_y_calculated为独立Y，square_y_calculated=X²，random_sign_y_calculated=sign×X。图形坐标有截断，数据未删除。
- `pair_summary.csv`：三种模型的经验Cov(X,Y)、相关系数及理论协方差0。
- `运行记录.txt`：问题、结果、解释、实际版本。

主要结果：正态/指数的经验Cov(X̄,S²)分别为0.01128634、1.62537230；理论分别为0、1.6。均值低于9.5时平均S²约3.994616、1.418619；高于10.5时约4.035819、8.046640。重复抽样拒绝比例为5.07%、9.97%，蒙特卡洛标准误分别约0.155、0.212个百分点。不要把9.97%推广成所有非正态数据的检验错误率。

## 理论核对：为什么指数模型的协方差是1.6？

对独立同分布、三阶绝对矩有限的样本，有Cov(X̄,S²)=μ₃/n，其中μ₃=E[(X−μ)³]。

一个简短推导：令Yᵢ=Xᵢ−μ，A=ΣYᵢ/n，则E[A]=0、S²=[ΣYᵢ²−nA²]/(n−1)。独立性及E[Yᵢ]=0使交叉项消失，得到E[AΣYᵢ²]=μ₃，以及E[A³]=μ₃/n²。因此Cov(A,S²)=[μ₃−n×μ₃/n²]/(n−1)=μ₃/n。

正态总体μ₃=0。尺度为θ的指数分布有E[W]=θ、E[W²]=2θ²、E[W³]=6θ³，故三阶中心矩为2θ³；θ=2时是16，平移8不改变中心矩。除以n=10得到1.6。这个协方差非零已足以证明此反例中均值与样本方差不独立；“协方差为零”的正态结论则还需正文给出的联合正态论证。

## 练习与科学边界

- 把筛选门槛改为9、11，先猜入选数量，再比较条件摘要的波动。
- 把所有人工数据整体加100，核对均值变化和方差不变。
- 区分经验不相关、理论不相关、相互独立、各自正态与联合正态。散点图和有限模拟均不能证明独立。
- 本课的独立性是在μ、σ²、n固定时的抽样分布性质，不是混合不同总体后的保证。
- 正文一般n的论证使用均值与整个偏差向量联合正态；偏差相加为0，所以这个联合分布可退化，仍可使用正态向量不同块间零协方差推出独立的性质。

## 阅读依据

- 正态样本与样本方差抽样分布：https://online.stat.psu.edu/stat414/Lesson26
- 联合正态定义与随机变号反例：https://www.probabilitycourse.com/chapter5/5_3_2_bivariate_normal_dist.php
- 联合正态的独立判定：https://statproofbook.github.io/P/mvn-ind.html
- t分布构造：https://stat.ethz.ch/R-manual/R-devel/library/stats/html/TDist.html

三语言实跑与数值核对见 `VALIDATION.md`；公众号稿为 `均值与方差的独立性_公众号稿_v2.md`。

## 封面与手机预览

`文章配图/00_第5次医学统计学课程封面_去姓名.png` 沿用第5次医学统计学课程原封面，从PPTX副本中删除姓名徽标组合后导出，未改动原课件。课程标题、校名、日期及图片布局保留；封面不是计算结果图。公开包仅含处理后的封面，不包含原始课件。

`文章配图/账号头像.png`用于手机预览的账号标识。预览和全部图片均为项目内相对路径，下载后可以本地打开HTML。重建预览：

```sh
python 代码/render_preview.py 均值与方差的独立性_公众号稿_v2.md 均值与方差的独立性_手机预览_v2.html --avatar 文章配图/账号头像.png --footer '统计课堂 · v2 · 本地审阅稿'
```

此脚本沿用公众号写作技能的渲染器，并为长目录名添加窄屏换行。只使用Python标准库。公众号上传/发布需在平台中另行完成。

当前稿件为v2：封面改为第5次课《抽样分布、χ² / t》（2026年9月28日），已去除姓名徽标；正文推导与计算不变。v1保留为历史稿件。
