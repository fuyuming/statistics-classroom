# 统计课堂10：秩和的抽样分布

两独立样本Wilcoxon秩和的教学补充。本文W是原始秩和，U=W-m(m+1)/2。不要与配对符号秩检验混用。

## 下载与运行

完整下载并解压仓库，进入 `stats_classroom/article10_rank_sum_distribution`，保留文件相对位置。

- R：打开 `统计课堂10.Rproj`，打开 `代码/rank_sum.R`，Source完整运行。只需基础R。
- Python：Spyder中“项目→新建项目→现有目录”，选择本篇目录，打开 `代码/rank_sum.py`，F5完整运行。主脚本自动调用 `plot_rank_sum.py`。在Spyder实际解释器中安装requirements.txt依赖。
- MATLAB：把Current Folder切到本篇目录，打开 `代码/rank_sum.m`，Run。只需基础MATLAB，不需统计工具箱。“设置路径”是函数搜索路径，不等于工作目录。

命令行备选：`Rscript 代码/rank_sum.R` 或 `python 代码/rank_sum.py`。

输出在 `运行结果/R`、`运行结果/Python`、`运行结果/MATLAB`。主脚本自动找到共享数据，每道题都在控制台打印数值与解释。

## 输入与输出

`数据/scores.csv` 为6行人工教学评分，没有真实实验或个人信息。group为A/B，value为无并列值例子，value_tied为将B的14改成12的例子。不含缺失值，不要自行把缺失当0。

各语言都会输出：

- `value_allocations.csv`：无并列值时20种记录分配及W。
- `value_pmf.csv`：精确概率表。
- `value_tied_allocations.csv`、`value_tied_pmf.csv`：平均秩下的条件置换分布。
- `theory.csv`：本篇8项理论结果。
- `large_simulation.csv`：m=n=20时20,000次名次分配模拟。
- `逐题结果.txt`：控制台解释。

R/MATLAB另有英文标签的正态近似核对图。Python生成正文三张中文图及标准正态曲线坐标；正文中的图1为示意，图2是穷举概率，图3是模拟密度。

## 核对与练习

无并列例：W=7、U=1、E(W)=10.5、Var(W)=5.25、双侧精确P=0.20。

并列例：W=7.5、Var(W)=5.10、双侧条件精确P=0.30。双侧按离零假设中心至少同样远定义。

全体分配就是完整概率空间，方差分母是20；模拟结果估计抽样方差时可用n-1分母，两者目的不同。

先逐行读20种分配，再修改模拟的m_large和n_large观察正态近似。穷举部分专门限制两组各3个，以防误把组合枚举用于很大样本造成卡顿。并列数据按观测记录选组，不把重复值记录删除。

R的wilcox.test两样本输出虽然标为W，数值已减去最小秩和，对应本文U。SciPy mannwhitneyu返回U；无并列精确算法不直接用于本篇并列例。脚本明确枚举并列记录的条件分布，避免软件版本差异。

## 中文

文件用UTF-8保存。代码乱码时按UTF-8重开；图片方框则是字体缺失。Python尝试微软雅黑、苹方、Noto Sans CJK SC等；找不到时先保存计算结果，再提示安装中文字体后重跑。R/MATLAB核对图用英文标签确保可移植。

`VALIDATION.md`记录实际运行结果。文章与HTML用于阅读审阅，未自动发布公众号。`课件阅读入口.md`提供可粘贴到课件的简短阅读说明，没有自动修改课件。
