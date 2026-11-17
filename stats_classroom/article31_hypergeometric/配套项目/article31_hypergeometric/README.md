# 统计课堂31配套项目：超几何分布与 HBV Fisher 精确检验

例题来源：颜艳/王彤《医学统计学》第5版例7-4，p.111。数据为预防组 4/22、对照组 5/11。

## 运行方式

- **RStudio**：打开 `article31_hypergeometric.Rproj`，打开脚本后可用 Source，或按带 `----` 的节标题顺序运行；先执行 `install.packages("ggplot2")`，完整运行：`Rscript article31_hypergeometric.R`。
- **Python/Spyder**：把项目根设为工作目录，在Spyder的IPython控制台执行 `%pip install -r requirements.txt`（或在对应解释器终端执行 `python -m pip install -r requirements.txt`），打开脚本按 `#%%` 顺序运行或按 F5 全部运行：`python article31_hypergeometric.py`。图用 seaborn + Matplotlib，默认不保存。
- **MATLAB R2025a**：把项目根设为 Current Folder，打开脚本按 `%%` 分节运行或运行整文件；需要 Statistics and Machine Learning Toolbox（`hygepdf`、`fishertest`）。

三种脚本都从项目根的相对路径 `data/hbv_table.csv` 读入数据；整文件从命令行运行时会定位脚本所在项目根；IDE逐节运行请先打开本项目并保持项目根为工作目录。脚本默认显示图。保存开关 `SAVE_OUTPUT` 默认为关闭，R改为`TRUE`、Python改为`True`、MATLAB改为`true` 后分别保存到 `results/r`、`results/python`、`results/matlab`。

## 例题计算链

`N=33, K=9, n=22`；`P(X=x)=choose(9,x) choose(24,22-x) / choose(33,22)`。脚本展示 x=0..9 的完整概率表及和、加权均值与方差、`P(X<=4)`、`P(X>=4)`，并按 R `fisher.test` 的“概率不大于观察表”规则核对双侧 P。双侧只报告 P，不讲 OR 区间。

## 练习

把抽样数改为 `n=11`，重算完整分布，并检查补集：抽取的11人与留下的22人阳性数之和为9，期望分别为3与6，方差均为1.5。这里始终在给定总阳性数9的条件分布中计算。


## 预期结果

观察表概率0.087627278172328；完整分布总概率1、均值6、方差1.5；左尾0.108147229114971；右尾0.979480049057357；概率排序双侧0.121044750577565。三语言默认的完整检验函数在本题给出同一双侧P。`results/`按语言保存本题的概率表、矩、函数核对及图形，可与自己的运行结果比较。

R Source后在Environment查看`probability_table`、`moments`和`fisher_check`；按节运行第7节可直接查看图。Python执行`plt.show()`显示图，MATLAB在图形窗口查看。Current Folder决定相对数据路径，加入搜索路径主要让MATLAB找到代码，二者作用不同。

## 大样本对照例入口

`data/large_sample_table.csv` 与 `large_sample_comparison.R/.py/.m` 对应本课 200 人例（新方案 50/100、对照 30/100）。三脚本均按组名读取 CSV，展示率差、合并率下期望表与 `npq`、手算 Z/正态双侧 P、近似计算及函数核对、超几何概率排序双侧 P 和 Fisher 精确 P；不含置信区间或新图。分别在 R/Python/MATLAB 的 `SAVE_OUTPUT` 开关打开后，将两份结果 CSV 写入对应 `results/<language>/`。

题干假设两组独立随机分组、各100人，4周用同一标准评价达标50/30；目标为总体达标率差，双侧 H0 为两率相等，显著性水平 alpha=0.05。脚本中的命名结果表逐项说明：合并率用于 H0 下期望数和 SE0，`npq` 用于近似条件检查，正态分布函数给出 Z 双侧尾概率；固定边际超几何 `A` 的概率按不大于观察表概率排序后求双侧 P，最后再用完整 Fisher 函数核对。

MATLAB对照脚本用`normcdf`计算Z的双侧尾概率；它是辅助分布函数。R用`prop.test(correct=FALSE)`、Python用`chi2_contingency(correction=False)`核对同一未校正近似P。三语均用完整Fisher检验函数核对条件精确P。
