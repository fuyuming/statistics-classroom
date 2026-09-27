# 统计课堂07：正态、χ²、t、F 的关系

公众号：明哥的微生物世界。本文从抽样分布解释小样本精确推断、大样本近似及常见检验用途。

## 下载与运行

本篇：https://github.com/fuyuming/statistics-classroom/tree/main/stats_classroom/article07_normal_chisq_t_f

全仓库 ZIP：https://github.com/fuyuming/statistics-classroom/archive/refs/heads/main.zip

解压后进入 `stats_classroom/article07_normal_chisq_t_f`，用 RStudio 打开 `代码/normal_chisq_t_f.R`，点击 Source 完整运行。只需基础 R，无额外包。命令行也可运行：

```sh
Rscript 代码/normal_chisq_t_f.R
```

脚本通过自身路径定位目录。逐行运行时须先将工作目录设为本篇目录。脚本及文字文件为 UTF-8；图中使用英文标签，正文提供中文解释，不需要额外中文字体。

## 文件与来源

- `正态卡方tF的关系_公众号稿.md`：完整文章及显式来源链接。
- `代码/normal_chisq_t_f.R`：全部计算、模拟、配图和恒等式检查。
- `数据/teaching_data.csv`：12 行人为编写的教学数据，**不是实验观测**。`sample_id` 为虚构编号，`group` 为 A/B 两组，`x` 为无量纲教学自变量，`response` 为任意单位的教学响应。用于分别演示组别回归和连续自变量回归，不将两个模型解释为同一真实实验的效应分析，不由这些数值判定正态假设成立。
- `运行结果/critical_values.csv`：理论分位点和错用 1.96 的理论拒绝概率；`n_one_sample` 仅对应单样本 df=n−1 的设定。
- `运行结果/t_normal_curves.csv`：正文图的完整横坐标与密度值。
- `运行结果/t_anova_regression_equivalence.csv`：合并方差双侧 t 检验、经典方差分析、组别回归的等价核对。
- `运行结果/simple_regression_equivalence.csv`：含截距一元回归斜率 t² 与整体 F 核对。
- `运行结果/construction_draws.csv.gz`：50,000 行模型模拟变量；可用 `read.csv(gzfile("运行结果/construction_draws.csv.gz"))` 读取。
- `运行结果/construction_simulation_quantiles.csv`：模拟分位点与理论值对照；允许抽样误差。
- `运行结果/逐项计算结果.txt`：完整运行输出与 R 环境记录。
- `文章配图/02-t分布与标准正态比较.png`：脚本生成的文章配图。

本文涉及的广义线性模型仅介绍分布与假设边界，没有实际 GLM 拟合案例；不存在遗漏的 GLM 原始数据或分析脚本。

## 核对结果

已在 R 4.6.1 完整运行：t(10) 的 0.975 分位点为 2.22813885；平方与 F(1,10) 的 0.95 分位点均为 4.96460274；错用 ±1.96 的理论第一类错误率为 7.843624%。两组教学数据的 t² 与 F 均为 19.7387606；连续自变量回归的斜率 t² 与模型 F 均为 119.394421。

脚本包含自动等价检查。模拟采用种子 20260927，理论计算不依赖随机数。图和 CSV 均为模型计算或教学结果，不是真实微生物实验数据。
