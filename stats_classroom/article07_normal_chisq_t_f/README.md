# 统计课堂07：正态、χ²、t、F 的关系

公众号：明哥的微生物世界。本文从抽样分布解释小样本精确推断、大样本近似及常见检验用途。

## 下载与运行

本篇：https://github.com/fuyuming/statistics-classroom/tree/main/stats_classroom/article07_normal_chisq_t_f

全仓库 ZIP：https://github.com/fuyuming/statistics-classroom/archive/refs/heads/main.zip

完整解压后进入 `stats_classroom/article07_normal_chisq_t_f`，不要只下载单个脚本。三套语言主文件同名，扩展名分别为 `.R`、`.py`、`.m`，共享 `数据/teaching_data.csv`。

### R / RStudio

双击本目录的 **`统计课堂07.Rproj`**，再打开 `代码/normal_chisq_t_f.R`，点击 **Source** 完整运行。只需基础 R，无额外包。也可在项目控制台运行：

```r
source("代码/normal_chisq_t_f.R", encoding = "UTF-8")
```

命令行：`Rscript 代码/normal_chisq_t_f.R`。项目设置不恢复旧工作空间，避免依赖上次会话的变量。

### Python / Spyder

Spyder 中选择 **Projects → New Project → Existing directory**，选择本篇目录；打开 `代码/normal_chisq_t_f.py`，按 F5 完整运行。工作目录应是含 `数据`、`代码`、`requirements.txt` 的这一层。脚本本身也通过 `__file__` 定位数据。

需要 NumPy、SciPy、Matplotlib。在实际运行脚本的 Python 环境中安装：

```sh
python -m pip install -r requirements.txt
python 代码/normal_chisq_t_f.py
```

若 Spyder 提示缺包，先用 `import sys; print(sys.executable)` 确认其解释器，不能仅凭另一个终端安装成功就判断环境已就绪。

### MATLAB

需要 **Statistics and Machine Learning Toolbox**。将 Current Folder 切到本篇目录，打开 `代码/normal_chisq_t_f.m`，点击 **Run**；或输入：

```matlab
run('代码/normal_chisq_t_f.m')
```

“当前文件夹”与“搜索路径”不同；无需为本课全局添加路径。完整运行脚本会根据自身位置定位项目目录。

### 输出与编码

R 的结果位于 `运行结果/`，正文配图位于 `文章配图/`；Python、MATLAB 各输出到这两个目录下同名的语言子目录，避免互相覆盖。重新运行会覆盖对应计算结果，不修改输入数据。

脚本及文字文件使用 UTF-8。图中使用英文标签，正文配中文图注，无须安装中文字体。源文件乱码应检查编码；更换字体不能修复读错的文本。三套脚本均计算分位点、生成三张分布图和曲线 CSV、验证 t²=F、保存模拟数据及运行记录。数值定义和参数相同，图形样式允许略有差异。

跨语言使用相同种子并不产生同一条随机数序列；模拟分位点只需接近理论值，不应逐行比较随机抽样。理论计算与共用教学数据的分析结果应在数值误差范围内一致。

## 文件与来源

- `正态卡方tF的关系_公众号稿.md`：完整文章及显式来源链接。
- `统计课堂07.Rproj`：RStudio 项目入口。
- `requirements.txt`：Python 依赖。
- `代码/normal_chisq_t_f.R`、`.py`、`.m`：三语言完整计算、模拟、配图和恒等式检查。
- `VALIDATION.md`：实际运行环境、跨语言核对方法与结果。
- `数据/teaching_data.csv`：12 行人为编写的教学数据，**不是实验观测**。`sample_id` 为虚构编号，`group` 为 A/B 两组，`x` 为无量纲教学自变量，`response` 为任意单位的教学响应。用于分别演示组别回归和连续自变量回归，不将两个模型解释为同一真实实验的效应分析，不由这些数值判定正态假设成立。
- `运行结果/critical_values.csv`：理论分位点和错用 1.96 的理论拒绝概率；`n_one_sample` 仅对应单样本 df=n−1 的设定。
- `运行结果/t_normal_curves.csv`：t 与正态比较图的完整横坐标与密度值。
- `运行结果/chisq_curves.csv`：χ² 图的完整坐标，df=1、2、3、5、10、30。
- `运行结果/f_curves.csv`：F 图的完整坐标，分别固定分子或分母自由度。
- `运行结果/t_anova_regression_equivalence.csv`：合并方差双侧 t 检验、经典方差分析、组别回归的等价核对。
- `运行结果/simple_regression_equivalence.csv`：含截距一元回归斜率 t² 与整体 F 核对。
- `运行结果/construction_draws.csv.gz`：50,000 行模型模拟变量；可用 `read.csv(gzfile("运行结果/construction_draws.csv.gz"))` 读取。
- `运行结果/construction_simulation_quantiles.csv`：模拟分位点与理论值对照；允许抽样误差。
- `运行结果/逐项计算结果.txt`：完整运行输出与 R 环境记录。
- `文章配图/02-t分布与标准正态比较.png`：脚本生成的 t 与正态比较图。
- `文章配图/03-卡方分布是一族曲线.png`：卡方分布的两面板图；df=1 零点附近高峰截断，图注明示。
- `文章配图/04-F分布的两个自由度.png`：F 分布的两面板图；竖线 1 不表示期望值。

本文涉及的广义线性模型仅介绍分布与假设边界，没有实际 GLM 拟合案例；不存在遗漏的 GLM 原始数据或分析脚本。

## 核对结果

已在 R 4.6.1 完整运行：t(10) 的 0.975 分位点为 2.22813885；平方与 F(1,10) 的 0.95 分位点均为 4.96460274；错用 ±1.96 的理论第一类错误率为 7.843624%。两组教学数据的 t² 与 F 均为 19.7387606；连续自变量回归的斜率 t² 与模型 F 均为 119.394421。

脚本包含自动等价检查。模拟采用种子 20260927，理论计算不依赖随机数。图和 CSV 均为模型计算或教学结果，不是真实微生物实验数据。
