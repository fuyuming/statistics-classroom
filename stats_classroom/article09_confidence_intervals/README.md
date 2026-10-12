# 统计课堂09｜置信区间的五个例题

**R、Python、MATLAB三种实现，使用同一套教材数据。** 按“题目—条件—输入—计算—解释”阅读；任选一种语言即可完成全部练习。

## 先下载完整项目

仓库ZIP：https://github.com/fuyuming/statistics-classroom/archive/refs/heads/main.zip

完整解压后进入 `stats_classroom/article09_confidence_intervals`。不要只下载一个脚本，否则会缺少共享CSV或R编号脚本。

## 怎样打开与运行

- **RStudio**：打开根目录 `统计课堂09.Rproj`，再打开 `代码/confidence_intervals.R`，点 **Source** 完整运行。只用R自带包。入口顺序调用01—06号脚本，编号脚本也可独立Source。Console显示计算，Plots显示图。
- **Spyder / Python**：打开 `代码/confidence_intervals.py`，按 **F5 / Run file** 整体运行。不要只选择语句或先运行中间单元格；脚本会自动根据自身位置找数据。需要numpy、pandas、scipy、matplotlib；若缺少包，在Spyder所用Python环境安装根目录requirements.txt。图保留在Plots/图形窗口。
- **MATLAB**：打开 `代码/confidence_intervals.m`，点 **Run** 完整运行。若提示更换当前文件夹，选择脚本所在“代码”目录即可，程序自动定位上级项目。需要 **Statistics and Machine Learning Toolbox**。Command Window显示过程，Figure窗口显示图。

结果自动保存到 `运行结果/R`、`运行结果/Python`、`运行结果/MATLAB`，互不覆盖。PNG是脚本生成的图，CSV是可检查的数值；无需先打开输出文件才能运行。

## 五个具体问题

| 例题 | 问题及方法 | 共享数据 | 95%区间（约） |
|---|---|---|---|
| 1 | 两批幼崽总体平均出生重差；独立正态、共同未知方差，合并t | birth_weight_g.csv：group为组别，weight_g为克 | −1.99635至6.09635 g |
| 2 | 总体客观缓解率；独立二项样本，Wald近似，并比较Wilson与精确二项 | response_counts.csv：events为缓解人数，n为总人数 | 0.352493至0.547507 |
| 3 | 低洼地减高坡地总体发病率差；独立二项、未合并标准误 | wheat_rust_counts.csv：group、events、n | 0.064531至0.164185 |
| 4 | 土壤铅浓度总体方差；独立正态模型，卡方反演 | lead_ug_g.csv：lead_ug_g为μg/g | 0.065495至0.620610 (μg/g)² |
| 5 | 品种1除以品种2的总体方差比；两组独立正态，F反演 | wheat_variance_summary.csv：n、sample_variance为样本量及样本方差 | 1.942118至31.479074 |

教材出处、单位转换及改编说明见 `数据/数据来源.md`。例题5仅有教材汇总值；没有虚构原始观测。CSV没有缺失值。样本方差均使用n−1分母；Python必须设置ddof=1。

例题2成功/失败为45、55；例题3四项为342、36、313、83，都满足本课采用的次数经验检查。总体率未知时，这是用样本率代入的检查，不是正态近似精度的保证。

## 正态性检查怎样读

例题4先画原始点图和正态Q-Q图。R与Python运行Shapiro–Wilk（W约0.95342，P约0.74561）。MATLAB自带工具箱没有同名检验，本项目明确使用 **Lilliefors** 辅助检查；两种检验的P值不是同一个量。MATLAB对超出查表范围的P值可能给出界限并警告，不能将界限当作精确P值。只有8个观测，未拒绝不等于证明总体正态；独立性还要根据采样设计判断。方差区间的精确性以正态模型为前提。

## 覆盖率模拟与改动练习

每轮从两个独立正态总体各抽10个观测，均数142、140 g，共同标准差2.5 g。重复1000轮，每轮重新计算合并t区间；CSV保存所有区间，图显示前40轮。

相同随机种子不保证三种语言得到同一批随机数，所以模拟覆盖比例不要求完全一致。五个确定性例题的端点应一致，核对见VALIDATION.md。

先核对例题1约[-2.00,6.10]g，再把置信水平改成99%，观察区间变宽。Python和MATLAB修改开头confidence；R编号脚本中将分位点0.975/0.025改成0.995/0.005，同时把内置函数conf.level和图表标签改为0.99/99%。避免只改标签而未改计算。

## 文件结构

- `代码/`：三种语言完整入口；R的六个编号脚本。
- `数据/`：共享CSV及教材来源。
- `文章配图/`：公众号正文用图；用于阅读，与脚本输出图区分。可运行 `代码/plot_article_figures.py` 复现公式图与教学图（先运行R主入口；需可用中文字体），封面除外。
- `运行结果/`：三种语言实际生成的数表、图形和控制台记录。
- `VALIDATION.md`：版本、确定性结果和模拟核对。
