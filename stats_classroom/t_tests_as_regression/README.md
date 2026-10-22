# 统计课堂｜t检验，原来也是回归

本项目把三个常见问题放在同一条线上：只有截距的模型对应单样本 t；0/1 组变量的普通 OLS 对应两独立样本的合并方差 t；配对设计先按 `id` 计算差值，再把差值交给只有截距的模型。三套脚本使用同一份 CSV 和同一套统计定义，任选 R、Python 或 MATLAB 即可运行。

## 先下载、解压、打开

仓库 ZIP：<https://github.com/fuyuming/statistics-classroom/archive/refs/heads/main.zip>

解压后进入 `stats_classroom/t_tests_as_regression`。打开 `统计课堂_t检验与回归.Rproj`，再运行 `代码/t_tests_as_regression.R`；或者打开 `代码/t_tests_as_regression.py` 用 Spyder 的 F5；也可以在 MATLAB 中打开并运行 `代码/t_tests_as_regression.m`。本地还没有把“运行结果”当作原始数据提交，脚本会在运行时重新生成它们。

## 运行方式

- **R / RStudio**：R 自带 `lm()`、`t.test()` 和绘图函数，模型计算无需额外包；报告演示需要先运行 `install.packages("report")`。打开本项目 Rproj 后 Source 完整的 `代码/t_tests_as_regression.R`。项目设置为 UTF-8，不恢复旧工作空间。脚本最后留下 `results` 结构对象，并把表格写入 `运行结果/R/`。
- **Python / Spyder**：按 F5 完整运行 `代码/t_tests_as_regression.py`；命令行是备选。需要 `numpy`、`pandas`、`scipy`、`statsmodels`、`matplotlib`，版本下限见 `requirements.txt`。脚本用 statsmodels 的 OLS 公式接口和 SciPy 的 t 检验核对，结果写入 `运行结果/Python/`。
- **MATLAB**：运行 `代码/t_tests_as_regression.m`，需要 Statistics and Machine Learning Toolbox（`fitlm`、`ttest`、`ttest2`）。脚本从自身位置定位项目，不要求手动把 `数据/` 加入搜索路径；表格写入 `运行结果/MATLAB/`。

三套脚本的主流程保留模型或结果对象，方便在 RStudio、Spyder 和 MATLAB 桌面端逐项查看；文字解释放在注释和本 README，CSV 是跨语言核对接口。代码没有用矩阵求逆来讲回归。R 中若希望进一步整理模型结果，可将 `lm` 对象交给 `broom::tidy()` 或 easystats 系列工具；整理结果仍需回看原始模型、设计和区间，不能把自动整理当成准确的中文草稿。

## 共享数据

`数据/rat_gain_long.csv` 是长表，一行一只大鼠：1月龄开始用高、低蛋白饲料饲养，3月龄测定增重量；A 高蛋白组 12 只，B 低蛋白组 7 只，`gain_g` 的单位为 g，`group01=1` 表示 A、`0` 表示 B。数值来自李春喜《生物统计学》第6版例4-8，pp.97–98。把 B 组均数与 100 g 比较是本项目提出的教学问题，不冒充原题。

`数据/paired_drinks.csv` 是按消费者排列的宽表，一行一人。`old`、`new` 是本次课件例2中旧款、新款饮料的评分；脚本按同一 `id` 计算 `delta = new - old`。它用于说明配对结构，不把两款评分当作独立观测。两张表均无缺失值，数据字典和来源见 `数据/数据来源.md`。

## 预期关键结果

以下结果使用普通同方差 OLS 推断，数值保留足够位数以便核对；不同语言的末位舍入可有微小差异。

| 问题 | 估计 / 自由度 | t / F | P值 | 95%区间 |
| --- | ---: | ---: | ---: | ---: |
| B组均数 − 100 g（只有截距） | 1.0000；df=6 | t=0.1283 | 0.9021 | [−18.0737, 20.0737] g |
| A−B（group01 的 OLS 系数） | 19.1667；df=17 | t=1.9157 | 0.07238 | [−1.9425, 40.2759] g |
| 两组模型整体检验 | df1=1，df2=17 | F=3.6698 | 0.07238 | F=t² |
| Welch（仅作口径对照） | df≈13.0159 | t=1.9319 | 0.07543 | [−2.2637, 40.5970] g |
| 新款−旧款（配对差值） | 1.3000；df=9 | t=2.7508 | 0.02245 | [0.2309, 2.3691] 分 |

`lm(y ~ 1)`、`smf.ols("y ~ 1", data=dat).fit()` 和 `fitlm(tbl,'y ~ 1')` 的截距都默认检验参考值 0。若目标是 `μ₀=100`，先用 `I(y - μ₀) ~ 1`（或等价的差值列）再看截距。两独立样本的 OLS 是合并方差推断；它不等于 Welch。两组时整体 F 与该单一组系数的 t² 相等，这个 F 检验均数模型，不是比较两个总体方差的 F 检验。

## 输出与图

每种语言分别写入自己的 `运行结果/<语言>/`，包括单样本、组系数、整体 F、合并方差与 Welch 对照、配对差值和 `all_key_results.csv`。两张图的生成参数在三套主脚本中相同：`文章配图/01_组均数与OLS截距.png` 画出 B 截距、A−B 组差与原始点；`文章配图/02_新旧款配对与差值.png` 画出每位消费者的新旧款连线及均数差。Python 入口按公众号阅读需要生成中文大字图；R、MATLAB 的图在缺少中文字体时仍可用英文标签生成，正文图注负责解释中文含义。

练习时可以先核对 CSV，再改一个明确的问题：把单样本参照值从 100 g 改为 110 g，或把配对数据中的一位消费者评分改回原始值后完整运行。不要通过挑选极端值来追求显著。

## 官方接口

- R `lm`：<https://stat.ethz.ch/R-manual/R-devel/library/stats/html/lm.html>
- R `t.test`：<https://stat.ethz.ch/R-manual/R-devel/library/stats/html/t.test.html>
- statsmodels `OLS`：<https://www.statsmodels.org/stable/generated/statsmodels.regression.linear_model.OLS.html>
- MATLAB `fitlm`：<https://www.mathworks.com/help/stats/fitlm.html>

## 中文与编码

脚本和数据使用 UTF-8。Python 绘图会查找系统中文字体（微软雅黑、黑体、Noto Sans CJK SC、苹方等）；若找不到会给出明确提示，避免生成方框字。安装任一上述字体后重启 Python。请保持整个项目解压后的目录结构，输出均留在本项目内。

## 2026-10-22 学生版补充

R第2.1节实际运行report并保存报告表和文字初稿；Python单样本和配对统一使用公式接口。步骤、参数含义、预期数字和思考题写在注释中。RStudio可选中一段运行；Spyder用单元格运行并在变量浏览器查看对象，F5完整运行时不会自动显示每个裸表达式。MATLAB注释使用 `%`。
