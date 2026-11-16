# 统计课堂29｜独立样本率差的置信区间

先解压本篇ZIP，进入这个含README的项目目录，选择一种语言即可。题目是假设的两批独立样本：A组9/20阳性、B组13/20阳性，每人检测一次。研究A减B的总体率差，双侧95%置信水平。无个体数据、无缺失值；两行汇总数据不表示配对。

## 第一次怎样逐节运行

- **R / RStudio**：双击`统计课堂29_独立率差.Rproj`。在Console运行`install.packages("DescTools")`安装依赖（首次需要联网）。打开`代码/independent_ci.R`，从“准备与参数”依次选择小节，用Run或Ctrl/Cmd+Enter执行。看Console表格与Plots图；也可Source完整文件。
- **Python / Spyder**：将项目根目录设为Current working directory，打开`代码/independent_ci.py`，从第一节按顺序用Run current cell执行，或F5运行整文件。无需`__file__`也可分节。依赖安装应在Spyder正在使用的解释器环境终端执行`python -m pip install -r requirements.txt`；Console里的`import sys; print(sys.executable)`可查解释器。结果见Console、Plots/图窗口。
- **MATLAB**：将项目根目录设为Current Folder，打开`代码/independent_ci.m`，从“准备与参数”依次Run Section，也可Run全文。本文件只用基础MATLAB，无需统计工具箱。Current Folder决定相对文件位置；仅把代码加入搜索路径不能代替把工作目录设正确。

三套脚本也支持Current Folder为`代码`目录。不要跳过准备/读取节；后面的小节会用前面建立的变量。首次运行先看字段`group,positive,total`及A/B行，再看阳性/阴性人数、样本率、标准误和区间。

## 应核对什么

原始点估计为−0.20，即−20个百分点。默认95%结果如下，端点按比例报告：

| 方法 | 下限 | 上限 |
| --- | --- | --- |
| Wald | −0.5020507 | 0.1020507 |
| Newcombe，无连续性校正 | −0.4555079 | 0.1006339 |
| Agresti–Caffo | −0.4711246 | 0.1074882 |

最后几位以实际完整精度输出为准。表中的`estimate`是原始样本率差，`lower/upper`是端点；`lower_pp/upper_pp/width_pp`单位为百分点。AC每组加1阳性、1阴性，计算中心为−0.181818，而非原始差−0.20。图中黑点是原始差，橙色菱形是AC中心；横线为区间，竖虚线是零差。近似区间包含0不证明相等或等效。

R三种区间调用DescTools::BinomDiffCI，方法参数依次`wald / score / ac`，输出`est, lwr.ci, upr.ci`。Python调用statsmodels的confint_proportions_2indep，方法依次`wald / newcomb / agresti-caffo`，`compare="diff"`表示第一组减第二组，返回下限、上限。Python的`score`不是R的同名Newcombe构造。MATLAB逐步实现Wald、Wilson/Newcombe和AC公式，是项目附带的教学代码，未调用R后端，也不声称存在同名内置函数。

## 第二次改参数

计数只改`数据/independent_counts.csv`，不要在脚本后面额外改几个中间变量。置信水平只改所选脚本顶部的`confidence`。每次只做一项，然后从准备/读取节起顺序运行：

1. A改为45/100、B改为65/100，原始点估计仍为−20个百分点；区间变窄。这是假设扩大样本的练习，复制原有人的记录不会增加独立信息。
2. 恢复9/20、13/20，把confidence从0.95改0.90，区间通常会缩窄；观察标题和输出同步。
3. 两组均改为0/20。Wald退化成[0,0]，Newcombe约为[−0.16112516,0.16112516]，AC约为[−0.12309452,0.12309452]。比较完恢复默认CSV与0.95。

## 保存开关与现有结果

默认`SAVE_RESULTS`为FALSE / False / false，首次只在IDE查看。需要保存时，把本语言脚本开关改为真，结果分别写入`运行结果/R`、`运行结果/Python`、`运行结果/MATLAB`中的`intervals.csv`、`intervals.png`，不同语言互不覆盖。文件夹中的现有结果是默认题目的实际运行示例；改参数运行但不打开保存时，旧文件不会自动更新，应以当前Console和图窗口为准。

## 环境与来源

本次实际运行：R 4.6.1、DescTools 0.99.60；Python 3.12.12、numpy 2.4.3、pandas 3.0.1、scipy 1.17.1、statsmodels 0.14.6、matplotlib 3.10.8；MATLAB R2025a。requirements列出所需包的兼容范围，不要求切换全局环境。

- R官方说明：https://andrisignorell.github.io/DescTools/reference/BinomDiffCI.html
- Python官方说明：https://www.statsmodels.org/stable/generated/statsmodels.stats.proportion.confint_proportions_2indep.html
- Newcombe1998独立样本方法10：https://doi.org/10.1002/%28SICI%291097-0258%2819980430%2917%3A8%3C873%3A%3AAID-SIM779%3E3.0.CO%3B2-I
- Agresti–Caffo2000：https://doi.org/10.1080/00031305.2000.10474560

本篇代码目录：https://github.com/fuyuming/statistics-classroom/tree/main/stats_classroom/article29_independent_proportion_ci

本篇ZIP：https://raw.githubusercontent.com/fuyuming/statistics-classroom/main/stats_classroom/article29_independent_proportion_ci/independent_proportion_ci_project.zip

## 对应文章

[公众号源稿](文章/独立样本率差_公众号稿.md) · [手机预览HTML](文章/独立样本率差_手机预览.html)
