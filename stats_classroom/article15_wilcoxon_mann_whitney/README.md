# 统计课堂15｜Wilcoxon 与 Mann–Whitney 配套练习

本篇项目：https://github.com/fuyuming/statistics-classroom/tree/main/stats_classroom/article15_wilcoxon_mann_whitney

完整下载：https://github.com/fuyuming/statistics-classroom/archive/refs/heads/main.zip

下载后完整解压，进入 `statistics-classroom-main/stats_classroom/article15_wilcoxon_mann_whitney/`。三种语言任选一种学习，输入数据和核对结果共用同一套。

| 语言 | 主脚本 |
| --- | --- |
| R | [code/wilcoxon_mannwhitney.R](code/wilcoxon_mannwhitney.R) |
| Python | [code/wilcoxon_mannwhitney.py](code/wilcoxon_mannwhitney.py) |
| MATLAB | [code/wilcoxon_mannwhitney.m](code/wilcoxon_mannwhitney.m) |

先打开 `Wilcoxon与MannWhitney.Rproj`；RStudio 运行 `code/wilcoxon_mannwhitney.R`，Spyder 打开并运行 `code/wilcoxon_mannwhitney.py`，MATLAB 打开 `code/wilcoxon_mannwhitney.m` 后运行。工作目录可设为项目根或 `code/`；三种运行方式均从共享 `data/` 读取。

三个脚本使用同一份 `data/textbook_example.csv`，复算教材例5-8，并用 `data/independent_RD.csv` 做实际课件 RD 数据的符号对照；第五节追加 `paired_vitaminA.csv` 的单样本/配对符号秩等价演示。脚本按“看数据→明问题→立假设”“选方法→查条件”“做检验→判显著→看效应”“讲结果→画图表→写结论”组织，末节仅可选导出。

Python 可整文件运行；也可在 Spyder 中从首节开始逐段运行，代码不依赖 `__file__`。R 需要基础 `stats`（随R安装）；Python 需要 `numpy`、`pandas`、`scipy`（见 `requirements.txt`）；MATLAB 需要 Statistics and Machine Learning Toolbox 的 `ranksum`/`tiedrank`。缺少数据时会给出清晰错误。三者保留各自软件对统计量的命名差异。

教材例5-8的方向为 A>B：A秩和37，B秩和18，UA=22，UB=3，胜出比例22/25=0.88；252个分配中 U_A≥22 的7个，单侧 P=7/252=0.0277778，双侧按两端计数为0.0555556。下限临界值 U=4，P(U≤4)=12/252，跨组差第5个有序值为2。R 的 W=22、HL=5、单侧位置区间[2, Inf)；旧版 R 可能显示95%标称置信水平，离散实际水平为20/21=0.952381。RD 对照只报告秩和、U、胜出比例与近似+连续性校正P，不追加置信区间重分析。

双侧枚举统一写作 `min(1, 2*min(P(U≥u), P(U≤u)))`；交换输入后改用等价的反方向检验，单侧P仍为0.0277778，双侧仍为0.0555556。可选末节分别写出 `results/R_results.csv`、`results/Python_results.csv`、`results/MATLAB_results.csv`，三者共享关键字段与实际CI置信水平。

配对数据的 `d=round(X-Y,2)` 有8个非零差、无绝对并列，正秩和35、负秩和1。R 的差值单样本和 `paired=TRUE` 都返回 V=35、P=0.015625、CI [0.275,1.200]，实际置信水平96.1%；二者是同一配对问题的等价写法，不是两个独立研究。Python/MATLAB只核对V和P，原样R输出与源代码在 `R运行输出/`。

脚本不使用 `print`、`cat`、`disp` 或 `fprintf` 叙述；输出以裸对象、Series、data.frame 或 table 为主，解释写在代码旁注。运行记录见 `VALIDATION.md`。


RD逐项编秩演示：`R运行输出/04-RD.R` 列出22个RD值的合并秩，再按组求秩和、换算U并执行含连续性校正的单侧近似检验。以 `R运行输出` 为工作目录运行；对应原始日志为 `04-RD.txt`。

RD的10×12表示A组10人与B组12人逐一跨组比较，共120对：A较大86对，B较大33对，相等1对；相等各计半分，因此UA=86.5、UB=33.5。
