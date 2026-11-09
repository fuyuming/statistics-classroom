# 配对样本率差：文章与三语言练习

本目录包括公众号文章、配图、共享数据和R/Python/MATLAB配套。

- 阅读文章：打开[手机预览](配对样本率差_手机预览_v1.html)，或阅读[Markdown正文](配对样本率差_公众号稿_v1.md)。下载ZIP后可在本地浏览器查看完整预览。
- 开始练习：先读[配套项目README](配套项目/README.md)，选择一种语言运行。
- 正文统计图：用R打开`配套项目/配对率差.Rproj`，Source `代码/paired_rate.R`，自动生成三张图；也可单独Source `代码/plot_figures.R`。
- 三套绘图入口：`plot_figures.R`、`plot_figures.py`、`plot_figures.m`；三套主入口均自动调用本语言绘图，Python和MATLAB结果各自保存。
- R精确区间：`配套项目/代码/paired_intervals.R`调用PairedCI；三语言另有共同的固定D排序精确算法，二者构造及结果分别说明。
- 完整下载：[本篇ZIP](https://github.com/fuyuming/statistics-classroom/raw/refs/heads/main/stats_classroom/paired_proportion_difference/%E9%85%8D%E5%AF%B9%E6%A0%B7%E6%9C%AC%E7%8E%87%E5%B7%AE_%E5%AE%8C%E6%95%B4%E9%A1%B9%E7%9B%AE.zip)。

完整脚本、数据定义、依赖、对应输出与运行结果见配套目录。数据均为教学构造。

- 逐项概率：[R计算的完整231组表](配套项目/运行结果/R/三类人数的231组概率.csv)，每行包含三类人数、人员安排种数、概率及尾部归属；三语言均有对应表。
